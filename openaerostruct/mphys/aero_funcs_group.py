import openmdao.api as om

from openaerostruct.aerodynamics.functionals import VLMFunctionals
from openaerostruct.functionals.total_aero_performance import TotalAeroPerformance
from openaerostruct.mphys.surface_contours import SurfaceContour
from openaerostruct.mphys.lift_distribution import LiftDistribution


class AeroFuncsGroup(om.Group):
    """
    Group to contain the total aerodynamic performance functions
    to be evaluated after the coupled states are solved.
    """

    def initialize(self):
        self.options.declare("surfaces", default=None, desc="oas surface dicts", recordable=False)
        self.options.declare("user_specified_Sref", types=bool, default=False)
        self.options.declare("write_solution", types=bool, default=True)
        self.options.declare("output_dir")
        self.options.declare("scenario_name", default=None)

    def setup(self):
        self.surfaces = self.options["surfaces"]
        self.user_specified_Sref = self.options["user_specified_Sref"]

        self.set_input_defaults("aoa", val=0.0, units="deg")
        self.set_input_defaults("yaw", val=0.0, units="deg")
        self.set_input_defaults("mach", val=0.0)

        proms_in = []
        for surface in self.surfaces:
            surf_name = surface["name"]
            self.add_subsystem(
                surf_name,
                VLMFunctionals(surface=surface),
                promotes_inputs=[
                    "v",
                    ("alpha", "aoa"),
                    ("beta", "yaw"),
                    ("Mach_number", "mach"),
                    ("re", "reynolds"),
                    "rho",
                ],
            )

            proms_in.append((surf_name + "_S_ref", surf_name + ".S_ref"))
            proms_in.append((surf_name + "_b_pts", surf_name + ".b_pts"))
            proms_in.append((surf_name + "_widths", surf_name + ".widths"))
            proms_in.append((surf_name + "_chords", surf_name + ".chords"))
            proms_in.append((surf_name + "_sec_forces", surf_name + ".sec_forces"))
            proms_in.append((surf_name + "_CL", surf_name + ".CL"))
            proms_in.append((surf_name + "_CD", surf_name + ".CD"))

        proms_out = ["CM", "CL", "CD", "L", "D"]
        if self.options["user_specified_Sref"]:
            proms_in.append("S_ref_total")
        else:
            proms_out.append("S_ref_total")

        # Add the total aero performance group to compute the CL, CD, and CM
        # of the total aircraft. This accounts for all lifting surfaces.
        self.add_subsystem(
            "total_perf",
            TotalAeroPerformance(surfaces=self.surfaces, user_specified_Sref=self.user_specified_Sref),
            promotes_inputs=proms_in + ["v", "rho", "cg"],
            promotes_outputs=proms_out,
        )

        proms_in = []
        for surface in self.surfaces:
            surf_name = surface["name"]
            proms_in.append((surf_name + "_sec_forces", surf_name + ".sec_forces"))

        if self.options["write_solution"]:
            sol_writer = self.add_subsystem("solution_writer", om.Group(), promotes=["*"])
            sol_writer.add_subsystem(
                "contour_writer",
                SurfaceContour(
                    surfaces=self.surfaces,
                    base_name=self.options["scenario_name"],
                    output_dir=self.options["output_dir"],
                ),
                promotes_inputs=proms_in + ["*"],
            )
            sol_writer.add_subsystem(
                "distribution_writer",
                LiftDistribution(
                    surfaces=self.surfaces,
                    base_name=self.options["scenario_name"],
                    output_dir=self.options["output_dir"],
                ),
                promotes_inputs=proms_in + ["*"],
            )
