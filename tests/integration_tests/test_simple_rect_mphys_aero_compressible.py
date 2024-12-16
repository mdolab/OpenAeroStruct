import unittest

import numpy as np

import openmdao.api as om
from openmdao.utils.assert_utils import assert_near_equal, assert_check_partials

from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.geometry.geometry_group import Geometry
from openaerostruct.aerodynamics.aero_groups import AeroPoint
from openaerostruct.mphys import AeroBuilder

# check if mphys is available
try:
    from mphys.core import Multipoint, MPhysVariables  # noqa: F401
    from mphys.scenarios import ScenarioAerodynamic  # noqa: F401

    mphys_flag = True
except ModuleNotFoundError:
    mphys_flag = False

# Global flight condition inputs to be used with both oas and mphys
mach = 0.84
aoa = 5.0
beta = 0.0
rho = 0.38
vel = 248.136
re = 1e6
cg = np.zeros((3))

FlowVars = MPhysVariables.Aerodynamics.FlowConditions


@unittest.skipUnless(mphys_flag, "MPhys is required.")
class Test(unittest.TestCase):
    def test(self):
        # Create a dictionary to store options about the surface
        mesh_dict = {"num_y": 5, "num_x": 2, "wing_type": "rect", "symmetry": True}

        mesh = generate_mesh(mesh_dict)

        surf_dict = {
            # Wing definition
            "name": "wing",  # name of the surface
            "type": "aero",
            "symmetry": True,  # if true, model one half of wing
            # reflected across the plane y = 0
            "S_ref_type": "wetted",  # how we compute the wing area,
            # can be 'wetted' or 'projected'
            "twist_cp": np.array([0.0]),
            "mesh": mesh,
            # Aerodynamic performance of the lifting surface at
            # an angle of attack of 0 (alpha=0).
            # These CL0 and CD0 values are added to the CL and CD
            # obtained from aerodynamic analysis of the surface to get
            # the total CL and CD.
            # These CL0 and CD0 values do not vary wrt alpha.
            "CL0": 0.0,  # CL of the surface at alpha=0
            "CD0": 0.015,  # CD of the surface at alpha=0
            # Airfoil properties for viscous drag calculation
            "k_lam": 0.05,  # percentage of chord with laminar
            # flow, used for viscous drag
            "t_over_c_cp": np.array([0.15]),  # thickness over chord ratio (NACA0015)
            "c_max_t": 0.303,  # chordwise location of maximum (NACA0015)
            # thickness
            "with_viscous": False,
            "with_wave": False,  # if true, compute wave drag
        }

        surfaces = [surf_dict]

        oas_sol = self.get_oas_solution(surfaces)

        mphys_sol = self.get_mphys_solution(surfaces)

        assert_near_equal(oas_sol["aero_point_0.wing_perf.CD"][0], mphys_sol["aero_point_0.wing.CD"][0], 1e-6)
        assert_near_equal(oas_sol["aero_point_0.wing_perf.CL"][0], mphys_sol["aero_point_0.wing.CL"][0], 1e-6)
        assert_near_equal(oas_sol["aero_point_0.CM"][0], mphys_sol["aero_point_0.CM"][0], 1e-6)
        assert_near_equal(oas_sol["aero_point_0.CM"][1], mphys_sol["aero_point_0.CM"][1], 1e-6)
        assert_near_equal(oas_sol["aero_point_0.CM"][2], mphys_sol["aero_point_0.CM"][2], 1e-6)

        # Check the partials at this point in the design space
        data = mphys_sol.check_partials(compact_print=True, out_stream=None, method="cs", step=1e-40)
        assert_check_partials(data, atol=1e20, rtol=1e-6)

    def get_mphys_solution(self, surfaces):
        class Top(Multipoint):
            def setup(self):
                dvs = self.add_subsystem("dvs", om.IndepVarComp(), promotes=["*"])
                dvs.add_output(FlowVars.ANGLE_OF_ATTACK, val=aoa, units="deg")
                dvs.add_output(FlowVars.YAW_ANGLE, val=beta, units="deg")
                dvs.add_output("rho", val=rho, units="kg/m**3")
                dvs.add_output(FlowVars.MACH_NUMBER, mach)
                dvs.add_output("v", vel, units="m/s")
                dvs.add_output(FlowVars.REYNOLDS_NUMBER, re, units="1/m")
                dvs.add_output("cg", val=cg, units="m")

                # Create mphys builder for aero solver
                options = {"compressible": True, "write_solution": False}
                aero_builder = AeroBuilder(surfaces, options)
                aero_builder.initialize(self.comm)

                # Create mesh component and connect with solver
                self.add_subsystem("mesh", aero_builder.get_mesh_coordinate_subsystem())
                self.mphys_add_scenario("aero_point_0", ScenarioAerodynamic(aero_builder=aero_builder))
                self.connect(
                    f"mesh.{MPhysVariables.Aerodynamics.Surface.Mesh.COORDINATES}",
                    f"aero_point_0.{MPhysVariables.Aerodynamics.Surface.COORDINATES}",
                )

                # Connect dv ivc's to solver
                for dv in [
                    FlowVars.ANGLE_OF_ATTACK,
                    FlowVars.YAW_ANGLE,
                    FlowVars.MACH_NUMBER,
                    FlowVars.REYNOLDS_NUMBER,
                    "rho",
                    "v",
                    "cg",
                ]:
                    self.connect(dv, f"aero_point_0.{dv}")

        prob = om.Problem()
        prob.model = Top()
        prob.setup()

        prob.run_model()

        return prob

    def get_oas_solution(self, surfaces):
        # Create the problem and the model group
        prob = om.Problem()

        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output("v", val=vel, units="m/s")
        indep_var_comp.add_output("alpha", val=aoa, units="deg")
        indep_var_comp.add_output("beta", val=beta, units="deg")
        indep_var_comp.add_output("Mach_number", val=mach)
        indep_var_comp.add_output("re", val=re, units="1/m")
        indep_var_comp.add_output("rho", val=rho, units="kg/m**3")
        indep_var_comp.add_output("cg", val=cg, units="m")

        prob.model.add_subsystem("prob_vars", indep_var_comp, promotes=["*"])

        # Loop over each surface in the surfaces list
        for surface in surfaces:
            geom_group = Geometry(surface=surface)

            # Add tmp_group to the problem as the name of the surface.
            # Note that is a group and performance group for each
            # individual surface.
            prob.model.add_subsystem(surface["name"], geom_group)

        # Loop through and add a certain number of aero points
        for i in range(1):
            # Create the aero point group and add it to the model
            aero_group = AeroPoint(surfaces=surfaces, compressible=True)
            point_name = "aero_point_{}".format(i)
            prob.model.add_subsystem(point_name, aero_group)

            # Connect flow properties to the analysis point
            prob.model.connect("v", point_name + ".v")
            prob.model.connect("alpha", point_name + ".alpha")
            prob.model.connect("Mach_number", point_name + ".Mach_number")
            prob.model.connect("re", point_name + ".re")
            prob.model.connect("rho", point_name + ".rho")
            prob.model.connect("cg", point_name + ".cg")

            # Connect the parameters within the model for each aero point
            for surface in surfaces:
                name = surface["name"]

                # Connect the mesh from the geometry component to the analysis point
                prob.model.connect(name + ".mesh", point_name + "." + name + ".def_mesh")

                # Perform the connections with the modified names within the
                # 'aero_states' group.
                prob.model.connect(name + ".mesh", point_name + ".aero_states." + name + "_def_mesh")

                prob.model.connect(name + ".t_over_c", point_name + "." + name + "_perf." + "t_over_c")

        # Set up the problem
        prob.setup()

        prob.run_model()

        return prob


if __name__ == "__main__":
    unittest.main()
