import openmdao.api as om
from openaerostruct.aerodynamics.lift_drag import LiftDrag
from openaerostruct.aerodynamics.coeffs import Coeffs
from openaerostruct.aerodynamics.total_lift import TotalLift
from openaerostruct.aerodynamics.total_drag import TotalDrag
from openaerostruct.aerodynamics.viscous_drag import ViscousDrag
from openaerostruct.aerodynamics.wave_drag import WaveDrag
from openaerostruct.aerodynamics.lift_coeff_2D import LiftCoeff2D


class VLMFunctionals(om.Group):
    """
    Group that contains the aerodynamic functionals used to evaluate
    performance. These are not included in the coupled aerostructural group,
    but are only used to compute aerodynamic performance. This includes
    computing lift, drag, CL, CD, viscous CD, and wave drag CD.
    """

    def initialize(self):
        self.options.declare("surface", types=dict)

    def setup(self):
        surface = self.options["surface"]

        self.add_subsystem(
            "liftcoeff",
            LiftCoeff2D(surface=surface),
            promotes_inputs=["v", "alpha", "rho", "widths", "chords", "sec_forces"],
            promotes_outputs=["Cl"],
        )

        self.add_subsystem(
            "liftdrag",
            LiftDrag(surface=surface),
            promotes_inputs=["alpha", "beta", "sec_forces"],
            promotes_outputs=["L", "D"],
        )

        self.add_subsystem(
            "coeffs", Coeffs(), promotes_inputs=["v", "rho", "S_ref", "L", "D"], promotes_outputs=["CL1", "CDi"]
        )

        self.add_subsystem("CL", TotalLift(surface=surface), promotes_inputs=["CL1"], promotes_outputs=["CL"])

        self.add_subsystem(
            "viscousdrag",
            ViscousDrag(surface=surface),
            promotes_inputs=["Mach_number", "re", "widths", "cos_sweep", "lengths", "S_ref", "t_over_c"],
            promotes_outputs=["CDv"],
        )

        self.add_subsystem(
            "wavedrag",
            WaveDrag(surface=surface),
            promotes_inputs=["Mach_number", "cos_sweep", "widths", "CL", "chords", "t_over_c"],
            promotes_outputs=["CDw"],
        )

        self.add_subsystem(
            "CD", TotalDrag(surface=surface), promotes_inputs=["CDv", "CDi", "CDw"], promotes_outputs=["CD"]
        )
