import openmdao.api as om
import numpy as np

from mphys.core import Multipoint, MPhysVariables
from mphys.scenarios import ScenarioAerodynamic

from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.mphys import AeroBuilder


class Top(Multipoint):
    def setup(self):
        # Create a dictionary to store options about the surface
        mesh_dict = {
            "num_y": 35,
            "num_x": 11,
            "wing_type": "rect",
            "symmetry": True,
            "span": 10.0,
            "root_chord": 1,
            "span_cos_spacing": 1.0,
            "chord_cos_spacing": 1.0,
        }

        # Generate half-wing mesh of rectangular wing
        mesh = generate_mesh(mesh_dict)

        surface = {
            # Wing definition
            "name": "wing",  # name of the surface
            "type": "aero",
            "symmetry": True,  # if true, model one half of wing
            # reflected across the plane y = 0
            "S_ref_type": "projected",  # how we compute the wing area,
            # can be 'wetted' or 'projected'
            "twist_cp": np.zeros(3),  # Define twist using 3 B-spline cp's
            # distributed along span
            "mesh": mesh,
            # Aerodynamic performance of the lifting surface at
            # an angle of attack of 0 (alpha=0).
            # These CL0 and CD0 values are added to the CL and CD
            # obtained from aerodynamic analysis of the surface to get
            # the total CL and CD.
            # These CL0 and CD0 values do not vary wrt alpha.
            "CL0": 0.0,  # CL of the surface at alpha=0
            "CD0": 0.0,  # CD of the surface at alpha=0
            # Airfoil properties for viscous drag calculation
            "k_lam": 0.05,  # percentage of chord with laminar
            # flow, used for viscous drag
            "t_over_c": 0.12,  # thickness over chord ratio (NACA0015)
            "c_max_t": 0.303,  # chordwise location of maximum (NACA0015)
            # thickness
            "with_viscous": False,  # if true, compute viscous drag,
            "with_wave": False,
        }  # end of surface dictionary

        mach = 0.0
        aoa = 5.0
        beta = 0.0
        rho = 0.38
        vel = 248.136
        re = 1e6

        dvs = self.add_subsystem("dvs", om.IndepVarComp(), promotes=["*"])
        dvs.add_output(MPhysVariables.Aerodynamics.FlowConditions.ANGLE_OF_ATTACK, val=aoa, units="deg")
        dvs.add_output(MPhysVariables.Aerodynamics.FlowConditions.YAW_ANGLE, val=beta, units="deg")
        dvs.add_output("rho", val=rho, units="kg/m**3")
        dvs.add_output(MPhysVariables.Aerodynamics.FlowConditions.MACH_NUMBER, mach)
        dvs.add_output("v", vel, units="m/s")
        dvs.add_output(MPhysVariables.Aerodynamics.FlowConditions.REYNOLDS_NUMBER, re, units="1/m")
        dvs.add_output("cg", val=np.zeros((3)), units="m")

        # Create mphys builder for aero solver
        aero_builder = AeroBuilder([surface])
        aero_builder.initialize(self.comm)

        # Create mesh component and connect with solver
        self.add_subsystem("mesh", aero_builder.get_mesh_coordinate_subsystem())
        self.mphys_add_scenario("cruise", ScenarioAerodynamic(aero_builder=aero_builder))
        self.connect(
            f"mesh.{MPhysVariables.Aerodynamics.Surface.Mesh.COORDINATES}",
            f"cruise.{MPhysVariables.Aerodynamics.Surface.COORDINATES}",
        )

        # Connect dv ivc's to solver
        for dv in [
            MPhysVariables.Aerodynamics.FlowConditions.ANGLE_OF_ATTACK,
            MPhysVariables.Aerodynamics.FlowConditions.YAW_ANGLE,
            MPhysVariables.Aerodynamics.FlowConditions.MACH_NUMBER,
            MPhysVariables.Aerodynamics.FlowConditions.REYNOLDS_NUMBER,
            "rho",
            "v",
            "cg",
        ]:
            self.connect(dv, f"cruise.{dv}")


prob = om.Problem()
prob.model = Top()
prob.setup()

# Create a n^2 diagram for user to view model connections
# om.n2(prob)

prob.run_model()

print("CL", prob["cruise.wing.CL"][0])
print("CD", prob["cruise.wing.CD"][0])
