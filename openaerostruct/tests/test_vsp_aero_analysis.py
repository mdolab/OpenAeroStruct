from openmdao.utils.assert_utils import assert_near_equal
import unittest
import numpy as np
import os

from openaerostruct.geometry.utils import generate_vsp_surfaces
from openaerostruct.geometry.geometry_group import Geometry
from openaerostruct.aerodynamics.aero_groups import AeroPoint

import openmdao.api as om

vsp_file = os.path.join(os.path.dirname(__file__), "vsp_model.vsp3")

# check if openvsp is available
try:
    generate_vsp_surfaces(vsp_file)

    openvsp_flag = True
except ImportError:
    openvsp_flag = False


@unittest.skipUnless(openvsp_flag, "OpenVSP is required.")
class Test(unittest.TestCase):
    def test(self):
        # Run a symmetric and full model aero case
        symm_prob = self.setup_prob(True)
        full_prob = self.setup_prob(False)

        # Run analysis
        symm_prob.run_model()
        full_prob.run_model()

        assert_near_equal(symm_prob["flight_condition_0.Wing_perf.CD"][0], 0.010722632534543076, 1e-6)
        assert_near_equal(symm_prob["flight_condition_0.Wing_perf.CL"][0], 0.5246182619241868, 1e-6)
        assert_near_equal(symm_prob["flight_condition_0.CM"][1], -0.591626120207946, 1e-6)

        assert_near_equal(
            symm_prob["flight_condition_0.Wing_perf.CD"][0], full_prob["flight_condition_0.Wing_perf.CD"][0], 1e-6
        )
        assert_near_equal(
            symm_prob["flight_condition_0.Wing_perf.CL"][0], full_prob["flight_condition_0.Wing_perf.CL"][0], 1e-6
        )
        assert_near_equal(symm_prob["flight_condition_0.CM"][1], full_prob["flight_condition_0.CM"][1], 1e-6)

    def setup_prob(self, symmetry):
        """
        Setup openMDAO problem for symmetric or full aerodynamic analysis.
        """
        # Generate half-body mesh, include only wing and tail surfaces
        surfaces = generate_vsp_surfaces(
            vsp_file, symmetry=symmetry, include=["Wing", "Horizontal_Tail", "Vertical_Tail"]
        )

        # Define input surface dictionary for our wing
        surf_options = {
            "type": "aero",
            "S_ref_type": "wetted",  # how we compute the wing area,
            # can be 'wetted' or 'projected'
            "twist_cp": np.zeros(3),  # Define twist using 3 B-spline cp's
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

        # Update each surface with default options
        for surface in surfaces:
            surface.update(surf_options)

        # -----------------------------------------------------------------------------
        # END SURFACES
        # -----------------------------------------------------------------------------
        # docs checkpoint 1

        # Instantiate the problem and the model group
        prob = om.Problem()

        # Define flight variables as independent variables of the model
        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output("v", val=248.136, units="m/s")  # Freestream Velocity
        indep_var_comp.add_output("alpha", val=5.0, units="deg")  # Angle of Attack
        indep_var_comp.add_output("beta", val=0.0, units="deg")  # Sideslip angle
        indep_var_comp.add_output("omega", val=np.zeros(3), units="deg/s")  # Rotation rate
        indep_var_comp.add_output("Mach_number", val=0.0)  # Freestream Mach number
        indep_var_comp.add_output("re", val=1.0e6, units="1/m")  # Freestream Reynolds number
        indep_var_comp.add_output("rho", val=0.38, units="kg/m**3")  # Freestream air density
        indep_var_comp.add_output("cg", val=np.zeros((3)), units="m")  # Aircraft center of gravity
        # Add vars to model, promoting is a quick way of automatically connecting inputs
        # and outputs of different OpenMDAO components
        prob.model.add_subsystem("flight_vars", indep_var_comp, promotes=["*"])

        # Add geometry group to the problem and add wing surface as a sub group.
        # These groups are responsible for manipulating the geometry of the mesh,
        # in this case spanwise twist.
        geom_group = om.Group()
        for surface in surfaces:
            geom_group.add_subsystem(surface["name"], Geometry(surface=surface))
        prob.model.add_subsystem("geom", geom_group, promotes=["*"])

        # Create the aero point group for this flight condition and add it to the model
        aero_group = AeroPoint(surfaces=surfaces, rotational=True, compressible=True)
        point_name = "flight_condition_0"
        prob.model.add_subsystem(
            point_name, aero_group, promotes_inputs=["v", "alpha", "beta", "omega", "Mach_number", "re", "rho", "cg"]
        )

        # Connect the mesh from the geometry component to the analysis point
        for surface in surfaces:
            name = surface["name"]
            prob.model.connect(name + ".mesh", point_name + "." + name + ".def_mesh")

            # Perform the connections with the modified names within the
            # 'aero_states' group.
            prob.model.connect(name + ".mesh", point_name + ".aero_states." + name + "_def_mesh")

        # Set up the problem
        prob.setup()

        return prob


if __name__ == "__main__":
    unittest.main()
