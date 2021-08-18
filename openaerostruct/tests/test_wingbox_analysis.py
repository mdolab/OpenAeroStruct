import numpy as np
import unittest
from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.structures.struct_groups import SpatialBeamAlone
from openmdao.utils.assert_utils import assert_check_partials, assert_near_equal
import openmdao.api as om

upper_x = np.array(
    [
        0.1,
        0.11,
        0.12,
        0.13,
        0.14,
        0.15,
        0.16,
        0.17,
        0.18,
        0.19,
        0.2,
        0.21,
        0.22,
        0.23,
        0.24,
        0.25,
        0.26,
        0.27,
        0.28,
        0.29,
        0.3,
        0.31,
        0.32,
        0.33,
        0.34,
        0.35,
        0.36,
        0.37,
        0.38,
        0.39,
        0.4,
        0.41,
        0.42,
        0.43,
        0.44,
        0.45,
        0.46,
        0.47,
        0.48,
        0.49,
        0.5,
        0.51,
        0.52,
        0.53,
        0.54,
        0.55,
        0.56,
        0.57,
        0.58,
        0.59,
        0.6,
    ],
    dtype="complex128",
)
lower_x = np.array(
    [
        0.1,
        0.11,
        0.12,
        0.13,
        0.14,
        0.15,
        0.16,
        0.17,
        0.18,
        0.19,
        0.2,
        0.21,
        0.22,
        0.23,
        0.24,
        0.25,
        0.26,
        0.27,
        0.28,
        0.29,
        0.3,
        0.31,
        0.32,
        0.33,
        0.34,
        0.35,
        0.36,
        0.37,
        0.38,
        0.39,
        0.4,
        0.41,
        0.42,
        0.43,
        0.44,
        0.45,
        0.46,
        0.47,
        0.48,
        0.49,
        0.5,
        0.51,
        0.52,
        0.53,
        0.54,
        0.55,
        0.56,
        0.57,
        0.58,
        0.59,
        0.6,
    ],
    dtype="complex128",
)
upper_y = np.array(
    [
        0.0447,
        0.046,
        0.0472,
        0.0484,
        0.0495,
        0.0505,
        0.0514,
        0.0523,
        0.0531,
        0.0538,
        0.0545,
        0.0551,
        0.0557,
        0.0563,
        0.0568,
        0.0573,
        0.0577,
        0.0581,
        0.0585,
        0.0588,
        0.0591,
        0.0593,
        0.0595,
        0.0597,
        0.0599,
        0.06,
        0.0601,
        0.0602,
        0.0602,
        0.0602,
        0.0602,
        0.0602,
        0.0601,
        0.06,
        0.0599,
        0.0598,
        0.0596,
        0.0594,
        0.0592,
        0.0589,
        0.0586,
        0.0583,
        0.058,
        0.0576,
        0.0572,
        0.0568,
        0.0563,
        0.0558,
        0.0553,
        0.0547,
        0.0541,
    ],
    dtype="complex128",
)
lower_y = np.array(
    [
        -0.0447,
        -0.046,
        -0.0473,
        -0.0485,
        -0.0496,
        -0.0506,
        -0.0515,
        -0.0524,
        -0.0532,
        -0.054,
        -0.0547,
        -0.0554,
        -0.056,
        -0.0565,
        -0.057,
        -0.0575,
        -0.0579,
        -0.0583,
        -0.0586,
        -0.0589,
        -0.0592,
        -0.0594,
        -0.0595,
        -0.0596,
        -0.0597,
        -0.0598,
        -0.0598,
        -0.0598,
        -0.0598,
        -0.0597,
        -0.0596,
        -0.0594,
        -0.0592,
        -0.0589,
        -0.0586,
        -0.0582,
        -0.0578,
        -0.0573,
        -0.0567,
        -0.0561,
        -0.0554,
        -0.0546,
        -0.0538,
        -0.0529,
        -0.0519,
        -0.0509,
        -0.0497,
        -0.0485,
        -0.0472,
        -0.0458,
        -0.0444,
    ],
    dtype="complex128",
)


class Test(unittest.TestCase):
    def test(self):

        # Create a dictionary to store options about the surface
        mesh_dict = {"num_y": 7, "wing_type": "uCRM_based", "symmetry": True, "num_twist_cp": 5}

        mesh, twist_cp = generate_mesh(mesh_dict)

        surf_dict = {
            # Wing definition
            "name": "wing",  # name of the surface
            # reflected across the plane y = 0
            "S_ref_type": "wetted",  # how we compute the wing area,
            # can be 'wetted' or 'projected'
            "fem_model_type": "wingbox",
            "symmetry": True,
            "spar_thickness_cp": np.array([0.004, 0.005, 0.005, 0.008, 0.008, 0.01]),  # [m]
            "skin_thickness_cp": np.array([0.005, 0.01, 0.015, 0.020, 0.025, 0.026]),
            "twist_cp": np.array([4.0, 5.0, 8.0, 8.0, 8.0, 9.0]),
            "mesh": mesh,
            "data_x_upper": upper_x,
            "data_x_lower": lower_x,
            "data_y_upper": upper_y,
            "data_y_lower": lower_y,
            "strength_factor_for_upper_skin": 1.0,
            # Aerodynamic performance of the lifting surface at
            # an angle of attack of 0 (alpha=0).
            # These CL0 and CD0 values are added to the CL and CD
            # obtained from aerodynamic analysis of the surface to get
            # the total CL and CD.
            # These CL0 and CD0 values do not vary wrt alpha.
            "CL0": 0.0,  # CL of the surface at alpha=0
            "CD0": 0.0078,  # CD of the surface at alpha=0
            # Airfoil properties for viscous drag calculation
            "k_lam": 0.05,  # percentage of chord with laminar
            # flow, used for viscous drag
            "t_over_c_cp": np.array([0.08, 0.08, 0.08, 0.10, 0.10, 0.08]),
            "original_wingbox_airfoil_t_over_c": 0.12,
            "c_max_t": 0.38,  # chordwise location of maximum thickness
            "with_viscous": True,
            "with_wave": True,  # if true, compute wave drag
            # Structural values are based on aluminum 7075
            "E": 73.1e9,  # [Pa] Young's modulus
            "G": (73.1e9 / 2 / 1.33),  # [Pa] shear modulus (calculated using E and the Poisson's ratio here)
            "yield": (420.0e6 / 1.5),  # [Pa] allowable yield stress
            "mrho": 2.78e3,  # [kg/m^3] material density
            "strength_factor_for_upper_skin": 1.0,  # the yield stress is multiplied by this factor for the upper skin
            # 'fem_origin' : 0.35,    # normalized chordwise location of the spar
            "wing_weight_ratio": 1.25,
            "struct_weight_relief": False,
            "distributed_fuel_weight": True,
            # Constraints
            "exact_failure_constraint": True,  # if false, use KS function
            "fuel_density": 803.0,  # [kg/m^3] fuel density (only needed if the fuel-in-wing volume constraint is used)
            "Wf_reserve": 15000.0,  # [kg] reserve fuel mass
        }

        # Create the problem and assign the model group
        prob = om.Problem()

        ny = surf_dict["mesh"].shape[1]

        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output("loads", val=np.ones((ny, 6)) * 2e5, units="N")
        indep_var_comp.add_output("load_factor", val=1.0)
        indep_var_comp.add_output("fuel_mass", val=10000.0, units="kg")
        struct_group = SpatialBeamAlone(surface=surf_dict)
        # Add indep_vars to the structural group
        struct_group.add_subsystem("indep_vars", indep_var_comp, promotes=["*"])
        prob.model.add_subsystem(surf_dict["name"], struct_group)
        if surf_dict["distributed_fuel_weight"]:
            prob.model.connect("wing.fuel_mass", "wing.struct_states.fuel_mass")
            prob.model.connect("wing.struct_setup.fuel_vols", "wing.struct_states.fuel_vols")

        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options["tol"] = 1e-9

        # Setup problem and add design variables, constraint, and objective
        prob.model.add_design_var("wing.spar_thickness_cp", lower=0.01, upper=0.5, ref=1e-1)
        prob.model.add_design_var("wing.skin_thickness_cp", lower=0.01, upper=0.5, ref=1e-1)
        prob.model.add_constraint("wing.failure", upper=0.0)
        # prob.model.add_constraint('wing.thickness_intersects', upper=0.)

        # Add design variables, constraisnt, and objective on the problem
        prob.model.add_objective("wing.structural_mass", scaler=1e-5)

        # Set up the problem
        prob.setup()

        prob.run_model()

        # the list of components that uses CS to compute derivatives. We exclude them from check_partials
        excludes_list = [
            "wing.wingbox_group.wingbox",
            "wing.struct_setup.fuel_vol",
            "wing.struct_states.fuel_loads",
            "wing.struct_funcs.vonmises",
        ]
        data = prob.check_partials(compact_print=True, out_stream=None, excludes=excludes_list, method="cs")
        assert_check_partials(data, atol=1e20, rtol=1e-6)

        prob.run_driver()
        assert_near_equal(prob["wing.structural_mass"], 16704.10113356, 1e-6)


if __name__ == "__main__":
    unittest.main()
