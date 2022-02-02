import unittest
import numpy as np

import openmdao.api as om
from openmdao.utils.assert_utils import assert_check_partials, assert_near_equal

from openaerostruct.aerodynamics.geometry import VLMGeometry
from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def setUp(self):
        # setup a random generator and fix the seed
        self.rng = np.random.default_rng(1)

    def test(self):
        surfaces = get_default_surfaces()

        group = om.Group()

        comp = VLMGeometry(surface=surfaces[0])

        indep_var_comp = om.IndepVarComp()

        indep_var_comp.add_output("def_mesh", val=surfaces[0]["mesh"], units="m")

        group.add_subsystem("geom", comp)
        group.add_subsystem("indep_var_comp", indep_var_comp)

        group.connect("indep_var_comp.def_mesh", "geom.def_mesh")

        run_test(self, group)

    def test_derivs_wetted(self):
        # This is a much richer test with the following attributes:
        # - 5x7 mesh so that each dimension of all variables is unique.
        # - Random values given as inputs.
        # - The mesh has been given a random z height at all points.
        # - S_ref_type option is "wetted"
        # Create a dictionary to store options about the mesh
        mesh_dict = {"num_y": 7, "num_x": 5, "wing_type": "CRM", "symmetry": True, "num_twist_cp": 5}

        # Generate the aerodynamic mesh based on the previous dictionary
        mesh, twist_cp = generate_mesh(mesh_dict)

        mesh[:, :, 2] = self.rng.random(mesh[:, :, 2].shape)

        # Create a dictionary with info and options about the aerodynamic
        # lifting surface
        surface = {
            # Wing definition
            "name": "wing",  # name of the surface
            "symmetry": True,  # if true, model one half of wing
            # reflected across the plane y = 0
            "S_ref_type": "wetted",  # how we compute the wing area,
            # can be 'wetted' or 'projected'
            "fem_model_type": "tube",
            "twist_cp": twist_cp,
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
            "with_viscous": True,  # if true, compute viscous drag
            "with_wave": False,  # if true, compute wave drag
        }

        # surfaces = get_default_surfaces()
        surfaces = [surface]

        prob = om.Problem()
        group = prob.model

        comp = VLMGeometry(surface=surfaces[0])

        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output("def_mesh", val=surfaces[0]["mesh"], units="m")

        group.add_subsystem("geom", comp)
        group.add_subsystem("indep_var_comp", indep_var_comp)

        group.connect("indep_var_comp.def_mesh", "geom.def_mesh")

        prob.setup()

        prob["geom.def_mesh"] = mesh + self.rng.random(prob["geom.def_mesh"].shape) * 0.1

        prob.run_model()

        check = prob.check_partials(compact_print=True, method="cs")
        assert_check_partials(check, atol=5e-5, rtol=1e-5)

    def test_derivs_projected(self):
        # This is a much richer test with the following attributes:
        # - 5x7 mesh so that each dimension of all variables is unique.
        # - Random values given as inputs.
        # - The mesh has been given a random z height at all points.
        # - S_ref_type option is "projected"

        # Create a dictionary to store options about the mesh
        mesh_dict = {"num_y": 7, "num_x": 5, "wing_type": "CRM", "symmetry": True, "num_twist_cp": 5}

        # Generate the aerodynamic mesh based on the previous dictionary
        mesh, twist_cp = generate_mesh(mesh_dict)

        mesh[:, :, 2] = self.rng.random(mesh[:, :, 2].shape)

        # Create a dictionary with info and options about the aerodynamic
        # lifting surface
        surface = {
            # Wing definition
            "name": "wing",  # name of the surface
            "symmetry": True,  # if true, model one half of wing
            # reflected across the plane y = 0
            "S_ref_type": "projected",  # how we compute the wing area,
            # can be 'wetted' or 'projected'
            "fem_model_type": "tube",
            "twist_cp": twist_cp,
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
            "with_viscous": True,  # if true, compute viscous drag
            "with_wave": False,  # if true, compute wave drag
        }

        # surfaces = get_default_surfaces()
        surfaces = [surface]

        prob = om.Problem()
        group = prob.model

        comp = VLMGeometry(surface=surfaces[0])

        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output("def_mesh", val=surfaces[0]["mesh"], units="m")

        group.add_subsystem("geom", comp)
        group.add_subsystem("indep_var_comp", indep_var_comp)

        group.connect("indep_var_comp.def_mesh", "geom.def_mesh")

        prob.setup()

        prob["geom.def_mesh"] = mesh + self.rng.random(prob["geom.def_mesh"].shape) * 0.1

        prob.run_model()

        check = prob.check_partials(compact_print=True, method="cs")

        assert_check_partials(check, atol=5e-5, rtol=1e-5)

    def test_outputs(self):
        surfaces = get_default_surfaces()

        group = om.Group()

        comp = VLMGeometry(surface=surfaces[0])

        indep_var_comp = om.IndepVarComp()

        indep_var_comp.add_output("def_mesh", val=surfaces[0]["mesh"], units="m")

        group.add_subsystem("geom", comp, promotes=["*"])
        group.add_subsystem("indep_var_comp", indep_var_comp, promotes=["*"])

        prob = om.Problem()
        prob.model.add_subsystem("group", group, promotes=["*"])
        prob.setup()
        prob.run_model()

        assert_near_equal(prob["widths"], np.array([11.95624787, 11.90425878, 11.44086572]), 1e-6)
        assert_near_equal(prob["cos_sweep"], np.array([9.7938336, 9.79384207, 9.79385053]), 1e-6)
        assert_near_equal(prob["S_ref"], np.array([415.02211208]), 1e-6)
        assert_near_equal(prob["chords"], np.array([2.72796, 5.1252628, 7.8891638, 13.6189974]), 1e-6)
        assert_near_equal(prob["lengths"], np.array([2.72796, 5.1252628, 7.8891638, 13.6189974]), 1e-6)


if __name__ == "__main__":
    unittest.main()
