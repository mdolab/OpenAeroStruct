""" Unit tests for each geometry mesh transformation component."""
import numpy as np

import unittest

import openmdao.api as om
from openmdao.utils.assert_utils import assert_check_partials, assert_near_equal

from openaerostruct.geometry.geometry_mesh_transformations import (
    Taper,
    ScaleX,
    Sweep,
    ShearX,
    Stretch,
    ShearY,
    Dihedral,
    ShearZ,
    Rotate,
)
from openaerostruct.geometry.utils import generate_mesh

# These have been chosen so that each dimension of the intermediate ndarrays is unique.
NY = 7
NX = 5


def get_mesh(symmetry):
    """
    Return a mesh for testing.
    """
    ny = (2 * NY - 1) if symmetry else NY

    # Create a dictionary to store options about the mesh
    mesh_dict = {
        "num_y": ny,
        "num_x": NX,
        "wing_type": "CRM",
        "symmetry": symmetry,
        "num_twist_cp": NY,
    }

    # Generate the aerodynamic mesh based on the previous dictionary
    mesh, twist_cp = generate_mesh(mesh_dict)

    surface = {}
    surface["symmetry"] = symmetry
    surface["type"] = "aero"

    # Random perturbations to the mesh so that we don't mask errors subtractively.
    rng = np.random.default_rng(0)
    mesh[:, :, 0] += 0.05 * rng.random(mesh[:, :, 2].shape)
    mesh[:, :, 1] += 0.05 * rng.random(mesh[:, :, 2].shape)
    mesh[:, :, 2] = rng.random(mesh[:, :, 2].shape)

    return mesh


class Test(unittest.TestCase):
    def setUp(self):
        # setup a random generator and fix the seed
        self.rng = np.random.default_rng(1)

    def test_taper(self):
        symmetry = False
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = self.rng.random(1)
        ref_axis_pos = self.rng.random(1)

        comp = Taper(val=val, mesh=mesh, symmetry=symmetry, ref_axis_pos=ref_axis_pos)
        group.add_subsystem("comp", comp)

        prob.setup()
        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_taper_symmetry(self):
        symmetry = True
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = self.rng.random(1)
        ref_axis_pos = self.rng.random(1)

        comp = Taper(val=val, mesh=mesh, symmetry=symmetry, ref_axis_pos=ref_axis_pos)
        group.add_subsystem("comp", comp)

        prob.setup()
        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_scalex(self):
        symmetry = False
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = self.rng.random(NY)
        ref_axis_pos = self.rng.random(1)

        comp = ScaleX(val=val, mesh_shape=mesh.shape, ref_axis_pos=ref_axis_pos)
        group.add_subsystem("comp", comp)

        prob.setup()

        prob["comp.in_mesh"] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_scalex_symmetry(self):
        symmetry = True
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = self.rng.random(NY)
        ref_axis_pos = self.rng.random(1)

        comp = ScaleX(val=val, mesh_shape=mesh.shape, ref_axis_pos=ref_axis_pos)
        group.add_subsystem("comp", comp)

        prob.setup()

        prob["comp.in_mesh"] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_scalex_ref_axis_trailing_edge(self):
        symmetry = True
        mesh = get_mesh(symmetry)

        # Test for chord_scaling_pos at trailing edge
        prob = om.Problem()
        group = prob.model

        val = self.rng.random(NY)

        comp = ScaleX(val=val, mesh_shape=mesh.shape, ref_axis_pos=1)
        group.add_subsystem("comp", comp)

        prob.setup()

        prob["comp.in_mesh"] = mesh

        prob.run_model()

        # If chord_scaling_pos = 1, TE should not move
        assert_near_equal(mesh[-1, :, :], prob["comp.mesh"][-1, :, :], tolerance=1e-10)

    def test_scalex_chord_value(self):
        # test actual chord value of a rectangular wing

        mesh_dict = {
            "num_y": 5,
            "num_x": 3,
            "wing_type": "rect",
            "symmetry": True,
            "span": 2.0,
            "root_chord": 0.2,
        }
        mesh_in = generate_mesh(mesh_dict)

        # initial chord
        chord_in = mesh_in[-1, 0, 0] - mesh_in[0, 0, 0]  # TE - LE

        prob = om.Problem()
        group = prob.model
        comp = ScaleX(val=np.ones(3), mesh_shape=mesh_in.shape)
        group.add_subsystem("comp", comp)

        prob.setup()
        chord_scaling = 1.3
        prob.set_val("comp.chord", val=chord_scaling, units=None)  # apply chord scaling factor
        prob.set_val("comp.in_mesh", val=mesh_in, units="m")

        prob.run_model()

        # chord after manipulation
        mesh_out = prob.get_val("comp.mesh", units="m")
        chord_out = mesh_out[-1, 0, 0] - mesh_out[0, 0, 0]  # TE - LE

        assert_near_equal(chord_in * chord_scaling, chord_out, tolerance=1e-10)

    def test_sweep(self):
        symmetry = False
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = self.rng.random(1)

        comp = Sweep(val=val, mesh_shape=mesh.shape, symmetry=symmetry)
        group.add_subsystem("comp", comp)

        prob.setup()

        prob["comp.in_mesh"] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_sweep_symmetry(self):
        symmetry = True
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = self.rng.random(1)

        comp = Sweep(val=val, mesh_shape=mesh.shape, symmetry=symmetry)
        group.add_subsystem("comp", comp)

        prob.setup()

        prob["comp.in_mesh"] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_shearx(self):
        symmetry = False
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = self.rng.random(NY)

        comp = ShearX(val=val, mesh_shape=mesh.shape)
        group.add_subsystem("comp", comp)

        prob.setup()

        prob["comp.in_mesh"] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_stretch(self):
        symmetry = False
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = self.rng.random(1)
        ref_axis_pos = self.rng.random(1)

        comp = Stretch(val=val, mesh_shape=mesh.shape, symmetry=symmetry, ref_axis_pos=ref_axis_pos)
        group.add_subsystem("comp", comp)

        prob.setup()

        prob["comp.in_mesh"] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_stretch_symmetry(self):
        symmetry = True
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = self.rng.random(1)
        ref_axis_pos = self.rng.random(1)

        comp = Stretch(val=val, mesh_shape=mesh.shape, symmetry=symmetry, ref_axis_pos=ref_axis_pos)
        group.add_subsystem("comp", comp)

        prob.setup()

        prob["comp.in_mesh"] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_sheary(self):
        symmetry = False
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = self.rng.random(NY)

        comp = ShearY(val=val, mesh_shape=mesh.shape)
        group.add_subsystem("comp", comp)

        prob.setup()

        prob["comp.in_mesh"] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_dihedral(self):
        symmetry = False
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = 15.0 * self.rng.random(1)

        comp = Dihedral(val=val, mesh_shape=mesh.shape, symmetry=symmetry)
        group.add_subsystem("comp", comp)

        prob.setup()

        prob["comp.in_mesh"] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_dihedral_symmetry(self):
        symmetry = True
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = self.rng.random(1)

        comp = Dihedral(val=val, mesh_shape=mesh.shape, symmetry=symmetry)
        group.add_subsystem("comp", comp)

        prob.setup()

        prob["comp.in_mesh"] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_shearz(self):
        symmetry = False
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = self.rng.random(NY)

        comp = ShearZ(val=val, mesh_shape=mesh.shape)
        group.add_subsystem("comp", comp)

        prob.setup()

        prob["comp.in_mesh"] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_rotate(self):
        symmetry = False
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = self.rng.random(NY)
        ref_axis_pos = self.rng.random(1)

        comp = Rotate(val=val, mesh_shape=mesh.shape, symmetry=symmetry, ref_axis_pos=ref_axis_pos)
        group.add_subsystem("comp", comp)

        prob.setup()

        prob["comp.in_mesh"] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_rotate_symmetry(self):
        symmetry = True
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = self.rng.random(NY)
        ref_axis_pos = self.rng.random(1)

        comp = Rotate(val=val, mesh_shape=mesh.shape, symmetry=symmetry, ref_axis_pos=ref_axis_pos)
        group.add_subsystem("comp", comp)

        prob.setup()

        prob["comp.in_mesh"] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_rotate_trailing_edge(self):
        symmetry = False
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = self.rng.random(NY)

        comp = Rotate(val=val, mesh_shape=mesh.shape, symmetry=symmetry, ref_axis_pos=1)
        group.add_subsystem("comp", comp)

        prob.setup()

        prob["comp.in_mesh"] = mesh

        prob.run_model()

        # If chord_scaling_pos = 1, TE should not move
        assert_near_equal(mesh[-1, :, :], prob["comp.mesh"][-1, :, :], tolerance=1e-10)


if __name__ == "__main__":
    unittest.main()
