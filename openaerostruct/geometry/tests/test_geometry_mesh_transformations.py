""" Unit tests for each geometry mesh transformation component."""
import numpy as np

import unittest

import openmdao.api as om
from openmdao.utils.assert_utils import assert_rel_error, assert_check_partials

from openaerostruct.geometry.geometry_mesh_transformations import \
     Taper, ScaleX, Sweep, ShearX, Stretch, ShearY, Dihedral, \
     ShearZ, Rotate
from openaerostruct.geometry.utils import generate_mesh

# These have been chosen so that each dimension of the intermediate ndarrays is unique.
NY = 7
NX = 5


def get_mesh(symmetry):
    """
    Return a mesh for testing.
    """
    ny = (2*NY - 1) if symmetry else NY

    # Create a dictionary to store options about the mesh
    mesh_dict = {'num_y' : ny,
                 'num_x' : NX,
                 'wing_type' : 'CRM',
                 'symmetry' : symmetry,
                 'num_twist_cp' : NY}

    # Generate the aerodynamic mesh based on the previous dictionary
    mesh, twist_cp = generate_mesh(mesh_dict)

    surface = {}
    surface['symmetry'] = symmetry
    surface['type'] = 'aero'

    # Random perturbations to the mesh so that we don't mask errors subtractively.
    mesh[:, :, 0] += 0.05*np.random.random(mesh[:, :, 2].shape)
    mesh[:, :, 1] += 0.05*np.random.random(mesh[:, :, 2].shape)
    mesh[:, :, 2] = np.random.random(mesh[:, :, 2].shape)

    return mesh


class Test(unittest.TestCase):

    def test_taper(self):
        symmetry = False
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = np.random.random(1)

        comp = Taper(val=val, mesh=mesh, symmetry=symmetry)
        group.add_subsystem('comp', comp)

        prob.setup()
        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_taper_symmetry(self):
        symmetry = True
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = np.random.random(1)

        comp = Taper(val=val, mesh=mesh, symmetry=symmetry)
        group.add_subsystem('comp', comp)

        prob.setup()
        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_scalex(self):
        symmetry = False
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = np.random.random(NY)

        comp = ScaleX(val=val, mesh_shape=mesh.shape)
        group.add_subsystem('comp', comp)

        prob.setup()

        prob['comp.in_mesh'] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_scalex_symmetry(self):
        symmetry = True
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = np.random.random(NY)

        comp = ScaleX(val=val, mesh_shape=mesh.shape)
        group.add_subsystem('comp', comp)

        prob.setup()

        prob['comp.in_mesh'] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_sweep(self):
        symmetry = False
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = np.random.random(1)

        comp = Sweep(val=val, mesh_shape=mesh.shape, symmetry=symmetry)
        group.add_subsystem('comp', comp)

        prob.setup()

        prob['comp.in_mesh'] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_sweep_symmetry(self):
        symmetry = True
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = np.random.random(1)

        comp = Sweep(val=val, mesh_shape=mesh.shape, symmetry=symmetry)
        group.add_subsystem('comp', comp)

        prob.setup()

        prob['comp.in_mesh'] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_shearx(self):
        symmetry = False
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = np.random.random(NY)

        comp = ShearX(val=val, mesh_shape=mesh.shape)
        group.add_subsystem('comp', comp)

        prob.setup()

        prob['comp.in_mesh'] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_stretch(self):
        symmetry = False
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = np.random.random(1)

        comp = Stretch(val=val, mesh_shape=mesh.shape, symmetry=symmetry)
        group.add_subsystem('comp', comp)

        prob.setup()

        prob['comp.in_mesh'] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_stretch_symmetry(self):
        symmetry = True
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = np.random.random(1)

        comp = Stretch(val=val, mesh_shape=mesh.shape, symmetry=symmetry)
        group.add_subsystem('comp', comp)

        prob.setup()

        prob['comp.in_mesh'] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_sheary(self):
        symmetry = False
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = np.random.random(NY)

        comp = ShearY(val=val, mesh_shape=mesh.shape)
        group.add_subsystem('comp', comp)

        prob.setup()

        prob['comp.in_mesh'] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_dihedral(self):
        symmetry = False
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = 15.0*np.random.random(1)

        comp = Dihedral(val=val, mesh_shape=mesh.shape, symmetry=symmetry)
        group.add_subsystem('comp', comp)

        prob.setup()

        prob['comp.in_mesh'] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_dihedral_symmetry(self):
        symmetry = True
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = np.random.random(1)

        comp = Dihedral(val=val, mesh_shape=mesh.shape, symmetry=symmetry)
        group.add_subsystem('comp', comp)

        prob.setup()

        prob['comp.in_mesh'] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_shearz(self):
        symmetry = False
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = np.random.random(NY)

        comp = ShearZ(val=val, mesh_shape=mesh.shape)
        group.add_subsystem('comp', comp)

        prob.setup()

        prob['comp.in_mesh'] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_rotate(self):
        symmetry = False
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = np.random.random(NY)

        comp = Rotate(val=val, mesh_shape=mesh.shape, symmetry=symmetry)
        group.add_subsystem('comp', comp)

        prob.setup()

        prob['comp.in_mesh'] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

    def test_rotate_symmetry(self):
        symmetry = True
        mesh = get_mesh(symmetry)

        prob = om.Problem()
        group = prob.model

        val = np.random.random(NY)

        comp = Rotate(val=val, mesh_shape=mesh.shape, symmetry=symmetry)
        group.add_subsystem('comp', comp)

        prob.setup()

        prob['comp.in_mesh'] = mesh

        prob.run_model()

        check = prob.check_partials(compact_print=True, abs_err_tol=1e-5, rel_err_tol=1e-5)
        assert_check_partials(check, atol=1e-6, rtol=1e-6)

if __name__ == '__main__':
    unittest.main()
