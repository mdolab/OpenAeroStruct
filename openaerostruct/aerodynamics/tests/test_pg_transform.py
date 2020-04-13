import unittest

import numpy as np

import openmdao.api as om
from openmdao.utils.assert_utils import assert_check_partials

from openaerostruct.aerodynamics.pg_transform import PGTransform
from openaerostruct.aerodynamics.pg_scale import ScaleFromPrandtlGlauert, ScaleToPrandtlGlauert
from openaerostruct.aerodynamics.pg_wind_rotation import RotateFromWindFrame, RotateToWindFrame
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):

    def test_to_wind_derivs(self):
        surfaces = get_default_surfaces()

        comp = RotateToWindFrame(surfaces=surfaces, rotational=True)

        prob = om.Problem()
        prob.model.add_subsystem('comp', comp)
        prob.setup(force_alloc_complex=True)

        prob['comp.alpha'] = np.random.random(prob['comp.alpha'].shape)
        prob['comp.beta'] = np.random.random(prob['comp.beta'].shape)
        prob['comp.coll_pts'] = np.random.random(prob['comp.coll_pts'].shape)
        prob['comp.bound_vecs'] = np.random.random(prob['comp.bound_vecs'].shape)
        prob['comp.rotational_velocities'] = np.random.random(prob['comp.rotational_velocities'].shape)
        prob['comp.wing_def_mesh'] = np.random.random(prob['comp.wing_def_mesh'].shape)
        prob['comp.tail_def_mesh'] = np.random.random(prob['comp.tail_def_mesh'].shape)
        prob['comp.tail_normals'] = np.random.random(prob['comp.tail_normals'].shape)
        prob['comp.wing_normals'] = np.random.random(prob['comp.wing_normals'].shape)

        prob.run_model()

        check = prob.check_partials(compact_print=True, method='cs', step=1e-40)

        assert_check_partials(check)

    def test_from_wind_derivs(self):
        surfaces = get_default_surfaces()

        comp = RotateFromWindFrame(surfaces=surfaces)

        prob = om.Problem()
        prob.model.add_subsystem('comp', comp)
        prob.setup(force_alloc_complex=True)

        prob['comp.alpha'] = np.random.random(prob['comp.alpha'].shape)
        prob['comp.beta'] = np.random.random(prob['comp.beta'].shape)
        prob['comp.wing_sec_forces_w_frame'] = np.random.random(prob['comp.wing_sec_forces_w_frame'].shape)
        prob['comp.tail_sec_forces_w_frame'] = np.random.random(prob['comp.tail_sec_forces_w_frame'].shape)

        prob.run_model()

        check = prob.check_partials(compact_print=True, method='cs', step=1e-40)

        assert_check_partials(check)

    def test_scale_to_pg(self):
        surfaces = get_default_surfaces()

        comp = ScaleToPrandtlGlauert(surfaces=surfaces, rotational=True)

        prob = om.Problem()
        prob.model.add_subsystem('comp', comp)
        prob.setup(force_alloc_complex=True)

        prob['comp.Mach_number'] = np.random.random(prob['comp.Mach_number'].shape)
        prob['comp.coll_pts_w_frame'] = np.random.random(prob['comp.coll_pts_w_frame'].shape)
        prob['comp.bound_vecs_w_frame'] = np.random.random(prob['comp.bound_vecs_w_frame'].shape)
        prob['comp.rotational_velocities_w_frame'] = np.random.random(prob['comp.rotational_velocities_w_frame'].shape)
        prob['comp.wing_def_mesh_w_frame'] = np.random.random(prob['comp.wing_def_mesh_w_frame'].shape)
        prob['comp.tail_def_mesh_w_frame'] = np.random.random(prob['comp.tail_def_mesh_w_frame'].shape)
        prob['comp.tail_normals_w_frame'] = np.random.random(prob['comp.tail_normals_w_frame'].shape)
        prob['comp.wing_normals_w_frame'] = np.random.random(prob['comp.wing_normals_w_frame'].shape)

        prob.run_model()

        check = prob.check_partials(compact_print=True, method='cs', step=1e-40)

        assert_check_partials(check)

    def test_scale_from_pg(self):
        surfaces = get_default_surfaces()

        comp = ScaleFromPrandtlGlauert(surfaces=surfaces)

        prob = om.Problem()
        prob.model.add_subsystem('comp', comp)
        prob.setup(force_alloc_complex=True)

        prob['comp.Mach_number'] = np.random.random(prob['comp.Mach_number'].shape)
        prob['comp.wing_sec_forces_pg'] = np.random.random(prob['comp.wing_sec_forces_pg'].shape)
        prob['comp.tail_sec_forces_pg'] = np.random.random(prob['comp.tail_sec_forces_pg'].shape)

        prob.run_model()

        check = prob.check_partials(compact_print=True, method='cs', step=1e-40)

        assert_check_partials(check)


if __name__ == '__main__':
    unittest.main()
