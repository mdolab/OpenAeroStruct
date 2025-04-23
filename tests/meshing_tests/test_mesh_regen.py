import unittest
from openmdao.utils.assert_utils import assert_near_equal


class Test(unittest.TestCase):
    def test_span_regen_rect(self):
        # Checks if a valid mesh is returned for various spanwise spacing factors and numbers of points
        import numpy as np
        from openaerostruct.meshing.mesh_generator import gen_rect_mesh
        from openaerostruct.meshing.utils import regen_spanwise_panels

        num_x = 2
        num_y = 11
        span = 10
        span_cos_spacing = [0.0, 0.25, 0.5, 0.75, 1.0, 2.0, 2.25, 2.5, 2.75, 3.0]

        mesh = gen_rect_mesh(num_x, num_y, span, chord=2, span_cos_spacing=0.0, chord_cos_spacing=0.0)
        rc = mesh[0, -1, 0] - mesh[-1, -1, 0]
        tc = mesh[0, 0, 0] - mesh[-1, 0, 0]
        b = mesh[0, 0, 1] - mesh[0, -1, 1]

        new_num_y = [7, 21, 31]

        for nny in new_num_y:
            for csp in span_cos_spacing:

                newMesh = regen_spanwise_panels(mesh, nny, span_cos_spacing=csp)

                new_rc = newMesh[0, -1, 0] - newMesh[-1, -1, 0]
                new_tc = newMesh[0, 0, 0] - newMesh[-1, 0, 0]
                new_b = newMesh[0, 0, 1] - newMesh[0, -1, 1]

                assert_near_equal(new_rc, rc, 1e-10)
                assert_near_equal(new_tc, tc, 1e-10)
                assert_near_equal(new_b, b, 1e-10)

    def test_span_regen_crm(self):
        # Checks if a valid mesh is returned for various spanwise spacing factors and numbers of points
        import numpy as np
        from openaerostruct.meshing.mesh_generator import gen_crm_mesh
        from openaerostruct.meshing.utils import regen_spanwise_panels

        num_x = 2
        num_y = 11

        span_cos_spacing = [0.0, 0.25, 0.5, 0.75, 1.0, 2.0, 2.25, 2.5, 2.75, 3.0]

        mesh, _, _ = gen_crm_mesh(num_x, num_y, span_cos_spacing=0.0, chord_cos_spacing=0.0)
        rc = mesh[0, -1, 0] - mesh[-1, -1, 0]
        tc = mesh[0, 0, 0] - mesh[-1, 0, 0]
        b = mesh[0, 0, 1] - mesh[0, -1, 1]

        new_num_y = [7, 21, 31]

        for nny in new_num_y:
            for csp in span_cos_spacing:

                newMesh = regen_spanwise_panels(mesh, nny, span_cos_spacing=csp)

                new_rc = newMesh[0, -1, 0] - newMesh[-1, -1, 0]
                new_tc = newMesh[0, 0, 0] - newMesh[-1, 0, 0]
                new_b = newMesh[0, 0, 1] - newMesh[0, -1, 1]

                assert_near_equal(new_rc, rc, 1e-10)
                assert_near_equal(new_tc, tc, 1e-10)
                assert_near_equal(new_b, b, 1e-10)


if __name__ == "__main__":
    unittest.main()
