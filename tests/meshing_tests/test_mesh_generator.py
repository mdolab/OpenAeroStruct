import unittest


class Test(unittest.TestCase):
    def test_span_cos_spacing_rect(self):
        # Checks if a valid mesh is returned for various spanwise spacing factors
        import numpy as np
        from openaerostruct.meshing.mesh_generator import gen_rect_mesh

        num_x = 2
        num_y = 11
        span = 10
        span_cos_spacing = [0.0, 0.25, 0.5, 0.75, 1.0, 2.0, 2.25, 2.5, 2.75, 3.0]

        for csp in span_cos_spacing:
            mesh = gen_rect_mesh(num_x, num_y, span, chord=2, span_cos_spacing=csp, chord_cos_spacing=0.0)

            self.assertTrue(np.all(np.diff(mesh[0, :, 1]) > 0.0))
            self.assertTrue(np.all(np.diff(mesh[:, 0, 0]) > 0.0))

    def test_chord_cos_spacing_rect(self):
        # Checks if a valid mesh is returned for various spanwise spacing factors
        import numpy as np
        from openaerostruct.meshing.mesh_generator import gen_rect_mesh

        num_x = 11
        num_y = 11
        span = 10
        chord_cos_spacing = [0.0, 0.25, 0.5, 0.75, 1.0]

        for csp in chord_cos_spacing:
            mesh = gen_rect_mesh(2, num_y, span, chord=2, span_cos_spacing=0.0, chord_cos_spacing=csp)

            self.assertTrue(np.all(np.diff(mesh[0, :, 1]) > 0.0))
            self.assertTrue(np.all(np.diff(mesh[:, 0, 0]) > 0.0))

    def test_span_cos_spacing_crm(self):
        # Checks if a valid mesh is returned for various spanwise spacing factors
        import numpy as np
        from openaerostruct.meshing.mesh_generator import gen_crm_mesh

        num_x = 2
        num_y = 6
        span_cos_spacing = [0.0, 0.25, 0.5, 0.75, 1.0, 2.0, 2.25, 2.5, 2.75, 3.0]

        for csp in span_cos_spacing:
            mesh, _, _ = gen_crm_mesh(num_x, num_y, span_cos_spacing=csp, chord_cos_spacing=0.0)

            self.assertTrue(np.all(np.diff(mesh[0, :, 1]) > 0.0))
            self.assertTrue(np.all(np.diff(mesh[:, 0, 0]) > 0.0))

    def test_chord_cos_spacing_crm(self):
        # Checks if a valid mesh is returned for various spanwise spacing factors
        import numpy as np
        from openaerostruct.meshing.mesh_generator import gen_crm_mesh

        num_x = 6
        num_y = 6
        chord_cos_spacing = [0.0, 0.25, 0.5, 0.75, 1.0]

        for csp in chord_cos_spacing:
            mesh, _, _ = gen_crm_mesh(num_x, num_y, span_cos_spacing=0.0, chord_cos_spacing=csp)

            self.assertTrue(np.all(np.diff(mesh[0, :, 1]) > 0.0))
            self.assertTrue(np.all(np.diff(mesh[:, 0, 0]) > 0.0))


if __name__ == "__main__":
    unittest.main()
