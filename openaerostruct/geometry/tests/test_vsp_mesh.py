import os
import unittest
from openmdao.utils.assert_utils import assert_near_equal

from openaerostruct.geometry.utils import generate_mesh, generate_vsp_surfaces

vsp_file = os.path.join(os.path.dirname(__file__), "rect_wing.vsp3")

# check if openvsp is available
try:
    generate_vsp_surfaces(vsp_file)

    openvsp_flag = True
except ImportError:
    openvsp_flag = False


@unittest.skipUnless(openvsp_flag, "OpenVSP is required.")
class Test(unittest.TestCase):
    def test_full(self):
        # Create a rectangular wing with uniform spacing
        mesh_dict = {
            "num_y": 9,
            "num_x": 3,
            "wing_type": "rect",
            "symmetry": False,
            "span": 10.0,
            "root_chord": 1,
            "span_cos_spacing": 0.0,
        }

        oas_mesh = generate_mesh(mesh_dict)

        # Read in equivilent wing from vsp
        vsp_surf_list = generate_vsp_surfaces(vsp_file, symmetry=False)

        assert_near_equal(vsp_surf_list[0]["mesh"], oas_mesh)

    def test_symm(self):
        # Create a rectangular wing with uniform spacing
        mesh_dict = {
            "num_y": 9,
            "num_x": 3,
            "wing_type": "rect",
            "symmetry": True,
            "span": 10.0,
            "root_chord": 1,
            "span_cos_spacing": 0.0,
        }

        oas_mesh = generate_mesh(mesh_dict)
        # Flip mesh to be right-handed
        oas_mesh = oas_mesh[:, ::-1]
        oas_mesh[:, :, 1] *= -1.0

        # Read in equivilent wing from vsp
        vsp_surf_list = generate_vsp_surfaces(vsp_file, symmetry=True)

        assert_near_equal(vsp_surf_list[0]["mesh"], oas_mesh)


if __name__ == "__main__":
    unittest.main()
