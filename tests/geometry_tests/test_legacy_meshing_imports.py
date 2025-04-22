import unittest


class Test(unittest.TestCase):
    def test_import(self):
        # Ensure meshing imports from old module works
        try:
            from openaerostruct.geometry.utils import generate_mesh
        except ImportError:
            self.fail("generate_mesh cannot be imported from utils module")

        try:
            from openaerostruct.geometry.utils import gen_rect_mesh
        except ImportError:
            self.fail("gen_rect_mesh cannot be imported from utils module")

        try:
            from openaerostruct.geometry.utils import gen_crm_mesh
        except ImportError:
            self.fail("gen_crm_mesh cannot be imported from utils module")

        try:
            from openaerostruct.geometry.utils import add_chordwise_panels
        except ImportError:
            self.fail("add_chordwise_panels cannot be imported from utils module")

        try:
            from openaerostruct.geometry.utils import get_default_geo_dict
        except ImportError:
            self.fail("get_default_geo_dict cannot be imported from utils module")

        try:
            from openaerostruct.geometry.utils import writeMesh
        except ImportError:
            self.fail("writeMesh cannot be imported from utils module")

        try:
            from openaerostruct.geometry.utils import getFullMesh
        except ImportError:
            self.fail("getFullMesh cannot be imported from utils module")


if __name__ == "__main__":
    unittest.main()
