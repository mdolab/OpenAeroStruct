import unittest


class Test(unittest.TestCase):
    def test_import(self):
        # Ensure meshing imports from old module works
        try:
            from openaerostruct.geometry.utils import generate_mesh  # noqa: F401
        except ImportError:
            self.fail("generate_mesh cannot be imported from utils module")

        try:
            from openaerostruct.geometry.utils import gen_rect_mesh  # noqa: F401
        except ImportError:
            self.fail("gen_rect_mesh cannot be imported from utils module")

        try:
            from openaerostruct.geometry.utils import gen_crm_mesh  # noqa: F401
        except ImportError:
            self.fail("gen_crm_mesh cannot be imported from utils module")

        try:
            from openaerostruct.geometry.utils import add_chordwise_panels  # noqa: F401
        except ImportError:
            self.fail("add_chordwise_panels cannot be imported from utils module")

        try:
            from openaerostruct.geometry.utils import get_default_geo_dict  # noqa: F401
        except ImportError:
            self.fail("get_default_geo_dict cannot be imported from utils module")

        try:
            from openaerostruct.geometry.utils import writeMesh  # noqa: F401
        except ImportError:
            self.fail("writeMesh cannot be imported from utils module")

        try:
            from openaerostruct.geometry.utils import getFullMesh  # noqa: F401
        except ImportError:
            self.fail("getFullMesh cannot be imported from utils module")


if __name__ == "__main__":
    unittest.main()
