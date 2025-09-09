import unittest
from openmdao.utils.assert_utils import assert_warning
from openmdao.utils.om_warnings import OMDeprecationWarning


class Test(unittest.TestCase):
    def test_import(self):
        # Ensure meshing imports from old module works
        try:
            from openaerostruct.geometry.utils import generate_mesh  # noqa: F401

            assert_warning(
                OMDeprecationWarning,
                msg="generate_mesh has been moved to mesh_generator.py. Importing from utils.py is deprecated and will be removed in a future release.",
            )
        except ImportError:
            self.fail("generate_mesh cannot be imported from utils module")

        try:
            from openaerostruct.geometry.utils import gen_rect_mesh  # noqa: F401

            assert_warning(
                OMDeprecationWarning,
                msg="gen_rect_mesh has been moved to mesh_generator.py. Importing from utils.py is deprecated and will be removed in a future release.",
            )
        except ImportError:
            self.fail("gen_rect_mesh cannot be imported from utils module")

        try:
            from openaerostruct.geometry.utils import gen_crm_mesh  # noqa: F401

            assert_warning(
                OMDeprecationWarning,
                msg="gen_crm_mesh has been moved to mesh_generator.py. Importing from utils.py is deprecated and will be removed in a future release.",
            )
        except ImportError:
            self.fail("gen_crm_mesh cannot be imported from utils module")

        try:
            from openaerostruct.geometry.utils import add_chordwise_panels  # noqa: F401

            assert_warning(
                OMDeprecationWarning,
                msg="add_chordwise_panels has been moved to mesh_generator.py and renamed to regen_chordwise_panels. Importing from utils.py is deprecated and will be removed in a future release.",
            )
        except ImportError:
            self.fail("add_chordwise_panels cannot be imported from utils module")

        try:
            from openaerostruct.geometry.utils import get_default_geo_dict  # noqa: F401

            assert_warning(
                OMDeprecationWarning,
                msg="get_default_geo_dict has been moved to mesh_generator.py and renamed to regen_chordwise_panels. Importing from utils.py is deprecated and will be removed in a future release.",
            )
        except ImportError:
            self.fail("get_default_geo_dict cannot be imported from utils module")

        try:
            from openaerostruct.geometry.utils import writeMesh  # noqa: F401

            assert_warning(
                OMDeprecationWarning,
                msg="writeMesh has been moved to mesh_generator.py and renamed to write_tecplot. Importing from utils.py is deprecated and will be removed in a future release.",
            )
        except ImportError:
            self.fail("writeMesh cannot be imported from utils module")

        try:
            from openaerostruct.geometry.utils import getFullMesh  # noqa: F401

            assert_warning(
                OMDeprecationWarning,
                msg="getFullMesh has been moved to mesh_generator.py. Importing from utils.py is deprecated and will be removed in a future release.",
            )
        except ImportError:
            self.fail("getFullMesh cannot be imported from utils module")


if __name__ == "__main__":
    unittest.main()
