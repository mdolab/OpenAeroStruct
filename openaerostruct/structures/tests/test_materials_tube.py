import unittest

from openaerostruct.structures.section_properties_tube import SectionPropertiesTube
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def test(self):
        surface = get_default_surfaces()[0]

        comp = SectionPropertiesTube(surface=surface)

        run_test(self, comp, complex_flag=True)


if __name__ == "__main__":
    unittest.main()
