import unittest

from openaerostruct.geometry.geometry_unification import GeomMultiUnification
from openaerostruct.geometry.geometry_group import build_sections
from openaerostruct.utils.testing import run_test, get_three_section_surface


class Test(unittest.TestCase):
    def test(self):
        (surface, chord_bspline) = get_three_section_surface()
        sec_dicts = build_sections(surface)

        comp = GeomMultiUnification(sections=sec_dicts, surface_name=surface["name"])

        run_test(self, comp, complex_flag=True, method="cs")


if __name__ == "__main__":
    unittest.main()
