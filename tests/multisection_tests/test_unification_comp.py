import unittest

import openmdao.api as om
from openaerostruct.geometry.geometry_unification import GeomMultiUnification
from openaerostruct.geometry.geometry_group import build_sections
from openaerostruct.geometry.utils import stretch, sweep, dihedral
from openaerostruct.utils.testing import run_test, get_three_section_surface


class Test(unittest.TestCase):
    def test_no_shift(self):
        (surface, chord_bspline) = get_three_section_surface()
        sec_dicts = build_sections(surface)

        comp = GeomMultiUnification(sections=sec_dicts, surface_name=surface["name"], shift_uni_mesh=False)

        group = om.Group()

        group.add_subsystem("comp", comp, promotes=["*"])

        group.set_input_defaults("sec0_def_mesh", sec_dicts[0]["mesh"])
        group.set_input_defaults("sec1_def_mesh", sec_dicts[1]["mesh"])
        group.set_input_defaults("sec2_def_mesh", sec_dicts[2]["mesh"])

        run_test(self, group, complex_flag=True, method="cs")

    def test_shift(self):
        (surface, chord_bspline) = get_three_section_surface()

        # Apply some scalar mesh transformations
        surface["span"] = [5.0, 5.0, 3.0]
        surface["sweep"] = [-10.0, 10.0, -20.0]
        surface["dihedral"] = [-10.0, 10.0, -20.0]

        sec_dicts = build_sections(surface)

        for i in range(surface["num_sections"]):
            sweep(sec_dicts[i]["mesh"], surface["sweep"][i], True)
            stretch(sec_dicts[i]["mesh"], surface["span"][i], True)
            dihedral(sec_dicts[i]["mesh"], surface["dihedral"][i], True)

        comp = GeomMultiUnification(sections=sec_dicts, surface_name=surface["name"], shift_uni_mesh=True)

        group = om.Group()

        group.add_subsystem("comp", comp, promotes=["*"])

        group.set_input_defaults("sec0_def_mesh", sec_dicts[0]["mesh"])
        group.set_input_defaults("sec1_def_mesh", sec_dicts[1]["mesh"])
        group.set_input_defaults("sec2_def_mesh", sec_dicts[2]["mesh"])

        run_test(self, group, complex_flag=True, method="cs")


if __name__ == "__main__":
    unittest.main()
