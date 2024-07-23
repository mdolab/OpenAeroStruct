import unittest

import openmdao.api as om

from openaerostruct.mphys.aero_funcs_group import AeroFuncsGroup
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def test(self):
        surfaces = get_default_surfaces()

        group = om.Group()

        ivc = group.add_subsystem("ivc", om.IndepVarComp())
        ivc.add_output("alpha", val=1.0)
        ivc.add_output("mach", val=0.6)
        ivc.add_output("wing_widths", val=[1.0, 1.0, 1.0])
        ivc.add_output("v", val=1.0)
        ivc.add_output("rho", val=1.0)

        group.add_subsystem("funcs", AeroFuncsGroup(surfaces=[surfaces[0]], write_solution=False))
        group.promotes("funcs", [(f"{surfaces[0]['name']}.widths", f"{surfaces[0]['name']}_widths")])

        group.connect("ivc.alpha", "funcs.aoa")
        group.connect("ivc.mach", "funcs.mach")
        group.connect("ivc.rho", "funcs.rho")
        group.connect("ivc.v", "funcs.v")
        group.connect("ivc.wing_widths", f"{surfaces[0]['name']}_widths")

        run_test(self, group, complex_flag=True, method="cs")


if __name__ == "__main__":
    unittest.main()
