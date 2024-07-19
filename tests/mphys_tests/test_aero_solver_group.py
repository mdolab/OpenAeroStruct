import unittest

import openmdao.api as om

from openaerostruct.mphys.aero_solver_group import AeroSolverGroup
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def test_incompressible(self):
        group = self.setup_solver(compressible=False)
        run_test(self, group, complex_flag=True, method="cs", atol=1e-4, rtol=1e-5)

    def test_compressible(self):
        group = self.setup_solver(compressible=True)
        run_test(self, group, complex_flag=True, method="cs", atol=1e-4, rtol=1e-5)

    def setup_solver(self, compressible=True):
        surfaces = get_default_surfaces()

        group = om.Group()

        ivc = group.add_subsystem("ivc", om.IndepVarComp())
        ivc.add_output(f"{surfaces[0]['name']}_def_mesh", val=surfaces[0]["mesh"])
        ivc.add_output("alpha", val=1.0)
        ivc.add_output("mach", val=0.6)

        group.add_subsystem("solver", AeroSolverGroup(surfaces=[surfaces[0]], compressible=compressible))
        group.promotes("solver", [(f"{surfaces[0]['name']}.normals", f"{surfaces[0]['name']}_normals")])

        group.connect(f"ivc.{surfaces[0]['name']}_def_mesh", f"solver.{surfaces[0]['name']}_def_mesh")
        group.connect("ivc.alpha", "solver.aoa")
        if compressible:
            group.connect("ivc.mach", "solver.mach")

        return group


if __name__ == "__main__":
    unittest.main()
