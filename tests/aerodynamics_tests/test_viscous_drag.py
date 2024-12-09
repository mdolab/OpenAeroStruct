import unittest

from openaerostruct.aerodynamics.viscous_drag import ViscousDrag
from openaerostruct.utils.testing import run_test, get_default_surfaces
import openmdao.api as om
import numpy as np


class Test(unittest.TestCase):
    def test(self):
        surface = get_default_surfaces()[0]
        surface["t_over_c_cp"] = np.array([0.1, 0.15, 0.2])

        ny = surface["mesh"].shape[1]
        n_cp = len(surface["t_over_c_cp"])

        group = om.Group()

        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output("t_over_c_cp", val=surface["t_over_c_cp"])
        group.add_subsystem("indep_var_comp", indep_var_comp, promotes=["*"])

        x_interp = np.linspace(0.0, 1.0, int(ny - 1))
        comp = group.add_subsystem(
            "t_over_c_bsp",
            om.SplineComp(
                method="bsplines", x_interp_val=x_interp, num_cp=n_cp, interp_options={"order": min(n_cp, 4)}
            ),
            promotes_inputs=["t_over_c_cp"],
            promotes_outputs=["t_over_c"],
        )
        comp.add_spline(y_cp_name="t_over_c_cp", y_interp_name="t_over_c", y_cp_val=np.zeros(n_cp))

        comp = ViscousDrag(surface=surface, with_viscous=True)
        group.add_subsystem("viscousdrag", comp, promotes=["*"])

        run_test(self, group, complex_flag=True)

    def test_2(self):
        surface = get_default_surfaces()[0]
        surface["k_lam"] = 0.5

        ny = surface["mesh"].shape[1]
        n_cp = len(surface["t_over_c_cp"])

        group = om.Group()

        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output("t_over_c_cp", val=surface["t_over_c_cp"])
        group.add_subsystem("indep_var_comp", indep_var_comp, promotes=["*"])

        x_interp = np.linspace(0.0, 1.0, int(ny - 1))
        comp = group.add_subsystem(
            "t_over_c_bsp",
            om.SplineComp(
                method="bsplines", x_interp_val=x_interp, num_cp=n_cp, interp_options={"order": min(n_cp, 4)}
            ),
            promotes_inputs=["t_over_c_cp"],
            promotes_outputs=["t_over_c"],
        )
        comp.add_spline(y_cp_name="t_over_c_cp", y_interp_name="t_over_c")

        comp = ViscousDrag(surface=surface, with_viscous=True)
        group.add_subsystem("viscousdrag", comp, promotes=["*"])

        run_test(self, group, complex_flag=True)


if __name__ == "__main__":
    unittest.main()
