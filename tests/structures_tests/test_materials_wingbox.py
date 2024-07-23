import unittest
import numpy as np
from openaerostruct.structures.section_properties_wingbox import SectionPropertiesWingbox
from openaerostruct.utils.testing import run_test, get_default_surfaces
from openmdao.utils.assert_utils import assert_near_equal
import openmdao.api as om


class Test(unittest.TestCase):
    @unittest.skip("Analytic derivatives not added yet.")
    def test(self):
        """
        This is for checking the partials.
        """

        surface = get_default_surfaces()[0]

        # turn down some of these properties, so the absolute deriv error isn't magnified
        surface["E"] = 7
        surface["G"] = 3
        surface["yield"] = 0.02

        surface["data_x_upper"] = np.array(
            [
                0.1,
                0.11,
                0.12,
                0.13,
                0.14,
                0.15,
                0.16,
                0.17,
                0.18,
                0.19,
                0.2,
                0.21,
                0.22,
                0.23,
                0.24,
                0.25,
                0.26,
                0.27,
                0.28,
                0.29,
                0.3,
                0.31,
                0.32,
                0.33,
                0.34,
                0.35,
                0.36,
                0.37,
                0.38,
                0.39,
                0.4,
                0.41,
                0.42,
                0.43,
                0.44,
                0.45,
                0.46,
                0.47,
                0.48,
                0.49,
                0.5,
                0.51,
                0.52,
                0.53,
                0.54,
                0.55,
                0.56,
                0.57,
                0.58,
                0.59,
                0.6,
            ],
            dtype="complex128",
        )

        surface["data_x_lower"] = np.array(
            [
                0.1,
                0.11,
                0.12,
                0.13,
                0.14,
                0.15,
                0.16,
                0.17,
                0.18,
                0.19,
                0.2,
                0.21,
                0.22,
                0.23,
                0.24,
                0.25,
                0.26,
                0.27,
                0.28,
                0.29,
                0.3,
                0.31,
                0.32,
                0.33,
                0.34,
                0.35,
                0.36,
                0.37,
                0.38,
                0.39,
                0.4,
                0.41,
                0.42,
                0.43,
                0.44,
                0.45,
                0.46,
                0.47,
                0.48,
                0.49,
                0.5,
                0.51,
                0.52,
                0.53,
                0.54,
                0.55,
                0.56,
                0.57,
                0.58,
                0.59,
                0.6,
            ],
            dtype="complex128",
        )
        surface["data_y_upper"] = np.array(
            [
                0.0447,
                0.046,
                0.0472,
                0.0484,
                0.0495,
                0.0505,
                0.0514,
                0.0523,
                0.0531,
                0.0538,
                0.0545,
                0.0551,
                0.0557,
                0.0563,
                0.0568,
                0.0573,
                0.0577,
                0.0581,
                0.0585,
                0.0588,
                0.0591,
                0.0593,
                0.0595,
                0.0597,
                0.0599,
                0.06,
                0.0601,
                0.0602,
                0.0602,
                0.0602,
                0.0602,
                0.0602,
                0.0601,
                0.06,
                0.0599,
                0.0598,
                0.0596,
                0.0594,
                0.0592,
                0.0589,
                0.0586,
                0.0583,
                0.058,
                0.0576,
                0.0572,
                0.0568,
                0.0563,
                0.0558,
                0.0553,
                0.0547,
                0.0541,
            ],
            dtype="complex128",
        )
        surface["data_y_lower"] = np.array(
            [
                -0.0447,
                -0.046,
                -0.0473,
                -0.0485,
                -0.0496,
                -0.0506,
                -0.0515,
                -0.0524,
                -0.0532,
                -0.054,
                -0.0547,
                -0.0554,
                -0.056,
                -0.0565,
                -0.057,
                -0.0575,
                -0.0579,
                -0.0583,
                -0.0586,
                -0.0589,
                -0.0592,
                -0.0594,
                -0.0595,
                -0.0596,
                -0.0597,
                -0.0598,
                -0.0598,
                -0.0598,
                -0.0598,
                -0.0597,
                -0.0596,
                -0.0594,
                -0.0592,
                -0.0589,
                -0.0586,
                -0.0582,
                -0.0578,
                -0.0573,
                -0.0567,
                -0.0561,
                -0.0554,
                -0.0546,
                -0.0538,
                -0.0529,
                -0.0519,
                -0.0509,
                -0.0497,
                -0.0485,
                -0.0472,
                -0.0458,
                -0.0444,
            ],
            dtype="complex128",
        )
        surface["original_wingbox_airfoil_t_over_c"] = 0.1
        comp = SectionPropertiesWingbox(surface=surface)

        run_test(self, comp, complex_flag=True, method="cs", step=1e-40)

    def test2(self):
        """
        This is for checking the computation.
        """

        surface = get_default_surfaces()[0]
        surface["t_over_c_cp"] = np.array([0.1, 0.15, 0.2])
        surface["spar_thickness_cp"] = np.array([0.004, 0.008, 0.02])
        surface["skin_thickness_cp"] = np.array([0.01, 0.015, 0.021])
        surface["fem_chords_cp"] = np.array([2.0, 3.0, 4.0])
        surface["streamwise_chords_cp"] = np.array([3.0, 4.0, 5.0])
        surface["fem_twists_cp"] = np.array([5.0, 3.0, 2.0]) / 180.0 * np.pi

        surface["data_x_upper"] = np.array(
            [
                0.1,
                0.11,
                0.12,
                0.13,
                0.14,
                0.15,
                0.16,
                0.17,
                0.18,
                0.19,
                0.2,
                0.21,
                0.22,
                0.23,
                0.24,
                0.25,
                0.26,
                0.27,
                0.28,
                0.29,
                0.3,
                0.31,
                0.32,
                0.33,
                0.34,
                0.35,
                0.36,
                0.37,
                0.38,
                0.39,
                0.4,
                0.41,
                0.42,
                0.43,
                0.44,
                0.45,
                0.46,
                0.47,
                0.48,
                0.49,
                0.5,
                0.51,
                0.52,
                0.53,
                0.54,
                0.55,
                0.56,
                0.57,
                0.58,
                0.59,
                0.6,
            ],
            dtype="complex128",
        )
        surface["data_x_lower"] = np.array(
            [
                0.1,
                0.11,
                0.12,
                0.13,
                0.14,
                0.15,
                0.16,
                0.17,
                0.18,
                0.19,
                0.2,
                0.21,
                0.22,
                0.23,
                0.24,
                0.25,
                0.26,
                0.27,
                0.28,
                0.29,
                0.3,
                0.31,
                0.32,
                0.33,
                0.34,
                0.35,
                0.36,
                0.37,
                0.38,
                0.39,
                0.4,
                0.41,
                0.42,
                0.43,
                0.44,
                0.45,
                0.46,
                0.47,
                0.48,
                0.49,
                0.5,
                0.51,
                0.52,
                0.53,
                0.54,
                0.55,
                0.56,
                0.57,
                0.58,
                0.59,
                0.6,
            ],
            dtype="complex128",
        )
        surface["data_y_upper"] = np.array(
            [
                0.0447,
                0.046,
                0.0472,
                0.0484,
                0.0495,
                0.0505,
                0.0514,
                0.0523,
                0.0531,
                0.0538,
                0.0545,
                0.0551,
                0.0557,
                0.0563,
                0.0568,
                0.0573,
                0.0577,
                0.0581,
                0.0585,
                0.0588,
                0.0591,
                0.0593,
                0.0595,
                0.0597,
                0.0599,
                0.06,
                0.0601,
                0.0602,
                0.0602,
                0.0602,
                0.0602,
                0.0602,
                0.0601,
                0.06,
                0.0599,
                0.0598,
                0.0596,
                0.0594,
                0.0592,
                0.0589,
                0.0586,
                0.0583,
                0.058,
                0.0576,
                0.0572,
                0.0568,
                0.0563,
                0.0558,
                0.0553,
                0.0547,
                0.0541,
            ],
            dtype="complex128",
        )
        surface["data_y_lower"] = np.array(
            [
                -0.0447,
                -0.046,
                -0.0473,
                -0.0485,
                -0.0496,
                -0.0506,
                -0.0515,
                -0.0524,
                -0.0532,
                -0.054,
                -0.0547,
                -0.0554,
                -0.056,
                -0.0565,
                -0.057,
                -0.0575,
                -0.0579,
                -0.0583,
                -0.0586,
                -0.0589,
                -0.0592,
                -0.0594,
                -0.0595,
                -0.0596,
                -0.0597,
                -0.0598,
                -0.0598,
                -0.0598,
                -0.0598,
                -0.0597,
                -0.0596,
                -0.0594,
                -0.0592,
                -0.0589,
                -0.0586,
                -0.0582,
                -0.0578,
                -0.0573,
                -0.0567,
                -0.0561,
                -0.0554,
                -0.0546,
                -0.0538,
                -0.0529,
                -0.0519,
                -0.0509,
                -0.0497,
                -0.0485,
                -0.0472,
                -0.0458,
                -0.0444,
            ],
            dtype="complex128",
        )
        surface["original_wingbox_airfoil_t_over_c"] = 0.1

        mesh = surface["mesh"]
        ny = mesh.shape[1]
        n_cp = len(surface["t_over_c_cp"])

        prob = om.Problem()

        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output("t_over_c_cp", val=surface["t_over_c_cp"])
        indep_var_comp.add_output("spar_thickness_cp", val=surface["spar_thickness_cp"], units="m")
        indep_var_comp.add_output("skin_thickness_cp", val=surface["skin_thickness_cp"], units="m")
        indep_var_comp.add_output("fem_chords_cp", val=surface["fem_chords_cp"], units="m")
        indep_var_comp.add_output("streamwise_chords_cp", val=surface["streamwise_chords_cp"], units="m")
        indep_var_comp.add_output("fem_twists_cp", val=surface["fem_twists_cp"], units="deg")
        prob.model.add_subsystem("indep_var_comp", indep_var_comp, promotes=["*"])

        x_interp = np.linspace(0.0, 1.0, int(ny - 1))
        comp = prob.model.add_subsystem(
            "bsplines_comp",
            om.SplineComp(
                method="bsplines", x_interp_val=x_interp, num_cp=n_cp, interp_options={"order": min(n_cp, 4)}
            ),
            promotes_inputs=["*"],
            promotes_outputs=["*"],
        )
        comp.add_spline(y_cp_name="t_over_c_cp", y_interp_name="t_over_c")
        comp.add_spline(y_cp_name="skin_thickness_cp", y_interp_name="skin_thickness", y_units="m")
        comp.add_spline(y_cp_name="spar_thickness_cp", y_interp_name="spar_thickness", y_units="m")
        comp.add_spline(y_cp_name="fem_chords_cp", y_interp_name="fem_chords", y_units="m")
        comp.add_spline(y_cp_name="streamwise_chords_cp", y_interp_name="streamwise_chords", y_units="m")
        comp.add_spline(y_cp_name="fem_twists_cp", y_interp_name="fem_twists", y_units="deg")

        comp = SectionPropertiesWingbox(surface=surface)
        prob.model.add_subsystem("sec_prop_wb", comp, promotes=["*"])

        prob.setup()
        # om.view_model(prob)

        prob.run_model()

        # print( prob['A'] )
        # print( prob['A_enc'] )
        # print( prob['A_int'] )
        # print( prob['Iy'] )
        # print( prob['Qz'] )
        # print( prob['Iz'] )
        # print( prob['J'] )
        # print( prob['htop'] )
        # print( prob['hbottom'] )
        # print( prob['hfront'] )
        # print( prob['hrear'] )

        assert_near_equal(prob["A"], np.array([0.02203548, 0.0563726, 0.11989703]), 1e-6)
        assert_near_equal(prob["A_enc"], np.array([0.3243776, 0.978003, 2.17591]), 1e-6)
        assert_near_equal(prob["A_int"], np.array([0.3132502, 0.949491, 2.11512]), 1e-6)
        assert_near_equal(prob["Iy"], np.array([0.00218612, 0.01455083, 0.06342765]), 1e-6)
        assert_near_equal(prob["Qz"], np.array([0.00169233, 0.00820558, 0.02707493]), 1e-6)
        assert_near_equal(prob["Iz"], np.array([0.00055292, 0.00520911, 0.02785168]), 1e-6)
        assert_near_equal(prob["J"], np.array([0.00124939, 0.01241967, 0.06649673]), 1e-6)
        assert_near_equal(prob["htop"], np.array([0.19106873, 0.36005945, 0.5907887]), 1e-6)
        assert_near_equal(prob["hbottom"], np.array([0.19906584, 0.37668887, 0.61850335]), 1e-6)
        assert_near_equal(prob["hfront"], np.array([0.52341176, 0.78649186, 1.04902676]), 1e-6)
        assert_near_equal(prob["hrear"], np.array([0.47524073, 0.71429312, 0.95303545]), 1e-6)

    def test3(self):
        """
        This is an extreme nonphysical case (large twist angles) for checking the computation.
        """

        surface = get_default_surfaces()[0]
        surface["t_over_c_cp"] = np.array([0.1, 0.15, 0.2])
        surface["spar_thickness_cp"] = np.array([0.004, 0.008, 0.02])
        surface["skin_thickness_cp"] = np.array([0.01, 0.015, 0.021])
        surface["fem_chords_cp"] = np.array([2.0, 3.0, 4.0])
        surface["streamwise_chords_cp"] = np.array([3.0, 4.0, 5.0])
        surface["fem_twists_cp"] = np.array([5.0, 3.0, 2.0])

        surface["data_x_upper"] = np.array(
            [
                0.1,
                0.11,
                0.12,
                0.13,
                0.14,
                0.15,
                0.16,
                0.17,
                0.18,
                0.19,
                0.2,
                0.21,
                0.22,
                0.23,
                0.24,
                0.25,
                0.26,
                0.27,
                0.28,
                0.29,
                0.3,
                0.31,
                0.32,
                0.33,
                0.34,
                0.35,
                0.36,
                0.37,
                0.38,
                0.39,
                0.4,
                0.41,
                0.42,
                0.43,
                0.44,
                0.45,
                0.46,
                0.47,
                0.48,
                0.49,
                0.5,
                0.51,
                0.52,
                0.53,
                0.54,
                0.55,
                0.56,
                0.57,
                0.58,
                0.59,
                0.6,
            ],
            dtype="complex128",
        )
        surface["data_x_lower"] = np.array(
            [
                0.1,
                0.11,
                0.12,
                0.13,
                0.14,
                0.15,
                0.16,
                0.17,
                0.18,
                0.19,
                0.2,
                0.21,
                0.22,
                0.23,
                0.24,
                0.25,
                0.26,
                0.27,
                0.28,
                0.29,
                0.3,
                0.31,
                0.32,
                0.33,
                0.34,
                0.35,
                0.36,
                0.37,
                0.38,
                0.39,
                0.4,
                0.41,
                0.42,
                0.43,
                0.44,
                0.45,
                0.46,
                0.47,
                0.48,
                0.49,
                0.5,
                0.51,
                0.52,
                0.53,
                0.54,
                0.55,
                0.56,
                0.57,
                0.58,
                0.59,
                0.6,
            ],
            dtype="complex128",
        )
        surface["data_y_upper"] = np.array(
            [
                0.0447,
                0.046,
                0.0472,
                0.0484,
                0.0495,
                0.0505,
                0.0514,
                0.0523,
                0.0531,
                0.0538,
                0.0545,
                0.0551,
                0.0557,
                0.0563,
                0.0568,
                0.0573,
                0.0577,
                0.0581,
                0.0585,
                0.0588,
                0.0591,
                0.0593,
                0.0595,
                0.0597,
                0.0599,
                0.06,
                0.0601,
                0.0602,
                0.0602,
                0.0602,
                0.0602,
                0.0602,
                0.0601,
                0.06,
                0.0599,
                0.0598,
                0.0596,
                0.0594,
                0.0592,
                0.0589,
                0.0586,
                0.0583,
                0.058,
                0.0576,
                0.0572,
                0.0568,
                0.0563,
                0.0558,
                0.0553,
                0.0547,
                0.0541,
            ],
            dtype="complex128",
        )
        surface["data_y_lower"] = np.array(
            [
                -0.0447,
                -0.046,
                -0.0473,
                -0.0485,
                -0.0496,
                -0.0506,
                -0.0515,
                -0.0524,
                -0.0532,
                -0.054,
                -0.0547,
                -0.0554,
                -0.056,
                -0.0565,
                -0.057,
                -0.0575,
                -0.0579,
                -0.0583,
                -0.0586,
                -0.0589,
                -0.0592,
                -0.0594,
                -0.0595,
                -0.0596,
                -0.0597,
                -0.0598,
                -0.0598,
                -0.0598,
                -0.0598,
                -0.0597,
                -0.0596,
                -0.0594,
                -0.0592,
                -0.0589,
                -0.0586,
                -0.0582,
                -0.0578,
                -0.0573,
                -0.0567,
                -0.0561,
                -0.0554,
                -0.0546,
                -0.0538,
                -0.0529,
                -0.0519,
                -0.0509,
                -0.0497,
                -0.0485,
                -0.0472,
                -0.0458,
                -0.0444,
            ],
            dtype="complex128",
        )
        surface["original_wingbox_airfoil_t_over_c"] = 0.1

        mesh = surface["mesh"]
        ny = mesh.shape[1]
        n_cp = len(surface["t_over_c_cp"])

        prob = om.Problem()

        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output("t_over_c_cp", val=surface["t_over_c_cp"])
        indep_var_comp.add_output("spar_thickness_cp", val=surface["spar_thickness_cp"], units="m")
        indep_var_comp.add_output("skin_thickness_cp", val=surface["skin_thickness_cp"], units="m")
        indep_var_comp.add_output("fem_chords_cp", val=surface["fem_chords_cp"], units="m")
        indep_var_comp.add_output("streamwise_chords_cp", val=surface["streamwise_chords_cp"], units="m")
        indep_var_comp.add_output("fem_twists_cp", val=surface["fem_twists_cp"], units="deg")
        prob.model.add_subsystem("indep_var_comp", indep_var_comp, promotes=["*"])

        x_interp = np.linspace(0.0, 1.0, int(ny - 1))
        comp = prob.model.add_subsystem(
            "bsplines_comp",
            om.SplineComp(
                method="bsplines", x_interp_val=x_interp, num_cp=n_cp, interp_options={"order": min(n_cp, 4)}
            ),
            promotes_inputs=["*"],
            promotes_outputs=["*"],
        )
        comp.add_spline(y_cp_name="t_over_c_cp", y_interp_name="t_over_c")
        comp.add_spline(y_cp_name="skin_thickness_cp", y_interp_name="skin_thickness", y_units="m")
        comp.add_spline(y_cp_name="spar_thickness_cp", y_interp_name="spar_thickness", y_units="m")
        comp.add_spline(y_cp_name="fem_chords_cp", y_interp_name="fem_chords", y_units="m")
        comp.add_spline(y_cp_name="streamwise_chords_cp", y_interp_name="streamwise_chords", y_units="m")
        comp.add_spline(y_cp_name="fem_twists_cp", y_interp_name="fem_twists", y_units="deg")

        comp = SectionPropertiesWingbox(surface=surface)
        prob.model.add_subsystem("sec_prop_wb", comp, promotes=["*"])

        prob.setup()
        # om.view_model(prob)

        prob.run_model()

        # print( prob['A'] )
        # print( prob['A_enc'] )
        # print( prob['A_int'] )
        # print( prob['Iy'] )
        # print( prob['Qz'] )
        # print( prob['Iz'] )
        # print( prob['J'] )
        # print( prob['htop'] )
        # print( prob['hbottom'] )
        # print( prob['hfront'] )
        # print( prob['hrear'] )

        assert_near_equal(prob["A"], np.array([0.0058738, -0.05739528, -0.05042289]), 1e-6)
        assert_near_equal(prob["A_enc"], np.array([0.3243776, 0.978003, 2.17591]), 1e-6)
        assert_near_equal(prob["A_int"], np.array([0.3132502, 0.949491, 2.11512]), 1e-6)
        assert_near_equal(prob["Iy"], np.array([3.59803239e-05, -1.52910019e-02, -4.01035510e-03]), 1e-6)
        assert_near_equal(prob["Qz"], np.array([0.00129261, 0.00870662, 0.02500053]), 1e-6)
        assert_near_equal(prob["Iz"], np.array([0.00056586, -0.00582207, -0.02877714]), 1e-6)
        assert_near_equal(prob["J"], np.array([0.00124939, 0.01241967, 0.06649673]), 1e-6)
        assert_near_equal(prob["htop"], np.array([0.53933652, -0.23509863, 0.71255343]), 1e-6)
        assert_near_equal(prob["hbottom"], np.array([0.50366564, -0.19185349, 0.73525459]), 1e-6)
        assert_near_equal(prob["hfront"], np.array([0.13442747, -0.78514756, -0.3919784]), 1e-6)
        assert_near_equal(prob["hrear"], np.array([0.12219305, -0.71214916, -0.35484131]), 1e-6)


if __name__ == "__main__":
    unittest.main()
