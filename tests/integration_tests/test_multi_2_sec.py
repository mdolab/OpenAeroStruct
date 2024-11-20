from openmdao.utils.assert_utils import assert_near_equal
import unittest


class Test(unittest.TestCase):
    def assert_check_results_sym(self, prob):
        """We put the assert_near_equal for the test_constraint and test construction cases here since they should both return the same results."""
        assert_near_equal(prob["aero_point_0.surface_perf.CD"][0], 0.02920058, 1e-3)
        assert_near_equal(prob["aero_point_0.surface_perf.CL"][0], 0.2999996, 1e-3)
        assert_near_equal(prob["aero_point_0.CM"][1], -0.07335945, 1e-3)

    def test_constraint(self):
        import numpy as np

        import openmdao.api as om

        from openaerostruct.geometry.geometry_group import MultiSecGeometry
        from openaerostruct.aerodynamics.aero_groups import AeroPoint
        from openaerostruct.geometry.geometry_group import build_sections
        from openaerostruct.geometry.geometry_unification import unify_mesh
        from openaerostruct.utils.testing import get_two_section_surface

        surface, sec_chord_cp = get_two_section_surface()

        # Create the OpenMDAO problem for the constrained version
        prob = om.Problem()

        # Create an independent variable component that will supply the flow
        # conditions to the problem.
        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output("v", val=1.0, units="m/s")
        indep_var_comp.add_output("alpha", val=10.0, units="deg")
        indep_var_comp.add_output("Mach_number", val=0.3)
        indep_var_comp.add_output("re", val=1.0e5, units="1/m")
        indep_var_comp.add_output("rho", val=0.38, units="kg/m**3")
        indep_var_comp.add_output("cg", val=np.zeros((3)), units="m")

        # Add this IndepVarComp to the problem model
        prob.model.add_subsystem("prob_vars", indep_var_comp, promotes=["*"])

        # Generate the sections and unified mesh here in addition to adding the components.
        # This has to also be done here since AeroPoint has to know the unified mesh size.
        section_surfaces = build_sections(surface)
        uni_mesh = unify_mesh(section_surfaces)
        surface["mesh"] = uni_mesh

        # Create and add a group that handles the geometry for the
        # aerodynamic lifting surface
        multi_geom_group = MultiSecGeometry(
            surface=surface, joining_comp=True, dim_constr=[np.array([1, 0, 0]), np.array([1, 0, 0])]
        )
        prob.model.add_subsystem(surface["name"], multi_geom_group)

        # Create the aero point group, which contains the actual aerodynamic
        # analyses
        aero_group = AeroPoint(surfaces=[surface])
        point_name = "aero_point_0"
        prob.model.add_subsystem(
            point_name, aero_group, promotes_inputs=["v", "alpha", "Mach_number", "re", "rho", "cg"]
        )

        # Get name of surface and construct unified mesh name
        name = surface["name"]
        unification_name = "{}_unification".format(surface["name"])

        # Connect the mesh from the mesh unification component to the analysis point
        prob.model.connect(
            name + "." + unification_name + "." + name + "_uni_mesh", point_name + "." + "surface" + ".def_mesh"
        )

        # Perform the connections with the modified names within the
        # 'aero_states' group.
        prob.model.connect(
            name + "." + unification_name + "." + name + "_uni_mesh",
            point_name + ".aero_states." + "surface" + "_def_mesh",
        )

        # Add DVs
        prob.model.add_design_var("surface.sec0.chord_cp", lower=0.1, upper=10.0, units=None)
        prob.model.add_design_var("surface.sec1.chord_cp", lower=0.1, upper=10.0, units=None)
        prob.model.add_design_var("alpha", lower=0.0, upper=10.0, units="deg")

        # Add joined mesh constraint
        prob.model.add_constraint("surface.surface_joining.section_separation", upper=0, lower=0)

        # Add CL constraint
        prob.model.add_constraint(point_name + ".CL", equals=0.3)

        # Add Area constraint
        prob.model.add_constraint(point_name + ".total_perf.S_ref_total", equals=2.0)

        # Add objective
        prob.model.add_objective(point_name + ".CD", scaler=1e4)

        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options["optimizer"] = "SLSQP"
        prob.driver.options["tol"] = 1e-7
        prob.driver.options["disp"] = True
        prob.driver.options["maxiter"] = 1000
        # prob.driver.options["debug_print"] = ["nl_cons", "objs", "desvars"]

        # Set up and run the optimization problem
        prob.setup()
        prob.run_driver()

        self.assert_check_results_sym(prob)

    def test_construction(self):
        import numpy as np

        import openmdao.api as om

        from openaerostruct.geometry.geometry_group import MultiSecGeometry
        from openaerostruct.aerodynamics.aero_groups import AeroPoint
        from openaerostruct.geometry.geometry_group import build_sections
        from openaerostruct.geometry.geometry_unification import unify_mesh
        from openaerostruct.geometry.multi_unified_bspline_utils import build_multi_spline, connect_multi_spline
        from openaerostruct.utils.testing import get_two_section_surface

        surface, sec_chord_cp = get_two_section_surface()

        # Create the OpenMDAO problem
        prob = om.Problem()

        # Create an independent variable component that will supply the flow
        # conditions to the problem.
        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output("v", val=1.0, units="m/s")
        indep_var_comp.add_output("alpha", val=10.0, units="deg")
        indep_var_comp.add_output("Mach_number", val=0.3)
        indep_var_comp.add_output("re", val=1.0e5, units="1/m")
        indep_var_comp.add_output("rho", val=0.38, units="kg/m**3")
        indep_var_comp.add_output("cg", val=np.zeros((3)), units="m")

        # Add this IndepVarComp to the problem model
        prob.model.add_subsystem("prob_vars", indep_var_comp, promotes=["*"])

        # Generate the sections and unified mesh here. It's needed to join the sections by construction.
        section_surfaces = build_sections(surface)
        uni_mesh = unify_mesh(section_surfaces)
        surface["mesh"] = uni_mesh

        # Build a component with B-spline control points that joins the sections by construction
        chord_comp = build_multi_spline("chord_cp", len(section_surfaces), sec_chord_cp)
        prob.model.add_subsystem("chord_bspline", chord_comp)

        # Connect the B-spline component to the section B-splines
        connect_multi_spline(prob, section_surfaces, sec_chord_cp, "chord_cp", "chord_bspline")

        # Create and add a group that handles the geometry for the
        # aerodynamic lifting surface
        multi_geom_group = MultiSecGeometry(surface=surface)
        prob.model.add_subsystem(surface["name"], multi_geom_group)

        # Create the aero point group, which contains the actual aerodynamic
        # analyses
        aero_group = AeroPoint(surfaces=[surface])
        point_name = "aero_point_0"
        prob.model.add_subsystem(
            point_name, aero_group, promotes_inputs=["v", "alpha", "Mach_number", "re", "rho", "cg"]
        )

        # Get name of surface and construct unified mesh name
        name = surface["name"]
        unification_name = "{}_unification".format(surface["name"])

        # Connect the mesh from the mesh unification component to the analysis point
        prob.model.connect(
            name + "." + unification_name + "." + name + "_uni_mesh", point_name + "." + "surface" + ".def_mesh"
        )

        # Perform the connections with the modified names within the
        # 'aero_states' group.
        prob.model.connect(
            name + "." + unification_name + "." + name + "_uni_mesh",
            point_name + ".aero_states." + "surface" + "_def_mesh",
        )

        # Add DVs
        prob.model.add_design_var("chord_bspline.chord_cp_spline", lower=0.1, upper=10.0, units=None)
        prob.model.add_design_var("alpha", lower=0.0, upper=10.0, units="deg")

        # Add CL constraint
        prob.model.add_constraint(point_name + ".CL", equals=0.3)

        # Add Area constraint
        prob.model.add_constraint(point_name + ".total_perf.S_ref_total", equals=2.0)

        # Add objective
        prob.model.add_objective(point_name + ".CD", scaler=1e4)

        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options["optimizer"] = "SLSQP"
        prob.driver.options["tol"] = 1e-7
        prob.driver.options["disp"] = True
        prob.driver.options["maxiter"] = 1000
        # prob.driver.options["debug_print"] = ["nl_cons", "objs", "desvars"]

        # Set up and run the optimization problem
        prob.setup()
        prob.run_driver()

        self.assert_check_results_sym(prob)

    def test_asymmetrical(self):
        import numpy as np

        import openmdao.api as om

        from openaerostruct.geometry.geometry_group import MultiSecGeometry
        from openaerostruct.aerodynamics.aero_groups import AeroPoint
        from openaerostruct.geometry.geometry_group import build_sections
        from openaerostruct.geometry.geometry_unification import unify_mesh
        from openaerostruct.utils.testing import get_two_section_surface

        surface, sec_chord_cp = get_two_section_surface(sym=False, visc=False)

        # Create the OpenMDAO problem
        prob = om.Problem()

        # Create an independent variable component that will supply the flow
        # conditions to the problem.
        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output("v", val=1.0, units="m/s")
        indep_var_comp.add_output("alpha", val=10.0, units="deg")
        indep_var_comp.add_output("Mach_number", val=0.3)
        indep_var_comp.add_output("re", val=1.0e5, units="1/m")
        indep_var_comp.add_output("rho", val=0.38, units="kg/m**3")
        indep_var_comp.add_output("cg", val=np.zeros((3)), units="m")

        # Add this IndepVarComp to the problem model
        prob.model.add_subsystem("prob_vars", indep_var_comp, promotes=["*"])

        # Create and add a group that handles the geometry for the
        # aerodynamic lifting surface
        multi_geom_group = MultiSecGeometry(
            surface=surface, joining_comp=True, dim_constr=[np.array([1, 0, 0]), np.array([1, 0, 0])]
        )
        prob.model.add_subsystem(surface["name"], multi_geom_group)

        # Generate the sections and unified mesh here in addition to adding the components.
        # This has to ALSO be done here since AeroPoint has to know the unified mesh size.
        section_surfaces = build_sections(surface)
        uni_mesh = unify_mesh(section_surfaces)
        surface["mesh"] = uni_mesh

        # Create the aero point group, which contains the actual aerodynamic
        # analyses
        aero_group = AeroPoint(surfaces=[surface])
        point_name = "aero_point_0"
        prob.model.add_subsystem(
            point_name, aero_group, promotes_inputs=["v", "alpha", "Mach_number", "re", "rho", "cg"]
        )

        # Get name of surface and construct unified mesh name
        name = surface["name"]
        unification_name = "{}_unification".format(surface["name"])

        # Connect the mesh from the mesh unification component to the analysis point
        prob.model.connect(
            name + "." + unification_name + "." + name + "_uni_mesh", point_name + "." + "surface" + ".def_mesh"
        )

        # Perform the connections with the modified names within the
        # 'aero_states' group.
        prob.model.connect(
            name + "." + unification_name + "." + name + "_uni_mesh",
            point_name + ".aero_states." + "surface" + "_def_mesh",
        )

        prob.model.add_design_var("surface.sec0.chord_cp", lower=0.1, upper=10.0, units=None)
        prob.model.add_design_var("surface.sec1.chord_cp", lower=0.1, upper=10.0, units=None)
        prob.model.add_design_var("alpha", lower=0.0, upper=10.0, units="deg")

        # Add joined mesh constraint
        prob.model.add_constraint("surface.surface_joining.section_separation", upper=0, lower=0)
        # prob.model.add_constraint('surface.surface_joining.section_separation',equals=0.0,scaler=1e-4)

        # Add CL constraint
        prob.model.add_constraint(point_name + ".CL", equals=0.3)

        # Add Area constraint
        prob.model.add_constraint(point_name + ".total_perf.S_ref_total", equals=2.0)

        # Add objective
        prob.model.add_objective(point_name + ".CD", scaler=1e4)

        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options["optimizer"] = "SLSQP"
        prob.driver.options["tol"] = 1e-3
        prob.driver.options["disp"] = True
        prob.driver.options["maxiter"] = 1000
        # prob.driver.options["debug_print"] = ["nl_cons", "objs", "desvars"]

        # Set up and run the optimization problem
        prob.setup()
        prob.run_driver()

        assert_near_equal(prob["aero_point_0.surface_perf.CD"][0], 0.02902683, 1e-6)
        assert_near_equal(prob["aero_point_0.surface_perf.CL"][0], 0.29999999, 1e-6)
        assert_near_equal(prob["aero_point_0.CM"][1], -0.07375071, 1e-6)

    def test_visc(self):
        import numpy as np

        import openmdao.api as om

        from openaerostruct.geometry.geometry_group import MultiSecGeometry
        from openaerostruct.aerodynamics.aero_groups import AeroPoint
        from openaerostruct.geometry.geometry_group import build_sections
        from openaerostruct.geometry.geometry_unification import unify_mesh
        from openaerostruct.geometry.multi_unified_bspline_utils import build_multi_spline, connect_multi_spline
        from openaerostruct.utils.testing import get_two_section_surface

        surface, sec_chord_cp = get_two_section_surface(sym=True, visc=True)

        # Create the OpenMDAO problem
        prob = om.Problem()

        # Create an independent variable component that will supply the flow
        # conditions to the problem.
        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output("v", val=1.0, units="m/s")
        indep_var_comp.add_output("alpha", val=10.0, units="deg")
        indep_var_comp.add_output("Mach_number", val=0.3)
        indep_var_comp.add_output("re", val=1.0e5, units="1/m")
        indep_var_comp.add_output("rho", val=0.38, units="kg/m**3")
        indep_var_comp.add_output("cg", val=np.zeros((3)), units="m")

        # Add this IndepVarComp to the problem model
        prob.model.add_subsystem("prob_vars", indep_var_comp, promotes=["*"])

        # Generate the sections and unified mesh here. It's needed to join the sections by construction.
        section_surfaces = build_sections(surface)
        uni_mesh = unify_mesh(section_surfaces)
        surface["mesh"] = uni_mesh

        # Build a component with B-spline control points that joins the sections by construction
        chord_comp = build_multi_spline("chord_cp", len(section_surfaces), sec_chord_cp)
        prob.model.add_subsystem("chord_bspline", chord_comp)

        # Connect the B-spline component to the section B-splines
        connect_multi_spline(prob, section_surfaces, sec_chord_cp, "chord_cp", "chord_bspline")

        # Create and add a group that handles the geometry for the
        # aerodynamic lifting surface
        multi_geom_group = MultiSecGeometry(surface=surface)
        prob.model.add_subsystem(surface["name"], multi_geom_group)

        # Create the aero point group, which contains the actual aerodynamic
        # analyses
        aero_group = AeroPoint(surfaces=[surface])
        point_name = "aero_point_0"
        prob.model.add_subsystem(
            point_name, aero_group, promotes_inputs=["v", "alpha", "Mach_number", "re", "rho", "cg"]
        )

        # Get name of surface and construct unified mesh name
        name = surface["name"]
        unification_name = "{}_unification".format(surface["name"])

        # Connect the mesh from the mesh unification component to the analysis point
        prob.model.connect(
            name + "." + unification_name + "." + name + "_uni_mesh", point_name + "." + "surface" + ".def_mesh"
        )

        # Perform the connections with the modified names within the
        # 'aero_states' group.
        prob.model.connect(
            name + "." + unification_name + "." + name + "_uni_mesh",
            point_name + ".aero_states." + "surface" + "_def_mesh",
        )

        # Connect t over c B-spline for viscous drag
        prob.model.connect(
            name + "." + unification_name + "." + name + "_uni_t_over_c",
            point_name + "." + name + "_perf." + "t_over_c",
        )

        # Add DVs
        prob.model.add_design_var("chord_bspline.chord_cp_spline", lower=0.1, upper=10.0, units=None)
        prob.model.add_design_var("alpha", lower=0.0, upper=10.0, units="deg")

        # Add CL constraint
        prob.model.add_constraint(point_name + ".CL", equals=0.3)

        # Add Area constraint
        prob.model.add_constraint(point_name + ".total_perf.S_ref_total", equals=2.0)

        # Add objective
        prob.model.add_objective(point_name + ".CD", scaler=1e4)

        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options["optimizer"] = "SLSQP"
        prob.driver.options["tol"] = 1e-7
        prob.driver.options["disp"] = True
        prob.driver.options["maxiter"] = 1000
        # prob.driver.options["debug_print"] = ["nl_cons", "objs", "desvars"]

        # Set up and run the optimization problem
        prob.setup()
        prob.run_driver()

        assert_near_equal(prob["aero_point_0.surface_perf.CD"][0], 0.04946752, 1e-6)
        assert_near_equal(prob["aero_point_0.surface_perf.CL"][0], 0.3, 1e-6)
        assert_near_equal(prob["aero_point_0.CM"][1], -0.05045758, 1e-6)


if __name__ == "__main__":
    unittest.main()
