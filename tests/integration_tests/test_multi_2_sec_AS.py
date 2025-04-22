from openmdao.utils.assert_utils import assert_near_equal
import unittest


class Test(unittest.TestCase):
    def assert_check_results_sym(self, prob):
        """We put the assert_near_equal for the test_constraint and test construction cases here since they should both return the same results."""
        assert_near_equal(prob["AS_point_0.fuelburn"][0], 109110.19599915, 1e-2)

    def test_constraint(self):
        import numpy as np

        import openmdao.api as om

        from openaerostruct.integration.aerostruct_groups import (
            AerostructPoint,
            MultiSecAerostructGeometry,
        )
        from openaerostruct.utils.constants import grav_constant
        from openaerostruct.geometry.utils import (
            build_section_dicts,
            unify_mesh,
        )
        from openaerostruct.utils.testing import get_two_section_surface_AS

        surface, sec_twist_cp = get_two_section_surface_AS()

        # Create the OpenMDAO problem
        prob = om.Problem(reports=False)

        # Create an independent variable component that will supply the flow
        # conditions to the problem.
        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output("v", val=248.136, units="m/s")
        indep_var_comp.add_output("alpha", val=9.0, units="deg")
        indep_var_comp.add_output("Mach_number", val=0.84)
        indep_var_comp.add_output("re", val=1.0e6, units="1/m")
        indep_var_comp.add_output("rho", val=0.38, units="kg/m**3")
        indep_var_comp.add_output("CT", val=grav_constant * 17.0e-6, units="1/s")
        indep_var_comp.add_output("R", val=11.165e6, units="m")
        indep_var_comp.add_output("W0", val=0.4 * 3e5, units="kg")
        indep_var_comp.add_output("speed_of_sound", val=295.4, units="m/s")
        indep_var_comp.add_output("load_factor", val=1.0)
        indep_var_comp.add_output("empty_cg", val=np.zeros((3)), units="m")

        # Add this IndepVarComp to the problem model
        prob.model.add_subsystem("prob_vars", indep_var_comp, promotes=["*"])

        # Generate the sections and unified mesh here. It's needed to join the sections by construction.
        section_surfaces = build_section_dicts(surface)
        uni_mesh = unify_mesh(section_surfaces)
        surface["mesh"] = uni_mesh

        # Create and add a group that handles the geometry for the
        # aerodynamic lifting surface
        multi_aerostruct_geom_group = MultiSecAerostructGeometry(
            surface=surface,
            joining_comp=True,
            dim_constr=[np.array([1, 0, 0]), np.array([1, 0, 0])],
            shift_uni_mesh=False,
        )
        prob.model.add_subsystem(surface["name"], multi_aerostruct_geom_group)

        name = surface["name"]

        # Create the aero point group, which contains the actual aerodynamic
        # analyses
        aerostruct_group = AerostructPoint(surfaces=[surface])
        point_name = "AS_point_0"
        prob.model.add_subsystem(point_name, aerostruct_group)

        # Connect Flight conditions
        prob.model.connect("v", point_name + ".v")
        prob.model.connect("alpha", point_name + ".alpha")
        prob.model.connect("Mach_number", point_name + ".Mach_number")
        prob.model.connect("re", point_name + ".re")
        prob.model.connect("rho", point_name + ".rho")
        prob.model.connect("CT", point_name + ".CT")
        prob.model.connect("R", point_name + ".R")
        prob.model.connect("W0", point_name + ".W0")
        prob.model.connect("speed_of_sound", point_name + ".speed_of_sound")
        prob.model.connect("empty_cg", point_name + ".empty_cg")
        prob.model.connect("load_factor", point_name + ".load_factor")

        # Connect aerostruct groups
        com_name = point_name + "." + name + "_perf"
        prob.model.connect(
            name + ".local_stiff_transformed", point_name + ".coupled." + name + ".local_stiff_transformed"
        )
        prob.model.connect(name + ".nodes", point_name + ".coupled." + name + ".nodes")

        # Connect aerodyamic mesh to coupled group mesh
        prob.model.connect(name + ".mesh", point_name + ".coupled." + name + ".mesh")

        # Connect performance calculation variables
        prob.model.connect(name + ".radius", com_name + ".radius")
        prob.model.connect(name + ".thickness", com_name + ".thickness")
        prob.model.connect(name + ".nodes", com_name + ".nodes")
        prob.model.connect(name + ".cg_location", point_name + "." + "total_perf." + name + "_cg_location")
        prob.model.connect(name + ".structural_mass", point_name + "." + "total_perf." + name + "_structural_mass")
        prob.model.connect(name + ".t_over_c", com_name + ".t_over_c")

        # Add DVs
        prob.model.add_design_var("alpha", lower=-10.0, upper=10.0)
        prob.model.add_design_var("surface.sec0.twist_cp", lower=-10.0, upper=15.0)
        prob.model.add_design_var("surface.sec1.twist_cp", lower=-10.0, upper=15.0)
        prob.model.add_design_var("surface.thickness_cp", lower=0.01, upper=0.5, scaler=1e2)

        # Add constraints
        prob.model.add_constraint(point_name + ".L_equals_W", equals=0.0)
        prob.model.add_constraint(point_name + ".surface_perf.failure", upper=0.0)
        prob.model.add_constraint(point_name + ".surface_perf.thickness_intersects", upper=0.0)
        prob.model.add_constraint("surface.surface_joining.section_separation", upper=0, lower=0, scaler=1e5)

        # Add objective
        prob.model.add_objective(point_name + ".fuelburn", scaler=1e-5)

        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options["optimizer"] = "SLSQP"
        prob.driver.options["tol"] = 1e-10
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

        from openaerostruct.integration.aerostruct_groups import (
            AerostructPoint,
            MultiSecAerostructGeometry,
        )
        from openaerostruct.utils.constants import grav_constant
        from openaerostruct.geometry.utils import (
            build_section_dicts,
            unify_mesh,
            connect_multi_spline,
            build_multi_spline,
        )
        from openaerostruct.utils.testing import get_two_section_surface_AS

        surface, sec_twist_cp = get_two_section_surface_AS()

        # Create the OpenMDAO problem
        prob = om.Problem(reports=False)

        # Create an independent variable component that will supply the flow
        # conditions to the problem.
        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output("v", val=248.136, units="m/s")
        indep_var_comp.add_output("alpha", val=9.0, units="deg")
        indep_var_comp.add_output("Mach_number", val=0.84)
        indep_var_comp.add_output("re", val=1.0e6, units="1/m")
        indep_var_comp.add_output("rho", val=0.38, units="kg/m**3")
        indep_var_comp.add_output("CT", val=grav_constant * 17.0e-6, units="1/s")
        indep_var_comp.add_output("R", val=11.165e6, units="m")
        indep_var_comp.add_output("W0", val=0.4 * 3e5, units="kg")
        indep_var_comp.add_output("speed_of_sound", val=295.4, units="m/s")
        indep_var_comp.add_output("load_factor", val=1.0)
        indep_var_comp.add_output("empty_cg", val=np.zeros((3)), units="m")

        # Add this IndepVarComp to the problem model
        prob.model.add_subsystem("prob_vars", indep_var_comp, promotes=["*"])

        # Generate the sections and unified mesh here. It's needed to join the sections by construction.
        section_surfaces = build_section_dicts(surface)
        uni_mesh = unify_mesh(section_surfaces)
        surface["mesh"] = uni_mesh

        # Build a component with B-spline control points that joins the sections by construction
        chord_comp = build_multi_spline("twist_cp", len(section_surfaces), sec_twist_cp)
        prob.model.add_subsystem("twist_bspline", chord_comp)

        # Connect the B-spline component to the section B-splines
        connect_multi_spline(prob, section_surfaces, sec_twist_cp, "twist_cp", "twist_bspline", surface["name"])

        # Create and add a group that handles the geometry for the
        # aerodynamic lifting surface
        multi_aerostruct_geom_group = MultiSecAerostructGeometry(surface=surface)
        prob.model.add_subsystem(surface["name"], multi_aerostruct_geom_group)

        name = surface["name"]

        # Create the aero point group, which contains the actual aerodynamic
        # analyses
        aerostruct_group = AerostructPoint(surfaces=[surface])
        point_name = "AS_point_0"
        prob.model.add_subsystem(point_name, aerostruct_group)

        # Connect Flight conditions
        prob.model.connect("v", point_name + ".v")
        prob.model.connect("alpha", point_name + ".alpha")
        prob.model.connect("Mach_number", point_name + ".Mach_number")
        prob.model.connect("re", point_name + ".re")
        prob.model.connect("rho", point_name + ".rho")
        prob.model.connect("CT", point_name + ".CT")
        prob.model.connect("R", point_name + ".R")
        prob.model.connect("W0", point_name + ".W0")
        prob.model.connect("speed_of_sound", point_name + ".speed_of_sound")
        prob.model.connect("empty_cg", point_name + ".empty_cg")
        prob.model.connect("load_factor", point_name + ".load_factor")

        # Connect aerostruct groups
        com_name = point_name + "." + name + "_perf"
        prob.model.connect(
            name + ".local_stiff_transformed", point_name + ".coupled." + name + ".local_stiff_transformed"
        )
        prob.model.connect(name + ".nodes", point_name + ".coupled." + name + ".nodes")

        # Connect aerodyamic mesh to coupled group mesh
        prob.model.connect(name + ".mesh", point_name + ".coupled." + name + ".mesh")

        # Connect performance calculation variables
        prob.model.connect(name + ".radius", com_name + ".radius")
        prob.model.connect(name + ".thickness", com_name + ".thickness")
        prob.model.connect(name + ".nodes", com_name + ".nodes")
        prob.model.connect(name + ".cg_location", point_name + "." + "total_perf." + name + "_cg_location")
        prob.model.connect(name + ".structural_mass", point_name + "." + "total_perf." + name + "_structural_mass")
        prob.model.connect(name + ".t_over_c", com_name + ".t_over_c")

        # Add DVs
        prob.model.add_design_var("alpha", lower=-10.0, upper=10.0)
        prob.model.add_design_var("twist_bspline.twist_cp_spline", lower=-10.0, upper=15.0)
        prob.model.add_design_var("surface.thickness_cp", lower=0.01, upper=0.5, scaler=1e2)

        # Add constraints
        prob.model.add_constraint(point_name + ".L_equals_W", equals=0.0)
        prob.model.add_constraint(point_name + ".surface_perf.failure", upper=0.0)
        prob.model.add_constraint(point_name + ".surface_perf.thickness_intersects", upper=0.0)

        # Add objective
        prob.model.add_objective(point_name + ".fuelburn", scaler=1e-5)

        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options["optimizer"] = "SLSQP"
        prob.driver.options["tol"] = 1e-10
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

        from openaerostruct.integration.aerostruct_groups import (
            AerostructPoint,
            MultiSecAerostructGeometry,
        )
        from openaerostruct.utils.constants import grav_constant
        from openaerostruct.geometry.utils import (
            build_section_dicts,
            unify_mesh,
            connect_multi_spline,
            build_multi_spline,
        )
        from openaerostruct.utils.testing import get_two_section_surface_AS

        surface, sec_twist_cp = get_two_section_surface_AS(sym=False)

        # Create the OpenMDAO problem
        prob = om.Problem(reports=False)

        # Create an independent variable component that will supply the flow
        # conditions to the problem.
        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output("v", val=248.136, units="m/s")
        indep_var_comp.add_output("alpha", val=9.0, units="deg")
        indep_var_comp.add_output("Mach_number", val=0.84)
        indep_var_comp.add_output("re", val=1.0e6, units="1/m")
        indep_var_comp.add_output("rho", val=0.38, units="kg/m**3")
        indep_var_comp.add_output("CT", val=grav_constant * 17.0e-6, units="1/s")
        indep_var_comp.add_output("R", val=11.165e6, units="m")
        indep_var_comp.add_output("W0", val=0.4 * 3e5, units="kg")
        indep_var_comp.add_output("speed_of_sound", val=295.4, units="m/s")
        indep_var_comp.add_output("load_factor", val=1.0)
        indep_var_comp.add_output("empty_cg", val=np.zeros((3)), units="m")

        # Add this IndepVarComp to the problem model
        prob.model.add_subsystem("prob_vars", indep_var_comp, promotes=["*"])

        # Generate the sections and unified mesh here. It's needed to join the sections by construction.
        section_surfaces = build_section_dicts(surface)
        uni_mesh = unify_mesh(section_surfaces)
        surface["mesh"] = uni_mesh

        # Build a component with B-spline control points that joins the sections by construction
        chord_comp = build_multi_spline("twist_cp", len(section_surfaces), sec_twist_cp)
        prob.model.add_subsystem("twist_bspline", chord_comp)

        # Connect the B-spline component to the section B-splines
        connect_multi_spline(prob, section_surfaces, sec_twist_cp, "twist_cp", "twist_bspline", surface["name"])

        # Create and add a group that handles the geometry for the
        # aerodynamic lifting surface
        multi_aerostruct_geom_group = MultiSecAerostructGeometry(surface=surface)
        prob.model.add_subsystem(surface["name"], multi_aerostruct_geom_group)

        name = surface["name"]

        # Create the aero point group, which contains the actual aerodynamic
        # analyses
        aerostruct_group = AerostructPoint(surfaces=[surface])
        point_name = "AS_point_0"
        prob.model.add_subsystem(point_name, aerostruct_group)

        # Connect Flight conditions
        prob.model.connect("v", point_name + ".v")
        prob.model.connect("alpha", point_name + ".alpha")
        prob.model.connect("Mach_number", point_name + ".Mach_number")
        prob.model.connect("re", point_name + ".re")
        prob.model.connect("rho", point_name + ".rho")
        prob.model.connect("CT", point_name + ".CT")
        prob.model.connect("R", point_name + ".R")
        prob.model.connect("W0", point_name + ".W0")
        prob.model.connect("speed_of_sound", point_name + ".speed_of_sound")
        prob.model.connect("empty_cg", point_name + ".empty_cg")
        prob.model.connect("load_factor", point_name + ".load_factor")

        # Connect aerostruct groups
        com_name = point_name + "." + name + "_perf"
        prob.model.connect(
            name + ".local_stiff_transformed", point_name + ".coupled." + name + ".local_stiff_transformed"
        )
        prob.model.connect(name + ".nodes", point_name + ".coupled." + name + ".nodes")

        # Connect aerodyamic mesh to coupled group mesh
        prob.model.connect(name + ".mesh", point_name + ".coupled." + name + ".mesh")

        # Connect performance calculation variables
        prob.model.connect(name + ".radius", com_name + ".radius")
        prob.model.connect(name + ".thickness", com_name + ".thickness")
        prob.model.connect(name + ".nodes", com_name + ".nodes")
        prob.model.connect(name + ".cg_location", point_name + "." + "total_perf." + name + "_cg_location")
        prob.model.connect(name + ".structural_mass", point_name + "." + "total_perf." + name + "_structural_mass")
        prob.model.connect(name + ".t_over_c", com_name + ".t_over_c")

        # Add DVs
        prob.model.add_design_var("alpha", lower=-10.0, upper=10.0)
        prob.model.add_design_var("twist_bspline.twist_cp_spline", lower=-10.0, upper=15.0)
        prob.model.add_design_var("surface.thickness_cp", lower=0.01, upper=0.5, scaler=1e2)

        # Add constraints
        prob.model.add_constraint(point_name + ".L_equals_W", equals=0.0)
        prob.model.add_constraint(point_name + ".surface_perf.failure", upper=0.0)
        prob.model.add_constraint(point_name + ".surface_perf.thickness_intersects", upper=0.0)

        # Add objective
        prob.model.add_objective(point_name + ".fuelburn", scaler=1e-5)

        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options["optimizer"] = "SLSQP"
        prob.driver.options["tol"] = 1e-10
        prob.driver.options["disp"] = True
        prob.driver.options["maxiter"] = 1000
        # prob.driver.options["debug_print"] = ["nl_cons", "objs", "desvars"]

        # Set up and run the optimization problem
        prob.setup()
        prob.run_driver()

        assert_near_equal(prob["AS_point_0.fuelburn"][0], 85719.42698562, 1e-3)


if __name__ == "__main__":
    unittest.main()
