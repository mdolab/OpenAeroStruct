from openmdao.utils.assert_utils import assert_near_equal
import unittest


class Test(unittest.TestCase):
    def test(self):
        import numpy as np

        import openmdao.api as om

        from openaerostruct.integration.aerostruct_groups import (
            AerostructGeometry,
            AerostructPoint,
            MultiSecAerostructGeometry,
        )
        from openaerostruct.aerodynamics.aero_groups import AeroPoint
        from openaerostruct.utils.constants import grav_constant
        from openaerostruct.geometry.utils import (
            build_section_dicts,
            unify_mesh,
            connect_multi_spline,
            build_multi_spline,
        )
        from openaerostruct.utils.testing import get_two_section_surface_AS, get_single_section_surface_AS

        # Setup and run the single-section aerostrucutral optimization
        """Create a dictionary with info and options about the aerodynamic
        single section lifting surface"""
        surfaceSingle = get_single_section_surface_AS()

        # Create the OpenMDAO problem
        prob = om.Problem()

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

        # Create and add a group that handles the geometry for the
        # aerodynamic single lifting surface
        aerostruct_geom_group = AerostructGeometry(surface=surfaceSingle)
        prob.model.add_subsystem(surfaceSingle["name"], aerostruct_geom_group)

        # Create the aero point group, which contains the actual aerodynamic
        # analyses for the single section surface
        aerostruct_group = AerostructPoint(surfaces=[surfaceSingle])
        point_name0 = "AS_point_0"
        prob.model.add_subsystem(point_name0, aerostruct_group)
        name = surfaceSingle["name"]

        # Connect Flight conditions
        prob.model.connect("v", point_name0 + ".v")
        prob.model.connect("alpha", point_name0 + ".alpha")
        prob.model.connect("Mach_number", point_name0 + ".Mach_number")
        prob.model.connect("re", point_name0 + ".re")
        prob.model.connect("rho", point_name0 + ".rho")
        prob.model.connect("CT", point_name0 + ".CT")
        prob.model.connect("R", point_name0 + ".R")
        prob.model.connect("W0", point_name0 + ".W0")
        prob.model.connect("speed_of_sound", point_name0 + ".speed_of_sound")
        prob.model.connect("empty_cg", point_name0 + ".empty_cg")
        prob.model.connect("load_factor", point_name0 + ".load_factor")

        # Connect geometry group to analysis point
        com_name = point_name0 + "." + name + "_perf"
        prob.model.connect(
            name + ".local_stiff_transformed", point_name0 + ".coupled." + name + ".local_stiff_transformed"
        )
        prob.model.connect(name + ".nodes", point_name0 + ".coupled." + name + ".nodes")

        # Connect aerodyamic mesh to coupled group mesh
        prob.model.connect(name + ".mesh", point_name0 + ".coupled." + name + ".mesh")

        # Connect performance calculation variables
        prob.model.connect(name + ".radius", com_name + ".radius")
        prob.model.connect(name + ".thickness", com_name + ".thickness")
        prob.model.connect(name + ".nodes", com_name + ".nodes")
        prob.model.connect(name + ".cg_location", point_name0 + "." + "total_perf." + name + "_cg_location")
        prob.model.connect(name + ".structural_mass", point_name0 + "." + "total_perf." + name + "_structural_mass")
        prob.model.connect(name + ".t_over_c", com_name + ".t_over_c")

        # Add DVs
        prob.model.add_design_var("alpha", lower=-20.0, upper=20.0)

        # Add L=W constraint
        prob.model.add_constraint("AS_point_0.L_equals_W", equals=0.0)

        # Add objective
        prob.model.add_objective(point_name0 + ".fuelburn", scaler=1e-5)

        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options["optimizer"] = "SLSQP"
        prob.driver.options["tol"] = 1e-6
        prob.driver.options["disp"] = True
        prob.driver.options["maxiter"] = 1000
        # prob.driver.options["debug_print"] = ["nl_cons", "objs", "desvars"]

        # Set up and run the optimization problem
        prob.setup()
        prob.run_driver()

        ## Setup and run the two-section aerostructural

        """Create a dictionary with info and options about the multi-section aerodynamic
        lifting surface"""
        surfaceMulti, sec_twist_cp = get_two_section_surface_AS()
        surfaceMulti["name"] = "surfaceMulti"

        # Create the other OpenMDAO problem
        prob2 = om.Problem(reports=False)

        # Create an independent variable component that will supply the flow
        # conditions to the problem.
        indep_var_comp2 = om.IndepVarComp()
        indep_var_comp2.add_output("v", val=248.136, units="m/s")
        indep_var_comp2.add_output("alpha", val=9.0, units="deg")
        indep_var_comp2.add_output("Mach_number", val=0.84)
        indep_var_comp2.add_output("re", val=1.0e6, units="1/m")
        indep_var_comp2.add_output("rho", val=0.38, units="kg/m**3")
        indep_var_comp2.add_output("CT", val=grav_constant * 17.0e-6, units="1/s")
        indep_var_comp2.add_output("R", val=11.165e6, units="m")
        indep_var_comp2.add_output("W0", val=0.4 * 3e5, units="kg")
        indep_var_comp2.add_output("speed_of_sound", val=295.4, units="m/s")
        indep_var_comp2.add_output("load_factor", val=1.0)
        indep_var_comp2.add_output("empty_cg", val=np.zeros((3)), units="m")

        # Add this IndepVarComp to the problem model
        prob2.model.add_subsystem("prob_vars", indep_var_comp2, promotes=["*"])

        # Generate the sections and unified mesh here in addition to adding the components.
        # This has to also be done here since AeroPoint has to know the unified mesh size.
        section_surfaces = build_section_dicts(surfaceMulti)
        uniMesh = unify_mesh(section_surfaces)
        surfaceMulti["mesh"] = uniMesh

        # Build a component with B-spline control points that joins the sections by construction
        twist_comp = build_multi_spline("twist_cp", len(section_surfaces), sec_twist_cp)
        prob2.model.add_subsystem("twist_bspline", twist_comp)

        # Connect the B-spline component to the section B-splines
        connect_multi_spline(prob2, section_surfaces, sec_twist_cp, "twist_cp", "twist_bspline", surfaceMulti["name"])

        # Create and add a group that handles the geometry for the
        # multi-section aerodynamic lifting surface
        multi_geom_group = MultiSecAerostructGeometry(surface=surfaceMulti)
        prob2.model.add_subsystem(surfaceMulti["name"], multi_geom_group)

        name = surfaceMulti["name"]

        # Create the aero point group, which contains the actual aerodynamic
        # analyses
        aerostruct_group = AerostructPoint(surfaces=[surfaceMulti])
        point_name1 = "AS_point_1"
        prob2.model.add_subsystem(point_name1, aerostruct_group)

        # Connect Flight conditions
        prob2.model.connect("v", point_name1 + ".v")
        prob2.model.connect("alpha", point_name1 + ".alpha")
        prob2.model.connect("Mach_number", point_name1 + ".Mach_number")
        prob2.model.connect("re", point_name1 + ".re")
        prob2.model.connect("rho", point_name1 + ".rho")
        prob2.model.connect("CT", point_name1 + ".CT")
        prob2.model.connect("R", point_name1 + ".R")
        prob2.model.connect("W0", point_name1 + ".W0")
        prob2.model.connect("speed_of_sound", point_name1 + ".speed_of_sound")
        prob2.model.connect("empty_cg", point_name1 + ".empty_cg")
        prob2.model.connect("load_factor", point_name1 + ".load_factor")

        com_name = point_name1 + "." + name + "_perf"
        prob2.model.connect(
            name + ".local_stiff_transformed", point_name1 + ".coupled." + name + ".local_stiff_transformed"
        )
        prob2.model.connect(name + ".nodes", point_name1 + ".coupled." + name + ".nodes")

        # Connect aerodyamic mesh to coupled group mesh
        prob2.model.connect(name + ".mesh", point_name1 + ".coupled." + name + ".mesh")

        # Connect performance calculation variables
        prob2.model.connect(name + ".radius", com_name + ".radius")
        prob2.model.connect(name + ".thickness", com_name + ".thickness")
        prob2.model.connect(name + ".nodes", com_name + ".nodes")
        prob2.model.connect(name + ".cg_location", point_name1 + "." + "total_perf." + name + "_cg_location")
        prob2.model.connect(name + ".structural_mass", point_name1 + "." + "total_perf." + name + "_structural_mass")
        prob2.model.connect(name + ".t_over_c", com_name + ".t_over_c")

        # Add DVs
        prob2.model.add_design_var("alpha", lower=-20.0, upper=20.0)

        # Add L=W constraint
        prob2.model.add_constraint("AS_point_1.L_equals_W", equals=0.0)

        # Add objective
        prob2.model.add_objective(point_name1 + ".fuelburn", scaler=1e-5)

        prob2.driver = om.ScipyOptimizeDriver()
        prob2.driver.options["optimizer"] = "SLSQP"
        prob2.driver.options["tol"] = 1e-6
        prob2.driver.options["disp"] = True
        prob2.driver.options["maxiter"] = 1000
        # prob.driver.options["debug_print"] = ["nl_cons", "objs", "desvars"]

        # Set up and run the optimization problem
        prob2.setup()
        prob2.run_driver()

        fuelburndiff = np.abs(prob[point_name0 + ".fuelburn"][0] - prob2[point_name1 + ".fuelburn"][0])

        assert_near_equal(fuelburndiff, 0.0, 1e-6)


if __name__ == "__main__":
    unittest.main()
