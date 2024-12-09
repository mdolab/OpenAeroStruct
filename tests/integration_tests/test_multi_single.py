from openmdao.utils.assert_utils import assert_near_equal
import unittest


class Test(unittest.TestCase):
    def test(self):
        import numpy as np

        import openmdao.api as om

        from openaerostruct.geometry.geometry_group import Geometry, MultiSecGeometry
        from openaerostruct.aerodynamics.aero_groups import AeroPoint
        from openaerostruct.geometry.geometry_group import build_sections
        from openaerostruct.geometry.geometry_unification import unify_mesh
        from openaerostruct.utils.testing import get_two_section_surface, get_single_section_surface

        """Create a dictionary with info and options about the aerodynamic
        single section lifting surface"""
        surfaceSingle = get_single_section_surface()

        """Create a dictionary with info and options about the multi-section aerodynamic
        lifting surface"""
        surfaceMulti, sec_chord_cp = get_two_section_surface()
        surfaceMulti["name"] = "surfaceMulti"

        # Create the OpenMDAO problem
        prob = om.Problem()

        # Create an independent variable component that will supply the flow
        # conditions to the problem.
        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output("v", val=1.0, units="m/s")
        indep_var_comp.add_output("alpha", val=10.0 * np.ones(2), units="deg")
        indep_var_comp.add_output("Mach_number", val=0.3)
        indep_var_comp.add_output("re", val=1.0e5, units="1/m")
        indep_var_comp.add_output("rho", val=0.38, units="kg/m**3")
        indep_var_comp.add_output("cg", val=np.zeros((3)), units="m")

        # Add this IndepVarComp to the problem model
        prob.model.add_subsystem("prob_vars", indep_var_comp, promotes=["*"])

        # Create and add a group that handles the geometry for the
        # aerodynamic single lifting surface
        geom_group = Geometry(surface=surfaceSingle)
        prob.model.add_subsystem(surfaceSingle["name"], geom_group)

        # Create the aero point group, which contains the actual aerodynamic
        # analyses for the single section surface
        aero_group = AeroPoint(surfaces=[surfaceSingle])
        point_name0 = "aero_point_0"
        prob.model.add_subsystem(point_name0, aero_group)
        name = surfaceSingle["name"]

        # Connect Flight conditions
        prob.model.connect("v", point_name0 + ".v")
        prob.model.connect("alpha", point_name0 + ".alpha", src_indices=[0])
        prob.model.connect("Mach_number", point_name0 + ".Mach_number")
        prob.model.connect("re", point_name0 + ".re")
        prob.model.connect("rho", point_name0 + ".rho")
        prob.model.connect("cg", point_name0 + ".cg")

        # Connect the mesh from the geometry component to the analysis point
        prob.model.connect(name + ".mesh", point_name0 + "." + name + ".def_mesh")
        # Perform the connections with the modified names within the
        # 'aero_states' group.
        prob.model.connect(name + ".mesh", point_name0 + ".aero_states." + name + "_def_mesh")

        # Create and add a group that handles the geometry for the
        # multi-section aerodynamic lifting surface
        multi_geom_group = MultiSecGeometry(
            surface=surfaceMulti, joining_comp=True, dim_constr=[np.array([1, 0, 0]), np.array([1, 0, 0])]
        )
        prob.model.add_subsystem(surfaceMulti["name"], multi_geom_group)

        # Generate the sections and unified mesh here in addition to adding the components.
        # This has to also be done here since AeroPoint has to know the unified mesh size.
        section_surfaces = build_sections(surfaceMulti)
        uni_mesh = unify_mesh(section_surfaces)
        surfaceMulti["mesh"] = uni_mesh

        # Create the aero point group, which contains the actual aerodynamic
        # analyses
        aero_group = AeroPoint(surfaces=[surfaceMulti])
        point_name1 = "aero_point_1"
        prob.model.add_subsystem(point_name1, aero_group)

        # Connect Flight conditions
        prob.model.connect("v", point_name1 + ".v")
        prob.model.connect("alpha", point_name1 + ".alpha", src_indices=[1])
        prob.model.connect("Mach_number", point_name1 + ".Mach_number")
        prob.model.connect("re", point_name1 + ".re")
        prob.model.connect("rho", point_name1 + ".rho")
        prob.model.connect("cg", point_name1 + ".cg")

        # Get name of surface and construct unified mesh name
        name = surfaceMulti["name"]
        unification_name = "{}_unification".format(surfaceMulti["name"])

        # Connect the mesh from the mesh unification component to the analysis point
        prob.model.connect(
            name + "." + unification_name + "." + name + "_uni_mesh", point_name1 + "." + name + ".def_mesh"
        )

        # Perform the connections with the modified names within the
        # 'aero_states' group.
        prob.model.connect(
            name + "." + unification_name + "." + name + "_uni_mesh",
            point_name1 + ".aero_states." + name + "_def_mesh",
        )

        # Add DVs
        prob.model.add_design_var("alpha", lower=0.0, upper=10.0, units="deg")

        # Add CL constraint
        prob.model.add_constraint(point_name0 + ".CL", equals=0.3)
        prob.model.add_constraint(point_name1 + ".CL", equals=0.3)

        # Add objective
        prob.model.add_objective(point_name0 + ".CD", scaler=1e4)

        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options["optimizer"] = "SLSQP"
        prob.driver.options["tol"] = 1e-6
        prob.driver.options["disp"] = True
        prob.driver.options["maxiter"] = 1000
        # prob.driver.options["debug_print"] = ["nl_cons", "objs", "desvars"]

        # Set up and run the optimization problem
        prob.setup()

        prob.run_driver()

        CDDiff = np.abs(prob["aero_point_1.surfaceMulti_perf.CD"][0] - prob["aero_point_0.surface_perf.CD"][0])
        CLDiff = np.abs(prob["aero_point_1.surfaceMulti_perf.CL"][0] - prob["aero_point_0.surface_perf.CL"][0])

        assert_near_equal(CDDiff, 0.0, 1e-6)
        assert_near_equal(CLDiff, 0.0, 1e-6)


if __name__ == "__main__":
    unittest.main()
