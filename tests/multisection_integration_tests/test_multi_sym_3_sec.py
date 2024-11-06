from openmdao.utils.assert_utils import assert_near_equal
import unittest


class Test(unittest.TestCase):
    def test(self):
        import numpy as np

        import openmdao.api as om

        from openaerostruct.geometry.geometry_group import MultiSecGeometry
        from openaerostruct.aerodynamics.aero_groups import AeroPoint
        from openaerostruct.geometry.geometry_group import build_sections
        from openaerostruct.geometry.geometry_unification import unify_mesh

        # Create a dictionary with info and options about the multi-section aerodynamic
        # lifting surface
        surface = {
            # Wing definition
            # Basic surface parameters
            "name": "surface",
            "isMultiSection": True,
            "num_sections": 3,  # The number of sections in the multi-section surface
            "sec_name": ["sec0", "sec1", "sec2"],  # names of the individual sections
            "symmetry": True,  # if true, model one half of wing. reflected across the midspan of the root section
            "S_ref_type": "wetted",  # how we compute the wing area,
            # can be 'wetted' or 'projected'
            # Geometry Parameters
            "taper": [1.0, 1.0, 1.0],  # Wing taper for each section
            "span": [1.0, 1.0, 1.0],  # Wing span for each section
            "sweep": [0.0, 0, 0.0],  # Wing sweep for each section
            "twist_cp": [np.zeros(2), np.zeros(2), np.zeros(2)],
            "chord_cp": [np.array([1, 1]), np.array([0.2, 1.0]), np.array([0.2, 1.0])],
            # "sec_chord_cp": [np.ones(1),2*np.ones(1),3*np.ones(1)], #Chord B-spline control points for each section
            "root_chord": 1.0,  # Wing root chord for each section
            # Mesh Parameters
            "meshes": "gen-meshes",  # Supply a mesh for each section or "gen-meshes" for automatic mesh generation
            "nx": 2,  # Number of chordwise points. Same for all sections
            "ny": [3, 3, 3],  # Number of spanwise points for each section
            # Aerodynamic Parameters
            "CL0": 0.0,  # CL of the surface at alpha=0
            "CD0": 0.015,  # CD of the surface at alpha=0
            # Airfoil properties for viscous drag calculation
            "k_lam": 0.05,  # percentage of chord with laminar
            # flow, used for viscous drag
            "t_over_c_cp": [
                np.array([0.15]),
                np.array([0.15]),
                np.array([0.15]),
            ],  # thickness over chord ratio (NACA0015)
            "c_max_t": 0.303,  # chordwise location of maximum (NACA0015)
            # thickness
            "with_viscous": False,  # if true, compute viscous drag
            "with_wave": False,  # if true, compute wave drag
            "groundplane": False,
        }

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
            surface=surface,
            joining_comp=True,
            dim_constr=[np.array([1, 0, 0]), np.array([1, 0, 0]), np.array([1, 0, 0])],
        )
        prob.model.add_subsystem(surface["name"], multi_geom_group)

        # Generate the sections and unified mesh here in addition to adding the components.
        # This has to ALSO be done here since AeroPoint has to know the unified mesh size.
        section_surfaces = build_sections(surface)
        uniMesh = unify_mesh(section_surfaces)
        surface["mesh"] = uniMesh

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
        prob.model.add_design_var("surface.sec2.chord_cp", lower=0.1, upper=10.0, units=None)
        prob.model.add_design_var("alpha", lower=0.0, upper=10.0, units="deg")

        # Add joined mesh constraint
        prob.model.add_constraint("surface.surface_joining.section_separation", upper=0, lower=0)
        # prob.model.add_constraint('surface.surface_joining.section_separation',equals=0.0)

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
        prob.driver.options["debug_print"] = ["nl_cons", "objs", "desvars"]

        # Set up and run the optimization problem
        prob.setup()

        # prob.run_model()
        prob.run_driver()
        # om.n2(prob)

        assert_near_equal(prob["aero_point_0.surface_perf.CD"][0], 0.02087887, 1e-6)
        assert_near_equal(prob["aero_point_0.surface_perf.CL"][0], 0.29999954, 1e-6)
        assert_near_equal(prob["aero_point_0.CM"][1], -0.10869776, 1e-6)


if __name__ == "__main__":
    unittest.main()
