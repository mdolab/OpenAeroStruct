"""
Parallel version of the multipoint aerostructural optimization example (Q400).
"""

import unittest
from openmdao.utils.assert_utils import assert_near_equal

try:
    from petsc4py import PETSc  # noqa: F401

    petsc_installed = True
except ModuleNotFoundError:
    petsc_installed = False

try:
    from mpi4py import MPI  # noqa: F401

    mpi_installed = True
except ModuleNotFoundError:
    mpi_installed = False


@unittest.skipUnless(petsc_installed and mpi_installed, "MPI and PETSc are required.")
class Test(unittest.TestCase):
    N_PROCS = 2

    def test_multipoint_MPI(self):
        import numpy as np
        import time

        from openaerostruct.integration.aerostruct_groups import AerostructGeometry, AerostructPoint
        from openaerostruct.structures.wingbox_fuel_vol_delta import WingboxFuelVolDelta
        import openmdao.api as om

        from mpi4py import MPI  # noqa: F811

        # Provide coordinates for a portion of an airfoil for the wingbox cross-section as an nparray with dtype=complex (to work with the complex-step approximation for derivatives).
        # These should be for an airfoil with the chord scaled to 1.
        # We use the 10% to 60% portion of the NASA SC2-0612 airfoil for this case
        # We use the coordinates available from airfoiltools.com. Using such a large number of coordinates is not necessary.
        # The first and last x-coordinates of the upper and lower surfaces must be the same

        # fmt: off
        upper_x = np.array([0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6], dtype="complex128")
        lower_x = np.array([0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6], dtype="complex128")
        upper_y = np.array([ 0.0447,  0.046,  0.0472,  0.0484,  0.0495,  0.0505,  0.0514,  0.0523,  0.0531,  0.0538, 0.0545,  0.0551,  0.0557, 0.0563,  0.0568, 0.0573,  0.0577,  0.0581,  0.0585,  0.0588,  0.0591,  0.0593,  0.0595,  0.0597,  0.0599,  0.06,    0.0601,  0.0602,  0.0602,  0.0602,  0.0602,  0.0602,  0.0601,  0.06,    0.0599,  0.0598,  0.0596,  0.0594,  0.0592,  0.0589,  0.0586,  0.0583,  0.058,   0.0576,  0.0572,  0.0568,  0.0563,  0.0558,  0.0553,  0.0547,  0.0541], dtype="complex128")  # noqa: E201, E241
        lower_y = np.array([-0.0447, -0.046, -0.0473, -0.0485, -0.0496, -0.0506, -0.0515, -0.0524, -0.0532, -0.054, -0.0547, -0.0554, -0.056, -0.0565, -0.057, -0.0575, -0.0579, -0.0583, -0.0586, -0.0589, -0.0592, -0.0594, -0.0595, -0.0596, -0.0597, -0.0598, -0.0598, -0.0598, -0.0598, -0.0597, -0.0596, -0.0594, -0.0592, -0.0589, -0.0586, -0.0582, -0.0578, -0.0573, -0.0567, -0.0561, -0.0554, -0.0546, -0.0538, -0.0529, -0.0519, -0.0509, -0.0497, -0.0485, -0.0472, -0.0458, -0.0444], dtype="complex128")
        # fmt: on

        # Here we create a custom mesh for the wing
        # It is evenly spaced with nx chordwise nodal points and ny spanwise nodal points for the half-span

        span = 28.42  # wing span in m
        root_chord = 3.34  # root chord in m

        # Let's use a finer mesh
        nx = 3  # number of chordwise nodal points (should be odd)
        ny = 11  # number of spanwise nodal points for the half-span

        # Initialize the 3-D mesh object. Chordwise, spanwise, then the 3D coordinates.
        mesh = np.zeros((nx, ny, 3))

        # Start away from the symmetry plane and approach the plane as the array indices increase.
        # The form of this 3-D array can be very confusing initially.
        # For each node we are providing the x, y, and z coordinates.
        # x is chordwise, y is spanwise, and z is up.
        # For example (for a mesh with 5 chordwise nodes and 15 spanwise nodes for the half wing), the node for the leading edge at the tip would be specified as mesh[0, 0, :] = np.array([1.1356, -14.21, 0.])
        # and the node at the trailing edge at the root would be mesh[4, 14, :] = np.array([3.34, 0., 0.]).
        # We only provide the left half of the wing because we use symmetry.
        # Print the following mesh and elements of the mesh to better understand the form.

        mesh[:, :, 1] = np.linspace(-span / 2, 0, ny)
        mesh[0, :, 0] = 0.34 * root_chord * np.linspace(1.0, 0.0, ny)
        mesh[2, :, 0] = root_chord * (np.linspace(0.4, 1.0, ny) + 0.34 * np.linspace(1.0, 0.0, ny))
        mesh[1, :, 0] = (mesh[2, :, 0] + mesh[0, :, 0]) / 2

        # print(mesh)

        surf_dict = {
            # Wing definition
            "name": "wing",  # name of the surface
            "symmetry": True,  # if true, model one half of wing
            "S_ref_type": "wetted",  # how we compute the wing area,
            # can be 'wetted' or 'projected'
            "mesh": mesh,
            "twist_cp": np.array([6.0, 7.0, 7.0, 7.0]),
            "fem_model_type": "wingbox",
            "data_x_upper": upper_x,
            "data_x_lower": lower_x,
            "data_y_upper": upper_y,
            "data_y_lower": lower_y,
            "spar_thickness_cp": np.array([0.004, 0.004, 0.004, 0.004]),  # [m]
            "skin_thickness_cp": np.array([0.003, 0.006, 0.010, 0.012]),  # [m]
            "original_wingbox_airfoil_t_over_c": 0.12,
            # Aerodynamic deltas.
            # These CL0 and CD0 values are added to the CL and CD
            # obtained from aerodynamic analysis of the surface to get
            # the total CL and CD.
            # These CL0 and CD0 values do not vary wrt alpha.
            # They can be used to account for things that are not included, such as contributions from the fuselage, nacelles, tail surfaces, etc.
            "CL0": 0.0,
            "CD0": 0.0142,
            "with_viscous": True,  # if true, compute viscous drag
            "with_wave": True,  # if true, compute wave drag
            # Airfoil properties for viscous drag calculation
            "k_lam": 0.05,  # percentage of chord with laminar
            # flow, used for viscous drag
            "c_max_t": 0.38,  # chordwise location of maximum thickness
            "t_over_c_cp": np.array([0.1, 0.1, 0.15, 0.15]),
            # Structural values are based on aluminum 7075
            "E": 73.1e9,  # [Pa] Young's modulus
            "G": (73.1e9 / 2 / 1.33),  # [Pa] shear modulus (calculated using E and the Poisson's ratio here)
            "yield": (420.0e6 / 1.5),  # [Pa] allowable yield stress
            "mrho": 2.78e3,  # [kg/m^3] material density
            "strength_factor_for_upper_skin": 1.0,  # the yield stress is multiplied by this factor for the upper skin
            "wing_weight_ratio": 1.25,
            "exact_failure_constraint": False,  # if false, use KS function
            "struct_weight_relief": True,
            "distributed_fuel_weight": True,
            "fuel_density": 803.0,  # [kg/m^3] fuel density (only needed if the fuel-in-wing volume constraint is used)
            "Wf_reserve": 500.0,  # [kg] reserve fuel mass
        }

        surfaces = [surf_dict]

        # Create the problem and assign the model group
        prob = om.Problem()

        # Add problem information as an independent variables component
        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output("v", val=np.array([0.5 * 310.95, 0.3 * 340.294]), units="m/s")
        indep_var_comp.add_output("alpha", val=0.0, units="deg")
        indep_var_comp.add_output("alpha_maneuver", val=0.0, units="deg")
        indep_var_comp.add_output("Mach_number", val=np.array([0.5, 0.3]))
        indep_var_comp.add_output(
            "re",
            val=np.array([0.569 * 310.95 * 0.5 * 1.0 / (1.56 * 1e-5), 1.225 * 340.294 * 0.3 * 1.0 / (1.81206 * 1e-5)]),
            units="1/m",
        )
        indep_var_comp.add_output("rho", val=np.array([0.569, 1.225]), units="kg/m**3")
        indep_var_comp.add_output("CT", val=0.43 / 3600, units="1/s")
        indep_var_comp.add_output("R", val=2e6, units="m")
        indep_var_comp.add_output("W0", val=25400 + surf_dict["Wf_reserve"], units="kg")
        indep_var_comp.add_output("speed_of_sound", val=np.array([310.95, 340.294]), units="m/s")
        indep_var_comp.add_output("load_factor", val=np.array([1.0, 2.5]))
        indep_var_comp.add_output("empty_cg", val=np.zeros((3)), units="m")
        indep_var_comp.add_output("fuel_mass", val=3000.0, units="kg")

        prob.model.add_subsystem("prob_vars", indep_var_comp, promotes=["*"])

        # Loop over each surface in the surfaces list
        for surface in surfaces:
            # Get the surface name and create a group to contain components
            # only for this surface
            name = surface["name"]

            aerostruct_group = AerostructGeometry(surface=surface)

            # Add group to the problem with the name of the surface.
            prob.model.add_subsystem(name, aerostruct_group)

        # ------ modification for parallel multipoint analysis -----
        # [rst Setup ParallelGroup (beg)]
        # Add a ParallelGroup to parallelize AS_point_0 and AS_point_1.
        # We now set max_procs=2, but you may want to increase it if you have more than 2 AS points.
        parallel = prob.model.add_subsystem("parallel", om.ParallelGroup(), promotes=["*"], min_procs=1, max_procs=2)

        # Loop through and add a certain number of aerostruct points.
        for i in range(2):
            point_name = "AS_point_{}".format(i)

            # Create the aero point group and add it to the model
            AS_point = AerostructPoint(surfaces=surfaces, internally_connect_fuelburn=False)

            # Now, add each point under the parallel group (instead of directly under prob.model)
            parallel.add_subsystem(point_name, AS_point)
            # [rst Setup ParallelGroup (end)]

            # Connect flow properties to the analysis point
            prob.model.connect("v", point_name + ".v", src_indices=[i])
            prob.model.connect("Mach_number", point_name + ".Mach_number", src_indices=[i])
            prob.model.connect("re", point_name + ".re", src_indices=[i])
            prob.model.connect("rho", point_name + ".rho", src_indices=[i])
            prob.model.connect("CT", point_name + ".CT")
            prob.model.connect("R", point_name + ".R")
            prob.model.connect("W0", point_name + ".W0")
            prob.model.connect("speed_of_sound", point_name + ".speed_of_sound", src_indices=[i])
            prob.model.connect("empty_cg", point_name + ".empty_cg")
            prob.model.connect("load_factor", point_name + ".load_factor", src_indices=[i])
            prob.model.connect("fuel_mass", point_name + ".total_perf.L_equals_W.fuelburn")
            prob.model.connect("fuel_mass", point_name + ".total_perf.CG.fuelburn")

            for surface in surfaces:
                name = surface["name"]

                if surf_dict["distributed_fuel_weight"]:
                    prob.model.connect("load_factor", point_name + ".coupled.load_factor", src_indices=[i])

                com_name = point_name + "." + name + "_perf."
                prob.model.connect(
                    name + ".local_stiff_transformed", point_name + ".coupled." + name + ".local_stiff_transformed"
                )
                prob.model.connect(name + ".nodes", point_name + ".coupled." + name + ".nodes")

                # Connect aerodyamic mesh to coupled group mesh
                prob.model.connect(name + ".mesh", point_name + ".coupled." + name + ".mesh")
                if surf_dict["struct_weight_relief"]:
                    prob.model.connect(name + ".element_mass", point_name + ".coupled." + name + ".element_mass")

                # Connect performance calculation variables
                prob.model.connect(name + ".nodes", com_name + "nodes")
                prob.model.connect(name + ".cg_location", point_name + "." + "total_perf." + name + "_cg_location")
                prob.model.connect(
                    name + ".structural_mass", point_name + "." + "total_perf." + name + "_structural_mass"
                )

                # Connect wingbox properties to von Mises stress calcs
                prob.model.connect(name + ".Qz", com_name + "Qz")
                prob.model.connect(name + ".J", com_name + "J")
                prob.model.connect(name + ".A_enc", com_name + "A_enc")
                prob.model.connect(name + ".htop", com_name + "htop")
                prob.model.connect(name + ".hbottom", com_name + "hbottom")
                prob.model.connect(name + ".hfront", com_name + "hfront")
                prob.model.connect(name + ".hrear", com_name + "hrear")

                prob.model.connect(name + ".spar_thickness", com_name + "spar_thickness")
                prob.model.connect(name + ".t_over_c", com_name + "t_over_c")

        prob.model.connect("alpha", "AS_point_0" + ".alpha")
        prob.model.connect("alpha_maneuver", "AS_point_1" + ".alpha")

        # Here we add the fuel volume constraint componenet to the model
        prob.model.add_subsystem("fuel_vol_delta", WingboxFuelVolDelta(surface=surface))
        prob.model.connect("wing.struct_setup.fuel_vols", "fuel_vol_delta.fuel_vols")
        prob.model.connect("AS_point_0.fuelburn", "fuel_vol_delta.fuelburn")

        if surf_dict["distributed_fuel_weight"]:
            prob.model.connect("wing.struct_setup.fuel_vols", "AS_point_0.coupled.wing.struct_states.fuel_vols")
            prob.model.connect("fuel_mass", "AS_point_0.coupled.wing.struct_states.fuel_mass")

            prob.model.connect("wing.struct_setup.fuel_vols", "AS_point_1.coupled.wing.struct_states.fuel_vols")
            prob.model.connect("fuel_mass", "AS_point_1.coupled.wing.struct_states.fuel_mass")

        comp = om.ExecComp("fuel_diff = (fuel_mass - fuelburn) / fuelburn", units="kg")
        prob.model.add_subsystem("fuel_diff", comp, promotes_inputs=["fuel_mass"], promotes_outputs=["fuel_diff"])
        prob.model.connect("AS_point_0.fuelburn", "fuel_diff.fuelburn")

        # Use these settings if you do not have pyOptSparse or SNOPT
        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options["optimizer"] = "SLSQP"
        prob.driver.options["tol"] = 1e-4

        # # The following are the optimizer settings used for the EngOpt conference paper
        # # Uncomment them if you can use SNOPT
        # prob.driver = om.pyOptSparseDriver()
        # prob.driver.options['optimizer'] = "SNOPT"
        # prob.driver.opt_settings['Major optimality tolerance'] = 5e-6
        # prob.driver.opt_settings['Major feasibility tolerance'] = 1e-8
        # prob.driver.opt_settings['Major iterations limit'] = 200

        recorder = om.SqliteRecorder("aerostruct.db")
        prob.driver.add_recorder(recorder)

        # We could also just use prob.driver.recording_options['includes']=['*'] here, but for large meshes the database file becomes extremely large. So we just select the variables we need.
        prob.driver.recording_options["includes"] = [
            "alpha",
            "rho",
            "v",
            "cg",
            "AS_point_1.cg",
            "AS_point_0.cg",
            "AS_point_0.coupled.wing_loads.loads",
            "AS_point_1.coupled.wing_loads.loads",
            "AS_point_0.coupled.wing.normals",
            "AS_point_1.coupled.wing.normals",
            "AS_point_0.coupled.wing.widths",
            "AS_point_1.coupled.wing.widths",
            "AS_point_0.coupled.aero_states.wing_sec_forces",
            "AS_point_1.coupled.aero_states.wing_sec_forces",
            "AS_point_0.wing_perf.CL1",
            "AS_point_1.wing_perf.CL1",
            "AS_point_0.coupled.wing.S_ref",
            "AS_point_1.coupled.wing.S_ref",
            "wing.geometry.twist",
            "wing.mesh",
            "wing.skin_thickness",
            "wing.spar_thickness",
            "wing.t_over_c",
            "wing.structural_mass",
            "AS_point_0.wing_perf.vonmises",
            "AS_point_1.wing_perf.vonmises",
            "AS_point_0.coupled.wing.def_mesh",
            "AS_point_1.coupled.wing.def_mesh",
        ]

        prob.driver.recording_options["record_objectives"] = True
        prob.driver.recording_options["record_constraints"] = True
        prob.driver.recording_options["record_desvars"] = True
        prob.driver.recording_options["record_inputs"] = True

        # --- define optimization problem ---
        # Design variables
        prob.model.add_design_var("wing.twist_cp", lower=-15.0, upper=15.0, scaler=0.1)
        prob.model.add_design_var("wing.spar_thickness_cp", lower=0.003, upper=0.1, scaler=1e2)
        prob.model.add_design_var("wing.skin_thickness_cp", lower=0.003, upper=0.1, scaler=1e2)
        prob.model.add_design_var("wing.geometry.t_over_c_cp", lower=0.07, upper=0.2, scaler=10.0)
        prob.model.add_design_var("fuel_mass", lower=0.0, upper=2e5, scaler=1e-5)
        prob.model.add_design_var("alpha_maneuver", lower=-15.0, upper=15)

        # ------ modification for parallel multipoint analysis -----
        # We have 6 functions (objective + constraints), which would require 6 linear solves for reverse-mode derivative computations in series.
        # Among these functions (outputs), 4 depend only on AS_point_0, and 2 depend only on AS_point_1.
        # Therefore we can form 2 pairs and perform linear solves in parallel. This is done by assigning the same parallel_deriv_color.
        # [rst Parallel deriv color setup 1 (beg)]
        if MPI.COMM_WORLD.size == 1:  # serial
            color1 = None
            color2 = None
        else:  # parallel
            # color1 will parallelize AS_point_0.fuelburn and AS_point_1.L_equals_W derivative computation
            color1 = "parcon1"
            # color2 will parallelize AS_point_0.CL and AS_point_1.wing_perf.failure derivative computation
            color2 = "parcon2"

        # AS_point_0.fuelburn, AS_point_0.CL, fuel_vol_delta.fuel_vol_delta, and fuel_diff depend only on AS_point_0
        # AS_point_1.L_equals_W and AS_point_1.wing_perf.failure depend only on AS_point_1
        prob.model.add_objective("AS_point_0.fuelburn", scaler=1e-5, parallel_deriv_color=color1)
        prob.model.add_constraint("AS_point_0.CL", equals=0.6, parallel_deriv_color=color2)
        prob.model.add_constraint("AS_point_1.L_equals_W", equals=0.0, parallel_deriv_color=color1)
        prob.model.add_constraint("AS_point_1.wing_perf.failure", upper=0.0, parallel_deriv_color=color2)
        prob.model.add_constraint("fuel_vol_delta.fuel_vol_delta", lower=0.0, parallel_deriv_color=None)
        prob.model.add_constraint("fuel_diff", equals=0.0, parallel_deriv_color=None)
        # [rst Parallel deriv color setup 1 (end)]

        # We will consider another dummy constraint on the sum of fuel burns from AS_point_0 and AS_point_1
        # (This constraint doesn't make any physical sense, but we do this to explain how parallel group works for reverse-mode derivatives)
        # fuel_sum depends on both AS_point_0 and AS_point_1. Linear solves for point_0 and point_1 will be parallelized.
        # [rst Parallel deriv color setup 2 (beg)]
        fuel_sum_comp = om.ExecComp("fuel_sum = fuelburn1 + fuelburn2", units="kg", shape=(1))
        prob.model.add_subsystem("fuel_sum", fuel_sum_comp, promotes_outputs=["fuel_sum"])
        prob.model.connect("AS_point_0.fuelburn", "fuel_sum.fuelburn1")
        prob.model.connect("AS_point_1.fuelburn", "fuel_sum.fuelburn2")
        prob.model.add_constraint("fuel_sum", lower=0, upper=1e10, scaler=1e-5)
        # [rst Parallel deriv color setup 2 (end)]

        # Set up the problem
        prob.setup()

        # change linear solvers.
        # [rst Change linear solver (beg)]
        if MPI.COMM_WORLD.size == 1:
            # serial
            prob.model.parallel.AS_point_1.coupled.linear_solver = om.LinearBlockGS(
                iprint=2, maxiter=30, use_aitken=True
            )
            prob.model.parallel.AS_point_0.coupled.linear_solver = om.LinearBlockGS(
                iprint=2, maxiter=30, use_aitken=True
            )
            prob.model.parallel.AS_point_1.coupled.linear_solver = om.LinearBlockGS(
                iprint=2, maxiter=30, use_aitken=True
            )

        elif MPI.COMM_WORLD.size == 2:
            # parallel. let each proc set the linear solver for their AS point
            if MPI.COMM_WORLD.rank == 0:
                prob.model.parallel.AS_point_0.coupled.linear_solver = om.LinearBlockGS(
                    iprint=2, maxiter=30, use_aitken=True
                )
            if MPI.COMM_WORLD.rank == 1:
                prob.model.parallel.AS_point_1.coupled.linear_solver = om.LinearBlockGS(
                    iprint=2, maxiter=30, use_aitken=True
                )
        # [rst Change linear solver (end)]

        # run aerostructural analysis
        start_time = time.time()
        prob.run_model()
        run_model_time = time.time() - start_time

        # compute derivatives
        start_time = time.time()
        totals = prob.compute_totals()
        derivs_time = time.time() - start_time

        print("I am processor", MPI.COMM_WORLD.rank, ", finished run_model and compute_totals")

        if MPI.COMM_WORLD.rank == 0:
            print("Analysis runtime: ", run_model_time, "[s]")
            print("Derivatives runtime: ", derivs_time, "[s]")

        assert_near_equal(MPI.COMM_WORLD.size, 2, 1e-8)
        assert_near_equal(prob.get_val("fuel_sum.fuel_sum", units="kg"), 5649.1290836, 1e-5)
        assert_near_equal(
            totals[("fuel_sum.fuel_sum", "wing.spar_thickness_cp")],
            np.array([[1712.12137573, 2237.99650867, 3036.45032547, 5065.16727605]]),
            1e-5,
        )


if __name__ == "__main__":
    unittest.main()
