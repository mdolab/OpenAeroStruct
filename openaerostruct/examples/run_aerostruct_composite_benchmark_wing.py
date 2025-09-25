"""
Aerostructural optimization example of the "Simple Composite Wing" from:
Gray and Martins, A Proposed Benchmark Model for Practical Aeroelastic Optimization of Aircraft Wings, AIAA 2024-2775.
https://www.researchgate.net/publication/377154425_A_Proposed_Benchmark_Model_for_Practical_Aeroelastic_Optimization_of_Aircraft_Wings
"""

import matplotlib.pyplot as plt
import numpy as np
from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.integration.aerostruct_groups import AerostructGeometry, AerostructPoint
from openaerostruct.structures.wingbox_fuel_vol_delta import WingboxFuelVolDelta
from openaerostruct.structures.utils import compute_composite_stiffness
import openmdao.api as om

# --- Case control ---
# case_name = "Analysis"   # trim analysis for an initial design.
case_name = "Case1"  # structural optimization
# case_name = "Case2"  # aerostructural optimization with fixed planform
# case_name = "Case3"  # aerostructural optimization with varying planform (span, chord, and sweep)

# --- Airfoil definition for wingbox structure model ---
# We use the 15% to 65% chord portion of the RAE2822 airfoil.
# fmt: off
upper_x = np.array([0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.27, 0.29, 0.31, 0.33, 0.35, 0.37, 0.39, 0.41, 0.43, 0.45, 0.47, 0.49, 0.51, 0.53, 0.55, 0.57, 0.59, 0.61, 0.63, 0.65], dtype="complex128")
lower_x = np.array([0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.27, 0.29, 0.31, 0.33, 0.35, 0.37, 0.39, 0.41, 0.43, 0.45, 0.47, 0.49, 0.51, 0.53, 0.55, 0.57, 0.59, 0.61, 0.63, 0.65], dtype="complex128")
upper_y = np.array([0.04591116991110611,   0.04839609476913142,   0.05062208326872383,  0.05261581485173187,   0.054397682405710426, 0.055983808645155245,  0.057386804014598544, 0.05861859395614192,  0.05968942721753942,  0.060589981875354845,  0.0613225097626915,   0.061904166920345005,  0.06233885550219253,   0.06262971153846154,  0.06277830920512399,   0.06277420344063064,  0.06258449250255363,  0.062233206407434585, 0.061713111844786825,  0.0610172702757916,  0.060139111717851856, 0.0590313764575464,   0.05772568052109181,   0.05624323825433285,  0.05459345628276678,   0.05278581315710423], dtype="complex128")  # noqa: E201, E241
lower_y = np.array([-0.04604238094970181, -0.048404951781293726, -0.05050769979841836, -0.052371999003239474, -0.05401536186939327, -0.055450351894347404, -0.05668039589416059, -0.057700571759671956, -0.05849340424509761, -0.05900613517054636, -0.05919639050765996, -0.059068267715203786, -0.058614112758404675, -0.05783551674937965, -0.056741090244703564, -0.05532706495319429, -0.05358094598569969, -0.05159666756338143, -0.049403792369772555, -0.0470318138100102, -0.04450996703071111, -0.04184532587452785, -0.039089564516129036, -0.03626824702443099, -0.033399043294435986, -0.03049965785521177], dtype="complex128")
# fmt: on

# --- Geometry and analysis options ---
# Create a dictionary to store options about the surface
mesh_dict = {
    "num_y": 51,
    "num_x": 7,
    "wing_type": "rect",
    "symmetry": True,
    "chord_cos_spacing": 0,
    "span_cos_spacing": 0,
    "root_chord": 1.0,  # apply chord later.
    "span": 28.0,
}

mesh = generate_mesh(mesh_dict)

surf_dict = {
    # Wing definition
    "name": "wing",  # name of the surface
    "symmetry": True,  # if true, model one half of wing
    "S_ref_type": "wetted",  # how we compute the wing area, can be 'wetted' or 'projected'
    "mesh": mesh,
    # planform defition
    "span": 28.0,
    "chord_cp": np.array([1.5, 5.0]),  # tip, root [m]
    "sweep": 28.18,  # leading edge sweep [deg]
    "ref_axis_pos": 0,  # set this to 0 to define the sweep for leading edge.
    # twist distribution
    "twist_cp": np.zeros(5),
    "fem_model_type": "wingbox",
    "data_x_upper": upper_x,
    "data_x_lower": lower_x,
    "data_y_upper": upper_y,
    "data_y_lower": lower_y,
    "spar_thickness_cp": np.linspace(0.003, 0.02, 10),  # initial guess, [m]
    "skin_thickness_cp": np.linspace(0.003, 0.02, 10),  # initial guess, [m]
    "original_wingbox_airfoil_t_over_c": 0.121,
    # Aerodynamic parameters and options
    "CL0": 0.0,
    "CD0": 0.01508,  # CD0 = 0.01508 for fuselage + nacelle + tails from Gray 2024 + wing drag calibration to match Case 1 L/D
    "with_viscous": True,  # if true, compute viscous drag
    "with_wave": True,  # if true, compute wave drag
    # Airfoil properties for viscous drag calculation
    "k_lam": 0.05,  # percentage of chord with laminar flow, used for viscous drag
    "c_max_t": 0.38,  # chordwise location of maximum thickness
    "t_over_c_cp": np.array([0.121, 0.121, 0.121, 0.121, 0.121]),
    "wing_weight_ratio": 3.03,  # from Case 1 results of Gray 2024 (they used nonlinear weight model)
    "exact_failure_constraint": False,  # if false, use KS function
    "struct_weight_relief": True,
    "distributed_fuel_weight": True,
    "fuel_density": 804.0,  # [kg/m^3] fuel density (only needed if the fuel-in-wing volume constraint is used)
    "Wf_reserve": 2000.0,  # [kg] reserve fuel mass
    # Composite material parameters
    "useComposite": True,
    "mrho": 1550,  # [kg/m^3]
    "safety_factor": 1.5,
    "ply_angles": [0, 45, -45, 90],
    "ply_fractions": [0.4441, 0.222, 0.222, 0.1119],  # skin layup from Gray 2024
    "E1": 117.7e9,
    "E2": 9.7e9,
    "nu12": 0.35,
    "G12": 4.8e9,
    "sigma_t1": 1648.0e6,  # longitudinal tensile strength
    "sigma_c1": 1034.0e6,  # longitudinal compressive strength
    "sigma_t2": 64.0e6,  # transverse tensile strength
    "sigma_c2": 228.0e6,  # transverse compressive strength
    "sigma_12max": 71.0e6,  # maximum shear strength
}

# Compute effective E and G for composite material
compute_composite_stiffness(surf_dict)

surfaces = [surf_dict]

# --- OpenMDAO problem setup ---
prob = om.Problem()

# flight conditions (cruise, 2.5g pull-up maneuver)
mach_number = np.array([0.77, 0.458])  # cruise, 2.5g maneuver
speed_of_sound = np.array([297.71, 340.294])
v = mach_number * speed_of_sound
rho_air = np.array([0.39263, 1.225])
mu_air = np.array([1.454e-5, 1.812e-5])

# Add problem information as an independent variables component
indep_var_comp = om.IndepVarComp()
indep_var_comp.add_output("v", val=v, units="m/s")
indep_var_comp.add_output("alpha", val=4.0, units="deg")  # initial guess
indep_var_comp.add_output("alpha_maneuver", val=9.0, units="deg")  # initial guess
indep_var_comp.add_output("Mach_number", val=mach_number)
indep_var_comp.add_output("re", val=rho_air * v / mu_air, units="1/m")
indep_var_comp.add_output("rho", val=rho_air, units="kg/m**3")
indep_var_comp.add_output("CT", val=0.64 / 3600, units="1/s")  # 0.64 lbm/lbf-h = 1.81e-5 kg/N-s
indep_var_comp.add_output("R", val=3815, units="km")
indep_var_comp.add_output("W0", val=39500 + surf_dict["Wf_reserve"], units="kg")
indep_var_comp.add_output("speed_of_sound", val=speed_of_sound, units="m/s")
indep_var_comp.add_output("load_factor", val=np.array([1.0, 2.5]))
indep_var_comp.add_output("empty_cg", val=np.zeros((3)), units="m")
indep_var_comp.add_output("fuel_mass", val=10000, units="kg")  # initial guess

prob.model.add_subsystem("prob_vars", indep_var_comp, promotes=["*"])

# Compute mid-cruise fuel mass
mid_cruise_fuel_mass_comp = om.ExecComp("fuel_mid_cruise = fuel_mass / 2", units="kg")
prob.model.add_subsystem("mid_cruise_fuel_mass", mid_cruise_fuel_mass_comp, promotes_inputs=["fuel_mass"], promotes_outputs=["fuel_mid_cruise"])

# Add geometry groups
for surface in surfaces:
    # Get the surface name and create a group to contain components
    # only for this surface
    name = surface["name"]

    aerostruct_group = AerostructGeometry(surface=surface)

    # Add group to the problem with the name of the surface.
    prob.model.add_subsystem(name, aerostruct_group)

# Add analysis flight conditions
for i in range(2):
    point_name = "AS_point_{}".format(i)
    # Connect the parameters within the model for each aerostruct point

    # Create the aerostruct point group and add it to the model
    AS_point = AerostructPoint(surfaces=surfaces, internally_connect_fuelburn=False)

    prob.model.add_subsystem(point_name, AS_point)

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
        prob.model.connect(name + ".structural_mass", point_name + "." + "total_perf." + name + "_structural_mass")

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
    prob.model.connect("fuel_mid_cruise", "AS_point_0.coupled.wing.struct_states.fuel_mass")

    prob.model.connect("wing.struct_setup.fuel_vols", "AS_point_1.coupled.wing.struct_states.fuel_vols")
    # set zero fuel for maneuver points
    prob.model.set_input_defaults("AS_point_1.coupled.wing.struct_states.fuel_mass", 0.0)

comp = om.ExecComp("fuel_diff = (fuel_mass - fuelburn) / fuelburn", units="kg")
prob.model.add_subsystem("fuel_diff", comp, promotes_inputs=["fuel_mass"], promotes_outputs=["fuel_diff"])
prob.model.connect("AS_point_0.fuelburn", "fuel_diff.fuelburn")

# compute wing loading
wing_loading_group = prob.model.add_subsystem("wing_loading", om.Group(), promotes=["*"])
mtow_comp = om.ExecComp("MTOW = W0 + fuel_mass + wing_str_mass", units="kg")
wing_loading_group.add_subsystem("MTOW", mtow_comp, promotes_inputs=["W0", "fuel_mass"], promotes_outputs=["MTOW"])
prob.model.connect("wing.structural_mass", "MTOW.wing_str_mass")

wing_loading_comp = om.ExecComp(
    "wing_loading = MTOW / wing_area",
    wing_loading={"units": "kg/m**2"},
    MTOW={"units": "kg"},
    wing_area={"units": "m**2"},
)
wing_loading_group.add_subsystem("wing_loading", wing_loading_comp, promotes_inputs=["MTOW"], promotes_outputs=["wing_loading"])
prob.model.connect("AS_point_0.coupled.wing.S_ref", "wing_loading.wing_area")

# --- Optimizer settings ---
## Use these settings if you do not have pyOptSparse or SNOPT
# prob.driver = om.ScipyOptimizeDriver()
# prob.driver.options["optimizer"] = "SLSQP"
# prob.driver.options["tol"] = 1e-8

prob.driver = om.pyOptSparseDriver()
prob.driver.options['print_results'] = True
prob.driver.options['optimizer'] = 'SNOPT'
prob.driver.opt_settings['Major iterations limit'] = 100
prob.driver.opt_settings['Major feasibility tolerance'] = 1e-6
prob.driver.opt_settings['Major optimality tolerance'] = 1e-6
prob.driver.opt_settings['Verify level'] = -1
prob.driver.opt_settings['Function precision'] = 1e-10
prob.driver.opt_settings['Hessian full memory'] = 1
prob.driver.opt_settings['Hessian frequency'] = 100

# --- Define optimization problem ---
if case_name == "Case1":
    prob.model.add_objective("wing.structural_mass", units="kg", ref=1000)
else:
    prob.model.add_objective("AS_point_0.fuelburn", units="kg", ref=10000)

# angle of attack and L=W constraint.
prob.model.add_design_var("alpha", lower=0.0, upper=20, ref=10)
prob.model.add_design_var("alpha_maneuver", lower=0.0, upper=30, ref=10)
prob.model.add_constraint("AS_point_0.L_equals_W", equals=0.0)
prob.model.add_constraint("AS_point_1.L_equals_W", equals=0.0)

# structural design variables and failure constraint
if case_name in ["Case1", "Case2", "Case3"]:
    prob.model.add_design_var("wing.spar_thickness_cp", lower=0.002, upper=0.1, ref=0.01)
    prob.model.add_design_var("wing.skin_thickness_cp", lower=0.002, upper=0.1, ref=0.01)
    prob.model.add_constraint("AS_point_1.wing_perf.failure", upper=0.0)

# twist and t/c variables
if case_name in ["Case2", "Case3"]:
    twist_indices = range(len(surf_dict["twist_cp"]) - 1)   # exclude root twist = last index
    prob.model.add_design_var("wing.twist_cp", lower=-15.0, upper=15.0, scaler=0.1, indices=twist_indices)
    prob.model.add_design_var("wing.geometry.t_over_c_cp", lower=0.07, upper=0.2, scaler=10.0)

# planform variable
if case_name in ["Case3"]:
    prob.model.add_design_var("wing.geometry.span", lower=25, upper=40, ref=28.0)
    prob.model.add_design_var("wing.geometry.chord_cp", lower=1.0, upper=8.0, ref=1.5)
    prob.model.add_design_var("wing.sweep", lower=15.0, upper=35.0, ref=28.18)

    # constraint wing loading
    prob.model.add_constraint("wing_loading", units="kg/m**2", upper=600, ref=600)

# fuel mass variable and consistency constraint
if case_name in ["Case2", "Case3"]:
    prob.model.add_design_var("fuel_mass", lower=5000, upper=20000, ref=10000)
    prob.model.add_constraint("fuel_diff", equals=0.0)
    prob.model.add_constraint("fuel_vol_delta.fuel_vol_delta", lower=0.0)

# Set up the problem
prob.setup()

# change linear solver for aerostructural coupled adjoint
prob.model.AS_point_0.coupled.linear_solver = om.PETScKrylov(assemble_jac=True, iprint=0, rhs_checking=True)
prob.model.AS_point_0.coupled.linear_solver.precon = om.LinearRunOnce(iprint=-1)
prob.model.AS_point_1.coupled.linear_solver = om.PETScKrylov(assemble_jac=True, iprint=0, rhs_checking=True)
prob.model.AS_point_1.coupled.linear_solver.precon = om.LinearRunOnce(iprint=-1)

# Use LinearBlockGS instead if you don't have PETSc installed
# prob.model.AS_point_0.coupled.linear_solver = om.LinearBlockGS(iprint=0, maxiter=30, use_aitken=True)
# prob.model.AS_point_1.coupled.linear_solver = om.LinearBlockGS(iprint=0, maxiter=30, use_aitken=True)

prob.run_driver()

om.n2(prob, show_browser=False, outfile="n2_simple_transonic_wing.html")

# --- Print results ---
fuel_burn = prob.get_val("AS_point_0.fuelburn", units="kg")[0]
wing_mass = prob.get_val("wing.structural_mass", units="kg")[0]
cruise_CD = prob.get_val("AS_point_0.CD")[0]
cruise_CL = prob.get_val("AS_point_0.CL")[0]
cruise_LbyD = cruise_CL / cruise_CD

print("\n--------------------------------")
print("Fuel burn =", fuel_burn, "[kg]")
print("Wingbox structure mass =", wing_mass / surf_dict["wing_weight_ratio"], "[kg]")
print("Wing total mass =", wing_mass, "[kg]")
print("Cruise CD =", cruise_CD)
print("Cruise CL =", cruise_CL)
print("Cruise L/D =", cruise_LbyD)

print("\nComposite effective stiffness:")
print('E =', surf_dict["E"] / 1e9)
print('G =', surf_dict["G"] / 1e9)

print("\nStructural variables")
print('Spar thickness cp =', prob.get_val("wing.spar_thickness_cp", units="mm"), "[mm] (tip to root)")
print('Skin thickness cp =', prob.get_val("wing.skin_thickness_cp", units="mm"), "[mm] (tip to root)")

print("\nAero/aerostructural variables")
print('Twist cp =', prob.get_val("wing.twist_cp", units="deg"), "[deg] (tip to root)")
print('t/c cp =', prob.get_val("wing.geometry.t_over_c_cp"), "(tip to root)")

print("\nPlanform variables")
print('Span =', prob.get_val("wing.geometry.span", units="m"), "[m]")
print('Chord cp =', prob.get_val("wing.geometry.chord_cp"), "[m] (tip, root)")
print('LE sweep =', prob.get_val("wing.sweep", units="deg"), "[deg]")

# --- plot planform ---
mesh = prob.get_val("wing.mesh")
mesh_x = mesh[:, :, 0]
mesh_y = mesh[:, ::-1, 1] * -1
le_root_x = mesh[0, -1, 0]
mesh_x -= le_root_x

plt.figure()
plt.plot(mesh_y[0, :], mesh_x[0, :], color='k')  # LE
plt.plot(mesh_y[-1, :], mesh_x[-1, :], color='k')  # TE
plt.plot(mesh_y[:, 0], mesh_x[:, 0], color='k')  # tip
plt.plot(mesh_y[:, -1], mesh_x[:, -1], color='k')  # root
plt.axis('equal')
plt.gca().invert_xaxis()
plt.grid()
plt.savefig("simple_transonic_wing_planform.pdf", bbox_inches='tight')

# --- plot structural thickness ---
y = mesh_y[0, :]
skin_thickness = prob.get_val("wing.skin_thickness", units="mm").ravel()[::-1]
spar_thickness = prob.get_val("wing.spar_thickness", units="mm").ravel()[::-1]
skin_t_step = np.concatenate([[skin_thickness[0]], skin_thickness])
spar_t_step = np.concatenate([[spar_thickness[0]], spar_thickness])

fig, axs = plt.subplots(2, 1, figsize=(6, 6))
axs[0].step(y, skin_t_step, where='pre', lw=2)
axs[0].set_ylabel("Skin Thickness (mm)")
axs[0].set_xticklabels([])

axs[1].step(y, spar_t_step, where='pre', lw=2)
axs[1].set_ylabel("Spar Thickness (mm)")
axs[1].set_xlabel("Spanwise (m)")
plt.savefig("simple_transonic_wing_thickness.pdf", bbox_inches='tight')

# --- plot twist and t/c---
# jig twist
jig_mesh = prob.get_val("wing.mesh", units="m")
jig_LE = jig_mesh[0, :, :]
jig_TE = jig_mesh[-1, :, :]
jig_chord = jig_TE[:, 0] - jig_LE[:, 0]
twist_jig = np.arctan2(jig_LE[:, 2] - jig_TE[:, 2], jig_chord) * 180 / np.pi

# cruise twist
cruise_mesh = prob.get_val("AS_point_0.coupled.wing.def_mesh", units="m")
cruise_LE = cruise_mesh[0, :, :]
cruise_TE = cruise_mesh[-1, :, :]
cruise_chord = cruise_TE[:, 0] - cruise_LE[:, 0]
twist_cruise = np.arctan2(cruise_LE[:, 2] - cruise_TE[:, 2], cruise_chord) * 180 / np.pi

# 2.5g twist
maneuver_mesh = prob.get_val("AS_point_1.coupled.wing.def_mesh", units="m")
maneuver_LE = maneuver_mesh[0, :, :]
maneuver_TE = maneuver_mesh[-1, :, :]
maneuver_chord = maneuver_TE[:, 0] - maneuver_LE[:, 0]
twist_maneuver = np.arctan2(maneuver_LE[:, 2] - maneuver_TE[:, 2], maneuver_chord) * 180 / np.pi

fig, axs = plt.subplots(2, 1, figsize=(6, 6))
axs[0].plot(y, twist_jig[::-1], color='darkgray', label='jig')
axs[0].plot(y, twist_cruise[::-1], color='C0', label='1g cruise')
axs[0].plot(y, twist_maneuver[::-1], color='C1', label='2.5g pull-up')
axs[0].set_xticklabels([])
axs[0].set_ylabel("Twist (deg)")
axs[0].legend()

# t/c
t_over_c = prob.get_val("wing.t_over_c").ravel()[::-1]
t_over_c_step = np.concatenate([[t_over_c[0]], t_over_c])
axs[1].step(y, t_over_c_step, where='pre', color='C0')
axs[1].set_xlabel("Spanwise (m)")
axs[1].set_ylabel("t/c")
axs[1].set_xlabel("Spanwise (m)")
plt.savefig("simple_transonic_wing_twist_tc.pdf", bbox_inches='tight')

plt.show()
