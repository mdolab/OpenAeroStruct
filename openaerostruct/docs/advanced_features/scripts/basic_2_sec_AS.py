"""Optimizes the section twist distribution of a two section symmetrical wing using the construction-based approach for section
joining and the aerostructural tube model. This example is referenced as part of the multi-section tutorial."""

# docs checkpoint 0
import numpy as np
import openmdao.api as om
from openaerostruct.integration.aerostruct_groups import MultiSecAerostructGeometry, AerostructPoint
from openaerostruct.utils.constants import grav_constant
from openaerostruct.geometry.utils import build_section_dicts, unify_mesh, build_multi_spline, connect_multi_spline
import matplotlib.pyplot as plt


# docs checkpoint 1

# The geometry parameterization used here is identical to the one in the two section contruction based example. However,
# instead of chord we apply the principle to the twist B-spline.

# Set-up B-splines for each section. Done here since this information will be needed multiple times.
sec_twist_cp = [np.zeros(2), np.zeros(2)]

# Note the additional of structural variables to the surface dictionary
surface = {
    # Wing definition
    # Basic surface parameters
    "name": "surface",
    "is_multi_section": True,
    "num_sections": 2,  # The number of sections in the multi-section surface
    "sec_name": ["sec0", "sec1"],  # names of the individual sections
    "symmetry": True,  # if true, model one half of wing. reflected across the midspan of the root section
    "S_ref_type": "wetted",  # how we compute the wing area, can be 'wetted' or 'projected'
    "root_section": 1,
    # Geometry Parameters
    "taper": [1.0, 1.0],  # Wing taper for each section
    "span": [10.0, 10.0],  # Wing span for each section
    "sweep": [0.0, 0.0],  # Wing sweep for each section
    "twist_cp": sec_twist_cp,
    "t_over_c_cp": [np.array([0.15]), np.array([0.15])],  # thickness over chord ratio (NACA0015)
    "root_chord": 5.0,  # Wing root chord for each section
    # Mesh Parameters
    "meshes": "gen-meshes",  # Supply a mesh for each section or "gen-meshes" for automatic mesh generation
    "nx": 2,  # Number of chordwise points. Same for all sections
    "ny": [3, 3],  # Number of spanwise points for each section
    # Aerodynamic Parameters
    "CL0": 0.0,  # CL of the surface at alpha=0
    "CD0": 0.015,  # CD of the surface at alpha=0
    # Airfoil properties for viscous drag calculation
    "k_lam": 0.05,  # percentage of chord with laminar
    # flow, used for viscous drag
    "c_max_t": 0.303,  # chordwise location of maximum (NACA0015)
    # thickness
    "with_viscous": True,  # if true, compute viscous drag
    "with_wave": False,  # if true, compute wave drag
    "groundplane": False,
    # Structural
    "fem_model_type": "tube",
    "thickness_cp": 0.1 * np.ones((2)),
    "E": 70.0e9,  # [Pa] Young's modulus of the spar
    "G": 30.0e9,  # [Pa] shear modulus of the spar
    "yield": 500.0e6 / 2.5,  # [Pa] yield stress divided by 2.5 for limiting case
    "mrho": 3.0e3,  # [kg/m^3] material density
    "fem_origin": 0.35,  # normalized chordwise location of the spar
    "wing_weight_ratio": 2.0,
    "struct_weight_relief": False,  # True to add the weight of the structure to the loads on the structure
    "distributed_fuel_weight": False,
    # Constraints
    "exact_failure_constraint": False,  # if false, use KS function
}

# docs checkpoint 2

# Create the problem and assign the model group
prob = om.Problem(reports=False)

# Add problem information as an independent variables component
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

prob.model.add_subsystem("prob_vars", indep_var_comp, promotes=["*"])


# Generate the sections and unified mesh here in addition to adding the components.
# This has to also be done here since AeroPoint has to know the unified mesh size.
section_surfaces = build_section_dicts(surface)
uniMesh = unify_mesh(section_surfaces)
surface["mesh"] = uniMesh


# Build a component with B-spline control points that joins the sections by construction
twist_comp = build_multi_spline("twist_cp", len(section_surfaces), sec_twist_cp)
prob.model.add_subsystem("twist_bspline", twist_comp)

# Connect the B-spline component to the section B-splines
connect_multi_spline(prob, section_surfaces, sec_twist_cp, "twist_cp", "twist_bspline", surface["name"])

# docs checkpoint 3

# Create and add a group that handles the geometry for the
# aerostructual multi-section lifting surface
multi_geom_group = MultiSecAerostructGeometry(surface=surface)
prob.model.add_subsystem(surface["name"], multi_geom_group)

name = surface["name"]

point_name = "AS_point_0"

# Create the aero point group and add it to the model
AS_point = AerostructPoint(surfaces=[surface])

prob.model.add_subsystem(
    point_name,
    AS_point,
    promotes_inputs=[
        "v",
        "alpha",
        "Mach_number",
        "re",
        "rho",
        "CT",
        "R",
        "W0",
        "speed_of_sound",
        "empty_cg",
        "load_factor",
    ],
)

com_name = point_name + "." + name + "_perf"
prob.model.connect(name + ".local_stiff_transformed", point_name + ".coupled." + name + ".local_stiff_transformed")
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


# docs checkpoint 4

prob.driver = om.ScipyOptimizeDriver()
prob.driver.options["optimizer"] = "SLSQP"
prob.driver.options["tol"] = 1e-6
prob.driver.options["disp"] = True
prob.driver.options["maxiter"] = 1000


# Setup problem and add design variables, constraint, and objective
prob.model.add_design_var("twist_bspline.twist_cp_spline", lower=-10.0, upper=15.0)
prob.model.add_design_var("surface.thickness_cp", lower=0.01, upper=0.5, scaler=1e2)
prob.model.add_constraint("AS_point_0.surface_perf.failure", upper=0.0)
prob.model.add_constraint("AS_point_0.surface_perf.thickness_intersects", upper=0.0)

# Add design variables, constraints, and objective on the problem
prob.model.add_design_var("alpha", lower=-10.0, upper=10.0)
prob.model.add_constraint("AS_point_0.L_equals_W", equals=0.0)
prob.model.add_objective("AS_point_0.fuelburn", scaler=1e-5)


# Set up the problem
prob.setup(check=True)

# Run the optimization
optResult = prob.run_driver()
# om.n2(prob, show_browser=False)


# docs checkpoint 5

# Get the unified mesh and plot the results
meshUni = prob.get_val(name + "." + "surface_unification" + "." + name + "_uni_mesh")


def plot_meshes(meshes):
    """this function plots to plot the mesh"""
    plt.figure(figsize=(8, 4))
    ax = plt.gca()
    for i, mesh in enumerate(meshes):
        mesh_x = mesh[:, :, 0]
        mesh_y = mesh[:, :, 1]
        color = "k"
        for i in range(mesh_x.shape[0]):
            plt.plot(mesh_y[i, :], mesh_x[i, :], color, lw=1)
            plt.plot(-mesh_y[i, :], mesh_x[i, :], color, lw=1)  # plots the other side of symmetric wing
        for j in range(mesh_x.shape[1]):
            plt.plot(mesh_y[:, j], mesh_x[:, j], color, lw=1)
            plt.plot(-mesh_y[:, j], mesh_x[:, j], color, lw=1)  # plots the other side of symmetric wing
    plt.axis("equal")
    plt.xlabel("y (m)")
    plt.ylabel("x (m)")
    plt.savefig("opt_planform_AS.pdf")


# plot_meshes([mesh1,mesh2])
plot_meshes([meshUni])

# Print fuelburn
print(prob.get_val("AS_point_0.fuelburn"))
# docs checkpoint 6
