"""Optimizes the section chord distribution of a two section symmetrical wing with thickenss and viscous effects accounted
for. This examples is similar to the inviscid case but explains how to connect the unified t/c B-spline to the OAS performance
component for correct viscous drag computation. This example is referenced as part of the multi-section tutorial."""

import numpy as np

import openmdao.api as om

from openaerostruct.geometry.geometry_group import MultiSecGeometry
from openaerostruct.aerodynamics.aero_groups import AeroPoint
from openaerostruct.geometry.geometry_group import build_sections
from openaerostruct.geometry.geometry_unification import unify_mesh
from openaerostruct.geometry.multi_unified_bspline_utils import build_multi_spline, connect_multi_spline
import matplotlib.pyplot as plt


# Set-up B-splines for each section. Done here since this information will be needed multiple times.
sec_chord_cp = [np.array([1.0, 1.0]), np.array([1.0, 1.0])]


# Create a dictionary with info and options about the multi-section aerodynamic
# lifting surface
surface = {
    # Wing definition
    # Basic surface parameters
    "name": "surface",
    "is_multi_section": True,  # This key must be present for the AeroPoint to correctly interpret this surface as multi-section
    "num_sections": 2,  # The number of sections in the multi-section surface
    "sec_name": [
        "sec0",
        "sec1",
    ],  # names of the individual sections. Each section must be named and the list length must match the specified number of sections.
    "symmetry": True,  # if true, model one half of wing. reflected across the midspan of the root section
    "S_ref_type": "wetted",  # how we compute the wing area, can be 'wetted' or 'projected'
    # Geometry Parameters
    "taper": [1.0, 1.0],  # Wing taper for each section. The list length must match the specified number of sections.
    "span": [2.0, 2.0],  # Wing span for each section. The list length must match the specified number of sections.
    "sweep": [0.0, 0.0],  # Wing sweep for each section. The list length must match the specified number of sections.
    "chord_cp": [
        np.array([1, 1]),
        np.array([1, 1]),
    ],  # The chord B-spline parameterization for EACH SECTION. The list length must match the specified number of sections.
    "twist_cp": [
        np.zeros(2),
        np.zeros(2),
    ],  # The twist B-spline parameterization for EACH SECTION. The list length must match the specified number of sections.
    "t_over_c_cp": [
        np.array([0.15]),
        np.array([0.15]),
    ],  # The thickenss to chord ratio B-spline parameterization for EACH SECTION. The list length must match the specified number of sections.
    "root_chord": 1.0,  # Root chord length of the section indicated as "root section"(required if using the built-in mesh generator)
    # Mesh Parameters
    "meshes": "gen-meshes",  # Supply a list of meshes for each section or "gen-meshes" for automatic mesh generation
    "nx": 2,  # Number of chordwise points. Same for all sections.(required if using the built-in mesh generator)
    "ny": [
        21,
        21,
    ],  # Number of spanwise points for each section. The list length must match the specified number of sections. (required if using the built-in mesh generator)
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

# Generate the sections and unified mesh here. It's needed to join the sections by construction.
section_surfaces = build_sections(surface)
uniMesh = unify_mesh(section_surfaces)
surface["mesh"] = uniMesh

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
prob.model.add_subsystem(point_name, aero_group, promotes_inputs=["v", "alpha", "Mach_number", "re", "rho", "cg"])

# Get name of surface and construct unified mesh name
name = surface["name"]
unification_name = "{}_unification".format(surface["name"])

# Connect the mesh from the mesh unification component to the analysis point
prob.model.connect(name + "." + unification_name + "." + name + "_uni_mesh", point_name + "." + "surface" + ".def_mesh")

# Perform the connections with the modified names within the
# 'aero_states' group.
prob.model.connect(
    name + "." + unification_name + "." + name + "_uni_mesh", point_name + ".aero_states." + "surface" + "_def_mesh"
)


# docs checkpoint 0

"""If the surface features a thickness to chord ratio B-spline, then the Muli-section geometry component will automatically
combine them together while enforcing C0 continuity. This process will occurs regardless of whether the constraint or construction-based
section joining method are used as otherwise OAS will not be able to calculate the viscous drag correctly. Pay close attention to how the
unified t_over_c B-spline(given a unique name by the multi-section geometry group) needs to be connect to the AeroPoint performance component."""
prob.model.connect(
    name + "." + unification_name + "." + name + "_uni_t_over_c", point_name + "." + name + "_perf." + "t_over_c"
)

# docs checkpoint 1

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
# om.n2(prob)


mesh1 = prob.get_val("surface.sec0.mesh", units="m")
mesh2 = prob.get_val("surface.sec1.mesh", units="m")

meshUni = prob.get_val(name + "." + unification_name + "." + name + "_uni_mesh")


def plot_meshes(meshes):
    """this function plots to plot the mesh"""
    plt.figure(figsize=(8, 4))
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
    # plt.savefig('opt_planform_visc_construction.png')


plot_meshes([meshUni])
# plt.show()
