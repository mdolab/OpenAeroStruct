"""Optimizes the section chord distribution of a two section symmetrical wing using the construction-based approach for section
joining. This example is referenced as part of the multi-section tutorial."""

# docs checkpoint 0
import numpy as np
import openmdao.api as om
from openaerostruct.geometry.geometry_group import MultiSecGeometry
from openaerostruct.aerodynamics.aero_groups import AeroPoint
from openaerostruct.geometry.utils import build_section_dicts
from openaerostruct.geometry.utils import unify_mesh
from openaerostruct.geometry.utils import build_multi_spline, connect_multi_spline
import matplotlib.pyplot as plt

# docs checkpoint 1

# The multi-section geometry parameterization number section from left to right starting with section #0. A two-section symmetric wing parameterization appears as follows.
# For a symmetrical wing the last section in the sequence will always be marked as the "root section" as it's adjacent to the geometric centerline of the wing.
# Geometeric parameters must be specified for each section using lists with values corresponding in order of the surface numbering. Section section supports all the
# standard OpenAeroStruct geometery transformations including B-splines.


"""

-----------------------------------------------  ^
|                      |                       | |
|                      |                       | |
|        sec 0         |         sec 1         | | root         symmetrical BC
|                      |     "root section"    | | chord
|______________________|_______________________| |
                                                 _
                                              y = 0 ------------------> + y

"""


# A multi-section surface dictionary is very similar to the standard one. However, it features some additional options and requires that the user specify
# parameters for each desired section. The multi-section geometery group also features an built in mesh generator so the wing mesh parameters can be specified right
# in the surface dictionary. Let's create a dictionary with info and options for a two-section aerodynamic lifting surface


# We will set up our chord_cp as seperate variable since we will need to use it several times in this example.
sec_chord_cp = [np.ones(2), np.ones(2)]


# Create a dictionary with info and options about the multi-section aerodynamic lifting surface
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
    "chord_cp": sec_chord_cp,  # The chord B-spline parameterization for each section. The list length must match the specified number of sections.
    "twist_cp": [
        np.zeros(2),
        np.zeros(2),
    ],  # The twist B-spline parameterization for each section. The list length must match the specified number of sections.
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
    "with_viscous": False,  # if true, compute viscous drag
    "with_wave": False,  # if true, compute wave drag
    "groundplane": False,
}

# docs checkpoint 2

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

# docs checkpoint 3

"""Instead of creating a standard geometery group, here we will create a multi-section geometry group that will accept our multi-section surface
dictionary. In this example we will constrain the sections into a C0 continuous surface with a construction approach that assigns the chord B-spline
control points at each section junction to an index of a global control vector. """

# In order to construct this global B-spline control vector we first need to generate the unified surface mesh.
# The unified surface mesh is simply all the individual section surface meshes combine into a single unified OAS mesh array.

# First we will call the utility function build_section_dicts which takes the surface dictionary and outputs a list of surface dictionaries corresponding to
# each section.
section_surfaces = build_section_dicts(surface)

# We can then call unify_mesh which outputs the unified mesh of all of the sections.
uniMesh = unify_mesh(section_surfaces)

# We can then assign the unified mesh as the mesh for the entire surface.
surface["mesh"] = uniMesh


"""This functions builds an OpenMDAO Independent Variable Component with the correct length input vector
corresponding to each section junction on the surface plus the wingtips. Refer to the functions documentions for input details. After
the compnent has been generated it needs to be added to the model."""
chord_comp = build_multi_spline("chord_cp", surface["num_sections"], sec_chord_cp)
prob.model.add_subsystem("chord_bspline", chord_comp)

"""In order to properly transform the surface geometry the surface's global input vector need to be connected
to the corresponding control points of the local B-spline component on each section. This function automates this
process as it can be tedious.

The figure below explains how the global control vector's outputs are connected to the control points of the local
section B-splines. In this example, each section features a two point B-spline with control points at the section tips
however the principle is the same for B-splines with more points.


                surface B-spline
0;;;;;;;;;;;;;;;;;;;;;;1;;;;;;;;;;;;;;;;;;;;;;;2
^                      ^                       ^
|                      |                       |
|                      |                       |
|                      |    sec 1 B-spline     |
     sec 0 B-spline    c:::::::::::::::::::::::d
a::::::::::::::::::::::b
-----------------------------------------------  ^
|                      |                       | |
|                      |                       | |
|        sec 0         |         sec 1         | | root
|                      |                       | | chord
|______________________|_______________________| |
                                                 _
                                              y = 0 ------------------> + y


An edge case in this process is when a section features a B-spline with a single control point. The same control point
cannot be assigned to two different control points on the surface B-spline. In these situations a constraint will need
to be used to maintain C0 continuity. See the connect_multi_spline documentation for details.
"""
connect_multi_spline(prob, section_surfaces, sec_chord_cp, "chord_cp", "chord_bspline", surface["name"])


""" With the surface B-spline connected we can add the multi-section geometry group."""
multi_geom_group = MultiSecGeometry(surface=surface)
prob.model.add_subsystem(surface["name"], multi_geom_group)

# docs checkpoint 4

# Create the aero point group, which contains the actual aerodynamic
# analyses. This step is exactly as it's normally done except the surface dictionary we pass in is the multi-surface one
aero_group = AeroPoint(surfaces=[surface])
point_name = "aero_point_0"
prob.model.add_subsystem(point_name, aero_group, promotes_inputs=["v", "alpha", "Mach_number", "re", "rho", "cg"])

# docs checkpoint 5

# The following steps are similar to a normal OAS surface script but note the differences in surface naming. Note that
# unified surface created by the multi-section geometry group needs to be connected to AeroPoint(be careful with the naming)

# Get name of surface and construct the name of the unified surface mesh
name = surface["name"]
unification_name = "{}_unification".format(surface["name"])

# Connect the mesh from the mesh unification component to the analysis point
prob.model.connect(name + "." + unification_name + "." + name + "_uni_mesh", point_name + "." + "surface" + ".def_mesh")

# Perform the connections with the modified names within the
# 'aero_states' group.
prob.model.connect(
    name + "." + unification_name + "." + name + "_uni_mesh", point_name + ".aero_states." + "surface" + "_def_mesh"
)

# docs checkpoint 6

# Next, we add the DVs to the OpenMDAO problem.
# Here we use the global independent variable component vector and the angle-of-attack as DVs.
prob.model.add_design_var("chord_bspline.chord_cp_spline", lower=0.1, upper=10.0, units=None)
prob.model.add_design_var("alpha", lower=0.0, upper=10.0, units="deg")

# Add CL constraint
prob.model.add_constraint(point_name + ".CL", equals=0.3)

# Add Wing total area constraint
prob.model.add_constraint(point_name + ".total_perf.S_ref_total", equals=2.0)

# Add objective
prob.model.add_objective(point_name + ".CD", scaler=1e4)

prob.driver = om.ScipyOptimizeDriver()
prob.driver.options["optimizer"] = "SLSQP"
prob.driver.options["tol"] = 1e-7
prob.driver.options["disp"] = True
prob.driver.options["maxiter"] = 1000

# Set up and run the optimization problem
prob.setup()
prob.run_driver()
# om.n2(prob)

# docs checkpoint 7

# Get each section mesh
mesh1 = prob.get_val("surface.sec0.mesh", units="m")
mesh2 = prob.get_val("surface.sec1.mesh", units="m")

# Get the unified mesh
meshUni = prob.get_val(name + "." + unification_name + "." + name + "_uni_mesh")


# Plot the results
def plot_meshes(meshes):
    """this function plots a list of meshes on the same plot."""
    plt.figure(figsize=(8, 4))
    for i, mesh in enumerate(meshes):
        mesh_x = mesh[:, :, 0]
        mesh_y = mesh[:, :, 1]
        color = "w"
        for i in range(mesh_x.shape[0]):
            plt.plot(mesh_y[i, :], 1 - mesh_x[i, :], color, lw=1)
            plt.plot(-mesh_y[i, :], 1 - mesh_x[i, :], color, lw=1)  # plots the other side of symmetric wing
        for j in range(mesh_x.shape[1]):
            plt.plot(mesh_y[:, j], 1 - mesh_x[:, j], color, lw=1)
            plt.plot(-mesh_y[:, j], 1 - mesh_x[:, j], color, lw=1)  # plots the other side of symmetric wing
    plt.axis("equal")
    plt.xlabel("y (m)")
    plt.ylabel("x (m)")
    plt.savefig("opt_planform_construction.pdf")


plot_meshes([meshUni])
# docs checkpoint 8
