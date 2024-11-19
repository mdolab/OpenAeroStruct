"""Optimizes the section chord distribution of a two section symmetrical wing using the constraint-based approach for section
joining. This example is referenced as part of the multi-section tutorial."""

# docs checkpoint 0
import numpy as np
import openmdao.api as om
from openaerostruct.geometry.geometry_group import MultiSecGeometry
from openaerostruct.aerodynamics.aero_groups import AeroPoint
from openaerostruct.geometry.geometry_group import build_sections
from openaerostruct.geometry.geometry_unification import unify_mesh
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

# Instead of creating a standard geometery group, here we will create a multi-section group that will accept our multi-section surface
# dictionary and allow us to specify any C0 continuity constraints between the sections. In this example we will constrain the sections
# into a C0 continuous surface using a component that the optimizer can use as a constraint. The joining constraint component returns the
# distance between the leading edge and trailing edge points at section interections. Any combination of the x,y, and z distances can be returned
# to constrain the surface in a particular direction.


"""
                        LE1              LE2   cLE = [LE2x-LE1x,LE2y-LE1y,LE2z-LE1z]
------------------------- <-------------> -------------------------
|                       |                 |                       |
|                       |                 |                       |
|        sec 0          |                 |         sec 1         |
|                       |                 |     "root section"    |
|______________________ | <-------------> |_______________________|
                        TE1              TE2    cTE = [TE2x-TE1x,TE2y-TE1y,TE2z-TE1z]



"""

# We pass in the multi-section surface dictionary to the MultiSecGeometry geometery group. We also enabled joining_comp and pass two array to dim_contr
# These two arrays should only consists of 1 and 0 and tell the joining component which of the x,y, and z distance constraints we wish to enforce at the LE and TE
# In this example, we only wish to constraint the x-distance between the sections at both the leading and trailing edge.

multi_geom_group = MultiSecGeometry(
    surface=surface, joining_comp=True, dim_constr=[np.array([1, 0, 0]), np.array([1, 0, 0])]
)
prob.model.add_subsystem(surface["name"], multi_geom_group)

# docs checkpoint 4

# In this next part, we will setup the aerodynamics group. First we use a utility function called build_sections which takes our multi-section surface dictionary and outputs a
# surface dictionary for each individual section. We then inputs these dictionaries into the mesh unification function unify_mesh to produce a single mesh array for the the entire surface.
# We then add this mesh to the multi-section surface dictionary
section_surfaces = build_sections(surface)
uniMesh = unify_mesh(section_surfaces)
surface["mesh"] = uniMesh

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

# Connect the mesh from the mesh unification component to the analysis point.
prob.model.connect(name + "." + unification_name + "." + name + "_uni_mesh", point_name + "." + "surface" + ".def_mesh")

# Perform the connections with the modified names within the
# 'aero_states' group.
prob.model.connect(
    name + "." + unification_name + "." + name + "_uni_mesh", point_name + ".aero_states." + "surface" + "_def_mesh"
)

# docs checkpoint 6

# Next, we add the DVs to the OpenMDAO problem. Note that each surface's geometeric parameters are under the given section names specified in the multi-surface dictionary earlier.
# Here we use the chord B-spline that we specified earlier for each section and the angle-of-attack as DVs.
prob.model.add_design_var("surface.sec0.chord_cp", lower=0.1, upper=10.0, units=None)
prob.model.add_design_var("surface.sec1.chord_cp", lower=0.1, upper=10.0, units=None)
prob.model.add_design_var("alpha", lower=0.0, upper=10.0, units="deg")


# Next, we add the C0 continuity constraint for this problem by constraining the x-distance between sections to 0.
# NOTE: SLSQP optimizer does not handle the joining equality constraint properly so the constraint needs to be specified as an inequality constraint
# All other optimizers like SNOPT can handle the equality constraint as is.

prob.model.add_constraint("surface.surface_joining.section_separation", upper=0, lower=0)  # FOR SLSQP
# prob.model.add_constraint('surface.surface_joining.section_separation',equals=0.0,scaler=1e-4) #FOR OTHER OPTIMIZERS

# Add CL constraint
prob.model.add_constraint(point_name + ".CL", equals=0.3)

# Add Wing total area constraint
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

# docs checkpoint 7

# Get each section mesh
mesh1 = prob.get_val("surface.sec0.mesh", units="m")
mesh2 = prob.get_val("surface.sec1.mesh", units="m")

# Get the unified mesh
meshUni = prob.get_val(name + "." + unification_name + "." + name + "_uni_mesh")


# Plot the results
def plot_meshes(meshes):
    """this function plots to plot the mesh"""
    plt.figure(figsize=(8, 4))
    for i, mesh in enumerate(meshes):
        mesh_x = mesh[:, :, 0]
        mesh_y = mesh[:, :, 1]
        color = "k"
        for i in range(mesh_x.shape[0]):
            plt.plot(mesh_y[i, :], 1 - mesh_x[i, :], color, lw=1)
            plt.plot(-mesh_y[i, :], 1 - mesh_x[i, :], color, lw=1)  # plots the other side of symmetric wing
        for j in range(mesh_x.shape[1]):
            plt.plot(mesh_y[:, j], 1 - mesh_x[:, j], color, lw=1)
            plt.plot(-mesh_y[:, j], 1 - mesh_x[:, j], color, lw=1)  # plots the other side of symmetric wing
    plt.axis("equal")
    plt.xlabel("y (m)")
    plt.ylabel("x (m)")
    plt.savefig("opt_planform_constraint.png")


# plot_meshes([mesh1,mesh2])
plot_meshes([meshUni])
# plt.show()
# docs checkpoint 8
