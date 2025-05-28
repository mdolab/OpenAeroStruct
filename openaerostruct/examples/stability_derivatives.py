"""
Example of aerodynamic analysis including stability derivatives (CL_alpha and CM_alpha).

We compute CL_alpha and CM_alpha by finite differencing with respect to alpha.
To do so, we instantiate AeroPoint at alpha and alpha + delta_alpha, and add an ExecComp to compute finite difference.
Finite differencing w.r.t. alpha allows us to compute the derivatives of CL_alpha and CM_alpha w.r.t. design variables.

Note that this example does not trim the aircraft (e.g. CM != 0).
"""

import numpy as np
import openmdao.api as om

from openaerostruct.meshing.mesh_generator import generate_mesh
from openaerostruct.geometry.geometry_group import Geometry
from openaerostruct.aerodynamics.aero_groups import AeroPoint

# Create a dictionary to store options about the surface.
# Here, we setup a simple rectangular mesh with chord = 1 m and span = 10 m.
# We then apply sweep after prob.setup()
mesh_dict = {
    "num_y": 15,
    "num_x": 3,
    "wing_type": "rect",
    "root_chord": 1.0,  # m
    "span": 10.0,  # m
    "symmetry": True,
    "num_twist_cp": 2,
    "span_cos_spacing": 0.0,
}

mesh = generate_mesh(mesh_dict)

surf_dict = {
    # Wing definition
    "name": "wing",  # name of the surface
    "symmetry": True,  # if true, model one half of wing
    # reflected across the plane y = 0
    "S_ref_type": "wetted",  # how we compute the wing area,
    # can be 'wetted' or 'projected'
    "mesh": mesh,
    "twist_cp": np.zeros(mesh_dict["num_twist_cp"]),  # twist control points
    "sweep": 0.0,  # wing sweep angle
    # Aerodynamic performance of the lifting surface at
    # an angle of attack of 0 (alpha=0).
    # These CL0 and CD0 values are added to the CL and CD
    # obtained from aerodynamic analysis of the surface to get
    # the total CL and CD.
    # These CL0 and CD0 values do not vary wrt alpha.
    "CL0": 0.0,  # CL of the surface at alpha=0
    "CD0": 0.015,  # CD of the surface at alpha=0
    # Airfoil properties for viscous drag calculation
    "k_lam": 0.05,  # percentage of chord with laminar
    # flow, used for viscous drag
    "t_over_c_cp": np.array([0.15]),  # thickness over chord ratio (NACA0015)
    "c_max_t": 0.303,  # chordwise location of maximum (NACA0015)
    # thickness
    "with_viscous": True,  # if true, compute viscous drag
    "with_wave": False,  # if true, compute wave drag
}

surfaces = [surf_dict]

# Although this example only considers cruise, hence single point,
# We still need to AeroPoint instances to finite difference CL and CM w.r.t. alpha
n_points = 2

# Create the problem and the model group
prob = om.Problem()

indep_var_comp = om.IndepVarComp()
indep_var_comp.add_output("v", val=248.136, units="m/s")
indep_var_comp.add_output("alpha", val=5.0, units="deg")
indep_var_comp.add_output("Mach_number", val=0.84)
indep_var_comp.add_output("re", val=1.0e6, units="1/m")
indep_var_comp.add_output("rho", val=0.38, units="kg/m**3")
# Set center of gravity to be 0.5 m from the leading edge.
# Note that the CG location does *not* move as a result of wing shape changes
indep_var_comp.add_output("cg", val=np.array([0.5, 0.0, 0.0]), units="m")

prob.model.add_subsystem("prob_vars", indep_var_comp, promotes=["*"])

# Compute alpha perturbation for finite difference
alpha_FD_stepsize = 1e-4  # deg
alpha_perturb_comp = om.ExecComp(
    "alpha_plus_delta = alpha + delta_alpha",
    units="deg",
    delta_alpha={"val": alpha_FD_stepsize, "constant": True},
)
prob.model.add_subsystem("alpha_for_FD", alpha_perturb_comp, promotes=["*"])

# Loop over each surface and create the geometry groups
for surface in surfaces:
    # Get the surface name and create a group to contain components only for this surface.
    name = surface["name"]
    geom_group = Geometry(surface=surface)

    # Add geom_group to the problem with the name of the surface.
    prob.model.add_subsystem(name + "_geom", geom_group)

# Create aero analysis point for alpha and alpha_plus_delta
point_names = ["aero_point", "aero_point_FD"]
for i in range(n_points):
    # Create the aero point group and add it to the model
    aero_group = AeroPoint(surfaces=surfaces)
    point_name = point_names[i]
    prob.model.add_subsystem(point_name, aero_group)

    # Connect flow properties to the analysis point
    prob.model.connect("v", point_name + ".v")
    prob.model.connect("Mach_number", point_name + ".Mach_number")
    prob.model.connect("re", point_name + ".re")
    prob.model.connect("rho", point_name + ".rho")
    prob.model.connect("cg", point_name + ".cg")

    # Connect angle of attack. Use perturbed alpha for the second point for finite difference
    alpha_name = "alpha" if i == 0 else "alpha_plus_delta"
    prob.model.connect(alpha_name, point_name + ".alpha")

    # Connect the parameters within the model for each aero point
    for surface in surfaces:
        name = surface["name"]

        # Connect the mesh from the geometry component to the analysis point
        prob.model.connect(name + "_geom.mesh", point_name + "." + name + ".def_mesh")

        # Perform the connections with the modified names within the 'aero_states' group.
        prob.model.connect(name + "_geom.mesh", point_name + ".aero_states." + name + "_def_mesh")
        prob.model.connect(name + "_geom.t_over_c", point_name + "." + name + "_perf." + "t_over_c")

# Compute stability derivatives by finite difference
stabibility_derivatives_comp = om.ExecComp(
    ["CL_alpha = (CL_FD - CL) / delta_alpha", "CM_alpha = (CM_FD - CM) / delta_alpha"],
    delta_alpha={"val": alpha_FD_stepsize, "constant": True},
    CL_alpha={"val": 0.0, "units": "1/deg"},
    CL_FD={"val": 0.0, "units": None},
    CL={"val": 0.0, "units": None},
    CM_alpha={"val": np.zeros(3), "units": "1/deg"},
    CM_FD={"val": np.zeros(3), "units": None},
    CM={"val": np.zeros(3), "units": None},
)
prob.model.add_subsystem("stability_derivs", stabibility_derivatives_comp, promotes_outputs=["*"])
# Connect CL and CM from aero points
prob.model.connect("aero_point.CL", "stability_derivs.CL")
prob.model.connect("aero_point.CM", "stability_derivs.CM")
prob.model.connect("aero_point_FD.CL", "stability_derivs.CL_FD")
prob.model.connect("aero_point_FD.CM", "stability_derivs.CM_FD")

# Compute static margin
static_margin_comp = om.ExecComp(
    "static_margin = -CM_alpha / CL_alpha",
    CM_alpha={"val": 0.0, "units": "1/deg"},  # for pitching moment
    CL_alpha={"val": 0.0, "units": "1/deg"},
    static_margin={"val": 0.0, "units": None},  # static margin
)
prob.model.add_subsystem("static_margin", static_margin_comp, promotes_outputs=["*"])
# Connect stability derivatives to static margin component
prob.model.connect("CL_alpha", "static_margin.CL_alpha")
prob.model.connect("CM_alpha", "static_margin.CM_alpha", src_indices=1)  # connect pitching moment deriv

# Set up the problem
prob.setup()

# Set sweep angle
prob.set_val("wing_geom.sweep", 10.0, units="deg")

# Run model and print outputs
prob.run_model()

print("Sweep angle   =", prob.get_val("wing_geom.sweep", units="deg"), "deg")
print("CL            =", prob.get_val("aero_point.CL"))
print("CD            =", prob.get_val("aero_point.CD"))
print("CM            =", prob.get_val("aero_point.CM"))
print("CL_alpha      =", prob.get_val("CL_alpha", units="1/deg"), "1/deg")
print("CM_alpha      =", prob.get_val("CM_alpha", units="1/deg"), "1/deg")
print("Static Margin =", prob.get_val("static_margin"))
