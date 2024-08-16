import numpy as np

import openmdao.api as om

from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.geometry.geometry_group import MultiSecGeometry
from openaerostruct.aerodynamics.aero_groups import AeroPoint
import openaerostruct.geometry.geometry_multi_sec_mesh as multiMesh
from openaerostruct.geometry.geometry_group import build_sections
from openaerostruct.geometry.geometry_unification import unify_mesh
import matplotlib.pyplot as plt
import niceplots

plt.style.use(niceplots.get_style("doumont-light"))
# Create a dictionary with info and options about the multi-section aerodynamic
# lifting surface
surface = {
    # Wing definition

    #Basic surface parameters
    "name":"surface",
    "num_sections": 3, #The number of sections in the multi-section surface
    "sec_name": ["sec0","sec1","sec2"],  # names of the individual sections
    "symmetry": True,  # if true, model one half of wing. reflected across the midspan of the root section
    "S_ref_type": "wetted",  # how we compute the wing area,
    # can be 'wetted' or 'projected'


    #Geometry Parameters
    "sec_taper": np.array([1.0,1.0,1.0]), #Wing taper for each section
    "sec_span":np.array([1.0,1.0,1.0]), #Wing span for each section
    "sec_sweep":np.array([0.0,0.0,0.0]), #Wing sweep for each section
    "sec_chord_cp": [np.array([0.5,0.1]),np.array([1,1]),np.array([1,0.5])],
    "sec_twist_cp": [np.zeros(2),np.zeros(2),np.zeros(2)],
    #"sec_chord_cp": [np.ones(1),2*np.ones(1),3*np.ones(1)], #Chord B-spline control points for each section
    "root_chord" : 1.0, #Wing root chord for each section

    #Mesh Parameters
    "meshes": "gen-meshes", #Supply a mesh for each section or "gen-meshes" for automatic mesh generation
    "nx" : 5, #Number of chordwise points. Same for all sections
    "sec_ny" : np.array([21,21,21]), #Number of spanwise points for each section
    
    #Aerodynamic Parameters
    "sec_CL0": np.array([0.0,0.0,0.0]),  # CL of the surface at alpha=0
    "sec_CD0": np.array([0.015,0.015,0.015]),  # CD of the surface at alpha=0
    # Airfoil properties for viscous drag calculation
    "k_lam": 0.05,  # percentage of chord with laminar
    # flow, used for viscous drag
    "sec_t_over_c_cp": [np.array([0.15]),np.array([0.15]),np.array([0.15])],  # thickness over chord ratio (NACA0015)
    "sec_c_max_t": [0.303,0.303,0.303],  # chordwise location of maximum (NACA0015)
    # thickness
    "with_viscous": False,  # if true, compute viscous drag
    "with_wave": False,  # if true, compute wave drag
    "groundplane":False,
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
multi_geom_group = MultiSecGeometry(surface=surface,joining_comp=True,dim_constr=[np.array([1,0,0]),np.array([1,0,0]),np.array([1,0,0])])
prob.model.add_subsystem(surface["name"], multi_geom_group)

#Generate the sections and unified mesh here in addition to adding the components. 
#This has to ALSO be done here since AeroPoint has to know the unified mesh size.
section_surfaces = build_sections(surface)
uniMesh = unify_mesh(section_surfaces)

# Create the aero point group, which contains the actual aerodynamic
# analyses
aero_group = AeroPoint(surfaces=section_surfaces,multiSection=True,unifiedMesh=uniMesh,msSurfName=surface["name"])
point_name = "aero_point_0"
prob.model.add_subsystem(
    point_name, aero_group, promotes_inputs=["v", "alpha", "Mach_number", "re", "rho", "cg"]
)

#Get name of surface and construct unified mesh name
name = surface["name"]
unification_name = '{}_unification'.format(surface["name"])

# Connect the mesh from the mesh unification component to the analysis point
prob.model.connect(name + "." + unification_name + "." + name + "_uni_mesh", point_name + "." + "surface" + ".def_mesh")

# Perform the connections with the modified names within the
# 'aero_states' group.
prob.model.connect(name + "." + unification_name + "." + name + "_uni_mesh", point_name + ".aero_states." + "surface" + "_def_mesh")


prob.model.add_design_var("surface.sec0.chord_cp", lower=0.1, upper=10.0, units=None)
prob.model.add_design_var("surface.sec1.chord_cp", lower=0.1, upper=10.0, units=None)
#prob.model.add_design_var("surface.sec2.chord_cp", lower=0.1, upper=10.0, units=None)
prob.model.add_design_var("alpha", lower=0.0, upper=10.0, units='deg')

#Add joined mesh constraint
#prob.model.add_constraint('surface.surface_joining.section_separation',upper=0,lower=0)
prob.model.add_constraint('surface.surface_joining.section_separation',equals=0.0)

#Add CL constraint
prob.model.add_constraint(point_name +'.CL',equals=0.3)

#Add Area constraint
prob.model.add_constraint(point_name + '.total_perf.S_ref_total',equals=2.0)

#Add objective
prob.model.add_objective(point_name + ".CD", scaler=1e4)

'''
prob.driver = om.ScipyOptimizeDriver()
prob.driver.options['optimizer'] = 'SLSQP'
prob.driver.options['tol'] = 1e-3
prob.driver.options['disp'] = True
prob.driver.options['maxiter'] = 1000
prob.driver.options["debug_print"] = ["nl_cons", "objs", "desvars"]
'''

prob.driver = om.pyOptSparseDriver()
prob.driver.options['optimizer'] = 'SNOPT'
prob.driver.opt_settings['Major feasibility tolerance'] = 1e-3
#prob.driver.opt_settings['ACC'] = 1e-3
#prob.driver.options['disp'] = True
#prob.driver.opt_settings['MAXIT'] = 1000
prob.driver.opt_settings['Major iterations limit'] = 1000
prob.driver.options["debug_print"] = ["nl_cons", "objs", "desvars"]

# Set up and run the optimization problem
prob.setup()
#prob.run_model()

prob.run_driver()
#om.n2(prob)


mesh1 = prob.get_val("surface.sec0.mesh", units="m")
mesh2 = prob.get_val("surface.sec1.mesh", units="m")
mesh3 = prob.get_val("surface.sec2.mesh", units="m")

meshUni = prob.get_val(name + "." + unification_name + "." + name + "_uni_mesh")

def plot_meshes(meshes):
    """ this function plots to plot the mesh """
    plt.figure(figsize=(8, 4))
    for i,mesh in enumerate(meshes):
        mesh_x = mesh[:, :, 0]
        mesh_y = mesh[:, :, 1]
        color = 'k'
        for i in range(mesh_x.shape[0]):
            plt.plot(mesh_y[i, :], mesh_x[i, :], color, lw=1)
            plt.plot(-mesh_y[i, :], mesh_x[i, :], color, lw=1)   # plots the other side of symmetric wing
        for j in range(mesh_x.shape[1]):
            plt.plot(mesh_y[:, j], mesh_x[:, j], color, lw=1)
            plt.plot(-mesh_y[:, j], mesh_x[:, j], color, lw=1)   # plots the other side of symmetric wing
    plt.axis('equal')
    plt.xlabel('y (m)')
    plt.ylabel('x (m)')
    plt.savefig('figure2.png')
    plt.savefig('figure2.pdf')
    #plt.legend()
  

plot_meshes([mesh1,mesh2,mesh3])



# ---------------------------------
# plot spanwise lift distribution
# ---------------------------------
# get spanwise coordinate, chord length, and sectional lift coefficient
mesh = prob.get_val(point_name + ".surface.def_mesh", units="m")
y_vertices = mesh[0, :, 1]   # spanwise coordinate of VLM panel vertices
y_center = 0.5 * (y_vertices[:-1] + y_vertices[1:])   # spanwise coordinate of VLM panel centers
Cl = prob.get_val(point_name + ".surface_perf.Cl")   # sectional lift coefficient
chord_edge = prob.get_val(point_name + ".surface.chords", units="m")   # chord length at panel edges
chord_center = 0.5 * (chord_edge[:-1] + chord_edge[1:])   # chord length at panel centers

# concatenate the other side of symmetric wing
x_center = np.concatenate((y_center, -y_center[::-1]))
chord_center = np.concatenate((chord_center, chord_center[::-1]))
Cl = np.concatenate((Cl, Cl[::-1]))

# compute reference elliptical lift distribution
CL = prob.get_val(point_name + ".surface_perf.CL")[0]   # CL of the total wing
Sref = prob.get_val(point_name + ".surface_perf.S_ref", units="m**2")[0]   # wing area
semi_span = mesh[0, -1, 1] - mesh[0, 0, 1]   # wing semi-span
cl_chord_max = 2 * CL * Sref / (np.pi * semi_span)   # Cl * chord at the center of elliptical lift distribution (based on the area of ellipse = 2 * CL * S)
# smooth line for elliptical lift distribution
y_ellipse = np.linspace(-semi_span, semi_span, 100)
Cl_chord_ellipse = cl_chord_max * np.sqrt(1 - (y_ellipse / semi_span) ** 2)

# plot lift distribution
plt.figure(figsize=(8, 4))
plt.plot(x_center, Cl * chord_center, color="C0", lw=2,label="Optimized Wing")
plt.plot(y_ellipse, Cl_chord_ellipse, color="C1", lw=1,label="Elliptical")
plt.xlabel("y (m)")
plt.ylabel("$C_{l}c(y)$")
#plt.grid()
plt.legend()

plt.savefig("lift_distribution_optimized.png", bbox_inches="tight")


#plot_meshes([meshUni])
plt.show()