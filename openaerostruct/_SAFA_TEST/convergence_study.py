import numpy as np

import openmdao.api as om

from openaerostruct.geometry.geometry_group import MultiSecGeometry
from openaerostruct.aerodynamics.aero_groups import AeroPoint
from openaerostruct.geometry.geometry_group import build_sections
from openaerostruct.geometry.geometry_unification import unify_mesh
import matplotlib.pyplot as plt


#Set-up B-splines for each section. Done here since this information will be needed multiple times.
sec_chord_cp = [np.array([1,1]),np.array([1.0,1.0])]

# Create a dictionary with info and options about the multi-section aerodynamic
# lifting surface
surface = {
    # Wing definition

    #Basic surface parameters
    "name":"surface",
    "num_sections": 2, #The number of sections in the multi-section surface
    "sec_name": ["sec0","sec1"],  # names of the individual sections
    "symmetry": True,  # if true, model one half of wing. reflected across the midspan of the root section
    "S_ref_type": "wetted",  # how we compute the wing area,
    # can be 'wetted' or 'projected'


    #Geometry Parameters
    "sec_taper": np.array([1.0,1.0]), #Wing taper for each section
    "sec_span":np.array([1.0,1.0]), #Wing span for each section
    "sec_sweep":np.array([0.0,0.0]), #Wing sweep for each section
    "sec_twist_cp":[np.array([0.0,0.0]), np.array([0.0,0.0])], 
    "sec_chord_cp": sec_chord_cp,
    #"sec_chord_cp": [np.ones(1),2*np.ones(1),3*np.ones(1)], #Chord B-spline control points for each section
    "root_chord" : 1.0, #Wing root chord for each section

    #Mesh Parameters
    "meshes": "gen-meshes", #Supply a mesh for each section or "gen-meshes" for automatic mesh generation
    "nx" : 2, #Number of chordwise points. Same for all sections
    "sec_ny" : np.array([2,2]), #Number of spanwise points for each section
    
    #Aerodynamic Parameters
    "sec_CL0": np.array([0.0,0.0]),  # CL of the surface at alpha=0
    "sec_CD0": np.array([0.015,0.015]),  # CD of the surface at alpha=0
    # Airfoil properties for viscous drag calculation
    "k_lam": 0.05,  # percentage of chord with laminar
    # flow, used for viscous drag
    #"sec_t_over_c_cp": [np.array([0.15]),np.array([0.15])],  # thickness over chord ratio (NACA0015)
    "sec_c_max_t": [0.303,0.303],  # chordwise location of maximum (NACA0015)
    # thickness
    "with_viscous": False,  # if true, compute viscous drag
    "with_wave": False,  # if true, compute wave drag
    "groundplane":False,
}


panel_counts = [2,6,11,16,21,26,31,36,41,46,51]
surfaces2 = []
import copy
for pc in panel_counts:
    surface["sec_ny"] = np.array([pc,pc])
    surfaces2.append(copy.deepcopy(surface))

CL = []
CD = []
for i,pc in enumerate(panel_counts):
    # Create the OpenMDAO problem
    prob = om.Problem()

    # Create an independent variable component that will supply the flow
    # conditions to the problem.
    indep_var_comp = om.IndepVarComp()
    indep_var_comp.add_output("v", val=1.0, units="m/s")
    indep_var_comp.add_output("alpha", val=5.0, units="deg")
    indep_var_comp.add_output("Mach_number", val=0.3)
    indep_var_comp.add_output("re", val=1.0e5, units="1/m")
    indep_var_comp.add_output("rho", val=0.38, units="kg/m**3")
    indep_var_comp.add_output("cg", val=np.zeros((3)), units="m")

    # Add this IndepVarComp to the problem model
    prob.model.add_subsystem("prob_vars", indep_var_comp, promotes=["*"])

    #Generate the sections and unified mesh here. It's needed to join the sections by construction.
    section_surfaces = build_sections(surfaces2[i])
    uniMesh = unify_mesh(section_surfaces)

    # Create and add a group that handles the geometry for the
    # aerodynamic lifting surface
    multi_geom_group = MultiSecGeometry(surface=surfaces2[i])
    prob.model.add_subsystem(surface["name"], multi_geom_group)

    # Create the aero point group, which contains the actual aerodynamic
    # analyses
    aero_group = AeroPoint(surfaces=section_surfaces,multiSection=True,unifiedMesh=uniMesh,MSSurfName='surface')
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


    # Set up and run the model
    prob.setup()
    #om.n2(prob)
    prob.run_model()
    CL.append(prob.get_val('aero_point_0.CL')[0])
    CD.append(prob.get_val('aero_point_0.CD')[0])


#Compute L2 relative error norms
CLerror = np.sqrt((np.diff(np.array(CL))/np.array(CL[:-1]))**2)
CDerror = np.sqrt((np.diff(np.array(CD))/np.array(CD[:-1]))**2)

dy = 1/(np.array(panel_counts[1:]))
ooA = []
for i in range(len(CLerror)-1):
    ooA.append(np.log(CLerror[i+1]/CLerror[i])/np.log(dy[i+1]/dy[i]))

plt.figure(1)
plt.plot(dy,CLerror,label='$C_L ||\\epsilon||_{L2}$',marker='o')
plt.plot(dy,CDerror,label='$C_D ||\\epsilon||_{L2}$',marker='o')
plt.yscale('log')
plt.xscale('log')
plt.grid()
plt.legend()
plt.xlabel('log(dy) (m)')
plt.ylabel('$log(||\\epsilon||_{L2})$')
plt.savefig('convergence_study.pdf')
plt.show()


print(ooA)

