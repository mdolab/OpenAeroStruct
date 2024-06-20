import numpy as np

import openmdao.api as om

from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.geometry.geometry_group import MultiSecGeometry
from openaerostruct.aerodynamics.aero_groups import AeroPoint
import openaerostruct.geometry.geometry_multi_sec_mesh as multiMesh
import matplotlib.pyplot as plt

# Create a dictionary to store options about the mesh
mesh_dict = {"num_y": 7, "num_x": 2, "wing_type": "rect", "symmetry": True}

# Generate the aerodynamic mesh based on the previous dictionary
mesh, twist_cp = generate_mesh(mesh_dict)


# Create a dictionary with info and options about the aerodynamic
# lifting surface
surface = {
    # Wing definition
    "name":"surface",
    "num_sections": 2,
    "sec_name": ["sec1","sec2"],  # name of the surface
    "symmetry": True,  # if true, model one half of wing
    # reflected across the plane y = 0
    "S_ref_type": "wetted",  # how we compute the wing area,
    # can be 'wetted' or 'projected'
    "sec_taper": np.array([1.0,1.0]),
    "sec_span":np.array([1.0,1.0]),
    "sec_sweep":np.array([0.0,0.0]),
    #"sec_chord_cp": [np.array([1,0.5,0.5]),np.array([1,0.2,0.5])],
    "sec_chord_cp": [np.ones(1),np.ones(1)],
    "meshes": "gen-meshes",
    "sec_CL0": np.array([0.0,0.0]),  # CL of the surface at alpha=0
    "sec_CD0": np.array([0.015,0.015]),  # CD of the surface at alpha=0
    "nx" : 2,
    "sec_ny" : np.array([2,2]),
    "sec_root_chord" : np.array([1.0,1.0]),
    # Airfoil properties for viscous drag calculation
    "k_lam": 0.05,  # percentage of chord with laminar
    # flow, used for viscous drag
    "sec_t_over_c_cp": [np.array([0.15]),np.array([0.15])],  # thickness over chord ratio (NACA0015)
    "sec_c_max_t": [0.303,0.303],  # chordwise location of maximum (NACA0015)
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
multi_geom_group = MultiSecGeometry(surface=surface)
prob.model.add_subsystem(surface["name"], multi_geom_group)

#Get Surfaces from Multi-section geometry group

def build_aero_sections(surface):
    #Get number of sections
    num_sections = surface["num_sections"]

    if surface["meshes"] == "gen-meshes":
        #Verify that all required inputs for automatic mesh generation are provided for each section
        if len(surface["sec_ny"]) != num_sections:
            raise ValueError("Number of spanwise points needs to be provided for each section")
        if len(surface["sec_taper"]) != num_sections:
            raise ValueError("Taper needs to be provided for each section")
        if len(surface["sec_root_chord"]) != num_sections:
            raise ValueError("Root chord length needs to be provided for each section")
        if len(surface["sec_span"]) != num_sections:
            raise ValueError("Span needs to be provided for each section")
        if len(surface["sec_sweep"]) != num_sections:
            raise ValueError("Sweep needs to be provided for each section")
        
        #Get required data for section mesh generation
        nx = surface["nx"]
        sec_ny = surface["sec_ny"]
        sec_taper = surface["sec_taper"]
        sec_span = surface["sec_span"]
        sec_root_chord = surface["sec_root_chord"]
        sec_sweep = surface["sec_sweep"]

        #Compute section aspect ratio
        sec_S = sec_span*(sec_root_chord*(np.ones(num_sections)+sec_taper))/2
        sec_AR = sec_span**2/sec_S

        #Create data array for mesh generator
        sectionData = np.hstack([np.ones(num_sections)[:,np.newaxis],sec_root_chord[:,np.newaxis],sec_AR[:,np.newaxis],np.zeros(num_sections)[:,np.newaxis]])

        symmetry  = surface["symmetry"]

        #Generate unified and individual section meshes
        mesh, sec_meshes = multiMesh.generateMesh(num_sections,sectionData,sec_ny-np.ones(num_sections,dtype=np.int32),nx-1,symmetry)
        
    else:
        #Allow user to provide mesh for each section
        if len(surface["meshes"]) != num_sections:
            raise ValueError("A mesh needs to be provided for each section.")
        sec_meshes = surface["meshes"]

    if len(surface["sec_name"]) != num_sections:
            raise ValueError("A name needs to be provided for each section.")
    section_surfaces = []
    for i in range(num_sections):
        section = {
            "name": surface["sec_name"][i],  # name of the surface
            "symmetry": surface["symmetry"],  
            "S_ref_type": surface["S_ref_type"], 
            "mesh": sec_meshes[i],
            "span":surface["sec_span"][i],
            "taper":surface["sec_taper"][i],
            "sweep":surface["sec_sweep"][i],
            "chord_cp": surface["sec_chord_cp"][i],  
            "ref_axis_pos": 0.25, 
            "CL0": surface["sec_CL0"][i], 
            "CD0": surface["sec_CD0"][i], 
            "k_lam": surface["k_lam"], 
            #"t_over_c_cp": surface["sec_t_over_c_cp"][i], 
            "c_max_t": surface["sec_c_max_t"][i],  
            "with_viscous": surface["with_viscous"], 
            "with_wave": surface["with_wave"],
            "groundplane": surface["groundplane"],
        }  # end of surface dictionary

        section_surfaces.append(section)
    return section_surfaces


section_surfaces = build_aero_sections(surface)


# Create the aero point group, which contains the actual aerodynamic
# analyses
aero_group = AeroPoint(surfaces=section_surfaces)
point_name = "aero_point_0"
prob.model.add_subsystem(
    point_name, aero_group, promotes_inputs=["v", "alpha", "Mach_number", "re", "rho", "cg"]
)


#name = surface["name"]

# Connect the mesh from the geometry component to the analysis point
prob.model.connect("surface.sec1" + ".mesh", point_name + "." + "sec1" + ".def_mesh")
prob.model.connect("surface.sec2" + ".mesh", point_name + "." + "sec2" + ".def_mesh")

# Perform the connections with the modified names within the
# 'aero_states' group.
prob.model.connect("surface.sec1" +  ".mesh", point_name + ".aero_states." + "sec1"  + "_def_mesh")
prob.model.connect("surface.sec2" +  ".mesh", point_name + ".aero_states." + "sec2"  + "_def_mesh")

#prob.model.connect("surface.sec2" + ".t_over_c", point_name + "." + name + "_perf." + "t_over_c")
'''
# Import the Scipy Optimizer and set the driver of the problem to use
# it, which defaults to an SLSQP optimization method
prob.driver = om.ScipyOptimizeDriver()
prob.driver.options["tol"] = 1e-9

# Setup problem and add design variables, constraint, and objective
prob.model.add_design_var("wing.twist_cp", lower=-10.0, upper=15.0)
prob.model.add_design_var("wing.chord_cp", lower=0.5, upper=1.5)
prob.model.add_design_var("wing.xshear_cp", lower=0.0, upper=1.0)
prob.model.add_design_var("wing.yshear_cp", lower=-1.0, upper=0.0)
prob.model.add_design_var("wing.zshear_cp", lower=-1.0, upper=1.0)
prob.model.add_constraint(point_name + ".wing_perf.CL", equals=0.5)
prob.model.add_objective(point_name + ".wing_perf.CD", scaler=1e4)
'''


#prob.model.add_design_var("surface.sec1.chord_cp", lower=0.1, upper=10.0, units=None)
prob.model.add_design_var("surface.sec2.chord_cp", lower=0.1, upper=10.0, units=None)
prob.model.add_design_var("alpha", lower=0.1, upper=10.0, units='deg')

testConnect = om.ExecComp('dist = sum((edge1 - edge2)**2)**(1/2)',edge1={'units':'m','shape':(2,3)},edge2={'units':'m','shape':(2,3)},dist={'units':'m'})
prob.model.add_subsystem('testConnect',testConnect)

#testObj = om.ExecComp('testOut = testIn',testIn={'units':'deg'})
#prob.model.add_subsystem('testObj',testObj)

#prob.model.connect('alpha','testObj.testIn')

prob.model.connect("surface.sec1.mesh",'testConnect.edge1',src_indices = om.slicer[:,-1,:])
prob.model.connect("surface.sec2.mesh",'testConnect.edge2',src_indices = om.slicer[:,-1,:])

prob.model.add_constraint('testConnect.dist',equals=0.0)

prob.model.add_objective(point_name + ".CD", scaler=1e4)


prob.driver = om.ScipyOptimizeDriver()
prob.driver.options['optimizer'] = 'SLSQP'
prob.driver.options['tol'] = 1e-3
prob.driver.options['disp'] = True
prob.driver.options['maxiter'] = 1000
#prob.driver.options["debug_print"] = ["nl_cons", "objs", "desvars"]

# Set up and run the optimization problem
prob.setup()

prob.run_model()

#prob.run_driver()
om.n2(prob)


mesh1 = prob.get_val("surface.sec1.mesh", units="m")
mesh2 = prob.get_val("surface.sec2.mesh", units="m")

def plot_meshes(meshes):
    """ this function plots to plot the mesh """
    plt.figure(figsize=(6, 3))
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
    plt.xlabel('span, m')
    plt.ylabel('chord, m')
    #plt.legend()
    plt.show()

plot_meshes([mesh1,mesh2])
#

#assert_near_equal(prob["aero_point_0.wing_perf.CD"][0], 0.02891508386825118, 1e-6)
#assert_near_equal(prob["aero_point_0.wing_perf.CL"][0], 0.5, 1e-6)
#assert_near_equal(prob["aero_point_0.CM"][1], -4.281635590978787, 1e-6)