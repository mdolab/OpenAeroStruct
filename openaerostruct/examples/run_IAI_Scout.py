# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 09:46:02 2021

@author: rcharayr
"""

################################################################################
# Script Description

################################################################################

import numpy as np

from openaerostruct.geometry.utils import generate_mesh

from openaerostruct.integration.aerostruct_groups import AerostructGeometry, AerostructPoint

import openmdao.api as om
from openaerostruct.utils.constants import grav_constant


###############################################################################

# Create a mesh dictionary for the wing  to feed to generate_mesh to actually create
# the mesh array of the wing.
mesh_dict_wing = {'num_y' : 7,
             'num_x' : 3,
             'span' : 3.60,
             'wing_type' : 'rect',
             'symmetry' : False,
             'root_chord' : 0.45,
             }

mesh_wing = generate_mesh(mesh_dict_wing)
twist_cp_wing = np.zeros(5)
twist_cp_wing[2] = 5
thickness_cp_wing = 0.3*np.ones(5)


#création d'un dictionnaires qui contient les options du maillage de l'empennage:
mesh_dict_tail = {'num_y' : 3,
             'num_x' : 2,
             'span' : 1.23,
             'wing_type' : 'rect',
             'symmetry' : False,
             'root_chord' : 0.27,
             'offset' : np.array([2.28,0,0]),
             }

# génération d'un maillage de l'empennage basé sur le dictionnaire grâce à la fonction generate_mesh() :
mesh_tail = generate_mesh(mesh_dict_tail)
twist_cp_tail = np.zeros(3)
thickness_cp_tail = 0.3*np.ones(3)

surf_dict_wing = {
            # Wing definition
            'name' : 'wing',        # name of the surface
            'symmetry' : False,     # if true, model one half of wing
                                    # reflected across the plane y = 0
            'S_ref_type' : 'wetted', # how we compute the wing area,
                                     # can be 'wetted' or 'projected'
            'fem_model_type' : 'tube',

            'thickness_cp' : thickness_cp_wing,
            'twist_cp' : twist_cp_wing,
            'mesh' : mesh_wing,

            # Aerodynamic performance of the lifting surface at
            # an angle of attack of 0 (alpha=0).
            # These CL0 and CD0 values are added to the CL and CD
            # obtained from aerodynamic analysis of the surface to get
            # the total CL and CD.
            # These CL0 and CD0 values do not vary wrt alpha.
            'CL0' : 0.0,            # CL of the surface at alpha=0
            'CD0' : 0.015,            # CD of the surface at alpha=0

            # Airfoil properties for viscous drag calculation
            'k_lam' : 0.05,         # percentage of chord with laminar
                                    # flow, used for viscous drag
            't_over_c_cp' : np.array([0.15]),      # thickness over chord ratio (NACA0015)
            'c_max_t' : .303,       # chordwise location of maximum (NACA0015)
                                    # thickness
            'with_viscous' : True,  # if true, compute viscous drag
            'with_wave' : False,     # if true, compute wave drag
            
            # Structural values are based on aluminum 7075
            'E' : 70.e9,            # [Pa] Young's modulus of the spar
            'G' : 30.e9,            # [Pa] shear modulus of the spar
            'yield' : 500.e6 / 2.5, # [Pa] yield stress divided by 2.5 for limiting case
            'mrho' : 3.e3,          # [kg/m^3] material density
            'fem_origin' : 0.35,    # normalized chordwise location of the spar
            'wing_weight_ratio' : 1., # multiplicative factor on the computed structural weight
            'struct_weight_relief' : False,    # True to add the weight of the structure to the loads on the structure
            'distributed_fuel_weight' : False,
            
            # Constraints
            'exact_failure_constraint' : False, # if false, use KS function

            }



surf_dict_tail = {
            # tail definition
            'name' : 'tail',        # name of the surface
            'symmetry' : False,     # if true, model one half of wing
                                    # reflected across the plane y = 0
            'S_ref_type' : 'wetted', # how we compute the wing area,
                                     # can be 'wetted' or 'projected'
            'fem_model_type' : 'tube',

            'thickness_cp' : thickness_cp_tail,
            'twist_cp' : twist_cp_tail,
            'mesh' : mesh_tail,

            # Aerodynamic performance of the lifting surface at
            # an angle of attack of 0 (alpha=0).
            # These CL0 and CD0 values are added to the CL and CD
            # obtained from aerodynamic analysis of the surface to get
            # the total CL and CD.
            # These CL0 and CD0 values do not vary wrt alpha.
            'CL0' : 0.0,            # CL of the surface at alpha=0
            'CD0' : 0.005,            # CD of the surface at alpha=0

            'fem_origin' : 0.35,

            # Airfoil properties for viscous drag calculation
            'k_lam' : 0.05,         # percentage of chord with laminar
                                    # flow, used for viscous drag
            't_over_c_cp' : np.array([0.15]),      # thickness over chord ratio (NACA0015)
            'c_max_t' : .303,       # chordwise location of maximum (NACA0015)
                                    # thickness
            'with_viscous' : True,  # if true, compute viscous drag
            'with_wave' : False,     # if true, compute wave drag
            
            # Structural values are based on aluminum 7075
            'E' : 70.e9,            # [Pa] Young's modulus of the spar
            'G' : 30.e9,            # [Pa] shear modulus of the spar
            'yield' : 500.e6 / 2.5, # [Pa] yield stress divided by 2.5 for limiting case
            'mrho' : 3.e3,          # [kg/m^3] material density
            'fem_origin' : 0.35,    # normalized chordwise location of the spar
            'wing_weight_ratio' : 1., # multiplicative factor on the computed structural weight
            'struct_weight_relief' : False,    # True to add the weight of the structure to the loads on the structure
            'distributed_fuel_weight' : False,
            
            # Constraints
            'exact_failure_constraint' : False, # if false, use KS function
            
            }

surfaces = [surf_dict_wing, surf_dict_tail]

# Create the problem and the model group
prob = om.Problem()

# Add problem information as an independent variables component
indep_var_comp = om.IndepVarComp()
indep_var_comp.add_output('v', val=22.876, units='m/s')
indep_var_comp.add_output('alpha', val=5., units='deg')
indep_var_comp.add_output('Mach_number', val=0.071)
indep_var_comp.add_output('re', val=1.e6, units='1/m')
indep_var_comp.add_output('rho', val=0.770816, units='kg/m**3')
indep_var_comp.add_output('CT', val=grav_constant * 8.6e-6, units='1/s')
indep_var_comp.add_output('R', val=1800e3, units='m')
indep_var_comp.add_output('W0', val=10.,  units='kg')
indep_var_comp.add_output('speed_of_sound', val=322.2, units='m/s')
indep_var_comp.add_output('load_factor', val=1.)
indep_var_comp.add_output('empty_cg', val=np.array([0.2, 0., 0.]), units='m')

prob.model.add_subsystem('prob_vars',
     indep_var_comp,
     promotes=['*'])

# Loop over each surface in the surfaces list
for surface in surfaces:
    aerostruct_group = AerostructGeometry(surface=surface)

    # Add tmp_group to the problem as the name of the surface.
    # Note that is a group and performance group for each
    # individual surface.
    prob.model.add_subsystem(surface['name'], aerostruct_group)

# Loop through and add a certain number of aero points
for i in range(1):
    # Create the aero struct point group and add it to the model
    AS_point = AerostructPoint(surfaces=surfaces)
    point_name = 'AS_point_{}'.format(i)
    prob.model.add_subsystem(point_name, AS_point, promotes_inputs=['v', 'alpha', 'Mach_number', 're', 'rho', 'CT', 'R', 'W0', 'speed_of_sound', 'empty_cg', 'load_factor'])

    # Connect the parameters within the model for each aero struct point
    for surface in surfaces:
        name = surface['name'] 
        
        # Connect aerodynamic mesh to coupled group mesh
        prob.model.connect(name + '.mesh', point_name + '.coupled.' + name + '.mesh')
        
        # Perform the connections with the modified names within the 'aero_states' group.
        # prob.model.connect(name + '.mesh', point_name + '.aero_states.' + name + '_def_mesh')
        # prob.model.connect(name + '.mesh', point_name + '.coupled.aero_states.' + name + '_def_mesh')

        # Issue quite a few connections within the model to make sure all of the
        # parameters are connected correctly.
        com_name = point_name + '.' + name + '_perf'
        prob.model.connect(name + '.local_stiff_transformed', point_name + '.coupled.' + name + '.local_stiff_transformed')
        prob.model.connect(name + '.nodes', point_name + '.coupled.' + name + '.nodes')
        
        # Connect performance calculation variables
        prob.model.connect(name + '.radius', com_name + '.radius')
        prob.model.connect(name + '.thickness', com_name + '.thickness')
        prob.model.connect(name + '.nodes', com_name + '.nodes')
        prob.model.connect(name + '.cg_location', point_name + '.' + 'total_perf.' + name + '_cg_location')
        prob.model.connect(name + '.structural_mass', point_name + '.' + 'total_perf.' + name + '_structural_mass')
        prob.model.connect(name + '.t_over_c', com_name + '.t_over_c')

        
prob.driver = om.ScipyOptimizeDriver()
# prob.driver.options['tol'] = 1e-9
prob.driver.options['tol'] = 1e-7

recorder = om.SqliteRecorder("aerostruct_IAI_Scout.db")
prob.driver.add_recorder(recorder)
prob.driver.recording_options['record_derivatives'] = True
prob.driver.recording_options['includes'] = ['*']


# Here we're varying twist, thickness and alpha.
prob.model.add_design_var('wing.twist_cp', lower=-5., upper=10.)
prob.model.add_design_var('wing.thickness_cp', lower=0.001, upper=0.01, scaler=1e3)
prob.model.add_design_var('tail.twist_cp', lower=-5., upper=10.)
prob.model.add_design_var('tail.thickness_cp', lower=0.001, upper=0.01, scaler=1e3)
prob.model.add_design_var('alpha', lower=-10., upper=10.)

# Make sure the spar doesn't fail, we meet the lift needs, and the aircraft
# is trimmed through CM=0.
prob.model.add_constraint('AS_point_0.wing_perf.failure', upper=0.)
prob.model.add_constraint('AS_point_0.wing_perf.thickness_intersects', upper=0.)
prob.model.add_constraint('AS_point_0.tail_perf.failure', upper=0.)
prob.model.add_constraint('AS_point_0.tail_perf.thickness_intersects', upper=0.)
prob.model.add_constraint('AS_point_0.L_equals_W', equals=0.)

# Instead of using an equality constraint here, we have to give it a little
# wiggle room to make SLSQP work correctly.
prob.model.add_constraint('AS_point_0.CM', lower=-0.001, upper=0.001)
prob.model.add_constraint('wing.twist_cp', lower=np.array([-1e20, -1e20, 5., -1e20, -1e20]), upper=np.array([1e20, 1e20, 5., 1e20, 1e20]))
prob.model.add_constraint('tail.twist_cp', lower=np.array([-1e20, -1e20, 5., -1e20, -1e20]), upper=np.array([1e20, 1e20, 5., 1e20, 1e20]))

# We're trying to minimize fuel burn
prob.model.add_objective('AS_point_0.fuelburn', scaler=.1)


# Set up the problem
# prob.setup(check=True)
prob.setup()

# display the n2 diagram
om.n2(prob)

# Run the model
# prob.run_model()

# Run the optimization
prob.run_driver()

print('wing.twist_cp', prob['wing.twist_cp'])
print('wing.thickness_cp', prob['wing.thickness_cp'])
print('tail.twist_cp', prob['tail.twist_cp'])
print('tail.thickness_cp', prob['tail.thickness_cp'])

print('AS_point_0.wing_perf.failure', prob['AS_point_0.wing_perf.failure'])
print('AS_point_0.wing_perf.thickness_intersects', prob['AS_point_0.wing_perf.thickness_intersects'])
print('AS_point_0.tail_perf.failure', prob['AS_point_0.tail_perf.failure'])
print('AS_point_0.tail_perf.thickness_intersects', prob['AS_point_0.tail_perf.thickness_intersects'])

print('alpha', prob['alpha'])

print('AS_point_0.L_equals_W', prob['AS_point_0.L_equals_W'])
print('AS_point_0.fuelburn', prob['AS_point_0.fuelburn'])
print('wing'+'.struct_setup.structural_mass.structural_mass', prob['wing'+'.struct_setup.structural_mass.structural_mass'])
print('tail'+'.struct_setup.structural_mass.structural_mass', prob['tail'+'.struct_setup.structural_mass.structural_mass'])













