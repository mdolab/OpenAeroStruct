# python function which runs aerostruct analysis based on dict input

from __future__ import print_function, division  # Python 2/3 compatability
from six import iteritems   

from os import sys, path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))

try:
    from OpenAeroStruct import OASProblem
except:
    from run_classes import OASProblem

import numpy as np

import warnings
warnings.filterwarnings("ignore")


iterable_vars = ['chord_cp','thickness_cp','radius_cp','twist_cp',
                'xshear_cp','yshear_cp','zshear_cp']

def OAS_setup(user_prob_dict={}, user_surf_list=[]):
    # default prob_dict and surf_dict
    prob_dict = {
        'type' : 'aerostruct',
        'optimize' : False,
        'with_viscous' : True,
        'cg' : np.array([30., 0., 5.]),
        # default design variables, applied to all surfaces
        'des_vars' : [
            'alpha',
            'wing.thickness_cp',
            'wing.twist_cp',
            'wing.sweep',
            'wing.dihedral',
            'wing.taper',
            'tail.thickness_cp',
            'tail.twist_cp',
            'tail.sweep',
            'tail.dihedral',
            'tail.taper'
        ],
        'output_vars' : [
            'fuelburn', 
            'CD', 
            'CL', 
            'weight'
        ]  
    }
    surf_list = [
        {
            'name' : 'wing',
            'num_y' : 7,
            'num_x' : 2,
            'wing_type' : 'CRM',
            'CD0' : 0.015,
            'symmetry' : True,
            'num_twist_cp' : 2,
            'num_thickness_cp': 2,
         },
         {
             'name' : 'tail',
             'num_y' : 7,
             'num_x' : 2,
             'span' : 20.,
             'root_chord' : 5.,
             'wing_type' : 'rect',
             'offset' : np.array([50., 0., 5.]),
             'twist_cp' : np.array([-9.5])
         }
    ]   
    prob_dict.update(user_prob_dict)

    # remove 'des_vars' key and value from prob_dict
    des_vars = prob_dict.pop('des_vars')

    if user_surf_list:
       #print('user surf')
       surf_list = user_surf_list

    # remove surface des_vars key/value from surface dicts
    for surf in surf_list:
        #print(surf)
        if 'des_vars' in surf:
            surf_vars = surf.pop('des_vars', None)
            for var in surf_vars:
                des_vars.append(surf['name']+'.'+var)

    # check that values in prob_dict and surf_list are the correct ones
    
    # when wrapping from Matlab, an array of a single value will always
    # be converted to a float in Python and not an iterable, which
    # causes problems.



    for surf in surf_list:
    	for key, val in iteritems(surf):
            if (key in iterable_vars) and (not hasattr(val,'__iter__')):
		surf[key] = np.array([val])  # make an ndarray from list

    #print('des_vars',des_vars)
    # Create OASProblem object 
    OASprob = OASProblem(prob_dict)

#    # Add design variables
#    # problem-specific design vars...
#    prob_des_vars = ['alpha']
#    # surface-specific design vars...
#    surf_des_vars = ['thickness_cp','twist_cp']

    # Add surfaces and surface design vars to OASProblem
    for surf in surf_list:
        #print(surf)
	#if 'twist_cp'
        OASprob.add_surface(surf)
    
    for var in des_vars:
        OASprob.add_desvar(var)

    # setup OpenMDAO components in OASProblem
    OASprob.setup()

    return OASprob

def OAS_run(user_des_vars={}, OASprob=None, *args, **kwargs):
    if not OASprob:
        #print('setup OAS')
        OASprob = OAS_setup()

    # set print option
    iprint = kwargs.get('iprint',0)  # set default to only print errors and convergence failures

    # set design variables
    if user_des_vars:
        for var, value in iteritems(user_des_vars):
            #print('$$$$$$      var=',var,'  value=',value)

            if not hasattr(value,'flat'):
		value = np.array([value])  # make an ndarray from list
            OASprob.prob[var] = value
    #print('run OAS')
    
    print('INPUT:')
    print(OASprob.prob.driver.desvars_of_interest())
    #for key, val in iteritems(OASprob.prob._desvars):
    #    print(key+'=',val)
    
    
    OASprob.run()
    #print('after run OAS') 
    
    output = {}
    # get overall output variables and constraints, return None if not there
    overall_vars = ['fuelburn','CD','CL','L_equals_W','CM','v','rho','cg','weighted_obj','total_weight']
    for item in overall_vars:
        output[item] = OASprob.prob[item] 
#        print('item=',item)
#        print('OASprob.prob[item]=',OASprob.prob[item])
        
    # get lifting surface specific variables and constraints, return None if not there
    surface_var_map = {
        'weight' : 'total_perf.<name>structural_weight',
        'CD' : 'total_perf.<name>CD',
        'CL' : 'total_perf.<name>CL',
        'failure' : '<name>perf.failure',
        'vonmises' : '<name>perf.vonmises',
        'thickness_intersects' : '<name>perf.thickness_intersects'       
    }

    # lifting surface coupling variables that need trailing "_" removed from surface name
    coupling_var_map = {
        'loads' : 'coupled.<name>.loads',
        'def_mesh' : 'coupled.<name>.def_mesh' 
    }
    
    for surf in OASprob.surfaces:
        for key, val in iteritems(surface_var_map):
            output.update({surf['name']+key:OASprob.prob[val.replace('<name>',surf['name'])]})
        for key, val in iteritems(coupling_var_map):
            output.update({surf['name']+key:OASprob.prob[val.replace('<name>',surf['name'][:-1])]})
            
    # pretty print output
    #print('OUTPUT:')
    #print(OASprob.prob.driver.outputs_of_interest())
    #for key, val in iteritems(OASoutput):
    #    print(key+' = ',val)
    
    return output

if __name__ == "__main__":
    print('--INIT--')
    OASobj = OAS_setup()
    desvars = {'alpha':0.25}
    out = OAS_run(desvars,OASobj)
    print('--END--')
    #print(out)


