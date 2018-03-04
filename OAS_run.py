# python function which runs aerostruct analysis based on dict input

from __future__ import print_function, division
from six import iteritems   # backwards compatability for python 2

from os import sys, path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))

try:
    from OpenAeroStruct import OASProblem
except:
    from run_classes import OASProblem

import numpy as np

import warnings
warnings.filterwarnings("ignore")

def OAS_setup(user_prob_dict={}, user_surf_list=[]):
    # default prob_dict and surf_dict
    prob_dict = {
        'type' : 'aerostruct',
        'optimize' : False,
        'with_viscous' : True,
        'cg' : np.array([30., 0., 5.]),
        # default design variables, applied to all surfaces
        'des_vars' : ['alpha']  
    }
    surf_list = [{
        'name' : 'wing',
        'num_y' : 7,
        'num_x' : 2,
        'wing_type' : 'CRM',
        'CD0' : 0.015,
        'symmetry' : True,
        'num_twist_cp' : 2,
        'num_thickness_cp': 2,
        'des_vars' : ['thickness_cp','twist_cp']
    }]   
    prob_dict.update(user_prob_dict)

    # remove 'des_vars' key and value from prob_dict
    des_vars = prob_dict.pop('des_vars')

    if user_surf_list:
       surf_list = user_surf_list

    # remove surface des_vars key/value from surface dicts
    for surf in surf_list:
        print(surf)
        if 'des_vars' in surf:
            surf_vars = surf.pop('des_vars', None)
            for var in surf_vars:
                des_vars.append(surf['name']+'.'+var)

    print('des_vars',des_vars)
    # Create OASProblem object 
    OASprob = OASProblem(prob_dict)

    # Add design variables
    # problem-specific design vars...
    prob_des_vars = ['alpha']
    # surface-specific design vars...
    surf_des_vars = ['thickness_cp','twist_cp']

    # Add surfaces and surface design vars to OASProblem
    for surf in surf_list:
        OASprob.add_surface(surf)
    
    for var in des_vars:
        OASprob.add_desvar(var)

    # setup OpenMDAO components in OASProblem
    OASprob.setup()

    return OASprob

def OAS_run(user_des_vars={}, OASprob=None):
    if not OASprob:
        print('setup OAS')
        OASprob = OAS_setup()

    # set design variables
    if user_des_vars:
        for var, value in iteritems(user_des_vars):
            OASprob[var] = value
    
    print('run OAS')
    OASprob.run()
    output = {'fuelburn': OASprob.prob['fuelburn']}
    return output

if __name__ == "__main__":
    print('init')
    out = OAS_run()
    print('end')
    print(out)


