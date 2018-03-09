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
from time import time

import warnings
warnings.filterwarnings("ignore")


def _get_default_prob_dict():
    ''' Settings match the multiple lifting surfaces option in run_aerostruct.py'''
    prob_dict = {
        'type' : 'aerostruct',
        'optimize' : True,
        'with_viscous' : True,
        'cg' : np.array([30., 0., 5.]),
        'des_vars' : [
            ('alpha', {'lower':-10., 'upper':10.}),
            ('wing.twist_cp', {'lower':-15., 'upper':15.}),
            ('wing.thickness_cp', {'lower':0.01, 'upper':0.5, 'scaler':1e2}),
            ('tail.twist_cp', {'lower':-15., 'upper':15.}),
            ('tail.thickness_cp', {'lower':0.01, 'upper':0.5, 'scaler':1e2}),
        ],
        'constraints' : [
            ('L_equals_W', {'equals':0.}),
            ('wing_perf.failure', {'upper':0.}),
            ('wing_perf.thickness_intersects', {'upper':0.}),
            ('tail_perf.failure', {'upper':0.}),
            ('tail_perf.thickness_intersects', {'upper':0.})
        ],
        'objectives' : [
            ('fuelburn', {'scaler':1e-5})
        ],
        'print_level' : 2
    }
    return prob_dict


def _get_default_surf_list():
    ''' Settings match the multiple lifting surfaces option in run_aerostruct.py'''
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
    return surf_list


def OAS_setup(user_prob_dict={}, user_surf_list=[]):

    # default prob_dict and surf_dict's from run_aerostruct.py
    default_prob_dict = _get_default_prob_dict()
    default_prob_dict.update(user_prob_dict)
    prob_dict = default_prob_dict

    # remove 'des_vars', 'contraints', and 'objectives' entries from prob_dict
    # so it doesn't potentially conflict with OASProblem object
    des_vars = prob_dict.pop('des_vars', [])
    constraints = prob_dict.pop('constraints', [])
    objectives = prob_dict.pop('objectives', [])

    if user_surf_list:
        # replace default_surf_list if user supplied one
        surf_list = user_surf_list
    else:
        surf_list = _get_default_surf_list()


    # when wrapping from Matlab, an array of a single value will always be
    # converted to a float in Python and not an iterable, which causes problems.
    iterable_vars = ['chord_cp','thickness_cp','radius_cp','twist_cp',
                    'xshear_cp','yshear_cp','zshear_cp']
    for surf in surf_list:
    	for key, val in iteritems(surf):
            if (key in iterable_vars) and (not hasattr(val,'__iter__')):
		        surf[key] = np.array([val])  # make an ndarray from list

    # Create OASProblem object
    OASprob = OASProblem(prob_dict)

    # Add surfaces to OASProblem
    for surf in surf_list:
        OASprob.add_surface(surf)

    # Add design variables to OASProblem
    for var_tuple in des_vars:
        OASprob.add_desvar(var_tuple[0], **var_tuple[1])

    # Add constraints to OASProblem
    for var_tuple in constraints:
        OASprob.add_constraint(var_tuple[0], **var_tuple[1])

    # Add objectives to OASProblem
    for var_tuple in objectives:
        OASprob.add_objective(var_tuple[0], **var_tuple[1])

    # setup OpenMDAO components in OASProblem
    OASprob.setup()

    return OASprob


def OAS_run(user_des_vars={}, OASprob=None, *args, **kwargs):

    if not OASprob:
        OASprob = OAS_setup()

    # set design variables
    if user_des_vars:
        for var, value in iteritems(user_des_vars):
            if not hasattr(value,'flat'):
                # when wrapping from Matlab, an array of a single value will always be
                # converted to a float in Python and not a numpy array, which causes problems.
		        value = np.array([value])  # make an ndarray from list
            OASprob.prob[var] = value

    # Run OpenAeroStruct
    st = time()
    OASprob.run()
    tend = time() - st

    output = {}

    output['runtime'] = tend   # store runtime

    # get overall output variables and constraints, return None if not there
    overall_vars = ['fuelburn','CD','CL','L_equals_W','CM','v','rho','cg','weighted_obj','total_weight']
    for item in overall_vars:
        output[item] = OASprob.prob[item]

    # get lifting surface specific variables and constraints, return None if not there
    # <name> will be replaced by lifting surface name without trailing "_"
    surface_var_map = {
        'weight' : 'total_perf.<name>_structural_weight',
        'CD' : 'total_perf.<name>_CD',
        'CL' : 'total_perf.<name>_CL',
        'failure' : '<name>_perf.failure',
        'vonmises' : '<name>_perf.vonmises',
        'thickness_intersects' : '<name>_perf.thickness_intersects',
        'radius' : '<name>.radius',
        'A' : '<name>.A',
        'Iy' : '<name>.Iy',
        'Iz' : '<name>.Iz',
        'loads' : 'coupled.<name>.loads',
        'def_mesh' : 'coupled.<name>.def_mesh'
    }

    for surf in OASprob.surfaces:
        for key, val in iteritems(surface_var_map):
            output.update({
                surf['name']+key : OASprob.prob[val.replace('<name>',surf['name'][:-1])]
            })

    return output


if __name__ == "__main__":
    print('--INIT--')

    # Settings to match analysis on multiple surfaces in run_aerostruct.py
    prob_dict = _get_default_prob_dict()
    surf_list = _get_default_surf_list()

    OASobj = OAS_setup(prob_dict, surf_list)

    desvars = {'alpha':0.25}

    # pretty print input
    print('INPUT:')
    for key, val in iteritems(desvars):
        print('{:>14} = {}'.format(key,val))

    print('\nOAS_run()...\n')
    out = OAS_run(desvars, OASobj)

    # pretty print output
    print('\nOUTPUT:')
    for key, val in iteritems(out):
        print('{:>28} = {}'.format(key,val))

    print("\nFuelburn: {}".format(out.get('fuelburn')))
    print("Time elapsed: {} secs".format(out.get('runtime')))

    print('--END--')
