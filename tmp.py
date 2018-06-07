""" Example runscript to perform aerostructural optimization.

Call as `python run_aerostruct.py 0` to run a single analysis, or
call as `python run_aerostruct.py 1` to perform optimization.

Call as `python run_aerostruct.py 0m` to run a single analysis with
multiple surfaces, or
call as `python run_aerostruct.py 1m` to perform optimization with
multiple surfaces.

"""

from __future__ import division, print_function
import sys
from time import time
import numpy as np
from six import iteritems


# Append the parent directory to the system path so we can call those Python
# files. If you have OpenAeroStruct in your PYTHONPATH, this is not necessary.
from os import sys, path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))

from OpenAeroStruct import OASProblem

# Suppress warnings
import warnings
warnings.filterwarnings("ignore")

def create_OAS_prob(optimize=False):
    # Set problem type
    prob_dict = {'type' : 'aerostruct',
                 'with_viscous' : True,
                 'alpha': 0.,
                 'print_level': 0,
                 'record_db': True,
                 'optimize': optimize
                 }
    # Instantiate problem and add default surface
    OAS_prob = OASProblem(prob_dict)

    # Create a dictionary to store options about the surface
    surf_dict = {'num_y' : 7,
                 'num_x' : 2,
                 'wing_type' : 'CRM',
                 'CD0' : 0.015,
                 'symmetry' : True,
                 'num_twist_cp' : 2,
                 'num_thickness_cp' : 2,
                 'num_chord_cp': 1,
                 'exact_failure_constraint': True,
                 'span_cos_spacing': 0.5}

    # Add the specified wing surface to the problem
    OAS_prob.add_surface(surf_dict)

    # Add design variables, constraint, and objective on the problem
    OAS_prob.add_constraint('L_equals_W', equals=0.)
    OAS_prob.add_objective('fuelburn', scaler=1e-5)

    # Add additional lifting surface
    surf_dict = {'name' : 'tail',
                 'num_y' : 7,
                 'num_x' : 2,
                 'span' : 20.,
                 'root_chord' : 5.,
                 'wing_type' : 'rect',
                 'offset' : np.array([50., 0., 5.]),
                 'twist_cp' : np.array([-9.5]),
                 'exact_failure_constraint': True}
    OAS_prob.add_surface(surf_dict)

    # Add design variables and constraints for both the wing and tail
    OAS_prob.add_desvar('wing.twist_cp', lower=-15., upper=15.)
    OAS_prob.add_desvar('wing.thickness_cp', lower=0.01, upper=0.5, scaler=1e2)
    OAS_prob.add_desvar('wing.chord_cp', lower=0.9, upper=1.1)
    OAS_prob.add_desvar('wing.taper', lower=0.2, upper=1.1)
    OAS_prob.add_constraint('wing_perf.failure', upper=0.)
    OAS_prob.add_constraint('wing_perf.thickness_intersects', upper=0.)

    return OAS_prob

if __name__ == "__main__":

    print('\nANALYSIS')
    OAS_prob_analysis = create_OAS_prob(optimize=False)
    OAS_prob_analysis.setup()
    input_dict = {
        'wing.thickness_cp': np.array([0.03777685, 0.07183272]),
        'wing.twist_cp': np.array([12.80374032, 14.73784563]),
        'wing.taper': 0.2,
        'wing.chord_cp': np.array([0.9])
    }
    for var, val in iteritems(input_dict):
        OAS_prob_analysis.setvar(var, val)
    print('initial values:')
    for var in OAS_prob_analysis.desvars:
        print(var,'=',OAS_prob_analysis.getvar(var))
    print('mesh =\n',OAS_prob_analysis.getvar('wing.mesh'))
    print('Run analysis...')
    # print('mesh =\n',OAS_prob_analysis.getvar('wing.mesh'))
    # out = OAS_prob_analysis.prob.run_once()
    out = OAS_prob_analysis.run()
    for var, val in iteritems(input_dict):
        print(var,'=',OAS_prob_analysis.getvar(var))
    print('mesh =\n',OAS_prob_analysis.getvar('wing.mesh'))
    print('fuelburn =',OAS_prob_analysis.getvar('fuelburn'))

    print('\nOPTIMIZE')
    OAS_prob_optimize = create_OAS_prob(optimize=True)
    OAS_prob_optimize.setup()
    print('initial values:')
    for var in OAS_prob_optimize.desvars:
        print(var,'=',OAS_prob_optimize.getvar(var))
    print('mesh =\n',OAS_prob_optimize.getvar('wing.mesh'))
    print('Run optimization...')
    out = OAS_prob_optimize.run()
    for var in OAS_prob_optimize.desvars:
        print(var,'=',OAS_prob_optimize.getvar(var))
    print('mesh =\n',OAS_prob_optimize.getvar('wing.mesh'))
    print('fuelburn =',OAS_prob_optimize.getvar('fuelburn'))

    # Compare all desvars
    print('\nCOMPARE DESVARS')
    params = OAS_prob_analysis.prob.root._params_dict
    for var, val in iteritems(params):
        top_var = params[var]['top_promoted_name']
        dif = OAS_prob_analysis.getvar(top_var)-OAS_prob_optimize.getvar(top_var)
        err = np.linalg.norm(dif)
        if np.isclose(err,0.0):
            print(top_var,' diff=',dif)


    # Compare all unknowns
    print('\nCOMPARE UNKNOWNS')
    unknowns = OAS_prob_analysis.prob.root._unknowns_dict
    for var, val in iteritems(unknowns):
        top_var = unknowns[var]['top_promoted_name']
        dif = OAS_prob_analysis.getvar(top_var)-OAS_prob_optimize.getvar(top_var)
        err = np.linalg.norm(dif)
        if np.isclose(err,0.0):
            print(top_var,' diff=',dif)
