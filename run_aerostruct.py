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

# Append the parent directory to the system path so we can call those Python
# files. If you have OpenAeroStruct in your PYTHONPATH, this is not necessary.
from os import sys, path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))

from OpenAeroStruct import OASProblem

# Suppress warnings
import warnings
warnings.filterwarnings("ignore")

if __name__ == "__main__":

    # Make sure that the user-supplied input is one of the valid options
    input_options = ['0', '0m', '1', '1m']
    print_str = ''.join(str(e) + ', ' for e in input_options)

    # Parse the user-supplied command-line input and store it as input_arg
    try:
        input_arg = sys.argv[1]
        if input_arg not in input_options:
            raise(IndexError)
    except IndexError:
        print('\n +---------------------------------------------------------------+')
        print(' | ERROR: Please supply a correct input argument to this script. |')
        print(' | Possible options are ' + print_str[:-2] + '                             |')
        print(' | See the docstring at the top of this file for more info.      |')
        print(' +---------------------------------------------------------------+\n')
        raise


    # Set problem type
    prob_dict = {'type' : 'aerostruct',
                 'with_viscous' : True,
                 'alpha': 0.,
                 'print_level': 0
                 # 'cg' : np.array([30., 0., 5.])
                 }

    if sys.argv[1].startswith('0'):  # run analysis once
        prob_dict.update({'optimize' : False})
    else:  # perform optimization
        prob_dict.update({'optimize' : True})

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
    # OAS_prob.add_desvar('alpha', lower=-10., upper=10.)
    OAS_prob.add_constraint('L_equals_W', equals=0.)
    OAS_prob.add_objective('fuelburn', scaler=1e-5)

    # Single lifting surface
    if not sys.argv[1].endswith('m'):

        # Setup problem and add design variables, constraint, and objective
        OAS_prob.add_desvar('wing.twist_cp', lower=-15., upper=15.)
        OAS_prob.add_desvar('wing.thickness_cp', lower=0.01, upper=0.5, scaler=1e2)
        OAS_prob.add_constraint('wing_perf.failure', upper=0.)
        OAS_prob.add_constraint('wing_perf.thickness_intersects', upper=0.)
        OAS_prob.setup()

    # Multiple lifting surfaces
    else:

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
        # OAS_prob.add_desvar('tail.twist_cp', lower=-15., upper=15.)
        # OAS_prob.add_desvar('tail.thickness_cp', lower=0.01, upper=0.5, scaler=1e2)
        # OAS_prob.add_constraint('tail_perf.failure', upper=0.)
        # OAS_prob.add_constraint('tail_perf.thickness_intersects', upper=0.)

        # Setup problem
        OAS_prob.setup()

    st = time()
    # Actually run the problem
    if not prob_dict['optimize']:
        thickness_cp = np.array([
            [0.03709672, 0.04913425],
            [0.03, 0.05],
            [0.02, 0.05],
            [0.03777685, 0.07183272]
        ])
        twist_cp = np.array([
            [5.65946593, -3.53034877],
            [5, -3],
            [0, 0],
            [12.80374032, 14.73784563]
        ])
        taper = np.array([
            [0.2],
            [1.0],
            [0.75],
            [0.2]
        ])
        chord_cp = np.array([
            [0.9],
            [1.0],
            [1.1],
            [0.9]
        ])
        out = []
        for i in range(thickness_cp.shape[0]):
            input_dict = {
                'wing.thickness_cp': thickness_cp[i],
                'wing.twist_cp': twist_cp[i],
                'wing.taper': taper[i],
                'wing.chord_cp': chord_cp[i]
            }
            out.append(OAS_prob.run(**input_dict))
            print(' ---  NEW RUN  ---')
            for key, val in out[i].iteritems():
                print(key, val)
    else:
        out = OAS_prob.run()
        for key, val in out.iteritems():
            print(key, val)

    # for key, val in OAS_prob.prob.root._unknowns_dict.iteritems():
    #     print(val['top_promoted_name'],' --> ',key)
    #
    #

    print("\nFuelburn:", OAS_prob.prob['fuelburn'])
    print("Time elapsed: {} secs".format(time() - st))
