""" Example runscript to perform aerodynamics-only optimization.

Call as `python run_vlm.py 0` to run a single analysis, or
call as `python run_vlm.py 1` to perform optimization.

To run with multiple lifting surfaces instead of a single one,
Call as `python run_vlm.py 0m` to run a single analysis, or
call as `python run_vlm.py 1m` to perform optimization.

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

prob_type = 'aerostruct'
mesh_level = 'L3'

# Set problem type
prob_dict = {'optimize' : True,
             'type' : prob_type,
             'cg' : np.array([.4, 0., 0.]),
             'optimizer' : 'SNOPT',
             'with_viscous' : True,
             'W0' : 12.,  # 14-18kg empty weight
             'a' : 322.2,  # m/s at 15,000 ft
             'rho' : 0.770816, # kg/m^3 at 15,000 ft
             'R' : 2500e3, # estimated range based on cruise speed and flight endurance
             'CT' : 9.80665 * 8.6e-6,
             'Re' : 4e5,
             'M' : .1,
             'compute_static_margin' : True,
             }

prob_dict.update({})

# Instantiate problem and add default surface
OAS_prob = OASProblem(prob_dict)

zshear_cp = np.zeros(10)
zshear_cp[0] = .3

xshear_cp = np.zeros(10)
xshear_cp[0] = .15

chord_cp = np.ones(10)
chord_cp[0] = .5

radius_cp = 0.02  * np.ones(10)
radius_cp[0] = 0.015

if mesh_level == 'L1':
    num_y = 101
    num_x = 5
if mesh_level == 'L1.5':
    num_y = 61
    num_x = 3
elif mesh_level == 'L2':
    num_y = 21
    num_x = 3
elif mesh_level == 'L2.5':
    num_y = 15
    num_x = 2
else:
    num_y = 7
    num_x = 2

# Create a dictionary to store options about the surface
surf_dict = {'num_y' : num_y,
             'num_x' : num_x,
             'wing_type' : 'rect',
             'symmetry' : True,
             'span_cos_spacing' : 1.,
             'span' : 3.11,
             'root_chord' : .3,  # estimate
             'sweep' : 20.,
             'taper' : .8,
             'zshear_cp' : zshear_cp,
             'xshear_cp' : xshear_cp,
             'chord_cp' : chord_cp,

             # Material properties taken from http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
             'E' : 85.e9,
             'G' : 25.e9,
             'yield' : 350.e6 / 2.5,
             'mrho' : 1.6e3,
             'CD0' : 0.02,

             }

# Add the specified wing surface to the problem
OAS_prob.add_surface(surf_dict)

if prob_type == 'aero':

    # Setup problem and add design variables, constraint, and objective
    # OAS_prob.add_desvar('alpha', lower=-10., upper=15.)
    OAS_prob.add_desvar('wing.twist_cp', lower=-10., upper=15.)

    # OAS_prob.add_desvar('wing.chord_cp', lower=0.5, upper=3.)
    # OAS_prob.add_desvar('wing.xshear_cp', lower=-10., upper=15.)
    OAS_prob.add_desvar('wing.sweep', lower=-60., upper=60.)
    # OAS_prob.add_desvar('wing.taper', lower=.5, upper=2.)
    OAS_prob.add_constraint('wing_perf.CL', equals=0.5)
    OAS_prob.add_constraint('CM', equals=0.)
    OAS_prob.add_objective('wing_perf.CD', scaler=1e4)

else:
    # Add design variables, constraint, and objective on the problem
    OAS_prob.add_desvar('alpha', lower=-10., upper=10.)
    OAS_prob.add_constraint('eq_con', equals=0.)
    OAS_prob.add_objective('fuelburn', scaler=0.1)

    # Setup problem and add design variables, constraint, and objective
    OAS_prob.add_desvar('wing.twist_cp', lower=-15., upper=15.)
    OAS_prob.add_desvar('wing.thickness_cp', lower=0.0001, upper=0.5, scaler=1e3)
    OAS_prob.add_desvar('wing.sweep', lower=-60., upper=60., scaler=1e-1)
    OAS_prob.add_constraint('wing_perf.failure', upper=0.)
    OAS_prob.add_constraint('wing_perf.thickness_intersects', upper=0.)
    OAS_prob.add_constraint('CM', equals=0.)

OAS_prob.setup()


st = time()

# Actually run the problem
OAS_prob.run()

print("\nWing CL:", OAS_prob.prob['wing_perf.CL'])
print("Wing CD:", OAS_prob.prob['wing_perf.CD'])
print("Time elapsed: {} secs".format(time() - st))
