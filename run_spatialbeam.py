""" Example runscript to perform structural-only optimization.

Call as `python run_spatialbeam.py 0` to run a single analysis, or
call as `python run_spatialbeam.py 1` to perform optimization.

To run with multiple structural components instead of a single one,
call as `python run_spatialbeam.py 0m` to run a single analysis, or
call as `python run_spatialbeam.py 1m` to perform optimization.

"""

from __future__ import division
from time import time
import sys
import numpy

from openmdao.api import IndepVarComp, Problem, Group, ScipyOptimizer, SqliteRecorder
from openmdao.devtools.partition_tree_n2 import view_tree
from geometry import GeometryMesh, Bspline, get_inds, gen_mesh
from spatialbeam import SpatialBeamStates, SpatialBeamFunctionals, radii
from materials import MaterialsTube
from b_spline import get_bspline_mtx

# Single structural component
if not sys.argv[1].endswith('m'):
    num_x = 2  # number of chordwise node points
    num_y = 11  # number of spanwise node points, can only be odd numbers
    span = 10.  # full wingspan
    chord = 5.  # root chord
    cosine_spacing = 0.  # spacing distribution; 0 is uniform, 1 is cosine
    mesh = gen_mesh(num_x, num_y, span, chord, cosine_spacing)

    r = radii(mesh)
    mesh = mesh.reshape(-1, mesh.shape[-1])
    aero_ind = numpy.atleast_2d(numpy.array([num_x, num_y]))
    fem_ind = [num_y]
    aero_ind, fem_ind = get_inds(aero_ind, fem_ind)

# Multiple structural components
else:
    # Wing mesh generation
    num_x = 2  # number of chordwise node points
    num_y = 5  # number of spanwise node points, can only be odd numbers
    span = 5.  # full wingspan
    chord = 5.  # root chord
    cosine_spacing = 0.  # spacing distribution; 0 is uniform, 1 is cosine
    mesh_wing = gen_mesh(num_x, num_y, span, chord, cosine_spacing)
    mesh_wing[:, :, 1] = mesh_wing[:, :, 1] - span/2
    num_twist = numpy.max([int((num_y - 1) / 5), 5])

    # Add wing indices to aero_ind and fem_ind
    r_wing = radii(mesh_wing)
    mesh_wing = mesh_wing.reshape(-1, mesh_wing.shape[-1])
    aero_ind = numpy.atleast_2d(numpy.array([num_x, num_y]))
    fem_ind = [num_y]

    # Tail mesh generation
    nx = 2  # number of chordwise node points
    ny = 5  # number of spanwise node points, can only be odd numbers
    span = 5.  # full tailspan
    chord = 5.  # tail root chord
    cosine_spacing = 0.  # spacing distribution; 0 is uniform, 1 is cosine
    mesh_tail = gen_mesh(nx, ny, span, chord, cosine_spacing)
    mesh_tail[:, :, 1] = mesh_tail[:, :, 1] + span/2

    r_tail = radii(mesh_tail)
    mesh_tail = mesh_tail.reshape(-1, mesh_tail.shape[-1])

    aero_ind = numpy.vstack((aero_ind, numpy.atleast_2d(numpy.array([nx, ny]))))
    mesh = numpy.vstack((mesh_wing, mesh_tail))
    r = numpy.hstack((r_wing, r_tail))

    fem_ind.append(ny)
    aero_ind, fem_ind = get_inds(aero_ind, fem_ind)

# Set the number of thickness control points and the initial thicknesses
num_twist = numpy.max([int((num_y - 1) / 5), 5])
num_thickness = num_twist
t = r / 20
r /= 5

# Set additional mesh parameters
dihedral = 0.  # dihedral angle in degrees
sweep = 0.  # shearing sweep angle in degrees
taper = 1.  # taper ratio

# Define the material properties
execfile('aluminum.py')

# Define the loads
tot_n_fem = numpy.sum(fem_ind[:, 0])
num_surf = fem_ind.shape[0]
loads = numpy.zeros((tot_n_fem, 6))
loads[0, 2] = loads[-1, 2] = 1e3  # tip load of 1 kN

# Create the top-level system
root = Group()

# Define Jacobians for b-spline controls
jac_twist = get_bspline_mtx(num_twist, num_y)
jac_thickness = get_bspline_mtx(num_thickness, tot_n_fem-num_surf)

# Define the independent variables
des_vars = [
    ('twist_cp', numpy.zeros(num_twist)),
    ('thickness_cp', numpy.ones(num_thickness)*numpy.max(t)),
    ('dihedral', dihedral),
    ('sweep', sweep),
    ('span', span),
    ('taper', taper),
    ('r', r),
    ('loads', loads),
    ('fem_ind', fem_ind),
    ('aero_ind', aero_ind),
]

# Add structural components to the top-level system
root.add('des_vars',
         IndepVarComp(des_vars),
         promotes=['*'])
root.add('twist_bsp',
         Bspline('twist_cp', 'twist', jac_twist),
         promotes=['*'])
root.add('thickness_bsp',
         Bspline('thickness_cp', 'thickness', jac_thickness),
         promotes=['*'])
root.add('mesh',
         GeometryMesh(mesh, aero_ind),
         promotes=['*'])
root.add('tube',
         MaterialsTube(fem_ind),
         promotes=['*'])
root.add('spatialbeamstates',
         SpatialBeamStates(aero_ind, fem_ind, E, G),
         promotes=['*'])
root.add('spatialbeamfuncs',
         SpatialBeamFunctionals(aero_ind, fem_ind, E, G, stress, mrho),
         promotes=['*'])

# Set the optimization problem settings
prob = Problem()
prob.root = root

try:  # Use SNOPT optimizer if installed
    from openmdao.api import pyOptSparseDriver
    prob.driver = pyOptSparseDriver()
    prob.driver.options['optimizer'] = "SNOPT"
    prob.driver.opt_settings = {'Major optimality tolerance': 1.0e-8,
                                'Major feasibility tolerance': 1.0e-8}
except:  # Use SLSQP optimizer if SNOPT not installed
    prob.driver = ScipyOptimizer()
    prob.driver.options['optimizer'] = 'SLSQP'
    prob.driver.options['disp'] = True

# Add design variables for the optimizer to control
# Note that the scaling is very important to get correct convergence
prob.driver.add_desvar('thickness_cp',
                       lower=numpy.ones((num_thickness)) * 0.0003,
                       upper=numpy.ones((num_thickness)) * 0.25,
                       scaler=1e5)

# Set objective (minimize weight)
prob.driver.add_objective('weight')

# Set constraint (no structural failure)
prob.driver.add_constraint('failure', upper=0.0)

# Record optimization history to a database
# Data saved here can be examined using `plot_all.py`
prob.driver.add_recorder(SqliteRecorder('spatialbeam.db'))

# Can finite difference over the entire model
# Generally faster than using component derivatives
# Note that for this case, you may need to loosen the optimizer tolerances
# prob.root.deriv_options['type'] = 'fd'

# Setup the problem and produce an N^2 diagram
prob.setup()
view_tree(prob, outfile="spatialbeam.html", show_browser=False)

st = time()
prob.run_once()
if sys.argv[1].startswith('0'):  # run analysis once
    # Uncomment the following line to check derivatives.
    # prob.check_partial_derivatives(compact_print=True)
    pass
elif sys.argv[1].startswith('1'):  # perform optimization
    prob.run()
print "weight", prob['weight']
print "run time", time() - st
