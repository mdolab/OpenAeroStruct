""" Example runscript to perform aerostructural optimization.

Call as `python run_aerostruct.py` to check derivatives, or
call as `python run_aerostruct.py 0` to run a single analysis, or
call as `python run_aerostruct.py 1` to perform optimization. """

from __future__ import division
import sys
from time import time
import numpy

from openmdao.api import IndepVarComp, Problem, Group, ScipyOptimizer, Newton, ScipyGMRES, LinearGaussSeidel, NLGaussSeidel, SqliteRecorder, profile
from openmdao.devtools.partition_tree_n2 import view_tree
from geometry import GeometryMesh, Bspline, gen_crm_mesh, gen_mesh, get_inds
from transfer import TransferDisplacements, TransferLoads
from vlm import VLMStates, VLMFunctionals
from spatialbeam import SpatialBeamStates, SpatialBeamFunctionals, radii
from materials import MaterialsTube
from functionals import FunctionalBreguetRange, FunctionalEquilibrium

from gs_newton import HybridGSNewton
from b_spline import get_bspline_mtx

# Single surface aerostructural optimization
if not sys.argv[1].endswith('m'):
    # Create the mesh with 2 inboard points and 3 outboard points.
    # This will be mirrored to produce a mesh with 7 spanwise points,
    # or 6 spanwise panels
    mesh = gen_crm_mesh(n_points_inboard=3, n_points_outboard=5, num_x=3)
    num_x, num_y = mesh.shape[:2]
    num_twist = numpy.max([int((num_y - 1) / 5), 5])

    r = radii(mesh)
    mesh = mesh.reshape(-1, mesh.shape[-1])
    aero_ind = numpy.atleast_2d(numpy.array([num_x, num_y]))
    fem_ind = [num_y]
    aero_ind, fem_ind = get_inds(aero_ind, fem_ind)

# Multiple surface aerostructural optimization (currently not working)
else:
    num_x = 2  # number of chordwise node points
    num_y = 11  # number of spanwise node points, can only be odd numbers
    span = 60.  # full wingspan
    chord = 15.  # root chord
    cosine_spacing = 0.  # spacing distribution; 0 is uniform, 1 is cosine
    mesh_wing = gen_mesh(num_x, num_y, span, chord, cosine_spacing)
    mesh_wing[:, :, 1] = mesh_wing[:, :, 1]
    num_twist = numpy.max([int((num_y - 1) / 5), 5])

    r_wing = radii(mesh_wing)
    mesh_wing = mesh_wing.reshape(-1, mesh_wing.shape[-1])
    aero_ind = numpy.atleast_2d(numpy.array([num_x, num_y]))
    fem_ind = [num_y]

    nx = 2  # number of chordwise node points
    ny = 5  # number of spanwise node points, can only be odd numbers
    span = 20.  # full tailspan
    chord = 5.  # root chord
    cosine_spacing = 0.  # spacing distribution; 0 is uniform, 1 is cosine
    mesh_tail = gen_mesh(nx, ny, span, chord, cosine_spacing)
    mesh_tail[:, :, 0] += 50

    r_tail = radii(mesh_tail)
    mesh_tail = mesh_tail.reshape(-1, mesh_tail.shape[-1])

    aero_ind = numpy.vstack((aero_ind, numpy.atleast_2d(numpy.array([nx, ny]))))
    mesh = numpy.vstack((mesh_wing, mesh_tail))
    r = numpy.hstack((r_wing, r_tail))

    r /= 5

    fem_ind.append(ny)
    aero_ind, fem_ind = get_inds(aero_ind, fem_ind)


# Set the number of thickness control points and the initial thicknesses
num_thickness = num_twist
t = r / 10

# Define the aircraft properties
execfile('CRM.py')

# Define the material properties
execfile('aluminum.py')

# Create the top-level system
root = Group()

# Define Jacobians for b-spline controls
tot_n_fem = numpy.sum(fem_ind[:, 0])
num_surf = fem_ind.shape[0]
jac_twist = get_bspline_mtx(num_twist, num_y)
jac_thickness = get_bspline_mtx(num_thickness, tot_n_fem-num_surf)

# Define the independent variables
indep_vars = [
    ('span', span),
    ('twist_cp', numpy.zeros(num_twist)),
    ('thickness_cp', numpy.ones(num_thickness)*numpy.max(t)),
    ('v', v),
    ('alpha', alpha),
    ('rho', rho),
    ('r', r),
    ('Re', 0.),  # set Re=0 if you don't want skin frcition drag added
    ('M', M),
    ('aero_ind', aero_ind),
    ('fem_ind', fem_ind)
]

# Add material components to the top-level system
root.add('indep_vars',
         IndepVarComp(indep_vars),
         promotes=['*'])
root.add('twist_bsp',
         Bspline('twist_cp', 'twist', jac_twist),
         promotes=['*'])
root.add('thickness_bsp',
         Bspline('thickness_cp', 'thickness', jac_thickness),
         promotes=['*'])
root.add('tube',
         MaterialsTube(fem_ind),
         promotes=['*'])

# Create a coupled group to contain the aero, sruct, and transfer components
coupled = Group()
coupled.add('mesh',
            GeometryMesh(mesh, aero_ind),
            promotes=['*'])
coupled.add('def_mesh',
            TransferDisplacements(aero_ind, fem_ind),
            promotes=['*'])
coupled.add('vlmstates',
            VLMStates(aero_ind),
            promotes=['*'])
coupled.add('loads',
            TransferLoads(aero_ind, fem_ind),
            promotes=['*'])
coupled.add('spatialbeamstates',
            SpatialBeamStates(aero_ind, fem_ind, E, G),
            promotes=['*'])

# Set solver properties
coupled.nl_solver = Newton()
coupled.nl_solver.options['iprint'] = 1
coupled.ln_solver = ScipyGMRES()
coupled.ln_solver.options['iprint'] = 1
coupled.ln_solver.preconditioner = LinearGaussSeidel()
coupled.vlmstates.ln_solver = LinearGaussSeidel()
coupled.spatialbeamstates.ln_solver = LinearGaussSeidel()

coupled.nl_solver = NLGaussSeidel()   ### Uncomment this out to use NLGS
coupled.nl_solver.options['iprint'] = 1
coupled.nl_solver.options['atol'] = 1e-5
coupled.nl_solver.options['rtol'] = 1e-12

coupled.nl_solver = HybridGSNewton()   ### Uncomment this out to use Hybrid GS Newton
coupled.nl_solver.nlgs.options['iprint'] = 1
coupled.nl_solver.nlgs.options['maxiter'] = 10
coupled.nl_solver.nlgs.options['atol'] = 1e-8
coupled.nl_solver.nlgs.options['rtol'] = 1e-12
coupled.nl_solver.newton.options['atol'] = 1e-7
coupled.nl_solver.newton.options['rtol'] = 1e-7
coupled.nl_solver.newton.options['maxiter'] = 5
coupled.nl_solver.newton.options['iprint'] = 1

# Add the coupled group and functional groups to compute performance
root.add('coupled',
         coupled,
         promotes=['*'])
root.add('vlmfuncs',
         VLMFunctionals(aero_ind, CL0, CD0),
         promotes=['*'])
root.add('spatialbeamfuncs',
         SpatialBeamFunctionals(aero_ind, fem_ind, E, G, stress, mrho),
         promotes=['*'])
root.add('fuelburn',
         FunctionalBreguetRange(W0, CT, a, R, M, aero_ind),
         promotes=['*'])
root.add('eq_con',
         FunctionalEquilibrium(W0, aero_ind),
         promotes=['*'])

# Set the optimization problem settings
prob = Problem()
prob.root = root
prob.print_all_convergence()

try:  # Use SNOPT optimizer if installed
    from openmdao.api import pyOptSparseDriver
    prob.driver = pyOptSparseDriver()
    prob.driver.options['optimizer'] = "SNOPT"
    prob.driver.opt_settings = {'Major optimality tolerance': 1.0e-7,
                                'Major feasibility tolerance': 1.0e-7}
except:  # Use SLSQP optimizer if SNOPT not installed
    prob.driver = ScipyOptimizer()
    prob.driver.options['optimizer'] = 'SLSQP'
    prob.driver.options['disp'] = True
    prob.driver.options['tol'] = 1.0e-8

# Add design variables for the optimizer to control
# Note that the scaling is very important to get correct convergence
prob.driver.add_desvar('twist_cp',lower= -10.,
                       upper=10., scaler=1e0)
#prob.driver.add_desvar('alpha', lower=-10., upper=10., scaler=1000)
prob.driver.add_desvar('thickness_cp',
                       lower= 0.01,
                       upper= 0.25, scaler=1e3)

# Set the objective (minimize fuelburn)
prob.driver.add_objective('fuelburn', scaler=1e-4)

# Set the constraints (no structural failure and lift = weight)
prob.driver.add_constraint('failure', upper=0.0)
prob.driver.add_constraint('eq_con', equals=0.0)

# Record optimization history to a database
# Data saved here can be examined using `plot_all.py`
prob.driver.add_recorder(SqliteRecorder('aerostruct.db'))

# Profile the problem to examine the time spent in individual functions
profile.setup(prob)
profile.start()

# Setup the problem and produce an N^2 diagram
prob.setup()
view_tree(prob, outfile="aerostruct.html", show_browser=False)

st = time()
prob.run_once()

# Run the problem as selected in the command line argument
if len(sys.argv) == 1:  # run analysis once
    pass
elif sys.argv[1].startswith('0'):  # check derivatives
    prob.check_partial_derivatives(compact_print=True)
    pass
elif sys.argv[1].startswith('1'):  # perform optimization
    prob.run()
print "runtime: ", time() - st
