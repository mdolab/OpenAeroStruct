"""
The OASProblem class contains all of the methods necessary to set up and run
aerostructural optimization using OpenAeroStruct.

Check the default dictionary functions to see the default options for the
lifting surfaces and for the entire problem.
Additionally, the setup() and run() methods for each type of analysis and
optimization are defined below.
The portions of the code concerning multiple surfaces may be confusing, but if
you are only interested in using one surface, you can gloss over some of the
details there.
"""

# =============================================================================
# Standard Python modules
# =============================================================================
from __future__ import division, print_function
import sys
from time import time
import numpy as np

# =============================================================================
# OpenMDAO modules
# =============================================================================
from openmdao.api import IndepVarComp, Problem, Group, ScipyOptimizer, Newton, ScipyGMRES, LinearGaussSeidel, NLGaussSeidel, SqliteRecorder, profile, CaseReader, DirectSolver
from openmdao.api import view_model
from six import iteritems

# =============================================================================
# OpenAeroStruct modules
# =============================================================================
from .geometry import GeometryMesh, Bspline, gen_crm_mesh, gen_rect_mesh, MonotonicConstraint
from .transfer import TransferDisplacements, TransferLoads
from .vlm import VLMStates, VLMFunctionals, VLMGeometry
from .spatialbeam import SpatialBeamStates, SpatialBeamFunctionals, SpatialBeamSetup, radii
from .materials import MaterialsTube
from .functionals import TotalPerformance, TotalAeroPerformance, FunctionalBreguetRange, FunctionalEquilibrium
from .gs_newton import HybridGSNewton

try:
    import OAS_API
    fortran_flag = True
    data_type = float
except:
    fortran_flag = False
    data_type = complex

class Error(Exception):
    """
    Format the error message in a box to make it clear this
    was a expliclty raised exception.
    """
    def __init__(self, message):
        msg = '\n+'+'-'*78+'+'+'\n' + '| OpenAeroStruct Error: '
        i = 23
        for word in message.split():
            if len(word) + i + 1 > 78: # Finish line and start new one
                msg += ' '*(78-i)+'|\n| ' + word + ' '
                i = 1 + len(word)+1
            else:
                msg += word + ' '
                i += len(word)+1
        msg += ' '*(78-i) + '|\n' + '+'+'-'*78+'+'+'\n'
        print(msg)
        Exception.__init__(self)

class OASWarning(object):
    """
    Format a warning message
    """
    def __init__(self, message):
        msg = '\n+'+'-'*78+'+'+'\n' + '| OpenAeroStruct Warning: '
        i = 25
        for word in message.split():
            if len(word) + i + 1 > 78: # Finish line and start new one
                msg += ' '*(78-i)+'|\n| ' + word + ' '
                i = 1 + len(word)+1
            else:
                msg += word + ' '
                i += len(word)+1
        msg += ' '*(78-i) + '|\n' + '+'+'-'*78+'+'+'\n'
        print(msg)


class OASProblem(object):
    """
    Contain surface and problem information for aerostructural optimization.

    Parameters
    ----------
    input_dict : dictionary
        The problem conditions and type of analysis desired. Note that there
        are default values defined by `get_default_prob_dict` that are overwritten
        based on the user-provided input_dict.
    """

    def __init__(self, input_dict={}):

        print('Fortran =', fortran_flag)

        # Update prob_dict with user-provided values after getting defaults
        self.prob_dict = self.get_default_prob_dict()
        self.prob_dict.update(input_dict)

        # Set the airspeed velocity based on the supplied Mach number
        # and speed of sound
        self.prob_dict['v'] = self.prob_dict['M'] * self.prob_dict['a']
        self.surfaces = []

        # Set the setup function depending on the problem type selected by the user
        if self.prob_dict['type'] == 'aero':
            self.setup = self.setup_aero
        if self.prob_dict['type'] == 'struct':
            self.setup = self.setup_struct
        if self.prob_dict['type'] == 'aerostruct':
            self.setup = self.setup_aerostruct

        # Set up dictionaries to hold user-supplied parameters for optimization
        self.desvars = {}
        self.constraints = {}
        self.objective = {}

    def get_default_prob_dict(self):
        """
        Obtain the default settings for the problem description. Note that
        these defaults are overwritten based on user input for the problem.

        Returns
        -------
        defaults : dict
            A python dict containing the default problem-level settings.
        """

        defaults = {
                    # Problem and solver options
                    'optimize' : False,      # flag for analysis or optimization
                    'optimizer' : 'SLSQP',   # default optimizer
                    'force_fd' : False,      # if true, we FD over the whole model
                    'with_viscous' : False,  # if true, compute viscous drag
                    'print_level' : 0,       # int to control output during optimization
                                             # 0 for no additional printing
                                             # 1 for nonlinear solver printing
                                             # 2 for nonlinear and linear solver printing
                    'previous_case_db' : None,  # name of the .db file for warm restart
                                                # example: 'aerostruct.db'
                    'record_db' : True,      # True to output .db file
                    'profile' : False,       # True to profile the problem's time costs
                                             # view results using `view_profile prof_raw.0`
                    'compute_static_margin' : False,  # if true, compute and print the
                                                      # static margin after the run is finished

                    # Flow/environment properties
                    'Re' : 1e6,              # Reynolds number
                    'reynolds_length' : 1.0, # characteristic Reynolds length
                    'alpha' : 5.,            # [degrees] angle of attack
                    'M' : 0.84,              # Mach number at cruise
                    'rho' : 0.38,            # [kg/m^3] air density at 35,000 ft
                    'a' : 295.4,             # [m/s] speed of sound at 35,000 ft
                    'g' : 9.80665,           # [m/s^2] acceleration due to gravity

                    # Aircraft properties
                    'CT' : 9.80665 * 17.e-6, # [1/s] (9.80665 N/kg * 17e-6 kg/N/s)
                                             # specific fuel consumption
                    'R' : 11.165e6,            # [m] maximum range (B777-300)
                    'cg' : np.zeros((3)), # Center of gravity for the
                                                 # entire aircraft. Used in trim
                                                 # and stability calculations.
                    'W0' : 0.4 * 3e5,       # [kg] weight of the airplane without
                                            # the wing structure and fuel.
                                            # The default is 40% of the MTOW of
                                            # B777-300 is 3e5 kg.
                    'beta' : 1.,            # weighting factor for mixed objective
                    }

        return defaults


    def get_default_surf_dict(self):
        """
        Obtain the default settings for the surface descriptions. Note that
        these defaults are overwritten based on user input for each surface.
        Each dictionary describes one surface.

        Returns
        -------
        defaults : dict
            A python dict containing the default surface-level settings.
        """

        defaults = {
                    # Wing definition
                    'name' : 'wing',        # name of the surface
                    'num_x' : 3,            # number of chordwise points
                    'num_y' : 5,            # number of spanwise points
                    'span_cos_spacing' : 1, # 0 for uniform spanwise panels
                                            # 1 for cosine-spaced panels
                                            # any value between 0 and 1 for
                                            # a mixed spacing
                    'chord_cos_spacing' : 0.,   # 0 for uniform chordwise panels
                                            # 1 for cosine-spaced panels
                                            # any value between 0 and 1 for
                                            # a mixed spacing
                    'wing_type' : 'rect',   # initial shape of the wing
                                            # either 'CRM' or 'rect'
                                            # 'CRM' can have different options
                                            # after it, such as 'CRM:alpha_2.75'
                                            # for the CRM shape at alpha=2.75
                    'offset' : np.zeros((3)), # coordinates to offset
                                    # the surface from its default location
                    'symmetry' : True,     # if true, model one half of wing
                                            # reflected across the plane y = 0
                    'S_ref_type' : 'wetted', # how we compute the wing area,
                                             # can be 'wetted' or 'projected'

                    # Simple Geometric Variables
                    'span' : 10.,           # full wingspan, even for symmetric cases
                    'root_chord' : 1.,      # root chord
                    'dihedral' : 0.,        # wing dihedral angle in degrees
                                            # positive is upward
                    'sweep' : 0.,           # wing sweep angle in degrees
                                            # positive sweeps back
                    'taper' : 1.,           # taper ratio; 1. is uniform chord
                    'S_ref' : None,         # [m^2] area of the lifting surface

                    # B-spline Geometric Variables. The number of control points
                    # for each of these variables can be specified in surf_dict
                    # by adding the prefix "num" to the variable (e.g. num_twist)
                    'twist_cp' : None,
                    'chord_cp' : None,
                    'xshear_cp' : None,
                    'yshear_cp' : None,
                    'zshear_cp' : None,
                    'thickness_cp' : None,
                    'radius_cp' : None,

                    # Geometric variables. The user generally does not need
                    # to change these geometry variables. This is simply
                    # a list of possible geometry variables that is later
                    # filtered down based on which are active.
                    'geo_vars' : ['sweep', 'dihedral', 'twist_cp', 'xshear_cp', 'yshear_cp',
                        'zshear_cp', 'span', 'chord_cp', 'taper', 'thickness_cp', 'radius_cp'],

                    # Aerodynamic performance of the lifting surface at
                    # an angle of attack of 0 (alpha=0).
                    # These CL0 and CD0 values are added to the CL and CD
                    # obtained from aerodynamic analysis of the surface to get
                    # the total CL and CD.
                    # These CL0 and CD0 values do not vary wrt alpha.
                    'CL0' : 0.0,            # CL of the surface at alpha=0
                    'CD0' : 0.0,            # CD of the surface at alpha=0

                    # Airfoil properties for viscous drag calculation
                    'k_lam' : 0.05,         # percentage of chord with laminar
                                            # flow, used for viscous drag
                    't_over_c' : 0.12,      # thickness over chord ratio (NACA0012)
                    'c_max_t' : .303,       # chordwise location of maximum (NACA0012)
                                            # thickness

                    # Structural values are based on aluminum 7075
                    'E' : 70.e9,            # [Pa] Young's modulus of the spar
                    'G' : 30.e9,            # [Pa] shear modulus of the spar
                    'yield' : 500.e6 / 2.5, # [Pa] yield stress divided by 2.5 for limiting case
                    'mrho' : 3.e3,          # [kg/m^3] material density
                    'fem_origin' : 0.35,    # normalized chordwise location of the spar
                    'loads' : None,         # [N] allow the user to input loads
                    'disp' : None,          # [m] nodal displacements of the FEM model

                    # Constraints
                    'exact_failure_constraint' : False, # if false, use KS function
                    'monotonic_con' : None, # add monotonic constraint to the given
                                            # distributed variable. Ex. 'chord_cp'
                    }
        return defaults


    def add_surface(self, input_dict={}):
        """
        Add a surface to the problem. One surface definition is needed for
        each planar lifting surface.

        Parameters
        ----------
        input_dict : dictionary
            Surface definition. Note that there are default values defined by
            `get_default_surface` that are overwritten based on the
            user-provided input_dict.
        """

        # Get defaults and update surface with the user-provided input
        surf_dict = self.get_default_surf_dict()
        surf_dict.update(input_dict)

        # Check to see if the user provides the mesh points. If they do,
        # get the chordwise and spanwise number of points
        if 'mesh' in surf_dict.keys():
            mesh = surf_dict['mesh']
            num_x, num_y = mesh.shape[:2]

        # If the user doesn't provide a mesh, obtain the values from surface
        # to create the mesh
        elif 'num_x' in surf_dict.keys():
            num_x = surf_dict['num_x']
            num_y = surf_dict['num_y']
            span = surf_dict['span']
            chord = surf_dict['root_chord']
            span_cos_spacing = surf_dict['span_cos_spacing']
            chord_cos_spacing = surf_dict['chord_cos_spacing']

            # Check to make sure that an odd number of spanwise points (num_y) was provided
            if not num_y % 2:
                Error('num_y must be an odd number.')

            # Generate rectangular mesh
            if surf_dict['wing_type'] == 'rect':
                mesh = gen_rect_mesh(num_x, num_y, span, chord,
                    span_cos_spacing, chord_cos_spacing)

            # Generate CRM mesh. Note that this outputs twist information
            # based on the data from the CRM definition paper, so we save
            # this twist information to the surf_dict.
            elif 'CRM' in surf_dict['wing_type']:
                mesh, eta, twist = gen_crm_mesh(num_x, num_y, span, chord,
                    span_cos_spacing, chord_cos_spacing, surf_dict['wing_type'])
                num_x, num_y = mesh.shape[:2]
                surf_dict['crm_twist'] = twist

            else:
                Error('wing_type option not understood. Must be either a type of ' +
                      '"CRM" or "rect".')

            # Chop the mesh in half if using symmetry during analysis.
            # Note that this means that the provided mesh should be the full mesh
            if surf_dict['symmetry']:
                num_y = int((num_y+1)/2)
                mesh = mesh[:, :num_y, :]

        else:
            Error("Please either provide a mesh or a valid set of parameters.")

        # Compute span. We need .real to make span to avoid OpenMDAO warnings.
        quarter_chord = 0.25 * mesh[-1] + 0.75 * mesh[0]
        surf_dict['span'] = max(quarter_chord[:, 1]).real - min(quarter_chord[:, 1]).real
        if surf_dict['symmetry']:
            surf_dict['span'] *= 2.

        # Apply the user-provided coordinate offset to position the mesh
        mesh = mesh + surf_dict['offset']

        # We need to initialize some variables to ones and some others to zeros.
        # Here we define the lists for each case.
        ones_list = ['chord_cp', 'thickness_cp', 'radius_cp']
        zeros_list = ['twist_cp', 'xshear_cp', 'yshear_cp', 'zshear_cp']
        surf_dict['bsp_vars'] = ones_list + zeros_list

        # Loop through bspline variables and set the number of control points if
        # the user hasn't initalized the array.
        for var in surf_dict['bsp_vars']:
            numkey = 'num_' + var
            if surf_dict[var] is None:
                if numkey not in input_dict:
                    surf_dict[numkey] = np.max([int((num_y - 1) / 5), min(5, num_y-1)])
            else:
                surf_dict[numkey] = len(surf_dict[var])

        # Interpolate the twist values from the CRM wing definition to the twist
        # control points
        if 'CRM' in surf_dict['wing_type']:
            num_twist = surf_dict['num_twist_cp']

            # If the surface is symmetric, simply interpolate the initial
            # twist_cp values based on the mesh data
            if surf_dict['symmetry']:
                twist = np.interp(np.linspace(0, 1, num_twist), eta, surf_dict['crm_twist'])
            else:

                # If num_twist is odd, create the twist vector and mirror it
                # then stack the two together, but remove the duplicated twist
                # value.
                if num_twist % 2:
                    twist = np.interp(np.linspace(0, 1, (num_twist+1)/2), eta, surf_dict['crm_twist'])
                    twist = np.hstack((twist[:-1], twist[::-1]))

                # If num_twist is even, mirror the twist vector and stack
                # them together
                else:
                    twist = np.interp(np.linspace(0, 1, num_twist/2), eta, surf_dict['crm_twist'])
                    twist = np.hstack((twist, twist[::-1]))

            # Continue to use the user-defined twist_cp if inputted to the
            # surface dictionary. Otherwise, use the prescribed CRM twist.
            if surf_dict['twist_cp'] is None:
                surf_dict['twist_cp'] = twist

        # Store updated values
        surf_dict['num_x'] = num_x
        surf_dict['num_y'] = num_y
        surf_dict['mesh'] = mesh

        radius = radii(mesh, surf_dict['t_over_c'])
        surf_dict['radius'] = radius

        # Set initial thicknesses
        surf_dict['thickness'] = radius / 10

        # We now loop through the possible bspline variables and populate
        # the 'initial_geo' list with the variables that the geometry
        # or user provided. For example, the CRM wing defines an initial twist.
        # We must treat this separately so we add a twist bspline component
        # even if it is not a desvar.
        surf_dict['initial_geo'] = []
        for var in surf_dict['geo_vars']:

            # Add the bspline variables when they're needed
            if var in surf_dict['bsp_vars']:
                numkey = 'num_' + var
                if surf_dict[var] is None:

                    # Add the intialized geometry variables to either ones or zeros.
                    # These initial values do not perturb the mesh.
                    if var in ones_list:
                        surf_dict[var] = np.ones(surf_dict[numkey], dtype=data_type)
                    elif var in zeros_list:
                        surf_dict[var] = np.zeros(surf_dict[numkey], dtype=data_type)
                else:
                    surf_dict['initial_geo'].append(var)

            # If the user provided a scalar variable (span, sweep, taper, etc),
            # then include that in the initial_geo list
            elif var in input_dict.keys():
                surf_dict['initial_geo'].append(var)

        if 'thickness_cp' not in surf_dict['initial_geo']:
            surf_dict['thickness_cp'] *= np.max(surf_dict['thickness'])

        if surf_dict['loads'] is None:
            # Set default loads at the tips
            loads = np.zeros((surf_dict['thickness'].shape[0] + 1, 6), dtype=data_type)
            loads[0, 2] = 1e4
            if not surf_dict['symmetry']:
                loads[-1, 2] = 1e4
            surf_dict['loads'] = loads

        if surf_dict['disp'] is None:
            # Set default disp if not provided
            surf_dict['disp'] = np.zeros((surf_dict['num_y'], 6), dtype=data_type)

        # Throw a warning if the user provides two surfaces with the same name
        name = surf_dict['name']
        for surface in self.surfaces:
            if name == surface['name']:
                OASWarning("Two surfaces have the same name.")

        # Append '_' to each repeated surface name
        if not name:
            surf_dict['name'] = name
        else:
            surf_dict['name'] = name + '_'

        # Add the individual surface description to the surface list
        self.surfaces.append(surf_dict)

    def setup_prob(self):
        """
        Short method to select the optimizer. Uses pyOptSparse if available,
        or Scipy's SLSQP otherwise.
        """

        try:  # Use pyOptSparse optimizer if installed
            from openmdao.api import pyOptSparseDriver
            self.prob.driver = pyOptSparseDriver()
            if self.prob_dict['optimizer'] == 'SNOPT':
                self.prob.driver.options['optimizer'] = "SNOPT"
                self.prob.driver.opt_settings = {'Major optimality tolerance': 1.0e-8,
                                                 'Major feasibility tolerance': 1.0e-8,
                                                 'Major iterations limit':400,
                                                 'Minor iterations limit':2000,
                                                 'Iterations limit':1000
                                                 }
            elif self.prob_dict['optimizer'] == 'ALPSO':
                self.prob.driver.options['optimizer'] = 'ALPSO'
                self.prob.driver.opt_settings = {'SwarmSize': 40,
                                                'maxOuterIter': 200,
                                                'maxInnerIter': 6,
                                                'rtol': 1e-5,
                                                'atol': 1e-5,
                                                'dtol': 1e-5,
                                                'printOuterIters': 1
                                                 }
            elif self.prob_dict['optimizer'] == 'NOMAD':
                self.prob.driver.options['optimizer'] = 'NOMAD'
                self.prob.driver.opt_settings = {'maxiter':1000,
                                                'minmeshsize':1e-12,
                                                'minpollsize':1e-12,
                                                'displaydegree':0,
                                                'printfile':1
                                                }
            elif self.prob_dict['optimizer'] == 'SLSQP':
                self.prob.driver.options['optimizer'] = 'SLSQP'
                self.prob.driver.opt_settings = {'ACC' : 1e-10
                                                }

        except:  # Use Scipy SLSQP optimizer if pyOptSparse not installed
            self.prob.driver = ScipyOptimizer()
            self.prob.driver.options['optimizer'] = 'SLSQP'
            self.prob.driver.options['disp'] = True
            self.prob.driver.options['tol'] = 1.0e-10

        # Actually call the OpenMDAO functions to add the design variables,
        # constraints, and objective.
        for desvar_name, desvar_data in iteritems(self.desvars):
            self.prob.driver.add_desvar(desvar_name, **desvar_data)
        for con_name, con_data in iteritems(self.constraints):
            self.prob.driver.add_constraint(con_name, **con_data)
        for obj_name, obj_data in iteritems(self.objective):
            self.prob.driver.add_objective(obj_name, **obj_data)

        # Use finite differences over the entire model if user selected it
        if self.prob_dict['force_fd']:
            self.prob.root.deriv_options['type'] = 'fd'

        # Record optimization history to a database.
        # Data saved here can be examined using `plot_all.py` or `OptView.py`
        if self.prob_dict['record_db']:
            recorder = SqliteRecorder(self.prob_dict['prob_name']+".db")
            recorder.options['record_params'] = True
            recorder.options['record_derivs'] = True
            self.prob.driver.add_recorder(recorder)

        # Profile (time) the problem
        if self.prob_dict['profile']:
            profile.setup(self.prob)
            profile.start()

        # Set up the problem
        self.prob.setup()

        # Use warm start from previous db file if desired.
        # Note that we only have access to the unknowns, not the gradient history.
        if self.prob_dict['previous_case_db'] is not None:

            # Open the previous case and start from the last iteration.
            # Change the -1 value in get_case() if you want to select a different iteration.
            cr = CaseReader(self.prob_dict['previous_case_db'])
            case = cr.get_case(-1)

            # Loop through the unknowns and set them for this problem.
            for param_name, param_data in iteritems(case.unknowns):
                self.prob[param_name] = param_data

    def add_desvar(self, *args, **kwargs):
        """
        Store the design variables and later add them to the OpenMDAO problem.
        """
        self.desvars[str(*args)] = dict(**kwargs)

    def add_constraint(self, *args, **kwargs):
        """
        Store the constraints and later add them to the OpenMDAO problem.
        """
        self.constraints[str(*args)] = dict(**kwargs)

    def add_objective(self, *args, **kwargs):
        """
        Store the objectives and later add them to the OpenMDAO problem.
        """
        self.objective[str(*args)] = dict(**kwargs)

    def run(self):
        """
        Method to actually run analysis or optimization. Also saves history in
        a .db file and creates an N2 diagram to view the problem hierarchy.
        """

        # Have more verbose output about optimization convergence
        if self.prob_dict['print_level']:
            self.prob.print_all_convergence()

        # Save an N2 diagram for the problem
        if self.prob_dict['record_db']:
            view_model(self.prob, outfile=self.prob_dict['prob_name']+".html", show_browser=False)

        # If `optimize` == True in prob_dict, perform optimization. Otherwise,
        # simply pass the problem since analysis has already been run.
        if not self.prob_dict['optimize']:
            # Run a single analysis loop. This shouldn't actually be
            # necessary, but sometimes the .db file is not complete unless we do this.
            self.prob.run_once()
        else:
            # Perform optimization
            self.prob.run()

        # If the problem type is aero or aerostruct, we can compute the static margin.
        # This is a naive tempoerary implementation that currently finite differences
        # over the entire model to obtain the static margin.
        if self.prob_dict['compute_static_margin'] and 'aero' in self.prob_dict['type']:

            # Turn off problem recording (so nothing for these computations
            # appears in the .db file) and get the current CL and CM.
            self.prob.driver.recorders._recorders = []
            CL = self.prob['wing_perf.CL']
            CM = self.prob['CM'][1]
            step = 1e-5

            # Perturb alpha and run an analysis loop to obtain the new CL and CM.
            self.prob['alpha'] += step
            self.prob.run_once()
            CL_new = self.prob['wing_perf.CL']
            CM_new = self.prob['CM'][1]

            # Un-perturb alpha and run a single analysis loop to get the problem
            # back to where it was before we finite differenced.
            self.prob['alpha'] -= step
            self.prob.run_once()

            # Compute, print, and save the static margin in metadata.
            static_margin = -(CM_new - CM) / (CL_new - CL)
            print("Static margin is:", static_margin)
            self.prob.root.add_metadata('static_margin', static_margin)

        # Uncomment this to check the partial derivatives of each component
        # self.prob.check_partial_derivatives(compact_print=True)


    def setup_struct(self):
        """
        Specific method to add the necessary components to the problem for a
        structural problem.
        """

        # Set the problem name if the user doesn't
        if 'prob_name' not in self.prob_dict.keys():
            self.prob_dict['prob_name'] = 'struct'

        # Create the base root-level group
        root = Group()

        # Create the problem and assign the root group
        self.prob = Problem()
        self.prob.root = root

        # Loop over each surface in the surfaces list
        for surface in self.surfaces:

            # Get the surface name and create a group to contain components
            # only for this surface.
            # This group's name is whatever the surface's name is.
            # The default is 'wing'.
            name = surface['name']
            tmp_group = Group()

            # Strip the surface names from the desvars list and save this
            # modified list as self.desvars
            desvar_names = []
            for desvar in self.desvars.keys():

                # Check to make sure that the surface's name is in the design
                # variable and only add the desvar to the list if it corresponds
                # to this surface.
                if name[:-1] in desvar:
                    desvar_names.append(''.join(desvar.split('.')[1:]))

            # Add independent variables that do not belong to a specific component.
            # Note that these are the only ones necessary for structual-only
            # analysis and optimization.
            # Here we check and only add the variables that are desvars or a
            # special var, radius, which is necessary to compute weight.
            indep_vars = [('loads', surface['loads'])]
            for var in surface['geo_vars']:
                if var in desvar_names or 'thickness' in var or var in surface['initial_geo']:
                    indep_vars.append((var, surface[var]))

            # Add structural components to the surface-specific group
            tmp_group.add('indep_vars',
                     IndepVarComp(indep_vars),
                     promotes=['*'])
            tmp_group.add('mesh',
                     GeometryMesh(surface, self.desvars),
                     promotes=['*'])
            tmp_group.add('tube',
                     MaterialsTube(surface),
                     promotes=['*'])
            tmp_group.add('struct_setup',
                     SpatialBeamSetup(surface),
                     promotes=['*'])
            tmp_group.add('struct_states',
                     SpatialBeamStates(surface),
                     promotes=['*'])
            tmp_group.add('struct_funcs',
                     SpatialBeamFunctionals(surface),
                     promotes=['*'])

            # Add bspline components for active bspline geometric variables.
            # We only add the component if the corresponding variable is a desvar
            # or special (radius).
            for var in surface['bsp_vars']:
                if var in desvar_names or var in surface['initial_geo'] or 'thickness' in var:
                    n_pts = surface['num_y']
                    if var in ['thickness_cp', 'radius_cp']:
                        n_pts -= 1
                    trunc_var = var.split('_')[0]
                    tmp_group.add(trunc_var + '_bsp',
                             Bspline(var, trunc_var, surface['num_'+var], n_pts),
                             promotes=['*'])

            # Add tmp_group to the problem with the name of the surface.
            # The default is 'wing'.
            root.add(name[:-1], tmp_group, promotes=[])

            root.add_metadata(surface['name'] + 'yield_stress', surface['yield'])
            root.add_metadata(surface['name'] + 'fem_origin', surface['fem_origin'])

        # Actually set up the problem
        self.setup_prob()

    def setup_aero(self):
        """
        Specific method to add the necessary components to the problem for an
        aerodynamic problem.
        """

        # Set the problem name if the user doesn't
        if 'prob_name' not in self.prob_dict.keys():
            self.prob_dict['prob_name'] = 'aero'

        # Create the base root-level group
        root = Group()

        # Create the problem and assign the root group
        self.prob = Problem()
        self.prob.root = root

        # Loop over each surface in the surfaces list
        for surface in self.surfaces:

            # Get the surface name and create a group to contain components
            # only for this surface
            name = surface['name']
            tmp_group = Group()

            # Strip the surface names from the desvars list and save this
            # modified list as self.desvars
            desvar_names = []
            for desvar in self.desvars.keys():

                # Check to make sure that the surface's name is in the design
                # variable and only add the desvar to the list if it corresponds
                # to this surface.
                if name[:-1] in desvar:
                    desvar_names.append(''.join(desvar.split('.')[1:]))

            # Add independent variables that do not belong to a specific component
            indep_vars = [('disp', surface['disp'])]
            for var in surface['geo_vars']:
                if var in desvar_names or var in surface['initial_geo']:
                    indep_vars.append((var, surface[var]))

            # Add aero components to the surface-specific group
            tmp_group.add('indep_vars',
                     IndepVarComp(indep_vars),
                     promotes=['*'])
            tmp_group.add('mesh',
                     GeometryMesh(surface, self.desvars),
                     promotes=['*'])
            tmp_group.add('def_mesh',
                     TransferDisplacements(surface),
                     promotes=['*'])
            tmp_group.add('vlmgeom',
                     VLMGeometry(surface),
                     promotes=['*'])

            # Add bspline components for active bspline geometric variables.
            # We only add the component if the corresponding variable is a desvar.
            for var in surface['bsp_vars']:
                if var in desvar_names or var in surface['initial_geo']:
                    n_pts = surface['num_y']
                    if var in ['thickness_cp', 'radius_cp']:
                        n_pts -= 1
                    trunc_var = var.split('_')[0]
                    tmp_group.add(trunc_var + '_bsp',
                             Bspline(var, trunc_var, surface['num_'+var], n_pts),
                             promotes=['*'])

            # Add monotonic constraints for selected variables
            if surface['monotonic_con'] is not None:
                if type(surface['monotonic_con']) is not list:
                    surface['monotonic_con'] = [surface['monotonic_con']]
                for var in surface['monotonic_con']:
                    tmp_group.add('monotonic_' + var,
                        MonotonicConstraint(var, surface), promotes=['*'])

            # Add tmp_group to the problem as the name of the surface.
            # Note that is a group and performance group for each
            # individual surface.
            name_orig = name.strip('_')
            root.add(name_orig, tmp_group, promotes=[])
            root.add(name_orig+'_perf', VLMFunctionals(surface, self.prob_dict),
                    promotes=["v", "alpha", "M", "re", "rho"])

        # Add problem information as an independent variables component
        if self.prob_dict['Re'] == 0:
            Error('Reynolds number must be greater than zero for viscous drag ' +
            'calculation. If only inviscid drag is desired, set with_viscous ' +
            'flag to False.')

        prob_vars = [('v', self.prob_dict['v']),
            ('alpha', self.prob_dict['alpha']),
            ('M', self.prob_dict['M']),
            ('re', self.prob_dict['Re']/self.prob_dict['reynolds_length']),
            ('rho', self.prob_dict['rho']),
            ('cg', self.prob_dict['cg'])]

        root.add('prob_vars',
                 IndepVarComp(prob_vars),
                 promotes=['*'])

        # Add a single 'aero_states' component that solves for the circulations
        # and forces from all the surfaces.
        # While other components only depends on a single surface,
        # this component requires information from all surfaces because
        # each surface interacts with the others.
        root.add('aero_states',
                 VLMStates(self.surfaces),
                 promotes=['circulations', 'v', 'alpha', 'rho'])

        # Explicitly connect parameters from each surface's group and the common
        # 'aero_states' group.
        # This is necessary because the VLMStates component requires information
        # from each surface, but this information is stored within each
        # surface's group.
        for surface in self.surfaces:
            name = surface['name']

            # Perform the connections with the modified names within the
            # 'aero_states' group.
            root.connect(name[:-1] + '.def_mesh', 'aero_states.' + name + 'def_mesh')
            root.connect(name[:-1] + '.b_pts', 'aero_states.' + name + 'b_pts')
            root.connect(name[:-1] + '.c_pts', 'aero_states.' + name + 'c_pts')
            root.connect(name[:-1] + '.normals', 'aero_states.' + name + 'normals')

            # Connect the results from 'aero_states' to the performance groups
            root.connect('aero_states.' + name + 'sec_forces', name + 'perf' + '.sec_forces')

            # Connect S_ref for performance calcs
            root.connect(name[:-1] + '.S_ref', name + 'perf' + '.S_ref')
            root.connect(name[:-1] + '.widths', name + 'perf' + '.widths')
            root.connect(name[:-1] + '.chords', name + 'perf' + '.chords')
            root.connect(name[:-1] + '.lengths', name + 'perf' + '.lengths')
            root.connect(name[:-1] + '.cos_sweep', name + 'perf' + '.cos_sweep')

            # Connect S_ref for performance calcs
            root.connect(name[:-1] + '.S_ref', 'total_perf.' + name + 'S_ref')
            root.connect(name[:-1] + '.widths', 'total_perf.' + name + 'widths')
            root.connect(name[:-1] + '.chords', 'total_perf.' + name + 'chords')
            root.connect(name[:-1] + '.b_pts', 'total_perf.' + name + 'b_pts')
            root.connect(name + 'perf' + '.CL', 'total_perf.' + name + 'CL')
            root.connect(name + 'perf' + '.CD', 'total_perf.' + name + 'CD')
            root.connect('aero_states.' + name + 'sec_forces', 'total_perf.' + name + 'sec_forces')

        root.add('total_perf',
                  TotalAeroPerformance(self.surfaces, self.prob_dict),
                  promotes=['CM', 'CL', 'CD', 'v', 'rho', 'cg'])

        # Actually set up the problem
        self.setup_prob()

    def setup_aerostruct(self):
        """
        Specific method to add the necessary components to the problem for an
        aerostructural problem.

        Because this code has been extended to work for multiple aerostructural
        surfaces, a good portion of it is spent doing the bookkeeping for parameter
        passing and ensuring that each component modifies the correct data.
        """

        # Set the problem name if the user doesn't
        if 'prob_name' not in self.prob_dict.keys():
            self.prob_dict['prob_name'] = 'aerostruct'

        # Create the base root-level group
        root = Group()
        coupled = Group()

        # Create the problem and assign the root group
        self.prob = Problem()
        self.prob.root = root

        # Loop over each surface in the surfaces list
        for surface in self.surfaces:

            # Get the surface name and create a group to contain components
            # only for this surface
            name = surface['name']
            tmp_group = Group()

            # Strip the surface names from the desvars list and save this
            # modified list as self.desvars
            desvar_names = []
            for desvar in self.desvars.keys():

                # Check to make sure that the surface's name is in the design
                # variable and only add the desvar to the list if it corresponds
                # to this surface.
                if name[:-1] in desvar:
                    desvar_names.append(''.join(desvar.split('.')[1:]))

            # Add independent variables that do not belong to a specific component
            indep_vars = []
            for var in surface['geo_vars']:
                if var in desvar_names or var in surface['initial_geo'] or 'thickness' in var:
                    indep_vars.append((var, surface[var]))

            # Add components to include in the surface's group
            tmp_group.add('indep_vars',
                     IndepVarComp(indep_vars),
                     promotes=['*'])
            tmp_group.add('tube',
                     MaterialsTube(surface),
                     promotes=['*'])
            tmp_group.add('mesh',
                     GeometryMesh(surface, self.desvars),
                     promotes=['*'])
            tmp_group.add('struct_setup',
                     SpatialBeamSetup(surface),
                     promotes=['*'])

            # Add bspline components for active bspline geometric variables.
            # We only add the component if the corresponding variable is a desvar,
            # a special parameter (radius), or if the user or geometry provided
            # an initial distribution.
            for var in surface['bsp_vars']:
                if var in desvar_names or var in surface['initial_geo'] or 'thickness' in var:
                    n_pts = surface['num_y']
                    if var in ['thickness_cp', 'radius_cp']:
                        n_pts -= 1
                    trunc_var = var.split('_')[0]
                    tmp_group.add(trunc_var + '_bsp',
                             Bspline(var, trunc_var, surface['num_'+var], n_pts),
                             promotes=['*'])

            # Add monotonic constraints for selected variables
            if surface['monotonic_con'] is not None:
                if type(surface['monotonic_con']) is not list:
                    surface['monotonic_con'] = [surface['monotonic_con']]
                for var in surface['monotonic_con']:
                    tmp_group.add('monotonic_' + var,
                        MonotonicConstraint(var, surface), promotes=['*'])

            # Add tmp_group to the problem with the name of the surface.
            name_orig = name
            name = name[:-1]
            root.add(name, tmp_group, promotes=[])

            # Add components to the 'coupled' group for each surface.
            # The 'coupled' group must contain all components and parameters
            # needed to converge the aerostructural system.
            tmp_group = Group()
            tmp_group.add('def_mesh',
                     TransferDisplacements(surface),
                     promotes=['*'])
            tmp_group.add('aero_geom',
                     VLMGeometry(surface),
                     promotes=['*'])
            tmp_group.add('struct_states',
                     SpatialBeamStates(surface),
                     promotes=['*'])
            tmp_group.struct_states.ln_solver = LinearGaussSeidel()
            tmp_group.struct_states.ln_solver.options['atol'] = 1e-20

            name = name_orig
            coupled.add(name[:-1], tmp_group, promotes=[])

            # Add a loads component to the coupled group
            coupled.add(name_orig + 'loads', TransferLoads(surface), promotes=[])

            # Add a performance group which evaluates the data after solving
            # the coupled system
            tmp_group = Group()

            tmp_group.add('struct_funcs',
                     SpatialBeamFunctionals(surface),
                     promotes=['*'])
            tmp_group.add('aero_funcs',
                     VLMFunctionals(surface, self.prob_dict),
                     promotes=['*'])

            root.add(name_orig + 'perf', tmp_group, promotes=["rho", "v", "alpha", "re", "M"])

            root.add_metadata(surface['name'] + 'yield_stress', surface['yield'])
            root.add_metadata(surface['name'] + 'fem_origin', surface['fem_origin'])

        # Add a single 'aero_states' component for the whole system within the
        # coupled group.
        coupled.add('aero_states',
                 VLMStates(self.surfaces),
                 promotes=['v', 'alpha', 'rho'])

        # Explicitly connect parameters from each surface's group and the common
        # 'aero_states' group.
        for surface in self.surfaces:
            name = surface['name']

            root.connect(name[:-1] + '.K', 'coupled.' + name[:-1] + '.K')

            # Perform the connections with the modified names within the
            # 'aero_states' group.
            root.connect('coupled.' + name[:-1] + '.def_mesh', 'coupled.aero_states.' + name + 'def_mesh')
            root.connect('coupled.' + name[:-1] + '.b_pts', 'coupled.aero_states.' + name + 'b_pts')
            root.connect('coupled.' + name[:-1] + '.c_pts', 'coupled.aero_states.' + name + 'c_pts')
            root.connect('coupled.' + name[:-1] + '.normals', 'coupled.aero_states.' + name + 'normals')

            # Connect the results from 'aero_states' to the performance groups
            root.connect('coupled.aero_states.' + name + 'sec_forces', name + 'perf' + '.sec_forces')

            # Connect the results from 'coupled' to the performance groups
            root.connect('coupled.' + name[:-1] + '.def_mesh', 'coupled.' + name + 'loads.def_mesh')
            root.connect('coupled.aero_states.' + name + 'sec_forces', 'coupled.' + name + 'loads.sec_forces')

            # Connect the output of the loads component with the FEM
            # displacement parameter. This links the coupling within the coupled
            # group that necessitates the subgroup solver.
            root.connect('coupled.' + name + 'loads.loads', 'coupled.' + name[:-1] + '.loads')

            # Connect aerodyamic mesh to coupled group mesh
            root.connect(name[:-1] + '.mesh', 'coupled.' + name[:-1] + '.mesh')

            # Connect performance calculation variables
            root.connect(name[:-1] + '.radius', name + 'perf.radius')
            root.connect(name[:-1] + '.A', name + 'perf.A')
            root.connect(name[:-1] + '.thickness', name + 'perf.thickness')

            # Connection performance functional variables
            root.connect(name + 'perf.structural_weight', 'total_perf.' + name + 'structural_weight')
            root.connect(name + 'perf.L', 'total_perf.' + name + 'L')
            root.connect(name + 'perf.CL', 'total_perf.' + name + 'CL')
            root.connect(name + 'perf.CD', 'total_perf.' + name + 'CD')
            root.connect('coupled.aero_states.' + name + 'sec_forces', 'total_perf.' + name + 'sec_forces')

            # Connect parameters from the 'coupled' group to the performance
            # groups for the individual surfaces.
            root.connect(name[:-1] + '.nodes', name + 'perf.nodes')
            root.connect('coupled.' + name[:-1] + '.disp', name + 'perf.disp')
            root.connect('coupled.' + name[:-1] + '.S_ref', name + 'perf.S_ref')
            root.connect('coupled.' + name[:-1] + '.widths', name + 'perf.widths')
            root.connect('coupled.' + name[:-1] + '.chords', name + 'perf.chords')
            root.connect('coupled.' + name[:-1] + '.lengths', name + 'perf.lengths')
            root.connect('coupled.' + name[:-1] + '.cos_sweep', name + 'perf.cos_sweep')

            # Connect parameters from the 'coupled' group to the total performance group.
            root.connect('coupled.' + name[:-1] + '.S_ref', 'total_perf.' + name + 'S_ref')
            root.connect('coupled.' + name[:-1] + '.widths', 'total_perf.' + name + 'widths')
            root.connect('coupled.' + name[:-1] + '.chords', 'total_perf.' + name + 'chords')
            root.connect('coupled.' + name[:-1] + '.b_pts', 'total_perf.' + name + 'b_pts')
            root.connect(name + 'perf.cg_location', 'total_perf.' + name + 'cg_location')

        # Set solver properties for the coupled group
        coupled.ln_solver = ScipyGMRES()
        coupled.ln_solver.preconditioner = LinearGaussSeidel()
        coupled.aero_states.ln_solver = LinearGaussSeidel()
        coupled.nl_solver = NLGaussSeidel()

        # This is only available in the most recent version of OpenMDAO.
        # It may help converge tightly coupled systems when using NLGS.
        try:
            coupled.nl_solver.options['use_aitken'] = True
            coupled.nl_solver.options['aitken_alpha_min'] = 0.01
            # coupled.nl_solver.options['aitken_alpha_max'] = 0.5
        except:
            pass

        if self.prob_dict['print_level'] == 2:
            coupled.ln_solver.options['iprint'] = 1
        if self.prob_dict['print_level']:
            coupled.nl_solver.options['iprint'] = 1

        # Add the coupled group to the root problem
        root.add('coupled', coupled, promotes=['v', 'alpha', 'rho'])

        # Add problem information as an independent variables component
        prob_vars = [('v', self.prob_dict['v']),
            ('alpha', self.prob_dict['alpha']),
            ('M', self.prob_dict['M']),
            ('re', self.prob_dict['Re']/self.prob_dict['reynolds_length']),
            ('rho', self.prob_dict['rho'])]

        root.add('prob_vars',
                 IndepVarComp(prob_vars),
                 promotes=['*'])

        # Add functionals to evaluate performance of the system.
        # Note that only the interesting results are promoted here; not all
        # of the parameters.
        root.add('total_perf',
                 TotalPerformance(self.surfaces, self.prob_dict),
                 promotes=['L_equals_W', 'fuelburn', 'CM', 'CL', 'CD', 'v', 'rho', 'cg', 'weighted_obj', 'total_weight'])

        # Actually set up the system
        self.setup_prob()
