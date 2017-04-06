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
from openmdao.api import IndepVarComp, Problem, Group, ScipyOptimizer, Newton, ScipyGMRES, LinearGaussSeidel, NLGaussSeidel, SqliteRecorder, profile
from openmdao.api import view_model

# =============================================================================
# OpenAeroStruct modules
# =============================================================================
from geometry import GeometryMesh, Bspline, gen_crm_mesh, gen_rect_mesh, MonotonicConstraint
from transfer import TransferDisplacements, TransferLoads
from vlm import VLMStates, VLMFunctionals, VLMGeometry
from spatialbeam import SpatialBeamStates, SpatialBeamFunctionals, radii
from materials import MaterialsTube
from functionals import FunctionalBreguetRange, FunctionalEquilibrium

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

        # Update prob_dict with user-provided values after getting defaults
        self.prob_dict = self.get_default_prob_dict()
        self.prob_dict.update(input_dict)

        # Set the airspeed velocity based on the inputted Mach number
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

    def get_default_surf_dict(self):
        """
        Obtain the default settings for the surface descriptions. Note that
        these defaults are overwritten based on user input for each surface.
        Each dictionary describes one surface.
        """

        defaults = {
                    # Wing definition
                    'name' : 'wing',        # name of the surface
                    'num_x' : 3,            # number of chordwise points
                    'num_y' : 5,            # number of spanwise points
                    'root_chord' : 1.,       # root chord
                    'span_cos_spacing' : 1,   # 0 for uniform spanwise panels
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
                    'offset' : np.array([0., 0., 0.]), # coordinates to offset
                                    # the surface from its default location
                    'symmetry' : True,     # if true, model one half of wing
                                            # reflected across the plane y = 0
                    'S_ref_type' : 'wetted',      # 'wetted' or 'projected'

                    # Simple Geometric Variables
                    'span' : 10.,           # full wingspan
                    'dihedral' : 0.,        # wing dihedral angle in degrees
                                            # positive is upward
                    'sweep' : 0.,           # wing sweep angle in degrees
                                            # positive sweeps back
                    'taper' : 1.,           # taper ratio; 1. is uniform chord

                    # B-spline Geometric Variables. The number of control points
                    # for each of these variables can be specified in surf_dict
                    # by adding the prefix "num" to the variable (e.g. num_twist)
                    'twist_cp' : None,
                    'chord_cp' : None,
                    'xshear_cp' : None,
                    'zshear_cp' : None,
                    'thickness_cp' : None,

                    # Active geometric variables. This list can be reduced to only the
                    # design variables if desired.
                    'active_geo_vars' : ['sweep', 'dihedral', 'twist_cp', 'xshear_cp',
                        'zshear_cp', 'span', 'chord_cp', 'taper', 'thickness_cp'],

                    # Zero-lift aerodynamic performance
                    'CL0' : 0.0,            # CL value at AoA (alpha) = 0
                    'CD0' : 0.0,            # CD value at AoA (alpha) = 0

                    # Airfoil properties for viscous drag calculation
                    'k_lam' : 0.05,         # percentage of chord with laminar
                                            # flow, used for viscous drag
                    't_over_c' : 0.12,      # thickness over chord ratio (NACA0012)
                    'c_max_t' : .303,       # chordwise location of maximum (NACA0012)
                                            # thickness

                    # Structural values are based on aluminum
                    'E' : 70.e9,            # [Pa] Young's modulus of the spar
                    'G' : 30.e9,            # [Pa] shear modulus of the spar
                    'stress' : 20.e6,       # [Pa] yield stress
                    'mrho' : 3.e3,          # [kg/m^3] material density
                    'fem_origin' : 0.35,    # chordwise location of the spar

                    'W0' : 0.4 * 3e5,       # [kg] MTOW of B777 is 3e5 kg with fuel

                    # Constraints
                    'exact_failure_constraint' : False, # if false, use KS function
                    'monotonic_con' : None, # add monotonic constraint to the given
                                                # distributed variable
                    }
        return defaults

    def get_default_prob_dict(self):
        """
        Obtain the default settings for the problem description. Note that
        these defaults are overwritten based on user input for the problem.
        """

        defaults = {'optimize' : False,      # flag for analysis or optimization
                    'optimizer' : 'SNOPT',   # default optimizer
                    'Re' : 1e6,              # Reynolds number
                    'reynolds_length' : 1.0, # characteristic Reynolds length
                    'alpha' : 5.,            # angle of attack
                    'CT' : 9.80665 * 17.e-6, # [1/s] (9.81 N/kg * 17e-6 kg/N/s)
                    'R' : 14.3e6,            # [m] maximum range
                    'M' : 0.84,              # Mach number at cruise
                    'rho' : 0.38,            # [kg/m^3] air density at 35,000 ft
                    'a' : 295.4,             # [m/s] speed of sound at 35,000 ft
                    'force_fd' : False,      # if true, we FD over the whole model
                    'with_viscous' : False,  # if true, compute viscous drag
                    'print_level' : 0,       # int to control output during optimization
                                             # 0 for no additional printing
                                             # 1 for nonlinear solver printing
                                             # 2 for nonlinear and linear solver printing
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
            num_x, num_y = mesh.shape

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

        # Compute span. Not why exactly, but we need .real to make span
        # only real and not complex to avoid OpenMDAO warnings.
        quarter_chord = 0.25 * mesh[-1] + 0.75 * mesh[0]
        surf_dict['span'] = max(quarter_chord[:, 1]).real - min(quarter_chord[:, 1]).real
        if surf_dict['symmetry']:
            surf_dict['span'] *= 2.

        # Apply the user-provided coordinate offset to position the mesh
        mesh = mesh + surf_dict['offset']

        # Get the spar radius
        r = radii(mesh)

        # Set the number of twist and thickness control points.
        # These b-spline control points are what the optimizer sees
        # and controls
        ones_list = ['chord_cp', 'thickness_cp']
        zeros_list = ['twist_cp', 'xshear_cp', 'zshear_cp']
        surf_dict['active_bsp_vars'] = list(set(surf_dict['active_geo_vars']) & set(ones_list + zeros_list))

        for var in surf_dict['active_bsp_vars']:
            numkey = 'num_' + var
            if surf_dict[var] is None:
                if numkey not in input_dict:
                    surf_dict[numkey] = np.max([int((num_y - 1) / 5), 5])
                if var in ones_list:
                    surf_dict[var] = np.ones(surf_dict[numkey], dtype=data_type)
                else:
                    surf_dict[var] = np.zeros(surf_dict[numkey], dtype=data_type)
            else:
                surf_dict[numkey] = len(surf_dict[var])

        # If the mesh generation provided an initial twist, set this within
        # the surf_dict object
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

            surf_dict['twist_cp'] = twist

        # Store updated values
        surf_dict['num_x'] = num_x
        surf_dict['num_y'] = num_y
        surf_dict['mesh'] = mesh
        surf_dict['r'] = r
        if 'CRM' in surf_dict['wing_type']:
            surf_dict['t'] = r / 10
        else:
            surf_dict['t'] = r / 20
        if 'thickness_cp' in surf_dict['active_geo_vars']:
            surf_dict['thickness_cp'] *= np.max(surf_dict['t'])

        # Set default loads at the tips
        loads = np.zeros((r.shape[0] + 1, 6), dtype='complex')
        loads[0, 2] = 1e3
        if not surf_dict['symmetry']:
            loads[-1, 2] = 1e3
        surf_dict['loads'] = loads

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
        Short method to select the optimizer. Uses SNOPT if available,
        or SLSQP otherwise.
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
                self.prob.driver.opt_settings = {
                                                }

        except:  # Use Scipy SLSQP optimizer if pyOptSparse not installed
            self.prob.driver = ScipyOptimizer()
            self.prob.driver.options['optimizer'] = 'SLSQP'
            self.prob.driver.options['disp'] = True
            self.prob.driver.options['tol'] = 1.0e-10

    def add_desvar(self, *args, **kwargs):
        """
        Helper function that calls the OpenMDAO method to add design variables.
        """
        self.prob.driver.add_desvar(*args, **kwargs)

    def add_constraint(self, *args, **kwargs):
        """
        Helper function that calls the OpenMDAO method to add constraints.
        """
        self.prob.driver.add_constraint(*args, **kwargs)

    def add_objective(self, *args, **kwargs):
        """
        Helper function that calls the OpenMDAO method to add objectives.
        """
        self.prob.driver.add_objective(*args, **kwargs)

    def run(self):
        """
        Method to actually run analysis or optimization. Also saves history in
        a .db file and creates an N2 diagram to view the problem hierarchy.
        """

        # Use finite differences over the entire model if user selected it
        if self.prob_dict['force_fd']:
            self.prob.root.deriv_options['type'] = 'fd'

        # Record optimization history to a database
        # Data saved here can be examined using `plot_all.py`
        recorder = SqliteRecorder(self.prob_dict['prob_name']+".db")
        recorder.options['record_params'] = True
        recorder.options['record_derivs'] = True
        self.prob.driver.add_recorder(recorder)

        # Profile (time) the problem
        # profile.setup(self.prob)
        # profile.start()

        # Set up the problem
        self.prob.setup()

        # Have more verbose output about optimization convergence
        if self.prob_dict['print_level']:
            self.prob.print_all_convergence()

        # Save an N2 diagram for the problem
        view_model(self.prob, outfile=self.prob_dict['prob_name']+".html", show_browser=False)

        self.prob.run_once()

        # If `optimize` == True in prob_dict, perform optimization. Otherwise,
        # simply pass the problem since analysis has already been run.
        if not self.prob_dict['optimize']:
            # Run a single analysis loop
            pass
        else:
            # Perform optimization
            self.prob.run()

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

            # Add independent variables that do not belong to a specific component.
            # Note that these are the only ones necessary for structual-only
            # analysis and optimization.
            indep_vars = [('r', surface['r']), ('loads', surface['loads'])]
            for var in surface['active_geo_vars']:
                indep_vars.append((var, surface[var]))

            # Add structural components to the surface-specific group
            tmp_group.add('indep_vars',
                     IndepVarComp(indep_vars),
                     promotes=['*'])
            tmp_group.add('mesh',
                     GeometryMesh(surface),
                     promotes=['*'])
            tmp_group.add('tube',
                     MaterialsTube(surface),
                     promotes=['*'])
            tmp_group.add('struct_states',
                     SpatialBeamStates(surface),
                     promotes=['*'])
            tmp_group.add('struct_funcs',
                     SpatialBeamFunctionals(surface),
                     promotes=['*'])
            # Add bspline components for active bspline geometric variables
            for var in surface['active_bsp_vars']:
                n_pts = surface['num_y']
                if var == 'thickness_cp':
                    n_pts -= 1
                trunc_var = var.split('_')[0]
                tmp_group.add(trunc_var + '_bsp',
                         Bspline(var, trunc_var, surface['num_'+var], n_pts),
                         promotes=['*'])

            # Add tmp_group to the problem with the name of the surface.
            # The default is 'wing'.
            root.add(name[:-1], tmp_group, promotes=[])

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

            # Add independent variables that do not belong to a specific component
            indep_vars = [('disp', np.zeros((surface['num_y'], 6), dtype=data_type))]
            for var in surface['active_geo_vars']:
                indep_vars.append((var, surface[var]))

            # Add aero components to the surface-specific group
            tmp_group.add('indep_vars',
                     IndepVarComp(indep_vars),
                     promotes=['*'])
            tmp_group.add('mesh',
                     GeometryMesh(surface),
                     promotes=['*'])
            tmp_group.add('def_mesh',
                     TransferDisplacements(surface),
                     promotes=['*'])
            tmp_group.add('vlmgeom',
                     VLMGeometry(surface),
                     promotes=['*'])
            # Add bspline components for active bspline geometric variables
            for var in surface['active_bsp_vars']:
                n_pts = surface['num_y']
                if var == 'thickness_cp':
                    n_pts -= 1
                trunc_var = var.split('_')[0]
                tmp_group.add(trunc_var + '_bsp',
                         Bspline(var, trunc_var, surface['num_'+var], n_pts),
                         promotes=['*'])
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
            ('rho', self.prob_dict['rho'])]
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
            root.connect(name[:-1] + '.lengths', name + 'perf' + '.lengths')
            root.connect(name[:-1] + '.cos_sweep', name + 'perf' + '.cos_sweep')

        # Actually set up the problem
        self.setup_prob()

    def setup_aerostruct(self):
        """
        Specific method to add the necessary components to the problem for an
        aerostructural problem.
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

            # Add independent variables that do not belong to a specific component
            indep_vars = [('r', surface['r'])]
            for var in surface['active_geo_vars']:
                indep_vars.append((var, surface[var]))

            # Add components to include in the surface's group
            tmp_group.add('indep_vars',
                     IndepVarComp(indep_vars),
                     promotes=['*'])
            tmp_group.add('tube',
                     MaterialsTube(surface),
                     promotes=['*'])
            tmp_group.add('mesh',
                     GeometryMesh(surface),
                     promotes=['*'])
            # Add bspline components for active bspline geometric variables
            for var in surface['active_bsp_vars']:
                n_pts = surface['num_y']
                if var == 'thickness_cp':
                    n_pts -= 1
                trunc_var = var.split('_')[0]
                tmp_group.add(trunc_var + '_bsp',
                         Bspline(var, trunc_var, surface['num_'+var], n_pts),
                         promotes=['*'])

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

        # Add a single 'aero_states' component for the whole system within the
        # coupled group.
        coupled.add('aero_states',
                 VLMStates(self.surfaces),
                 promotes=['v', 'alpha', 'rho'])

        # Explicitly connect parameters from each surface's group and the common
        # 'aero_states' group.
        for surface in self.surfaces:
            name = surface['name']

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
            root.connect('coupled.' + name + 'loads.loads', name + 'perf.loads')

            # Connect the output of the loads component with the FEM
            # displacement parameter. This links the coupling within the coupled
            # group that necessitates the subgroup solver.
            root.connect('coupled.' + name + 'loads.loads', 'coupled.' + name[:-1] + '.loads')

            # Connect aerodyamic mesh to coupled group mesh
            root.connect(name[:-1] + '.mesh', 'coupled.' + name[:-1] + '.mesh')

            # Connect structural design variables
            root.connect(name[:-1] + '.A', 'coupled.' + name[:-1] + '.A')
            root.connect(name[:-1] + '.Iy', 'coupled.' + name[:-1] + '.Iy')
            root.connect(name[:-1] + '.Iz', 'coupled.' + name[:-1] + '.Iz')
            root.connect(name[:-1] + '.J', 'coupled.' + name[:-1] + '.J')

            # Connect performance calculation variables
            root.connect(name[:-1] + '.r', name + 'perf.r')
            root.connect(name[:-1] + '.A', name + 'perf.A')

            # Connection performance functional variables
            root.connect(name + 'perf.weight', 'fuelburn.' + name + 'weight')
            root.connect(name + 'perf.weight', 'eq_con.' + name + 'weight')
            root.connect(name + 'perf.L', 'eq_con.' + name + 'L')
            root.connect(name + 'perf.CL', 'fuelburn.' + name + 'CL')
            root.connect(name + 'perf.CD', 'fuelburn.' + name + 'CD')

            # Connect paramters from the 'coupled' group to the performance
            # group.
            root.connect('coupled.' + name[:-1] + '.nodes', name + 'perf.nodes')
            root.connect('coupled.' + name[:-1] + '.disp', name + 'perf.disp')
            root.connect('coupled.' + name[:-1] + '.S_ref', name + 'perf.S_ref')
            root.connect('coupled.' + name[:-1] + '.widths', name + 'perf.widths')
            root.connect('coupled.' + name[:-1] + '.lengths', name + 'perf.lengths')
            root.connect('coupled.' + name[:-1] + '.cos_sweep', name + 'perf.cos_sweep')

        # Set solver properties for the coupled group
        coupled.ln_solver = ScipyGMRES()
        coupled.ln_solver.preconditioner = LinearGaussSeidel()
        coupled.aero_states.ln_solver = LinearGaussSeidel()

        coupled.nl_solver = NLGaussSeidel()

        if self.prob_dict['print_level'] == 2:
            coupled.ln_solver.options['iprint'] = 1
        if self.prob_dict['print_level']:
            coupled.nl_solver.options['iprint'] = 1

        # Ensure that the groups are ordered correctly within the coupled group
        # so that a system with multiple surfaces is solved corretly.
        order_list = []
        for surface in self.surfaces:
            order_list.append(surface['name'][:-1])
        order_list.append('aero_states')
        for surface in self.surfaces:
            order_list.append(surface['name']+'loads')
        coupled.set_order(order_list)

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
        root.add('fuelburn',
                 FunctionalBreguetRange(self.surfaces, self.prob_dict),
                 promotes=['fuelburn'])
        root.add('eq_con',
                 FunctionalEquilibrium(self.surfaces, self.prob_dict),
                 promotes=['eq_con', 'fuelburn'])

        # Actually set up the system
        self.setup_prob()
