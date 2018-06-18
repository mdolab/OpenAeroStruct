# Univ. Michigan Aerostructural model.
# Based on OpenAeroStruct by John Hwang, and John Jasa (github.com/mdolab/OpenAeroStruct)
# author: Sam Friedman  (samfriedman@tamu.edu)
# date:   4/12/2017

<<<<<<< HEAD
# make compatible Python 2.x to 3.x
from __future__ import print_function, division
# from future.builtins import range  # make compatible Python 2.x to 3.x
# __all__ = ['setup','aerodynamics','structures']
=======
"""
analysis.py

This module contains wrapper functions for each part of the multidisciplinary
analysis of the OpenAeroStruct model. Specifically, this is the
solve_nonlinear() method to each OpenMDAO component in OpenAeroStruct. To use
them, first call the setup() function, which returns an OASProblem object. This
object contains the following attributes:

    OASProblem.prob_dict :   Dictionary of problem parameters
    OASProblem.surfaces  :   List of surface dictionaries defining properties of
                                each lifting surface
    OASProblem.comp_dict :   Dictionary of OpenAeroStruct component objects
                                which contain the analysis of each with a
                                dictionary of problem parameters

For each wrapper function, optionally pass in the necessary component object
from the comp_dict dictionary. Using pre-initialized components drastically
reduces the computation time for a full multidisciplinary analysis. Without
pre-initialization of the component, another argument must be given to initialize
the component within the function. This extra argument is usually the surface
dictionary, but can be other problem or surface parameters. An example with
pre-initiazation is shown in aerodynamics() and structures(). A example without
pre-initialization is shown in aerodynamics2() and structures2().

An example of the multidisciplinary analysis of the coupled system is in the
if __name__=="__main__" function. It uses fixed point iteration to converge the
coupled system of loads and displacements.

Current list of function wrappers available:
    vlm_geometry
    assemble_aic
    aero_circulations
    vlm_forces
    compute_nodes
    assemble_k
    spatial_beam_fem
    spatial_beam_disp
    materials_tube
    geometry_mesh
    transfer_displacements
    transfer_loads

For now, these functions only support a single lifting surface, and does not
support B-spline customization of the lifting surface.

Future work required:
    - Extend functions to be used with multiple lifting surfaces
    - Write wrappers for remaining components in functionals.py, VLMFunctionals,
        SpatialBeamFunctionals
    - Fix BSpline surface customization
    - Complete example of full multidisciplinary analysis in
        if __name__=="__main__" function

"""

# make compatible Python 2.x to 3.x
from __future__ import print_function, division
# from future.builtins import range  # make compatible Python 2.x to 3.x
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3

import numpy as np
import math

from materials import MaterialsTube
from spatialbeam import ComputeNodes, AssembleK, SpatialBeamFEM, SpatialBeamDisp#, SpatialBeamEnergy, SpatialBeamWeight, SpatialBeamVonMisesTube, SpatialBeamFailureKS
from transfer import TransferDisplacements, TransferLoads
from vlm import VLMGeometry, AssembleAIC, AeroCirculations, VLMForces#, VLMLiftDrag, VLMCoeffs, TotalLift, TotalDrag
<<<<<<< HEAD
from geometry import GeometryMesh, Bspline#, gen_crm_mesh, gen_rect_mesh, MonotonicConstraint
from run_classes import OASProblem
# from functionals import FunctionalBreguetRange, FunctionalEquilibrium

# to disable OpenMDAO warnings which will create an error in Matlab
# import warnings
# warnings.filterwarnings("ignore")
=======
from geometry import GeometryMesh#, Bspline, MonotonicConstraint
from run_classes import OASProblem
from openmdao.api import Component
# from functionals import FunctionalBreguetRange, FunctionalEquilibrium

# to disable OpenMDAO warnings which will create an error in Matlab
import warnings
warnings.filterwarnings("ignore")
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3

try:
    import OAS_API
    fortran_flag = True
    data_type = float
except:
    fortran_flag = False
    data_type = complex

"""
<<<<<<< HEAD
--------------------------------------------------------------------------------

                            GEOMETRY / SETUP

--------------------------------------------------------------------------------
=======
================================================================================
                            GEOMETRY / SETUP
================================================================================
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
From run_classes.py: Manipulate geometry mesh based on high-level design parameters """


def setup(prob_dict={}, surfaces=[{}]):
    ''' Setup the aerostruct mesh

    Default wing mesh (single lifting surface):
    -------------------------------------------
    name = 'wing'            # name of the surface
    num_x = 3                # number of chordwise points
    num_y = 5                # number of spanwise points
    root_chord = 1.          # root chord
    span_cos_spacing = 1     # 0 for uniform spanwise panels
                             # 1 for cosine-spaced panels
                             # any value between 0 and 1 for a mixed spacing
    chord_cos_spacing = 0.   # 0 for uniform chordwise panels
                             # 1 for cosine-spaced panels
                             # any value between 0 and 1 for a mixed spacing
    wing_type = 'rect'       # initial shape of the wing either 'CRM' or 'rect'
                             # 'CRM' can have different options after it, such as 'CRM:alpha_2.75' for the CRM shape at alpha=2.75
    offset = np.array([0., 0., 0.]) # coordinates to offset the surface from its default location
    symmetry = True          # if true, model one half of wing reflected across the plane y = 0
    S_ref_type = 'wetted'    # 'wetted' or 'projected'

    # Simple Geometric Variables
    span = 10.               # full wingspan
    dihedral = 0.            # wing dihedral angle in degrees positive is upward
    sweep = 0.               # wing sweep angle in degrees positive sweeps back
    taper = 1.               # taper ratio; 1. is uniform chord

    # B-spline Geometric Variables. The number of control points for each of these variables can be specified in surf_dict
    # by adding the prefix "num" to the variable (e.g. num_twist)
    twist_cp = None
    chord_cp = None
    xshear_cp = None
    zshear_cp = None
    thickness_cp = None

    Default wing parameters:
    ------------------------
    Zero-lift aerodynamic performance
        CL0 = 0.0            # CL value at AoA (alpha) = 0
        CD0 = 0.0            # CD value at AoA (alpha) = 0
    Airfoil properties for viscous drag calculation
        k_lam = 0.05         # percentage of chord with laminar flow, used for viscous drag
        t_over_c = 0.12      # thickness over chord ratio (NACA0012)
        c_max_t = .303       # chordwise location of maximum (NACA0012) thickness
    Structural values are based on aluminum
        E = 70.e9            # [Pa] Young's modulus of the spar
        G = 30.e9            # [Pa] shear modulus of the spar
<<<<<<< HEAD
        stress = 20.e6       # [Pa] yield stress
=======
        yield = 20.e6        # [Pa] yield stress
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
        mrho = 3.e3          # [kg/m^3] material density
        fem_origin = 0.35    # chordwise location of the spar
    Other
        W0 = 0.4 * 3e5       # [kg] MTOW of B777 is 3e5 kg with fuel

    Default problem parameters:
    ---------------------------
    Re = 1e6                 # Reynolds number
    reynolds_length = 1.0    # characteristic Reynolds length
    alpha = 5.               # angle of attack
    CT = 9.80665 * 17.e-6    # [1/s] (9.81 N/kg * 17e-6 kg/N/s)
    R = 14.3e6               # [m] maximum range
    M = 0.84                 # Mach number at cruise
    rho = 0.38               # [kg/m^3] air density at 35,000 ft
    a = 295.4                # [m/s] speed of sound at 35,000 ft
    with_viscous = False     # if true, compute viscous drag

    '''
    # Use steps in run_aerostruct.py to add wing surface to problem

    # Set problem type
    prob_dict.update({'type' : 'aerostruct'})  # this doesn't really matter since we aren't calling OASProblem.setup()

    # Instantiate problem
    OAS_prob = OASProblem(prob_dict)

    for surface in surfaces:
        # Add SpatialBeamFEM size
        FEMsize = 6 * surface['num_y'] + 6
        surface.update({'FEMsize': FEMsize})
        # Add the specified wing surface to the problem.
        OAS_prob.add_surface(surface)

    # Add materials properties for the wing surface to the surface dict in OAS_prob
    for idx, surface in enumerate(OAS_prob.surfaces):
<<<<<<< HEAD
        A, Iy, Iz, J = materials_tube(surface['r'], surface['t'], surface)
=======
        A, Iy, Iz, J = materials_tube(surface['radius'], surface['thickness'], surface)
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
        OAS_prob.surfaces[idx].update({
            'A': A,
            'Iy': Iy,
            'Iz': Iz,
            'J': J
        })

    # Get total panels and save in prob_dict
    tot_panels = 0
    for surface in OAS_prob.surfaces:
        ny = surface['num_y']
        nx = surface['num_x']
        tot_panels += (nx - 1) * (ny - 1)
    OAS_prob.prob_dict.update({'tot_panels': tot_panels})

    # Assume we are only using a single lifting surface for now
    surface = OAS_prob.surfaces[0]

    # Initialize the OpenAeroStruct components and save them in a component dictionary
    comp_dict = {}
    comp_dict['MaterialsTube'] = MaterialsTube(surface)
    comp_dict['GeometryMesh'] = GeometryMesh(surface)
    comp_dict['TransferDisplacements'] = TransferDisplacements(surface)
    comp_dict['VLMGeometry'] = VLMGeometry(surface)
    comp_dict['AssembleAIC'] = AssembleAIC([surface])
    comp_dict['AeroCirculations'] = AeroCirculations(OAS_prob.prob_dict['tot_panels'])
    comp_dict['VLMForces'] = VLMForces([surface])
    comp_dict['TransferLoads'] = TransferLoads(surface)
    comp_dict['ComputeNodes'] = ComputeNodes(surface)
    comp_dict['AssembleK'] = AssembleK(surface)
    comp_dict['SpatialBeamFEM'] = SpatialBeamFEM(surface['FEMsize'])
    comp_dict['SpatialBeamDisp'] = SpatialBeamDisp(surface)
    OAS_prob.comp_dict = comp_dict

<<<<<<< HEAD
    # return the surfaces list, problem dict, and component dict
    surfaces = [surface]
    prob_dict = OAS_prob.prob_dict
    # return surfaces, prob_dict, comp_dict
=======
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    return OAS_prob



def gen_init_mesh(surface, comp_dict=None):
    ''' Generate initial def_mesh '''
    if comp_dict:
<<<<<<< HEAD
        mesh = geometry_mesh(surface, comp=comp_dict['GeometryMesh'])
=======
        mesh = geometry_mesh(surface, comp_dict['GeometryMesh'])
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
        disp = np.zeros((surface['num_y'], 6), dtype=data_type)  # zero displacement
        def_mesh = transfer_displacements(mesh, disp, comp=comp_dict['TransferDisplacements'])
    else:
        mesh = geometry_mesh(surface)
        disp = np.zeros((surface['num_y'], 6), dtype=data_type)  # zero displacement
        def_mesh = transfer_displacements(mesh, disp, surface)
    return def_mesh


<<<<<<< HEAD
def aerodynamics(def_mesh, surface, prob_dict, comp_dict=None):
=======
def aerodynamics(def_mesh, surface, prob_dict, comp_dict):
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    ''' Use pre-initialized components '''

    # Unpack variables
    v = prob_dict.get('v')
    alpha = prob_dict.get('alpha')
    size = prob_dict.get('tot_panels')
    rho = prob_dict.get('rho')

<<<<<<< HEAD
    b_pts, c_pts, widths, cos_sweep, lengths, normals, S_ref = vlm_geometry(def_mesh, comp=comp_dict['VLMGeometry'])
    AIC, rhs= assemble_aic(surface, def_mesh, b_pts, c_pts, normals, v, alpha, comp=comp_dict['AssembleAIC'])
    circulations = aero_circulations(AIC, rhs, comp=comp_dict['AeroCirculations'])
    sec_forces = vlm_forces(surface, def_mesh, b_pts, circulations, alpha, v, rho, comp=comp_dict['VLMForces'])
    loads = transfer_loads(def_mesh, sec_forces, comp=comp_dict['TransferLoads'])
=======
    b_pts, c_pts, widths, cos_sweep, lengths, normals, S_ref = vlm_geometry(def_mesh, comp_dict['VLMGeometry'])
    AIC, rhs= assemble_aic(surface, def_mesh, b_pts, c_pts, normals, v, alpha, comp_dict['AssembleAIC'])
    circulations = aero_circulations(AIC, rhs, comp_dict['AeroCirculations'])
    sec_forces = vlm_forces(surface, def_mesh, b_pts, circulations, alpha, v, rho, comp_dict['VLMForces'])
    loads = transfer_loads(def_mesh, sec_forces, comp_dict['TransferLoads'])
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3

    return loads


def aerodynamics2(def_mesh, surface, prob_dict):
    ''' Don't use pre-initialized components '''

    # Unpack variables
    v = prob_dict.get('v')
    alpha = prob_dict.get('alpha')
    size = prob_dict.get('tot_panels')
    rho = prob_dict.get('rho')

    b_pts, c_pts, widths, cos_sweep, lengths, normals, S_ref = vlm_geometry(def_mesh, surface)
    AIC, rhs= assemble_aic(surface, def_mesh, b_pts, c_pts, normals, v, alpha)
    circulations = aero_circulations(AIC, rhs, size)
    sec_forces = vlm_forces(surface, def_mesh, b_pts, circulations, alpha, v, rho)
    loads = transfer_loads(def_mesh, sec_forces, surface)

    return loads


<<<<<<< HEAD
def structures(loads, surface, prob_dict, comp_dict=None):
=======
def structures(loads, surface, prob_dict, comp_dict):
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    ''' Use pre-initialized components '''

    # Unpack variables
    A = surface.get('A')
    Iy = surface.get('Iy')
    Iz = surface.get('Iz')
    J = surface.get('J')
    mesh = surface.get('mesh')
    v = prob_dict.get('v')
    alpha = prob_dict.get('alpha')
    size =  prob_dict.get('tot_panels')

<<<<<<< HEAD
    nodes = compute_nodes(mesh, comp=comp_dict['ComputeNodes'])
    K, forces = assemble_k(A, Iy, Iz, J, nodes, loads, comp=comp_dict['AssembleK'])
    disp_aug = spatial_beam_fem(K, forces, comp=comp_dict['SpatialBeamFEM'])
    disp = spatial_beam_disp(disp_aug, comp=comp_dict['SpatialBeamDisp'])
    def_mesh = transfer_displacements(mesh, disp, comp=comp_dict['TransferDisplacements'])
=======
    nodes = compute_nodes(mesh, comp_dict['ComputeNodes'])
    K, forces = assemble_k(A, Iy, Iz, J, nodes, loads, comp_dict['AssembleK'])
    disp_aug = spatial_beam_fem(K, forces, comp_dict['SpatialBeamFEM'])
    disp = spatial_beam_disp(disp_aug, comp_dict['SpatialBeamDisp'])
    def_mesh = transfer_displacements(mesh, disp, comp_dict['TransferDisplacements'])
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3

    return def_mesh  # Output the def_mesh matrix


def structures2(loads, surface, prob_dict):
    ''' Don't use pre-initialized components '''

    # Unpack variables
    A = surface.get('A')
    Iy = surface.get('Iy')
    Iz = surface.get('Iz')
    J = surface.get('J')
    mesh = surface.get('mesh')
    FEMsize =  surface.get('FEMsize')
    v = prob_dict.get('v')
    alpha = prob_dict.get('alpha')
     # Add the specified wing surface to the problem.

    nodes = compute_nodes(mesh, surface)
    K, forces = assemble_k(A, Iy, Iz, J, nodes, loads, surface)
    disp_aug = spatial_beam_fem(K, forces, FEMsize)
    disp = spatial_beam_disp(disp_aug, surface)
    def_mesh = transfer_displacements(mesh, disp, surface)

    return def_mesh  # Output the def_mesh matrix


# def cp2pt(cp, jac):
#     """
#     General function to translate from control points to actual points
#     using a b-spline representation.
#     """
#     pt = np.zeros(jac.shape[0])
#     pt = jac.dot(cp)
#     return pt


def geometry_mesh(surface, comp=None):
    """
    OpenMDAO component that performs mesh manipulation functions. It reads in
    the initial mesh from the surface dictionary and outputs the altered
    mesh based on the geometric design variables.

    Parameters
    ----------
    sweep : float
        Shearing sweep angle in degrees.
    dihedral : float
        Dihedral angle in degrees.
    twist[ny] : numpy array
        1-D array of rotation angles for each wing slice in degrees.
    chord_dist[ny] : numpy array
        Chord length for each panel edge.
    taper : float
        Taper ratio for the wing; 1 is untapered, 0 goes to a point at the tip.
<<<<<<< HEAD
=======
    comp : (optional) OpenAeroStruct component object.
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Modified mesh based on the initial mesh in the surface dictionary and
        the geometric design variables.
    """
    if not comp:
        comp = GeometryMesh(surface)
    params = {}
    #
    # The following is copied from the __init__() method of GeometryMesh()
    #
    ny = surface['num_y']
    ones_list = ['taper', 'chord_cp']     # Variables that should be initialized to one
    zeros_list = ['sweep', 'dihedral', 'twist_cp', 'xshear_cp', 'zshear_cp']     # Variables that should be initialized to zero
    set_list = ['span']     # Variables that should be initialized to given value
    all_geo_vars = ones_list + zeros_list + set_list
    geo_params = {}
    for var in all_geo_vars:
        if len(var.split('_')) > 1:
            param = var.split('_')[0]
            if var in ones_list:
                val = np.ones(ny)
            elif var in zeros_list:
                val = np.zeros(ny)
            else:
                val = surface[var]
        else:
            param = var
            if var in ones_list:
                val = 1.0
            elif var in zeros_list:
                val = 0.0
            else:
                val = surface[var]
        geo_params[param] = val
<<<<<<< HEAD
        if var in surface['active_geo_vars']:
=======
        if var in surface['geo_vars']:
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
            params.update({param: val})
    unknowns = {
        'mesh': comp.mesh
    }
    resids = None
    comp.solve_nonlinear(params, unknowns, resids)
    mesh = unknowns.get('mesh')
    return mesh


# def b_spline_surface(surface):
#     """
#     General function to translate from control points to actual points
#     using a b-spline representation.

#     Parameters
#     ----------
#     cpname : string
#         Name of the OpenMDAO component containing the control point values.
#     ptname : string
#         Name of the OpenMDAO component that will contain the interpolated
#         b-spline values.
#     n_input : int
#         Number of input control points.
#     n_output : int
#         Number of outputted interpolated b-spline points.
#     """
#     comp = Bspline(cpname, ptname, n_input, n_output)
#     params = {
#         cpname: cpname
#     }
#     unknowns = {
#         ptname: np.zeros(n_output)
#     }
#     resids = None
#     comp.solve_nonlinear(params, unknowns, resids)
#     ptname_out = unknowns.get(ptname)
#     return ptname_out


<<<<<<< HEAD
def transfer_displacements(mesh, disp, surface=None, comp=None):
=======
def transfer_displacements(mesh, disp, comp):
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    """
    Perform displacement transfer.

    Apply the computed displacements on the original mesh to obtain
    the deformed mesh.

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Flattened array defining the lifting surfaces.
    disp[ny, 6] : numpy array
        Flattened array containing displacements on the FEM component.
        Contains displacements for all six degrees of freedom, including
        displacements in the x, y, and z directions, and rotations about the
        x, y, and z axes.
<<<<<<< HEAD
=======
    comp : Either OpenAeroStruct component object (better), or surface dict.
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3

    Returns
    -------
    def_mesh[nx, ny, 3] : numpy array
        Flattened array defining the lifting surfaces after deformation.
    """
<<<<<<< HEAD
    if not comp:
=======
    if not isinstance(comp, Component):
        surface = comp
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
        comp = TransferDisplacements(surface)
    params = {
        'mesh': mesh,
        'disp': disp
    }
    unknowns = {
        'def_mesh': np.zeros((comp.nx, comp.ny, 3), dtype=data_type)
    }
    resids = None
    comp.solve_nonlinear(params, unknowns, resids)
    def_mesh = unknowns.get('def_mesh')
    return def_mesh


"""
<<<<<<< HEAD
--------------------------------------------------------------------------------

                                AERODYNAMICS

--------------------------------------------------------------------------------
From vlm.py: """


def vlm_geometry(def_mesh, surface=None, comp=None):
=======
================================================================================
                                AERODYNAMICS
================================================================================
From vlm.py: """


def vlm_geometry(def_mesh, comp):
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    """ Compute various geometric properties for VLM analysis.

    Parameters
    ----------
    def_mesh[nx, ny, 3] : numpy array
        Array defining the nodal coordinates of the lifting surface.
<<<<<<< HEAD
=======
    comp : Either OpenAeroStruct component object (better), or surface dict.
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3

    Returns
    -------
    b_pts[nx-1, ny, 3] : numpy array
        Bound points for the horseshoe vortices, found along the 1/4 chord.
    c_pts[nx-1, ny-1, 3] : numpy array
        Collocation points on the 3/4 chord line where the flow tangency
        condition is satisfed. Used to set up the linear system.
    widths[nx-1, ny-1] : numpy array
        The spanwise widths of each individual panel.
    lengths[ny] : numpy array
        The chordwise length of the entire airfoil following the camber line.
    normals[nx-1, ny-1, 3] : numpy array
        The normal vector for each panel, computed as the cross of the two
        diagonals from the mesh points.
    S_ref : float
        The reference area of the lifting surface.
    """
<<<<<<< HEAD
    if not comp:
=======
    if not isinstance(comp, Component):
        surface = comp
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
        comp = VLMGeometry(surface)
    params = {
        'def_mesh': def_mesh
    }
    unknowns = {
        'b_pts': np.zeros((comp.nx-1, comp.ny, 3), dtype=data_type),
        'c_pts': np.zeros((comp.nx-1, comp.ny-1, 3)),
        'widths': np.zeros((comp.ny-1)),
        'cos_sweep': np.zeros((comp.ny-1)),
        'lengths': np.zeros((comp.ny)),
        'normals': np.zeros((comp.nx-1, comp.ny-1, 3)),
        'S_ref': 0.
    }
    resids=None
    comp.solve_nonlinear(params, unknowns, resids)
    b_pts=unknowns.get('b_pts')
    c_pts=unknowns.get('c_pts')
    widths=unknowns.get('widths')
    cos_sweep=unknowns.get('cos_sweep')
    lengths=unknowns.get('lengths')
    normals=unknowns.get('normals')
    S_ref=unknowns.get('S_ref')
    return b_pts, c_pts, widths, cos_sweep, lengths, normals, S_ref


def assemble_aic(surface, def_mesh, b_pts, c_pts, normals, v, alpha, comp=None):
    """
    Compute the circulations based on the AIC matrix and the panel velocities.
    Note that the flow tangency condition is enforced at the 3/4 chord point.
    There are multiple versions of the first four parameters with one
    for each surface defined.
    Each of these parameters has the name of the surface prepended on the
    actual parameter name.

    Parameters
    ----------
    def_mesh[nx, ny, 3] : numpy array
        Array defining the nodal coordinates of the lifting surface.
    b_pts[nx-1, ny, 3] : numpy array
        Bound points for the horseshoe vortices, found along the 1/4 chord.
    c_pts[nx-1, ny-1, 3] : numpy array
        Collocation points on the 3/4 chord line where the flow tangency
        condition is satisfed. Used to set up the linear system.
    normals[nx-1, ny-1, 3] : numpy array
        The normal vector for each panel, computed as the cross of the two
        diagonals from the mesh points.

    v : float
        Freestream air velocity in m/s.
    alpha : float
        Angle of attack in degrees.
<<<<<<< HEAD
=======
    comp : (Optional) OpenAeroStruct component object.
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3

    Returns
    -------
    AIC[(nx-1)*(ny-1), (nx-1)*(ny-1)] : numpy array
        The aerodynamic influence coefficient matrix. Solving the linear system
        of AIC * circulations = n * v gives us the circulations for each of the
        horseshoe vortices.
    rhs[(nx-1)*(ny-1)] : numpy array
        The right-hand-side of the linear system that yields the circulations.
    """
    surfaces = [surface]
    if not comp:
        comp=AssembleAIC(surfaces)
<<<<<<< HEAD

=======
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    params = {}
    ny=surface['num_y']
    nx=surface['num_x']
    name=surface['name']
    params.update({
        name + 'def_mesh': def_mesh,
        name + 'b_pts': b_pts,
        name + 'c_pts': c_pts,
        name + 'normals': normals
    })
    params.update({
        'v': v,
        'alpha': alpha
    })
    unknowns={
        'AIC': np.zeros((comp.tot_panels, comp.tot_panels), dtype = data_type),
        'rhs': np.zeros((comp.tot_panels), dtype = data_type)
    }
    resids=None
    comp.solve_nonlinear(params, unknowns, resids)
    AIC=unknowns.get('AIC')
    rhs=unknowns.get('rhs')
    return AIC, rhs


<<<<<<< HEAD
def aero_circulations(AIC, rhs, size=None, comp=None):
=======
def aero_circulations(AIC, rhs, comp):
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    """
    Compute the circulation strengths of the horseshoe vortices by solving the
    linear system AIC * circulations = n * v.
    This component is copied from OpenMDAO's LinearSystem component with the
    names of the parameters and outputs changed to match our problem formulation.

    Parameters
    ----------
    AIC[(nx-1)*(ny-1), (nx-1)*(ny-1)] : numpy array
        The aerodynamic influence coefficient matrix. Solving the linear system
        of AIC * circulations = n * v gives us the circulations for each of the
        horseshoe vortices.
    rhs[(nx-1)*(ny-1)] : numpy array
        The right-hand-side of the linear system that yields the circulations.
<<<<<<< HEAD
=======
    comp : Either OpenAeroStruct component object (better), or tot_panels.

>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3

    Returns
    -------
    circulations[(nx-1)*(ny-1)] : numpy array
        Augmented displacement array. Obtained by solving the system
        AIC * circulations = n * v.
    """
<<<<<<< HEAD
    if not comp:
        comp = AeroCirculations(size)
    if not size:
        size = comp.size
=======
    if not isinstance(comp, Component):
        tot_panels = comp
        comp = AeroCirculations(tot_panels)
    size = comp.size
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    params = {
        'AIC': AIC,
        'rhs': rhs
    }
    unknowns = {
        'circulations': np.zeros((size), dtype=data_type)
    }
    resids = {
        'circulations': np.zeros((size), dtype=data_type)
    }
    comp.solve_nonlinear(params, unknowns, resids)
    circulations = unknowns.get('circulations')
    return circulations


def vlm_forces(surface, def_mesh, b_pts, circulations, alpha, v, rho, comp=None):
    """ Compute aerodynamic forces acting on each section.

    Note that the first two parameters and the unknown have the surface name
    prepended on it. E.g., 'def_mesh' on a surface called 'wing' would be
    'wing.def_mesh', etc.

    Parameters
    ----------
    def_mesh[nx, ny, 3] : numpy array
        Array defining the nodal coordinates of the lifting surface.
    b_pts[nx-1, ny, 3] : numpy array
        Bound points for the horseshoe vortices, found along the 1/4 chord.

    circulations : numpy array
        Flattened vector of horseshoe vortex strengths calculated by solving
        the linear system of AIC_mtx * circulations = rhs, where rhs is
        based on the air velocity at each collocation point.
    alpha : float
        Angle of attack in degrees.
    v : float
        Freestream air velocity in m/s.
    rho : float
        Air density in kg/m^3.
<<<<<<< HEAD
=======
    comp : (optional) OpenAeroStruct component object.
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3

    Returns
    -------
    sec_forces[nx-1, ny-1, 3] : numpy array
        Flattened array containing the sectional forces acting on each panel.
        Stored in Fortran order (only relevant with more than one chordwise
        panel).

    """
    surfaces = [surface]
    if not comp:
        comp=VLMForces(surfaces)
    params = {}
    unknowns = {}
    tot_panels = 0
    name = surface['name']
    ny = surface['num_y']
    nx = surface['num_x']
    tot_panels += (nx - 1) * (ny - 1)
    params.update({
        name+'def_mesh': def_mesh,
        name+'b_pts': b_pts
    })
    unknowns.update({
        name+'sec_forces': np.zeros((nx-1, ny-1, 3), dtype=data_type)
    })
    params.update({
        'circulations': circulations,
        'alpha': alpha,
        'v': v,
        'rho': rho
    })
    resids=None
    comp.solve_nonlinear(params, unknowns, resids)
    sec_forces=unknowns.get(name+'sec_forces')
    return sec_forces


<<<<<<< HEAD
def transfer_loads(def_mesh, sec_forces, surface=None, comp=None):
=======
def transfer_loads(def_mesh, sec_forces, comp):
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    """
    Perform aerodynamic load transfer.

    Apply the computed sectional forces on the aerodynamic surfaces to
    obtain the deformed mesh FEM loads.

    Parameters
    ----------
    def_mesh[nx, ny, 3] : numpy array
        Flattened array defining the lifting surfaces after deformation.
    sec_forces[nx-1, ny-1, 3] : numpy array
        Flattened array containing the sectional forces acting on each panel.
        Stored in Fortran order (only relevant when more than one chordwise
        panel).
<<<<<<< HEAD
=======
    comp : Either OpenAeroStruct component object (better), or surface dict.
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3

    Returns
    -------
    loads[ny, 6] : numpy array
        Flattened array containing the loads applied on the FEM component,
        computed from the sectional forces.
    """
<<<<<<< HEAD
    if not comp:
=======
    if not isinstance(comp, Component):
        surface = comp
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
        comp=TransferLoads(surface)
    params={
        'def_mesh': def_mesh,
        'sec_forces': sec_forces
    }
    unknowns={
        'loads': np.zeros((comp.ny, 6), dtype=complex)
    }
    resids=None
    comp.solve_nonlinear(params, unknowns, resids)
    loads=unknowns.get('loads')
    return loads


"""
<<<<<<< HEAD
--------------------------------------------------------------------------------

                                   STRUCTURES

--------------------------------------------------------------------------------
From spatialbeam.py: Define the structural analysis component using spatial beam theory. """

def spatial_beam_fem(K, forces, size=None, comp=None):
=======
================================================================================
                                   STRUCTURES
================================================================================
From spatialbeam.py: Define the structural analysis component using spatial beam theory. """

def spatial_beam_fem(K, forces, comp):
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    """
    Compute the displacements and rotations by solving the linear system
    using the structural stiffness matrix.
    This component is copied from OpenMDAO's LinearSystem component with the
    names of the parameters and outputs changed to match our problem formulation.

    Parameters
    ----------
    K[6*(ny+1), 6*(ny+1)] : numpy array
        Stiffness matrix for the entire FEM system. Used to solve the linear
        system K * u = f to obtain the displacements, u.
    forces[6*(ny+1)] : numpy array
        Right-hand-side of the linear system. The loads from the aerodynamic
        analysis or the user-defined loads.
<<<<<<< HEAD

    Optional Parameters
    -------------------
    size : int
        The total number of panels on the surface mesh. Equivalent to pro
=======
    comp : Either OpenAeroStruct component object (better), or FEMsize of surface.


>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    Returns
    -------
    disp_aug[6*(ny+1)] : numpy array
        Augmented displacement array. Obtained by solving the system
        K * u = f, where f is a flattened version of loads.

    """
<<<<<<< HEAD
    if not comp:
        comp=SpatialBeamFEM(size)
    if not size:
        size = comp.size
=======
    if not isinstance(comp, Component):
        FEMsize = comp
        comp=SpatialBeamFEM(FEMsize)
    else:
        FEMsize = comp.size
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    params={
        'K': K,
        'forces': forces
    }
    unknowns={
<<<<<<< HEAD
        'disp_aug': np.zeros((size), dtype=data_type)
    }
    resids={
        'disp_aug': np.zeros((size), dtype=data_type)
=======
        'disp_aug': np.zeros((FEMsize), dtype=data_type)
    }
    resids={
        'disp_aug': np.zeros((FEMsize), dtype=data_type)
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    }
    comp.solve_nonlinear(params, unknowns, resids)
    disp_aug=unknowns.get('disp_aug')
    return disp_aug


<<<<<<< HEAD
def spatial_beam_disp(disp_aug, surface=None, comp=None):
=======
def spatial_beam_disp(disp_aug, comp):
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    """
    Reshape the flattened displacements from the linear system solution into
    a 2D array so we can more easily use the results.

    The solution to the linear system has additional results due to the
    constraints on the FEM model. The displacements from this portion of
    the linear system are not needed, so we select only the relevant
    portion of the displacements for further calculations.

    Parameters
    ----------
    disp_aug[6*(ny+1)] : numpy array
        Augmented displacement array. Obtained by solving the system
        K * disp_aug = forces, where forces is a flattened version of loads.
<<<<<<< HEAD
=======
    comp : Either OpenAeroStruct component object (better), or surface dict.
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3

    Returns
    -------
    disp[6*ny] : numpy array
        Actual displacement array formed by truncating disp_aug.

    """
<<<<<<< HEAD
    if not comp:
=======
    if not isinstance(comp, Component):
        surface = comp
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
        comp=SpatialBeamDisp(surface)
    params={
        'disp_aug': disp_aug
    }
    unknowns={
        'disp': np.zeros((comp.ny, 6), dtype=data_type)
    }
    resids=None
    comp.solve_nonlinear(params, unknowns, resids)
    disp=unknowns.get('disp')
    return disp


<<<<<<< HEAD
def compute_nodes(mesh, surface=None, comp=None):
=======
def compute_nodes(mesh, comp):
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    """
    Compute FEM nodes based on aerodynamic mesh.

    The FEM nodes are placed at fem_origin * chord,
    with the default fem_origin = 0.35.

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Array defining the nodal points of the lifting surface.
<<<<<<< HEAD
=======
    comp : Either OpenAeroStruct component object (better), or surface dict.
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3

    Returns
    -------
    nodes[ny, 3] : numpy array
        Flattened array with coordinates for each FEM node.

    """
<<<<<<< HEAD
    if not comp:
=======
    if not isinstance(comp, Component):
        surface = comp
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
        comp=ComputeNodes(surface)
    params={
        'mesh': mesh
    }
    unknowns={
        'nodes': np.zeros((comp.ny, 3), dtype=data_type)
    }
    resids=None
    comp.solve_nonlinear(params, unknowns, resids)
    nodes=unknowns.get('nodes')
    return nodes


<<<<<<< HEAD
def assemble_k(A, Iy, Iz, J, nodes, loads, surface=None, comp=None):
=======
def assemble_k(A, Iy, Iz, J, nodes, loads, comp):
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    """
    Compute the displacements and rotations by solving the linear system
    using the structural stiffness matrix.

    Parameters
    ----------
    A[ny-1] : numpy array
        Areas for each FEM element.
    Iy[ny-1] : numpy array
<<<<<<< HEAD
        Mass moment of inertia around the y-axis for each FEM element.
    Iz[ny-1] : numpy array
        Mass moment of inertia around the z-axis for each FEM element.
=======
        Area moment of inertia around the y-axis for each FEM element.
    Iz[ny-1] : numpy array
        Area moment of inertia around the z-axis for each FEM element.
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    J[ny-1] : numpy array
        Polar moment of inertia for each FEM element.
    nodes[ny, 3] : numpy array
        Flattened array with coordinates for each FEM node.
    loads[ny, 6] : numpy array
        Flattened array containing the loads applied on the FEM component,
        computed from the sectional forces.
<<<<<<< HEAD
=======
    comp : Either OpenAeroStruct component object (better), or surface dict.
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3

    Returns
    -------
    K[(nx-1)*(ny-1), (nx-1)*(ny-1)] : numpy array
        Stiffness matrix for the entire FEM system. Used to solve the linear
        system K * u = f to obtain the displacements, u.
    forces[(nx-1)*(ny-1)] : numpy array
        Right-hand-side of the linear system. The loads from the aerodynamic
        analysis or the user-defined loads.
    """
<<<<<<< HEAD
    if not comp:
=======
    if not isinstance(comp, Component):
        surface = comp
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
        comp = AssembleK(surface)  # if component is not passed in, surface must be
    params = {
        'A': A,
        'Iy': Iy,
        'Iz': Iz,
        'J': J,
        'nodes': nodes,
        'loads': loads
    }
    unknowns = {
        'K': np.zeros((comp.size, comp.size), dtype=data_type),
        'forces': np.zeros((comp.size), dtype=data_type)
    }
    resids = None
    comp.solve_nonlinear(params, unknowns, resids)
    K = unknowns.get('K')
    forces = unknowns.get('forces')
    return K, forces


"""
<<<<<<< HEAD
--------------------------------------------------------------------------------

                                MATERIALS

--------------------------------------------------------------------------------
From materials.py: """


def materials_tube(r, thickness, surface=None, comp=None):
=======
================================================================================
                                MATERIALS
================================================================================
From materials.py: """


def materials_tube(r, thickness, comp):
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    """ Compute geometric properties for a tube element.

    Parameters
    ----------
    r : array_like
        Radii for each FEM element.
    thickness : array_like
        Tube thickness for each FEM element.
<<<<<<< HEAD
=======
    comp : Either OpenAeroStruct component object (better), or surface dict.
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3

    Returns
    -------
    A : array_like
        Areas for each FEM element.
    Iy : array_like
<<<<<<< HEAD
        Mass moment of inertia around the y-axis for each FEM element.
    Iz : array_like
        Mass moment of inertia around the z-axis for each FEM element.
=======
        Area moment of inertia around the y-axis for each FEM element.
    Iz : array_like
        Area moment of inertia around the z-axis for each FEM element.
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    J : array_like
        Polar moment of inertia for each FEM element.

    """
<<<<<<< HEAD
    if not comp:
        comp=MaterialsTube(surface)
    # if not r:
    #     r = surface['r']  # this is already contained in surface dict
    # if not thickness:
    #     thickness = surface['t']  # this is already contained in surface dict
    params={
        'r': r,
=======
    if not isinstance(comp, Component):
        surface = comp
        comp=MaterialsTube(surface)
    # if not r:
    #     r = surface['radius']  # this is already contained in surface dict
    # if not thickness:
    #     thickness = surface['thickness']  # this is already contained in surface dict
    params={
        'radius': r,
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
        'thickness': thickness
    }
    unknowns={
        'A': np.zeros((comp.ny - 1)),
        'Iy': np.zeros((comp.ny - 1)),
        'Iz': np.zeros((comp.ny - 1)),
        'J': np.zeros((comp.ny - 1))
    }
    resids = None
    comp.solve_nonlinear(params, unknowns, resids)
    A=unknowns.get('A')
    Iy=unknowns.get('Iy')
    Iz=unknowns.get('Iz')
    J=unknowns.get('J')
    return A, Iy, Iz, J

    """
<<<<<<< HEAD
    --------------------------------------------------------------------------------

                                    FUNCTIONALS

    --------------------------------------------------------------------------------
=======
================================================================================
                                FUNCTIONALS
================================================================================
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
    From functionals.py:

        to be added here...

        """


if __name__ == "__main__":
    ''' Test the coupled system with default parameters

     To change problem parameters, input the prob_dict dictionary, e.g.
     prob_dict = {
        'rho' : 0.35,
<<<<<<< HEAD
        'R': 14.0e6
=======
        'thickness': 14.0e6
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
     }
    '''
    print('Fortran Flag = {0}'.format(fortran_flag))
    # Define parameters
    prob_dict = {} # use default

    # Define surface
    surface = {
        'wing_type' : 'CRM',
        'num_x': 2,   # number of chordwise points
        'num_y': 9    # number of spanwise points
    }

    # Define fixed point iteration options
    # default options from OpenMDAO nonlinear solver NLGaussSeidel
    fpi_opt = {
        # 'atol': float(1e-06),   # Absolute convergence tolerance (unused)
        # 'err_on_maxiter': bool(False),  # raise AnalysisError if not converged at maxiter (unused)
        # 'print': int(0),        # Print option (unused)
        'maxiter': int(100),    # Maximum number of iterations
        # 'rtol': float(1e-06),   # Relative convergence tolerance (unused)
        'utol': float(1e-12)    # Convergence tolerance on the change in the unknowns
    }

    print('Run analysis.setup()')
    OAS_prob = setup(prob_dict=prob_dict, surfaces=[surface])
    # print('OAS_prob.surfaces = ')
    # print(OAS_prob.surfaces)
    # print('OAS_prob.prob_dict = ')
    # print(OAS_prob.prob_dict)
    # print('OAS_prob.comp_dict = ')
    # print(OAS_prob.comp_dict)

    print('Run coupled system analysis with fixed point iteration')

    # Make local functions for coupled system analysis
    def f_aero(def_mesh):
<<<<<<< HEAD
        loads = aerodynamics(def_mesh, OAS_prob.surfaces[0], OAS_prob.prob_dict, OAS_prob.comp_dict)
        return loads
    def f_struct(loads):
        def_mesh = structures(loads, OAS_prob.surfaces[0], OAS_prob.prob_dict, OAS_prob.comp_dict)
=======
        # loads = aerodynamics(def_mesh, OAS_prob.surfaces[0], OAS_prob.prob_dict, OAS_prob.comp_dict)
        loads = aerodynamics2(def_mesh, OAS_prob.surfaces[0], OAS_prob.prob_dict)
        return loads
    def f_struct(loads):
        # def_mesh = structures(loads, OAS_prob.surfaces[0], OAS_prob.prob_dict, OAS_prob.comp_dict)
        def_mesh = structures2(loads, OAS_prob.surfaces[0], OAS_prob.prob_dict)
>>>>>>> 7eefd15e6c26c95cbbd9f303d23ecb4716c2fea3
        return def_mesh

    # Define FPI parameters
    utol = fpi_opt['utol']
    maxiter = fpi_opt['maxiter']
    # Generate initial mesh with zero deformation
    def_mesh = gen_init_mesh(OAS_prob.surfaces[0], OAS_prob.comp_dict)
    x0 = def_mesh
    u_norm = 1.0e99
    iter_count = 0
    # Run fixed point iteration on coupled aerodynamics-structures system
    while (iter_count < maxiter) and (u_norm > utol):
          # Update iteration counter
          iter_count += 1
          # Run iteration and evaluate norm of residual
          loads = f_aero(x0)
          def_mesh = x = f_struct(loads)
          u_norm = np.linalg.norm(x - x0)
          x0 = x

    if iter_count >= maxiter or math.isnan(u_norm):
        msg = 'FAILED to converge after {0:d} iterations'.format(iter_count)
    else:
        msg = 'Converged in {0:d} iterations'.format(iter_count)

    print(msg)
    print('def_mesh=\n',def_mesh)
    print('loads=\n',loads)
