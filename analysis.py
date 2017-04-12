# Univ. Michigan Aerostructural model.
# Based on OpenAeroStruct by John Hwang, and John Jasa (github.com/mdolab/OpenAeroStruct)
# author: Sam Friedman  (samfriedman@tamu.edu)
# date:   4/12/2017

# make compatible Python 2.x to 3.x
from __future__ import print_function, division
# from future.builtins import range  # make compatible Python 2.x to 3.x
# __all__ = ['setup','aerodynamics','structures']

import numpy as np

from materials import MaterialsTube
from spatialbeam import ComputeNodes, AssembleK, SpatialBeamFEM, SpatialBeamDisp#, SpatialBeamEnergy, SpatialBeamWeight, SpatialBeamVonMisesTube, SpatialBeamFailureKS
from transfer import TransferDisplacements, TransferLoads
from vlm import VLMGeometry, AssembleAIC, AeroCirculations, VLMForces#, VLMLiftDrag, VLMCoeffs, TotalLift, TotalDrag
from geometry import GeometryMesh, Bspline#, gen_crm_mesh, gen_rect_mesh, MonotonicConstraint
from run_classes import OASProblem
# from functionals import FunctionalBreguetRange, FunctionalEquilibrium

# to disable OpenMDAO warnings which will create an error in Matlab
import warnings
warnings.filterwarnings("ignore")

try:
    import OAS_API
    fortran_flag = True
    data_type = float
except:
    fortran_flag = False
    data_type = complex

"""
--------------------------------------------------------------------------------

                            GEOMETRY / SETUP

--------------------------------------------------------------------------------
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
        stress = 20.e6       # [Pa] yield stress
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
    '''
    Output after calling OAS_prob.prob_dict:
    In [12]: OAS_prob.prob_dict
    Out[12]:
    {'CT': 0.00016671305,
     'M': 0.84,
     'R': 14300000.0,
     'Re': 1000000.0,
     'a': 295.4,
     'alpha': 5.0,
     'force_fd': False,
     'optimize': False,
     'optimizer': 'SNOPT',
     'print_level': 0,
     'reynolds_length': 1.0,
     'rho': 0.38,
     'type': 'aerostruct',
     'v': 248.13599999999997,
     'with_viscous': False}
    '''

    for surface in surfaces:
        OAS_prob.add_surface(surface)    # Add the specified wing surface to the problem.
    
    # Add materials properties for the wing surface to the surface dict in OAS_prob
    for idx, surface in enumerate(OAS_prob.surfaces):
        A, Iy, Iz, J = materials_tube(surface['r'], surface['t'], surface)
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

    '''
    Extract parameters and variables from OAS_prob to pass through to
    other discipline functions. For now, assume we are only using one lifting
    surface.

    Output after calling OAS_prob.add_surface(surf_dict):
    In [8]: OAS_prob.surfaces
    Out[8]:
    [{'CD0': 0.015,
    'CL0': 0.2,
    'E': 70000000000.0,
    'G': 30000000000.0,
    'S_ref_type': 'wetted',
    'W0': 120000.0,
    'active_bsp_vars': ['thickness_cp',
     'twist_cp',
     'xshear_cp',
     'chord_cp',
     'zshear_cp'],
    'active_geo_vars': ['sweep',
     'dihedral',
     'twist_cp',
     'xshear_cp',
     'zshear_cp',
     'span',
     'chord_cp',
     'taper',
     'thickness_cp'],
    'c_max_t': 0.303,
    'chord_cos_spacing': 0.0,
    'chord_cp': array([ 1.+0.j,  1.+0.j,  1.+0.j,  1.+0.j,  1.+0.j]),
    'crm_twist': array([-3.75  , -3.1248, -2.5773, -2.2772, -2.0301, -1.8158, -1.635 ,
         -1.4526, -1.2067, -0.9436, -0.6782, -0.2621,  0.4285,  0.9379,
          1.5252,  2.2419,  2.2419,  3.6063,  4.4402,  6.7166]),
    'dihedral': 0.0,
    'exact_failure_constraint': False,
    'fem_origin': 0.35,
    'k_lam': 0.05,
    'loads': array([[    0.+0.j,     0.+0.j,  1000.+0.j,     0.+0.j,     0.+0.j,
              0.+0.j],
         [    0.+0.j,     0.+0.j,     0.+0.j,     0.+0.j,     0.+0.j,
              0.+0.j],
         [    0.+0.j,     0.+0.j,     0.+0.j,     0.+0.j,     0.+0.j,
              0.+0.j],
         [    0.+0.j,     0.+0.j,     0.+0.j,     0.+0.j,     0.+0.j,
              0.+0.j]]),
    'mesh': array([[[  4.52307198e+01+0.j,  -2.93815262e+01+0.j,   6.70120580e+00+0.j],
          [  4.22333963e+01+0.j,  -2.54451497e+01+0.j,   6.02337146e+00+0.j],
          [  3.40443566e+01+0.j,  -1.46907758e+01+0.j,   4.78508060e+00+0.j],
          [  2.29690676e+01+0.j,  -1.79909493e-15+0.j,   4.42280040e+00+0.j]],

         [[  4.79586798e+01+0.j,  -2.93815262e+01+0.j,   6.70120580e+00+0.j],
          [  4.59248861e+01+0.j,  -2.54451497e+01+0.j,   6.02337146e+00+0.j],
          [  4.03682708e+01+0.j,  -1.46907758e+01+0.j,   4.78508060e+00+0.j],
          [  3.65880650e+01+0.j,  -1.79909493e-15+0.j,   4.42280040e+00+0.j]]]),
    'monotonic_con': None,
    'mrho': 3000.0,
    'name': 'wing_',
    'num_chord_cp': 5,
    'num_thickness_cp': 5,
    'num_twist_cp': 5,
    'num_x': 2L,
    'num_xshear_cp': 5,
    'num_y': 4,
    'num_zshear_cp': 5,
    'offset': array([ 0.,  0.,  0.]),
    'r': array([ 0.48145874+0.j,  0.75115530+0.j,  1.49571837+0.j]),
    'root_chord': 1.0,
    'span': 58.763052399999985,
    'span_cos_spacing': 1,
    'stress': 20000000.0,
    'sweep': 0.0,
    'symmetry': True,
    't': array([ 0.04814587+0.j,  0.07511553+0.j,  0.14957184+0.j]),
    't_over_c': 0.12,
    'taper': 1.0,
    'thickness_cp': array([ 0.14957184+0.j,  0.14957184+0.j,  0.14957184+0.j,  0.14957184+0.j,
          0.14957184+0.j]),
    'twist_cp': array([-3.75  , -2.0301, -0.9436,  1.5252,  6.7166]),
    'wing_type': 'CRM',
    'xshear_cp': array([ 0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j]),
    'zshear_cp': array([ 0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j])}]
    '''
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
    comp_dict['SpatialBeamFEM'] = SpatialBeamFEM(OAS_prob.prob_dict['tot_panels'])
    comp_dict['SpatialBeamDisp'] = SpatialBeamDisp(surface)
    OAS_prob.comp_dict = comp_dict

    # return the surfaces list, problem dict, and component dict
    surfaces = [surface]
    prob_dict = OAS_prob.prob_dict
    # return surfaces, prob_dict, comp_dict
    return OAS_prob



def gen_init_mesh(surface, comp_dict=None):
    ''' Generate initial def_mesh '''
    if comp_dict:
        mesh = geometry_mesh(surface, comp=comp_dict['GeometryMesh'])
        disp = np.zeros((surface['num_y'], 6), dtype=data_type)  # zero displacement
        def_mesh = transfer_displacements(mesh, disp, comp=comp_dict['TransferDisplacements'])
    else:
        mesh = geometry_mesh(surface)
        disp = np.zeros((surface['num_y'], 6), dtype=data_type)  # zero displacement
        def_mesh = transfer_displacements(mesh, disp, surface)
    return def_mesh


def aerodynamics(def_mesh, surface, prob_dict, comp_dict=None):
    ''' Use pre-initialized components '''

    # Unpack variables
    v = prob_dict.get('v')
    alpha = prob_dict.get('alpha')
    size = prob_dict.get('tot_panels')
    rho = prob_dict.get('rho')

    b_pts, c_pts, widths, cos_sweep, lengths, normals, S_ref = vlm_geometry(def_mesh, comp=comp_dict['VLMGeometry'])
    AIC, rhs= assemble_aic(surface, def_mesh, b_pts, c_pts, normals, v, alpha, comp=comp_dict['AssembleAIC'])
    circulations = aero_circulations(AIC, rhs, comp=comp_dict['AeroCirculations'])
    sec_forces = vlm_forces(surface, def_mesh, b_pts, circulations, alpha, v, rho, comp=comp_dict['VLMForces'])
    loads = transfer_loads(def_mesh, sec_forces, comp=comp_dict['TransferLoads'])

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


def structures(loads, surface, prob_dict, comp_dict=None):
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

    nodes = compute_nodes(mesh, comp=comp_dict['ComputeNodes'])
    K, forces = assemble_k(A, Iy, Iz, J, nodes, loads, comp=comp_dict['AssembleK'])
    disp_aug = spatial_beam_fem(K, forces, comp=comp_dict['SpatialBeamFEM'])
    disp = spatial_beam_disp(disp_aug, comp=comp_dict['SpatialBeamDisp'])
    def_mesh = transfer_displacements(mesh, disp, comp=comp_dict['TransferDisplacements'])

    return def_mesh  # Output the def_mesh matrix
    
    
def structures2(loads, surface, prob_dict):
    ''' Don't use pre-initialized components '''

    # Unpack variables
    A = surface.get('A')
    Iy = surface.get('Iy')
    Iz = surface.get('Iz')
    J = surface.get('J')
    mesh = surface.get('mesh')
    v = prob_dict.get('v')
    alpha = prob_dict.get('alpha')
    size =  prob_dict.get('tot_panels')

    nodes = compute_nodes(mesh, surface)
    K, forces = assemble_k(A, Iy, Iz, J, nodes, loads, surface)
    disp_aug = spatial_beam_fem(K, forces, size)
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
        if var in surface['active_geo_vars']:
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


def transfer_displacements(mesh, disp, surface=None, comp=None):
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

    Returns
    -------
    def_mesh[nx, ny, 3] : numpy array
        Flattened array defining the lifting surfaces after deformation.
    """
    if not comp:
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
--------------------------------------------------------------------------------

                                AERODYNAMICS

--------------------------------------------------------------------------------
From vlm.py: """


def vlm_geometry(def_mesh, surface=None, comp=None):
    """ Compute various geometric properties for VLM analysis.

    Parameters
    ----------
    def_mesh[nx, ny, 3] : numpy array
        Array defining the nodal coordinates of the lifting surface.

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
    if not comp:
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


def aero_circulations(AIC, rhs, size=None, comp=None):
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

    Returns
    -------
    circulations[(nx-1)*(ny-1)] : numpy array
        Augmented displacement array. Obtained by solving the system
        AIC * circulations = n * v.
    """
    if not comp:
        comp = AeroCirculations(size)
    if not size:
        size = comp.size
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


def transfer_loads(def_mesh, sec_forces, surface=None, comp=None):
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

    Returns
    -------
    loads[ny, 6] : numpy array
        Flattened array containing the loads applied on the FEM component,
        computed from the sectional forces.
    """
    if not comp:
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
--------------------------------------------------------------------------------

                                   STRUCTURES

--------------------------------------------------------------------------------
From spatialbeam.py: Define the structural analysis component using spatial beam theory. """

def spatial_beam_fem(K, forces, size=None, comp=None):
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

    Returns
    -------
    disp_aug[6*(ny+1)] : numpy array
        Augmented displacement array. Obtained by solving the system
        K * u = f, where f is a flattened version of loads.

    """
    if not comp:
        comp=SpatialBeamFEM(size)
    if not size:
        size = comp.size
    params={
        'K': K,
        'forces': forces
    }
    unknowns={
        'disp_aug': np.zeros((size), dtype=data_type)
    }
    resids={
        'disp_aug': np.zeros((size), dtype=data_type)
    }
    comp.solve_nonlinear(params, unknowns, resids)
    disp_aug=unknowns.get('disp_aug')
    return disp_aug


def spatial_beam_disp(disp_aug, surface=None, comp=None):
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

    Returns
    -------
    disp[6*ny] : numpy array
        Actual displacement array formed by truncating disp_aug.

    """
    if not comp:
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


def compute_nodes(mesh, surface=None, comp=None):
    """
    Compute FEM nodes based on aerodynamic mesh.

    The FEM nodes are placed at fem_origin * chord,
    with the default fem_origin = 0.35.

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Array defining the nodal points of the lifting surface.

    Returns
    -------
    nodes[ny, 3] : numpy array
        Flattened array with coordinates for each FEM node.

    """
    if not comp:
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


def assemble_k(A, Iy, Iz, J, nodes, loads, surface=None, comp=None):
    """
    Compute the displacements and rotations by solving the linear system
    using the structural stiffness matrix.

    Parameters
    ----------
    A[ny-1] : numpy array
        Areas for each FEM element.
    Iy[ny-1] : numpy array
        Mass moment of inertia around the y-axis for each FEM element.
    Iz[ny-1] : numpy array
        Mass moment of inertia around the z-axis for each FEM element.
    J[ny-1] : numpy array
        Polar moment of inertia for each FEM element.
    nodes[ny, 3] : numpy array
        Flattened array with coordinates for each FEM node.
    loads[ny, 6] : numpy array
        Flattened array containing the loads applied on the FEM component,
        computed from the sectional forces.

    Returns
    -------
    K[(nx-1)*(ny-1), (nx-1)*(ny-1)] : numpy array
        Stiffness matrix for the entire FEM system. Used to solve the linear
        system K * u = f to obtain the displacements, u.
    forces[(nx-1)*(ny-1)] : numpy array
        Right-hand-side of the linear system. The loads from the aerodynamic
        analysis or the user-defined loads.
    """
    if not comp:
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
--------------------------------------------------------------------------------

                                MATERIALS

--------------------------------------------------------------------------------
From materials.py: """


def materials_tube(r, thickness, surface=None, comp=None):
    """ Compute geometric properties for a tube element.

    Parameters
    ----------
    r : array_like
        Radii for each FEM element.
    thickness : array_like
        Tube thickness for each FEM element.

    Returns
    -------
    A : array_like
        Areas for each FEM element.
    Iy : array_like
        Mass moment of inertia around the y-axis for each FEM element.
    Iz : array_like
        Mass moment of inertia around the z-axis for each FEM element.
    J : array_like
        Polar moment of inertia for each FEM element.

    """
    if not comp:
        comp=MaterialsTube(surface)
    # if not r:
    #     r = surface['r']  # this is already contained in surface dict
    # if not thickness:
    #     thickness = surface['t']  # this is already contained in surface dict
    params={
        'r': r,
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
    --------------------------------------------------------------------------------

                                    FUNCTIONALS

    --------------------------------------------------------------------------------
    From functionals.py: 
        
        to be added here...
        
        """


if __name__ == "__main__":
    ''' Test the coupled system with default parameters 
    
     To change problem parameters, input the prob_dict dictionary, e.g.
     prob_dict = {
        'rho' : 0.35,
        'R': 14.0e6
     }
    '''
    print('Fortran Flag = {0}'.format(fortran_flag))
    
    print('Run analysis.setup()...')
    surface = {
        'wing_type' : 'CRM',
        'num_x': 2,   # number of chordwise points
        'num_y': 7    # number of spanwise points
    }
    OAS_prob = setup(prob_dict={}, surfaces=[surface])
    # print('OAS_prob.surfaces = ')
    # print(OAS_prob.surfaces)
    # print('OAS_prob.prob_dict = ')
    # print(OAS_prob.prob_dict)
    # print('OAS_prob.comp_dict = ')
    # print(OAS_prob.comp_dict)
    
    def_mesh = gen_init_mesh(OAS_prob.surfaces[0], OAS_prob.comp_dict)
    print('def_mesh = ')
    print(def_mesh)

    print('\nRun analysis.aerodynamics()...')
    loads = aerodynamics(def_mesh, OAS_prob.surfaces[0], OAS_prob.prob_dict, OAS_prob.comp_dict)
    print('loads = ')
    print(loads)
    #
    print('\nRun analysis.structures()...')
    def_mesh = structures(loads, OAS_prob.surfaces[0], OAS_prob.prob_dict, OAS_prob.comp_dict)
    print('def_mesh = ')
    print(def_mesh)