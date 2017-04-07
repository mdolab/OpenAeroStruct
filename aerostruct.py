# Univ. Michigan Aerostructural model.
# Based on OpenAeroStruct by John Hwang, and John Jasa (github.com/johnjasa/OpenAeroStruct)
# author: Sam Friedman  (samfriedman@tamu.edu)
# date:   4/6/2017

# make compatible Python 2.x to 3.x
from __future__ import print_function, division
# from future.builtins import range  # make compatible Python 2.x to 3.x
# __all__ = ['setup','aerodynamics','structures']
import warnings
import sys
import os
import numpy as np
import scipy.sparse
from scipy.linalg import lu_factor, lu_solve

from materials import MaterialsTube
from spatialbeam import ComputeNodes, SpatialBeamFEM, SpatialBeamDisp, SpatialBeamEnergy, SpatialBeamWeight, SpatialBeamVonMisesTube, SpatialBeamFailureKS
from transfer import TransferDisplacements, TransferLoads
from vlm import VLMGeometry, AssembleAIC, AeroCirculations, VLMForces, VLMLiftDrag, VLMCoeffs, TotalLift, TotalDrag
from b_spline import get_bspline_mtx
# from geometry import get_inds, rotate, sweep, dihedral, stretch, taper, mirror
from geometry import GeometryMesh, Bspline, gen_crm_mesh, gen_rect_mesh, MonotonicConstraint
from functionals import FunctionalBreguetRange, FunctionalEquilibrium
from OpenAeroStruct import OASProblem

try:
    import OAS_API
    fortran_flag = True
    data_type = float
except:
    fortran_flag = False
    data_type = complex

# to disable OpenMDAO warnings which will create an error in Matlab
warnings.filterwarnings("ignore")

# In Matlab code, add this before calling Python functions:
#   if count(py.sys.path,'') == 0
#       insert(py.sys.path,int32(0),'');
#   end

"""
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

                                GEOMETRY / SETUP

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
From run_classes.py: Manipulate geometry mesh based on high-level design parameters """


def setup(num_x=2, num_y=7):
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
    prob_dict = {'type' : 'aerostruct'}

    # To update problem parameters, update the prob_dict dictionary.
    # The dictionary key is a string, e.g.
    #   prob_dict.update({'rho' : 0.35,
    #                     'R': 14.0e6
    #   })

    # Instantiate problem
    OAS_prob = OASProblem(prob_dict)

    # Create a dictionary to store options about the wing surface
    surf_dict = {'name' : 'wing',
                 'symmetry' : True,
                 'num_y' : num_y,       # from input parameters
                 'num_x' : num_x,       # from input parameters
                 'wing_type' : 'CRM',
                 'CL0' : 0.2,
                 'CD0' : 0.015,
                #  'span_cos_spacing' : 1.,
                #  'chord_cos_spacing' : .8
                 }
    # Add the specified wing surface to the problem.
    OAS_prob.add_surface(surf_dict)

    '''
    Extract parameters and variables from OAS_prob to pass through to
    other discipline functions. For now, assume we are only using one lifting
    surface and hardcode the variable names for the wing. Later, I will create
    a class object to hold all of the surface variables.

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
    wing = surfaces[0]

    # wing_CD0 = wing.get('CD0')
    # wing_CL0 = wing.get('CL0')
    # wing_E = wing.get('E')
    # wing_G = wing.get('G')
    # wing_S_ref_type = wing.get('S_ref_type')
    # wing_W0 = wing.get('W0')
    # wing_active_bsp_vars = wing.get('active_bsp_vars')
    # wing_active_geo_vars = wing.get('active_geo_vars')
    # wing_c_max_t = wing.get('c_max_t')
    # wing_chord_cos_spacing = wing.get('chord_cos_spacing')
    # wing_chord_cp = wing.get('chord_cp')
    # wing_crm_twist = wing.get('crm_twist')
    # wing_dihedral = wing.get('dihedral')
    # wing_exact_failure_constraint = wing.get('exact_failure_constraint')
    # wing_fem_origin = wing.get('fem_origin')
    # wing_k_lam = wing.get('k_lam')
    # wing_loads = wing.get('loads')
    # wing_mesh = wing.get('mesh')
    # wing_monotonic_con = wing.get('monotonic_con')
    # wing_mrho = wing.get('mrho')
    # wing_name = wing.get('name')
    # wing_num_chord_cp = wing.get('num_chord_cp')
    # wing_num_thickness_cp = wing.get('num_thickness_cp')
    # wing_num_twist_cp = wing.get('num_twist_cp')
    # wing_num_x=wing.get('num_x')
    # wing_num_xshear_cp = wing.get('num_xshear_cp')
    # wing_num_y = wing.get('num_y')
    # wing_num_zshear_cp = wing.get('num_zshear_cp')
    # wing_offset = wing.get('offset')
    # wing_r = wing.get('r')
    # wing_root_chord = wing.get('root_chord')
    # wing_span = wing.get('span')
    # wing_span_cos_spacing = wing.get('span_cos_spacing')
    # wing_stress = wing.get('stress')
    # wing_sweep = wing.get('sweep')
    # wing_symmetry = wing.get('symmetry')
    # wing_t = wing.get('t')
    # wing_t_over_c = wing.get('t_over_c')
    # wing_taper = wing.get('taper')
    # wing_thickness_cp = wing.get('thickness_cp')
    # wing_twist_cp = wing.get('twist_cp')
    # wing_wing_type = wing.get('wing_type')
    # wing_xshear_cp = wing.get('xshear_cp')
    # wing_zshear_cp = wing.get('zshear_cp')

    # or just return the surface dict
    return wing


def aerodynamics(def_mesh=None, params=None):
    # Unpack variables
    aero_ind = params.get('aero_ind')
    alpha = params.get('alpha')
    v = params.get('v')
    rho = params.get('rho')
    fem_ind = params.get('fem_ind')
    fem_origin = params.get('fem_origin', 0.35)

    # num_y = params.get('num_y')
    # span = params.get('span')
    # twist_cp = params.get('twist_cp')
    # thickness_cp = params.get('thickness_cp')
    # num_thickness = params.get('num_thickness')
    # num_twist = params.get('num_twist')
    # sweep = params.get('sweep')
    # taper = params.get('taper')
    # disp = params.get('disp')
    # dihedral = params.get('dihedral')

    # # Define Jacobians for b-spline controls
    # tot_n_fem = np.sum(fem_ind[:, 0])
    # num_surf = fem_ind.shape[0]
    # jac_twist = get_bspline_mtx(num_twist, num_y)
    # jac_thickness = get_bspline_mtx(num_thickness, tot_n_fem - num_surf)

    b_pts, mid_b, c_pts, widths, normals, S_ref = vlm_geometry(
        aero_ind, def_mesh)
    circulations = vlm_circulations(
        aero_ind, def_mesh, b_pts, c_pts, normals, v, alpha)

    sec_forces = vlm_forces(def_mesh, aero_ind, b_pts,
                            mid_b, circulations, alpha, v, rho)
    loads = transfer_loads(def_mesh, sec_forces, aero_ind, fem_ind, fem_origin)
    # for i, j in enumerate([circulations, sec_forces, loads]):
    # print('shape of ',i,' = ',j.shape)
    return loads


def structures(loads, params):
    # Unpack variables
    mesh = params.get('mesh')
    twist_cp = params.get('twist_cp')
    thickness_cp = params.get('thickness_cp')
    r = params.get('r')
    aero_ind = params.get('aero_ind')
    fem_ind = params.get('fem_ind')
    E = params.get('E')
    G = params.get('G')
    jac_twist = params.get('jac_twist')
    jac_thickness = params.get('jac_thickness')
    fem_origin = params.get('fem_origin', 0.35)
    cg_x = params.get('cg_x', 5)

    twist = cp2pt(twist_cp, jac_twist)
    thickness = cp2pt(thickness_cp, jac_thickness)
    geometry_mesh(mesh, aero_ind, twist)
    A, Iy, Iz, J = materials_tube(r, thickness, fem_ind)
    nodes = compute_nodes(mesh, fem_ind, aero_ind, fem_origin)
    disp_aug = spatial_beam_FEM(
        A, Iy, Iz, J, nodes, loads, aero_ind, fem_ind, E, G, cg_x)
    disp = spatial_beam_disp(disp_aug, fem_ind)
    def_mesh = transfer_displacements(
        mesh, disp, aero_ind, fem_ind, fem_origin)

    return def_mesh  # Output the def_mesh matrix


def cp2pt(cp, jac):
    """
    General function to translate from control points to actual points
    using a b-spline representation.
    """
    pt = np.zeros(jac.shape[0])
    pt = jac.dot(cp)
    return pt


def gen_crm_mesh(n_points_inboard=3, n_points_outboard=4,
                 num_x=2):
    """
    Build the right hand side of the CRM wing with specified number
    of inboard and outboard panels.

    n_points_inboard : int
        Number of spanwise points between the wing root and yehudi break per
        wing side.
    n_points_outboard : int
        Number of spanwise points between the yehudi break and wingtip per
        wing side.
    num_x : int
        Number of chordwise points.

    """
    #   crm base mesh from crm_data.py
    # eta, xle, yle, zle, twist, chord
    raw_crm_points = np.array([
        [0., 904.294, 0.0, 174.126, 6.7166, 536.181],  # 0
        [.1, 989.505, 115.675, 175.722, 4.4402, 468.511],
        [.15, 1032.133, 173.513, 176.834, 3.6063, 434.764],
        [.2, 1076.030, 231.351, 177.912, 2.2419, 400.835],
        [.25, 1120.128, 289.188, 177.912, 2.2419, 366.996],
        [.3, 1164.153, 347.026, 178.886, 1.5252, 333.157],
        [.35, 1208.203, 404.864, 180.359, .9379, 299.317],  # 6 yehudi break
        [.4, 1252.246, 462.701, 182.289, .4285, 277.288],
        [.45, 1296.289, 520.539, 184.904, -.2621, 263],
        [.5, 1340.329, 578.377, 188.389, -.6782, 248.973],
        [.55, 1384.375, 636.214, 192.736, -.9436, 234.816],
        [.60, 1428.416, 694.052, 197.689, -1.2067, 220.658],
        [.65, 1472.458, 751.890, 203.294, -1.4526, 206.501],
        [.7, 1516.504, 809.727, 209.794, -1.6350, 192.344],
        [.75, 1560.544, 867.565, 217.084, -1.8158, 178.186],
        [.8, 1604.576, 925.402, 225.188, -2.0301, 164.029],
        [.85, 1648.616, 983.240, 234.082, -2.2772, 149.872],
        [.9, 1692.659, 1041.078, 243.625, -2.5773, 135.714],
        [.95, 1736.710, 1098.915, 253.691, -3.1248, 121.557],
        [1., 1780.737, 1156.753, 263.827, -3.75, 107.4]  # 19
    ])
    # le = np.vstack((raw_crm_points[:,1],
    #                 raw_crm_points[:,2],
    #                 raw_crm_points[:,3]))
    # te = np.vstack((raw_crm_points[:,1]+raw_crm_points[:,5],
    #                 raw_crm_points[:,2],
    #                 raw_crm_points[:,3]))
    # mesh = np.empty((2,20,3))
    # mesh[0,:,:] = le.T
    # mesh[1,:,:] = te.T
    # mesh *= 0.0254 # convert to meters
    # pull out the 3 key y-locations to define the two linear regions of the
    # wing
    crm_base_points = raw_crm_points[(0, 6, 19), :]
    le_base = np.vstack((crm_base_points[:, 1],
                         crm_base_points[:, 2],
                         crm_base_points[:, 3]))
    te_base = np.vstack((crm_base_points[:, 1] + crm_base_points[:, 5],
                         crm_base_points[:, 2],
                         crm_base_points[:, 3]))
    mesh = np.empty((2, 3, 3))
    mesh[0, :, :] = le_base.T
    mesh[1, :, :] = te_base.T
    mesh[:, :, 2] = 0  # get rid of the z deflection
    mesh *= 0.0254  # convert to meters
    # LE pre-yehudi
    s1 = (mesh[0, 1, 0] - mesh[0, 0, 0]) / (mesh[0, 1, 1] - mesh[0, 0, 1])
    o1 = mesh[0, 0, 0]
    # TE pre-yehudi
    s2 = (mesh[1, 1, 0] - mesh[1, 0, 0]) / (mesh[1, 1, 1] - mesh[1, 0, 1])
    o2 = mesh[1, 0, 0]
    # LE post-yehudi
    s3 = (mesh[0, 2, 0] - mesh[0, 1, 0]) / (mesh[0, 2, 1] - mesh[0, 1, 1])
    o3 = mesh[0, 2, 0] - s3 * mesh[0, 2, 1]
    # TE post-yehudi
    s4 = (mesh[1, 2, 0] - mesh[1, 1, 0]) / (mesh[1, 2, 1] - mesh[1, 1, 1])
    o4 = mesh[1, 2, 0] - s4 * mesh[1, 2, 1]
    n_points_total = n_points_inboard + n_points_outboard - 1
    half_mesh = np.zeros((2, n_points_total, 3))
    # generate inboard points
    dy = (mesh[0, 1, 1] - mesh[0, 0, 1]) / (n_points_inboard - 1)
    for i in range(n_points_inboard):
        y = half_mesh[0, i, 1] = i * dy
        half_mesh[0, i, 0] = s1 * y + o1  # le point
        half_mesh[1, i, 1] = y
        half_mesh[1, i, 0] = s2 * y + o2  # te point
    yehudi_break = mesh[0, 1, 1]
    # generate outboard points
    dy = (mesh[0, 2, 1] - mesh[0, 1, 1]) / (n_points_outboard - 1)
    for j in range(n_points_outboard):
        i = j + n_points_inboard - 1
        y = half_mesh[0, i, 1] = j * dy + yehudi_break
        half_mesh[0, i, 0] = s3 * y + o3  # le point
        half_mesh[1, i, 1] = y
        half_mesh[1, i, 0] = s4 * y + o4  # te point
    full_mesh = mirror(half_mesh)
    full_mesh = add_chordwise_panels(full_mesh, num_x)
    full_mesh[:, :, 1] -= np.mean(full_mesh[:, :, 1])
    return full_mesh


def add_chordwise_panels(mesh, num_x):
    """ Divide the wing into multiple chordwise panels. """
    le = mesh[0, :, :]
    te = mesh[-1, :, :]
    new_mesh = np.zeros((num_x, mesh.shape[1], 3))
    new_mesh[0, :, :] = le
    new_mesh[-1, :, :] = te
    for i in range(1, num_x - 1):
        w = float(i) / (num_x - 1)
        new_mesh[i, :, :] = (1 - w) * le + w * te
    return new_mesh


def gen_mesh(num_x, num_y, span, chord, cosine_spacing=0.):
    """ Generate simple rectangular wing mesh. """
    mesh = np.zeros((num_x, num_y, 3))
    ny2 = (num_y + 1) / 2
    beta = np.linspace(0, np.pi / 2, ny2)
    # mixed spacing with w as a weighting factor
    cosine = .5 * np.cos(beta)  # cosine spacing
    uniform = np.linspace(0, .5, ny2)[::-1]  # uniform spacing
    half_wing = cosine * cosine_spacing + (1 - cosine_spacing) * uniform
    full_wing = np.hstack((-half_wing[:-1], half_wing[::-1])) * span
    for ind_x in range(num_x):
        for ind_y in range(num_y):
            mesh[ind_x, ind_y, :] = [ind_x / (num_x - 1) * chord,
                                     full_wing[ind_y], 0]
    return mesh


def geometry_mesh(mesh, aero_ind, twist, sweep_angle=0, dihedral_angle=0, taper_ratio=1, span=58.7630524):
    ny = aero_ind[0, 1]
    nx = aero_ind[0, 0]
    n = nx * ny
    wing_mesh = mesh[: n, :].reshape(nx, ny, 3).astype('complex')
    twist = np.zeros(ny)  # default option
    wing_mesh = mesh[: n, :].reshape(nx, ny, 3).astype('complex')
    # stretch(wing_mesh, params['span'])
    sweep(wing_mesh, sweep_angle)
    rotate(wing_mesh, twist)
    dihedral(wing_mesh, dihedral_angle)
    taper(wing_mesh, taper_ratio)
    mesh[: n, :] = wing_mesh.reshape(n, 3).astype('complex')
    return mesh


def transfer_displacements(mesh, disp, aero_ind, fem_ind, fem_origin=0.35):
    """
    Perform displacement transfer.

    Apply the computed displacements on the original mesh to obtain
    the deformed mesh.

    Parameters
    ----------
    mesh : array_like
        Flattened array defining the lifting surfaces.
    disp : array_like
        Flattened array containing displacements on the FEM component.
        Contains displacements for all six degrees of freedom, including
        displacements in the x, y, and z directions, and rotations about the
        x, y, and z axes.

    Returns
    -------
    def_mesh : array_like
        Flattened array defining the lifting surfaces after deformation.

    """
    _Component = TransferDisplacements(aero_ind, fem_ind, fem_origin)
    params = {
        'mesh': mesh,
        'disp': disp
    }
    unknowns = {
        'def_mesh': np.zeros((np.sum(aero_ind[:, 2]), 3), dtype="complex")
    }
    resids = None
    _Component.solve_nonlinear(params, unknowns, resids)
    def_mesh = unknowns.get('def_mesh')
    return def_mesh


"""
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

                                    AERODYNAMICS

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
From vlm.py: """


def vlm_geometry(def_mesh, surface):
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
    _Component = VLMGeometry(surface)
    params = {
        'def_mesh': def_mesh
    }
    unknowns = {
        'b_pts': np.zeros((_Component.nx, _Component.ny, 3), dtype=data_type),
        'mid_b': np.zeros((_Component.nx - 1, _Component.ny, 3), dtype=data_type)),
        'c_pts': np.zeros((_Component.nx - 1, _Component.ny - 1, 3)),
        'widths': np.zeros((_Component.ny - 1))),
        'cos_sweep': np.zeros((_Component.ny - 1)),
        'lengths': np.zeros((_Component.ny)),
        'normals': np.zeros((_Component.nx - 1, _Component.ny - 1, 3)),
        'S_ref': val=0.
    }
    resids=None
    _Component.solve_nonlinear(params, unknowns, resids)
    b_pts=unknowns.get('b_pts')
    mid_b=unknowns.get('mid_b')
    c_pts=unknowns.get('c_pts')
    widths=unknowns.get('widths')
    cos_sweep=unknowns.get('cos_sweep')
    lengths=unknowns.get('lengths')
    normals=unknowns.get('normals')
    S_ref=unknowns.get('S_ref')
    return b_pts, mid_b, c_pts, widths, cos_sweep, lengths, normals, S_ref


def vlm_circulations(aero_ind, def_mesh, b_pts, c_pts, normals, v, alpha):
    """
    Compute the circulations based on the AIC matrix and the panel velocities.
    Note that the flow tangency condition is enforced at the 3/4 chord point.

    Parameters
    ----------
    def_mesh : array_like
        Flattened array defining the lifting surfaces.
    b_pts : array_like
        Bound points for the horseshoe vortices, found along the 1/4 chord.
    c_pts : array_like
        Collocation points on the 3/4 chord line where the flow tangency
        condition is satisfed. Used to set up the linear system.
    normals : array_like
        The normal vector for each panel, computed as the cross of the two
        diagonals from the mesh points.
    v : float
        Freestream air velocity in m/s.
    alpha : float
        Angle of attack in degrees.

    Returns
    -------
    circulations : array_like
        Flattened vector of horseshoe vortex strengths calculated by solving
        the linear system of AIC_mtx * circulations = rhs, where rhs is
        based on the air velocity at each collocation point.

    """
    _Component=VLMCirculations(aero_ind)
    params={
        'def_mesh': def_mesh,
        'b_pts': b_pts,
        'c_pts': c_pts,
        'normals': normals,
        'v': v,
        'alpha': alpha
    }
    unknowns={
        'circulations': np.zeros((np.sum(aero_ind[:, 4])), dtype = "complex")
    }
    resids=None
    _Component.solve_nonlinear(params, unknowns, resids)
    circulations=unknowns.get('circulations')
    return circulations


def assemble_aic(aero_ind, def_mesh, b_pts, c_pts, normals, v, alpha, surfaces):
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
    _Component=AssembleAIC(def_mesh, b_pts, c_pts, normals, surfaces)
    params={}
    for surface in surfaces:
        _Component.surface=surface
        ny=surface['num_y']
        nx=surface['num_x']
        name=surface['name']
        params.update({
            name + 'def_mesh': def_mesh,  # this can't be right
            name + 'b_pts': b_pts,
            name + 'c_pts': c_pts,
            name + 'normals': normals
        })
    unknowns={
        'AIC': np.zeros((_Component.tot_panels, _Component.tot_panels), dtype = data_type),
        'rhs': np.zeros((_Component.tot_panels), dtype = data_type)
    }
    resids=None
    _Component.solve_nonlinear(params, unknowns, resids)
    AIC=unknowns.get('AIC')
    rhs=unknowns.get('rhs')
    return AIC, rhs


def vlm_forces(def_mesh, aero_ind, b_pts, mid_b, circulations, alpha = 3, v = 10, rho = 3):
    """ Compute aerodynamic forces acting on each section.

    Parameters
    ----------
    def_mesh : array_like
        Flattened array defining the lifting surfaces.
    b_pts : array_like
        Bound points for the horseshoe vortices, found along the 1/4 chord.
    mid_b : array_like
        Midpoints of the bound vortex segments, used as the collocation
        points to compute drag.
    circ : array_like   (circulations)
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
    sec_forces : array_like
        Flattened array containing the sectional forces acting on each panel.
        Stored in Fortran order (only relevant when more than one chordwise
        panel).

    """
    _Component=VLMForces(aero_ind)
    params={
        'def_mesh': def_mesh,
        'b_pts': b_pts,
        'mid_b': mid_b,
        'circulations': circulations,
        'alpha': alpha,
        'v': v,
        'rho': rho
    }
    unknowns={
        'sec_forces': np.zeros((np.sum(aero_ind[:, 4]), 3))
    }
    resids=None
    _Component.solve_nonlinear(params, unknowns, resids)
    sec_forces=unknowns.get('sec_forces')
    return sec_forces


def transfer_loads(def_mesh, sec_forces, aero_ind, fem_ind, fem_origin = 0.35):
    """
    Perform aerodynamic load transfer.

    Apply the computed sectional forces on the aerodynamic surfaces to
    obtain the deformed mesh FEM loads.

    Parameters
    ----------
    def_mesh : array_like
        Flattened array defining the lifting surfaces after deformation.
    sec_forces : array_like
        Flattened array containing the sectional forces acting on each panel.
        Stored in Fortran order (only relevant when more than one chordwise
        panel).

    Returns
    -------
    loads : array_like
        Flattened array containing the loads applied on the FEM component,
        computed from the sectional forces.

    """
    _Component=TransferLoads(aero_ind, fem_ind, fem_origin)
    params={
        'def_mesh': def_mesh,
        'sec_forces': sec_forces
    }
    unknowns={
        'loads': np.zeros((np.sum(fem_ind[:, 0]), 6))
    }
    resids=None
    _Component.solve_nonlinear(params, unknowns, resids)
    loads=unknowns.get('loads')
    return loads


def vlm_lift_drag(sec_forces, alpha, aero_ind):
    """
    Calculate total lift and drag in force units based on section forces.

    Parameters
    ----------
    sec_forces : array_like
        Flattened array containing the sectional forces acting on each panel.
        Stored in Fortran order (only relevant when more than one chordwise
        panel).
    alpha : float
        Angle of attack in degrees.

    Returns
    -------
    L : array_like
        Total lift for each lifting surface.
    D : array_like
        Total drag for each lifting surface.

    """
    _Component=VLMLiftDrag(aero_ind)
    num_surf=aero_ind.shape[0]
    params={
        'sec_forces': sec_forces,
        'alpha': alpha
    }
    unknowns={
        'L': np.zeros((num_surf)),
        'D': np.zeros((num_surf))
    }
    resids=None
    _Component.solve_nonlinear(params, unknowns, resids)
    L=unknowns.get('L')
    D=unknowns.get('D')
    return L, D


def vlm_coeffs(S_ref, L, D, v, rho, aero_ind):
    """ Compute lift and drag coefficients.

    Parameters
    ----------
    S_ref : array_like
        The reference areas of each lifting surface.
    L : array_like
        Total lift for each lifting surface.
    D : array_like
        Total drag for each lifting surface.
    v : float
        Freestream air velocity in m/s.
    rho : float
        Air density in kg/m^3.

    Returns
    -------
    CL1 : array_like
        Induced coefficient of lift (CL) for each lifting surface.
    CDi : array_like
        Induced coefficient of drag (CD) for each lifting surface.

    """
    _Component=VLMCoeffs(aero_ind)
    num_surf=aero_ind.shape[0]
    params={
        'S_ref': S_ref,
        'L': L,
        'D': D,
        'v': v,
        'rho': rho
    }
    unknowns={
        'CL1': np.zeros((num_surf)),
        'CDi': np.zeros((num_surf))
    }
    resids=None
    _Component.solve_nonlinear(params, unknowns, resids)
    CL1=unknowns.get('CL1')
    CDi=unknowns.get('CDi')
    return CL1, CDi


def total_lift(CL1, CL0, aero_ind):
    """ Calculate total lift in force units.

    Parameters
    ----------
    CL1 : array_like
        Induced coefficient of lift (CL) for each lifting surface.

    Returns
    -------
    CL : array_like
        Total coefficient of lift (CL) for each lifting surface.
    CL_wing : float
        CL of the main wing, used for CL constrained optimization.

    """
    _Component=TotalLift(CL0, aero_ind)
    params={
        'CL1': CL1
    }
    unknowns={
        'CL': np.zeros((_Component.num_surf)),
        'CL_wing': 0.
    }
    resids=None
    _Component.solve_nonlinear(params, unknowns, resids)
    CL=unknowns.get('CL')
    CL_wing=unknowns.get('CL_wing')
    return CL, CL_wing


def total_drag(CL0, aero_ind):
    """ Calculate total drag in force units.

    Parameters
    ----------
    CDi : array_like
        Induced coefficient of drag (CD) for each lifting surface.

    Returns
    -------
    CD : array_like
        Total coefficient of drag (CD) for each lifting surface.
    CD_wing : float
        CD of the main wing, used for CD minimization.

    """
    _Component=TotalDrag(CDi, CL0, aero_ind)
    params={
        'CDi': CDi
    }
    unknowns={
        'CD': np.zeros((_Component.num_surf)),
        'CD_wing': 0.
    }
    resids=None
    _Component.solve_nonlinear(params, unknowns, resids)
    CD=unknowns.get('CD')
    CD_wing=unknowns.get('CD_wing')
    return CD, CD_wing


"""
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

                                    STRUCTURES

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
From spatialbeam.py: Define the structural analysis component using spatial beam theory. """


def norm(vec):
    return np.sqrt(np.sum(vec**2))


def unit(vec):
    return vec / norm(vec)


def radii(mesh, t_c = 0.15):
    """ Obtain the radii of the FEM component based on chord. """
    vectors=mesh[-1, :, :] - mesh[0, :, :]
    # print('sss mesh.shape',mesh.shape)
    # print('vectors.shape',vectors.shape)
    chords=np.sqrt(np.sum(vectors**2, axis=1))
    chords=0.5 * chords[: -1] + 0.5 * chords[1:]
    return t_c * chords


def spatial_beam_FEM(A, Iy, Iz, J, nodes, loads, aero_ind, fem_ind, E, G, cg_x = 5):
    """
    Compute the displacements and rotations by solving the linear system
    using the structural stiffness matrix.

    Parameters
    ----------
    A : array_like
        Areas for each FEM element.
    Iy : array_like
        Mass moment of inertia around the y-axis for each FEM element.
    Iz : array_like
        Mass moment of inertia around the z-axis for each FEM element.
    J : array_like
        Polar moment of inertia for each FEM element.
    nodes : array_like
        Flattened array with coordinates for each FEM node.
    loads : array_like
        Flattened array containing the loads applied on the FEM component,
        computed from the sectional forces.

    Returns
    -------
    disp_aug : array_like
        Augmented displacement array. Obtained by solving the system
        mtx * disp_aug = rhs, where rhs is a flattened version of loads.

    """
    _Component=SpatialBeamFEM(aero_ind, fem_ind, E, G, cg_x)
    params={
        'A': A,
        'Iy': Iy,
        'Iz': Iz,
        'J': J,
        'nodes': nodes,
        'loads': loads
    }
    unknowns={
        'disp_aug': np.zeros((_Component.size), dtype = "complex")
    }
    resids=None
    _Component.solve_nonlinear(params, unknowns, resids)
    disp_aug=unknowns.get('disp_aug')
    return disp_aug


def spatial_beam_disp(disp_aug, fem_ind):
    """
    Select displacements from augmented vector.

    The solution to the linear system has additional results due to the
    constraints on the FEM model. The displacements from this portion of
    the linear system is not needed, so we select only the relevant
    portion of the displacements for further calculations.

    Parameters
    ----------
    disp_aug : array_like
        Augmented displacement array. Obtained by solving the system
        mtx * disp_aug = rhs, where rhs is a flattened version of loads.

    Returns
    -------
    disp : array_like
        Actual displacement array formed by truncating disp_aug.

    """
    _Component=SpatialBeamDisp(fem_ind)
    params={
        'disp_aug': disp_aug
    }
    unknowns={
        'disp': np.zeros((_Component.tot_n_fem, 6))
    }
    resids=None
    _Component.solve_nonlinear(params, unknowns, resids)
    disp=unknowns.get('disp')
    return disp


def compute_nodes(mesh, fem_ind, aero_ind, fem_origin = 0.35):
    """
    Compute FEM nodes based on aerodynamic mesh.

    The FEM nodes are placed at 0.35*chord, or based on the fem_origin value.

    Parameters
    ----------
    mesh : array_like
        Flattened array defining the lifting surfaces.

    Returns
    -------
    nodes : array_like
        Flattened array with coordinates for each FEM node.

    """
    ComputeNodes_comp=ComputeNodes(fem_ind, aero_ind, fem_origin)
    params={
        'mesh': mesh
    }
    unknowns={
        'nodes': np.zeros((ComputeNodes_comp.tot_n_fem, 3))
    }
    resids=None
    ComputeNodes_comp.solve_nonlinear(params, unknowns, resids)
    nodes=unknowns.get('nodes')
    return nodes


def spatial_beam_energy(disp, loads, aero_ind, fem_ind):
    """ Compute strain energy.

    Parameters
    ----------
    disp : array_like
        Actual displacement array formed by truncating disp_aug.
    loads : array_like
        Flattened array containing the loads applied on the FEM component,
        computed from the sectional forces.

    Returns
    -------
    energy : float
        Total strain energy of the structural component.

    """
    _Component.SpatialBeamEnergy(aero_ind, fem_ind)
    params={
        'disp': disp,
        'loads': loads
    }
    unknowns={
        'energy': np.zeros((_Component.n, 6))
    }
    resids=None
    _Component.solve_nonlinear(params, unknowns, resids)
    energy=unknowns.get('energy')
    return energy


def spatial_beam_weight(A, nodes, aero_ind, fem_ind, mrho):
    """ Compute total weight.

    Parameters
    ----------
    A : array_like
        Areas for each FEM element.
    nodes : array_like
        Flattened array with coordinates for each FEM node.

    Returns
    -------
    weight : float
        Total weight of the structural component."""
    _Component=SpatialBeamWeight(aero_ind, fem_ind, mrho)
    params={
        'A': A,
        'nodes': nodes
    }
    unknowns={
        'weight': 0.
    }
    resids=None
    _Component.solve_nonlinear(params, unknowns, resids)
    weight=unknowns.get('weight')
    return weight


def spatial_beam_vonmises_tube(nodes, r, disp, aero_ind, fem_ind, E, G):
    """ Compute the max von Mises stress in each element.

    Parameters
    ----------
    r : array_like
        Radii for each FEM element.
    nodes : array_like
        Flattened array with coordinates for each FEM node.
    disp : array_like
        Displacements of each FEM node.

    Returns
    -------
    vonmises : array_like
        von Mises stress magnitudes for each FEM element.

    """
    _Component=SpatialBeamVonMisesTube(aero_ind, fem_ind, E, G)
    num_surf=fem_ind.shape[0]
    params={
        'nodes': nodes,
        'r': r,
        'disp': disp
    }
    unknowns={
        'vonmises': np.zeros((_Component.tot_n_fem - num_surf, 2), dtype = "complex")
    }
    resids=None
    _Component.solve_nonlinear(params, unknowns, resids)
    vonmises=unknowns.get('vonmises')
    return vonmises


def spatial_beam_failure_ks(vonmises, fem_ind, sigma, rho = 10):
    """
    Aggregate failure constraints from the structure.

    To simplify the optimization problem, we aggregate the individual
    elemental failure constraints using a Kreisselmeier-Steinhauser (KS)
    function.

    Parameters
    ----------
    vonmises : array_like
        von Mises stress magnitudes for each FEM element.

    Returns
    -------
    failure : float
        KS aggregation quantity obtained by combining the failure criteria
        for each FEM node. Used to simplify the optimization problem by
        reducing the number of constraints.

    """
    _Component=SpatialBeamFailureKS(fem_ind, sigma, rho)
    params={
        'vonmises': vonmises
    }
    unknowns={
        'failure'=0.
    }
    resids=None
    _Component.solve_nonlinear(params, unknowns, resids)
    failure=unknowns.get('failure')
    return failure


"""
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

                                    MATERIALS

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
From materials.py: """


def materials_tube(r, thickness, fem_ind):
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
    MaterialsTube_comp=MaterialsTube(fem_ind)
    params={
        'r': r,
        'thickness': thickness
    }
    unknowns={
        'A': np.zeros((np.sum(fem_ind[:, 0] - fem_ind.shape[0]))),
        'Iy': np.zeros((np.sum(fem_ind[:, 0] - fem_ind.shape[0]))),
        'Iz': np.zeros((np.sum(fem_ind[:, 0] - fem_ind.shape[0]))),
        'J': np.zeros((np.sum(fem_ind[:, 0] - fem_ind.shape[0])))
    }
    resids=None
    MaterialsTube_comp.solve_nonlinear(params, unknowns, resids)
    A=unknowns.get('A', None)
    Iy=unknowns.get('Iy', None)
    Iz=unknowns.get('Iz', None)
    J=unknowns.get('J', None)
    return A, Iy, Iz, J

    """
    --------------------------------------------------------------------------------
    --------------------------------------------------------------------------------

                                        FUNCTIONALS

    --------------------------------------------------------------------------------
    --------------------------------------------------------------------------------
    From functionals.py: """

    def functional_breguet_range(CL, CD, weight, W0, CT, a, R, M, aero_ind):
        """ Computes the fuel burn using the Breguet range equation """
        _Component=FunctionalBreguetRange(W0, CT, a, R, M, aero_ind)
        n_surf=aero_ind.shape[0]
        params={
            'CL': CL,
            'CD': CD,
            'weight': weight
        }
        unknowns={
            'fuelburn': 0.
        }
        resids=None
        _Component.solve_nonlinear(params, unknowns, resids)
        fuelburn=unknowns.get('fuelburn')
        return fuelburn

    def functional_equilibrium(L, weight, fuelburn, W0, aero_ind):
        """ L = W constraint """
        _Component=FunctionalEquilibrium(W0, aero_ind)
        params={
            'L': L,
            'weight': weight,
            'fuelburn': fuelburn
        }
        unknowns={
            'eq_con': 0.
        }
        resids=None
        _Component.solve_nonlinear(params, unknowns, resids)
        eq_con=unknowns.get('eq_con')
        return eq_con
