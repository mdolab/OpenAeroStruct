# Univ. Michigan Aerostructural model.
# Based on OpenAeroStruct by John Hwang, and John Jasa (github.com/johnjasa/OpenAeroStruct)
#
# Cython implementation by Sam Friedman

# Main file for coupled system components

# make compatible Python 2.x to 3.x
from __future__ import print_function, division
__all__ = ['setup','aerodynamics','structures']
# from future.builtins import range  # make compatible Python 2.x to 3.x
import warnings
import sys
import os
import scipy.sparse
from scipy.linalg import lu_factor, lu_solve

from materials import MaterialsTube
from spatialbeam import ComputeNodes, SpatialBeamFEM, SpatialBeamDisp
from transfer import TransferDisplacements, TransferLoads
from vlm import VLMGeometry, VLMCirculations, VLMForces

try:
    import lib
    fortran_flag = True
except:
    fortran_flag = False
sparse_flag = False  # don't use sparse functions

import numpy as np
DTYPE = np.float64  # double precision float
# import cython
# cimport cython
# import numpy
# cimport numpy
# DTYPE = numpy.float64
# ctypedef numpy.float64_t DTYPE_t
# @cython.boundscheck(False)
# @cython.wraparound(False)
# @cython.nonecheck(False)
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
From geometry.py: Manipulate geometry mesh based on high-level design parameters """


def setup(num_inboard=3, num_outboard=4):
    ''' Setup the aerostruct mesh '''
    # Define the aircraft properties, from CRM.py
    span = 58.7630524  # [m] baseline CRM
    # W0 = 0.5 * 2.5e6 # [N] (MTOW of B777 is 3e5 kg with fuel)
    # CT = 9.81 * 17.e-6 # [1/s] (9.81 N/kg * 17e-6 kg/N/s)
    # R = 14.3e6 # [m] maximum range
    M = 0.84  # at cruise
    alpha = 3.  # [deg.]
    rho = 0.38  # [kg/m^3] at 35,000 ft
    a = 295.4  # [m/s] at 35,000 ft
    v = a * M
    # CL0 = 0.2
    # CD0 = 0.015
    # Define spatialbeam properties, from aluminum.py
    E = 200.e9  # [Pa]
    G = 30.e9  # [Pa]
    stress = 20.e6  # [Pa]
    mrho = 3.e3  # [kg/m^3]
    # Create the mesh with 3 inboard points and 4 outboard points.
    # This will be mirrored to produce a mesh with ... spanwise points,
    # or ...-1 spanwise panels
    mesh = gen_crm_mesh(int(num_inboard), int(num_outboard), num_x=2)
    num_x, num_y = mesh.shape[: 2]
    num_twist = np.max([int((num_y - 1) / 5), 5])
    # print('234mesh.shape',mesh.shape)
    r = radii(mesh)
    # Set the number of thickness control points and the initial thicknesses
    num_thickness = num_twist
    t = r / 10
    mesh = mesh.reshape(-1, mesh.shape[-1])
    aero_ind = np.atleast_2d(np.array([num_x, num_y]))
    # print('..... aero_ind.shape',aero_ind.shape)
    # print(aero_ind)
    fem_ind = [num_y]
    aero_ind, fem_ind = get_inds(aero_ind, fem_ind)
    # Set additional mesh parameters
    dihedral = 0.  # dihedral angle in degrees
    sweep = 0.  # shearing sweep angle in degrees
    taper = 1.  # taper ratio
    fem_origin = 0.35
    # Initial displacements of zero
    tot_n_fem = np.sum(fem_ind[:, 0])
    disp = np.zeros((tot_n_fem, 6))
    # # Define Jacobians for b-spline controls
    tot_n_fem = np.sum(fem_ind[:, 0])
    num_surf = fem_ind.shape[0]
    jac_twist = get_bspline_mtx(num_twist, num_y)
    jac_thickness = get_bspline_mtx(num_thickness, tot_n_fem - num_surf)
    # # Define ...
    twist_cp = np.zeros(num_twist)
    thickness_cp = np.ones(num_thickness) * np.max(t)
    twist = cp2pt(twist_cp, jac_twist)
    thickness = cp2pt(thickness_cp, jac_thickness)
    mesh = geometry_mesh(mesh, aero_ind, twist, 0, 0, 1, span=58.7630524)
    #print('mesh.shape',mesh.shape)
    def_mesh = transfer_displacements(
        mesh, disp, aero_ind, fem_ind, fem_origin=0.35)
    # Output the def_mesh for the aero modules
    # Other variables needed for aero and struct modules
    params = {
        'mesh': mesh,
        'num_x': num_x,
        'num_y': num_y,
        'span': span,
        'twist_cp': twist_cp,
        'thickness_cp': thickness_cp,
        'v': v,
        'alpha': alpha,
        'rho': rho,
        'r': r,
        't': t,
        'aero_ind': aero_ind,
        'fem_ind': fem_ind,
        'num_thickness': num_thickness,
        'num_twist': num_twist,
        'sweep': sweep,
        'taper': taper,
        'dihedral': dihedral,
        'E': E,
        'G': G,
        'stress': stress,
        'mrho': mrho,
        'tot_n_fem': tot_n_fem,
        'num_surf': num_surf,
        'jac_twist': jac_twist,
        'jac_thickness': jac_thickness,
        'fem_origin': fem_origin
    }

    return (def_mesh, params)


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
        #print('shape of ',i,' = ',j.shape)
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


def get_bspline_mtx(num_cp, num_pt, order=4):
    """ Create Jacobian to fit a bspline to a set of data.
    ...from b_spline.py

    Parameters
    ----------
    num_cp : int
        Number of control points.
    num_pt : int
        Number of points.
    order : int, optional
        Order of b-spline fit.

    Returns
    -------
    out : CSR sparse matrix
        Matrix that gives the points vector when multiplied by the control
        points vector.

    """
    knots = np.zeros(num_cp + order)
    knots[order - 1:num_cp + 1] = np.linspace(0, 1, num_cp - order + 2)
    knots[num_cp + 1:] = 1.0
    t_vec = np.linspace(0, 1, num_pt)
    basis = np.zeros(order)
    arange = np.arange(order)
    data = np.zeros((num_pt, order))
    rows = np.zeros((num_pt, order), int)
    cols = np.zeros((num_pt, order), int)
    for ipt in range(num_pt):
        t = t_vec[ipt]
        i0 = -1
        for ind in range(order, num_cp + 1):
            if (knots[ind - 1] <= t) and (t < knots[ind]):
                i0 = ind - order
        if t == knots[-1]:
            i0 = num_cp - order
        basis[:] = 0.
        basis[-1] = 1.
        for i in range(2, order + 1):
            l = i - 1
            j1 = order - l
            j2 = order
            n = i0 + j1
            if knots[n + l] != knots[n]:
                basis[j1 - 1] = (knots[n + l] - t) / \
                    (knots[n + l] - knots[n]) * basis[j1]
            else:
                basis[j1 - 1] = 0.
            for j in range(j1 + 1, j2):
                n = i0 + j
                if knots[n + l - 1] != knots[n - 1]:
                    basis[j - 1] = (t - knots[n - 1]) / \
                        (knots[n + l - 1] - knots[n - 1]) * basis[j - 1]
                else:
                    basis[j - 1] = 0.
                if knots[n + l] != knots[n]:
                    basis[j - 1] += (knots[n + l] - t) / \
                        (knots[n + l] - knots[n]) * basis[j]
            n = i0 + j2
            if knots[n + l - 1] != knots[n - 1]:
                basis[j2 - 1] = (t - knots[n - 1]) / \
                    (knots[n + l - 1] - knots[n - 1]) * basis[j2 - 1]
            else:
                basis[j2 - 1] = 0.
        data[ipt, :] = basis
        rows[ipt, :] = ipt
        cols[ipt, :] = i0 + arange
    data, rows, cols = data.flatten(), rows.flatten(), cols.flatten()
    return scipy.sparse.csr_matrix((data, (rows, cols)), shape=(num_pt, num_cp))


def get_inds(aero_ind, fem_ind):
    """
    Calculate and store indices to describe panels for aero and
    structural analysis.

    Takes in aero_ind with each row containing [nx, ny] and fem_ind with
    each row containing [n_fem].

    Each outputted row has information for each individually defined surface,
    stored in the order [nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels]
    with the indices    [ 0,  1, 2,      3,        4, 5,      6,        7]

    nx : number of nodes in the chordwise direction
    ny : number of nodes in the spanwise direction
    n : total number of nodes
    n_bpts : total number of b_pts nodes
    n_panels : total number of panels
    i : current index of nodes when considering all surfaces
    i_bpts: current index of b_pts nodes when considering all surfaces
    i_panels : current index of panels when considering all surfaces

    Simpler than the aero case, the fem_ind array contains:
    [n_fem, i_fem]

    n_fem : number of fem nodes per surface
    i_fem : current index of fem nodes when considering all fem nodes

    """
    new_aero_ind = np.zeros((aero_ind.shape[0], 8), dtype=int)
    new_aero_ind[:, 0:2] = aero_ind
    for i, row in enumerate(aero_ind):
        nx, ny = aero_ind[i, :]
        new_aero_ind[i, 2] = nx * ny
        new_aero_ind[i, 3] = (nx - 1) * ny
        new_aero_ind[i, 4] = (nx - 1) * (ny - 1)
        new_aero_ind[i, 5] = np.sum(np.product(aero_ind[:i], axis=1))
        new_aero_ind[i, 6] = np.sum((aero_ind[:i, 0] - 1) * aero_ind[:i, 1])
        new_aero_ind[i, 7] = np.sum(np.product(aero_ind[:i] - 1, axis=1))
    new_fem_ind = np.zeros((len(fem_ind), 2), dtype=int)
    new_fem_ind[:, 0] = fem_ind
    for i, row in enumerate(fem_ind):
        new_fem_ind[i, 1] = np.sum(fem_ind[:i])
    return new_aero_ind, new_fem_ind


def rotate(mesh, thetas):
    """
    Compute rotation matrices given mesh and rotation angles in degrees.

    """
    te = mesh[-1]
    le = mesh[0]
    quarter_chord = 0.25 * te + 0.75 * le
    ny = mesh.shape[1]
    nx = mesh.shape[0]
    rad_thetas = thetas * np.pi / 180.
    mats = np.zeros((ny, 3, 3), dtype=DTYPE)
    mats[:, 0, 0] = np.cos(rad_thetas)
    mats[:, 0, 2] = np.sin(rad_thetas)
    mats[:, 1, 1] = 1
    mats[:, 2, 0] = -np.sin(rad_thetas)
    mats[:, 2, 2] = np.cos(rad_thetas)
    for ix in range(nx):
        row = mesh[ix]
        row[:] = np.einsum("ikj, ij -> ik", mats, row - quarter_chord)
        row += quarter_chord
    return mesh


def sweep(mesh, angle):
    """ Apply shearing sweep. Positive sweeps back. """
    num_x, num_y, _ = mesh.shape
    ny2 = int((num_y - 1) / 2)
    le = mesh[0]
    y0 = le[ny2, 1]
    p180 = np.pi / 180
    tan_theta = np.tan(p180 * angle)
    dx_right = (le[ny2:, 1] - y0) * tan_theta
    dx_left = -(le[:ny2, 1] - y0) * tan_theta
    dx = np.hstack((dx_left, dx_right))
    for i in range(num_x):
        mesh[i, :, 0] += dx
    return mesh


def dihedral(mesh, angle):
    """ Apply dihedral angle. Positive bends up. """
    num_x, num_y, _ = mesh.shape
    ny2 = int((num_y - 1) / 2)
    le = mesh[0]
    y0 = le[ny2, 1]
    p180 = np.pi / 180
    tan_theta = np.tan(p180 * angle)
    dx_right = (le[ny2:, 1] - y0) * tan_theta
    dx_left = -(le[:ny2, 1] - y0) * tan_theta
    dx = np.hstack((dx_left, dx_right))
    for i in range(num_x):
        mesh[i, :, 2] += dx
    return mesh


def stretch(mesh, length):
    """ Stretch mesh in spanwise direction to reach specified length. """
    le = mesh[0]
    num_x, num_y, _ = mesh.shape
    span = le[-1, 1] - le[0, 1]
    dy = (length - span) / (num_y - 1) * np.arange(1, num_y)
    for i in range(num_x):
        mesh[i, 1:, 1] += dy
    return mesh


def taper(mesh, taper_ratio):
    """ Alter the spanwise chord to produce a tapered wing. """
    le = mesh[0]
    te = mesh[-1]
    num_x, num_y, _ = mesh.shape
    ny2 = int((num_y + 1) / 2)
    center_chord = .5 * te + .5 * le
    taper = np.linspace(1, taper_ratio, ny2)[::-1]
    jac = get_bspline_mtx(ny2, ny2, order=2)
    taper = jac.dot(taper)
    dx = np.hstack((taper, taper[::-1][1:]))
    for i in range(num_x):
        for ind in range(3):
            mesh[i, :, ind] = (
                mesh[i, :, ind] - center_chord[:, ind]) * dx + center_chord[:, ind]
    return mesh


def mirror(mesh, right_side=True):
    """
    Take a half geometry and mirror it across the symmetry plane.

    If right_side==True, it mirrors from right to left,
    assuming that the first point is on the symmetry plane. Else
    it mirrors from left to right, assuming the last point is on the
    symmetry plane.

    """
    num_x, num_y, _ = mesh.shape
    new_mesh = np.empty((num_x, 2 * num_y - 1, 3))
    mirror_y = np.ones(mesh.shape)
    mirror_y[:, :, 1] *= -1.0
    if right_side:
        new_mesh[:, :num_y, :] = mesh[:, ::-1, :] * mirror_y
        new_mesh[:, num_y:, :] = mesh[:,   1:, :]
    else:
        new_mesh[:, :num_y, :] = mesh[:, ::-1, :]
        new_mesh[:, num_y:, :] = mesh[:,   1:, :] * mirror_y[:, 1:, :]
    # shift so 0 is at the left wing tip (structures wants it that way)
    y0 = new_mesh[0, 0, 1]
    new_mesh[:, :, 1] -= y0
    return new_mesh


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


def calc_vorticity(A, B, P):
    """ Calculates the influence coefficient for a vortex filament.

    Parameters
    ----------
    A[3] : array_like
        Coordinates for the start point of the filament.
    B[3] : array_like
        Coordinates for the end point of the filament.
    P[3] : array_like
        Coordinates for the collocation point where the influence coefficient
        is computed.

    Returns
    -------
    out[3] : array_like
        Influence coefficient contribution for the described filament.

    """
    r1 = P - A
    r2 = P - B
    r1_mag = norm(r1)
    r2_mag = norm(r2)
    return (r1_mag + r2_mag) * np.cross(r1, r2) / (r1_mag * r2_mag * (r1_mag * r2_mag + r1.dot(r2)))


def get_lengths(A, B, axis):
    return np.sqrt(np.sum((B - A)**2, axis=axis))


def vlm_geometry(aero_ind, def_mesh):
    """ Compute various geometric properties for VLM analysis.

    Parameters
    ----------
    def_mesh : array_like
        Flattened array defining the lifting surfaces.

    Returns
    -------
    b_pts : array_like
        Bound points for the horseshoe vortices, found along the 1/4 chord.
    mid_b : array_like
        Midpoints of the bound vortex segments, used as the collocation
        points to compute drag.
    c_pts : array_like
        Collocation points on the 3/4 chord line where the flow tangency
        condition is satisfed. Used to set up the linear system.
    widths : array_like
        The spanwise widths of each individual panel.
    normals : array_like
        The normal vector for each panel, computed as the cross of the two
        diagonals from the mesh points.
    S_ref : array_like
        The reference areas of each lifting surface.

    """
    _Component = VLMGeometry(aero_ind)
    params = {
        'def_mesh': def_mesh
    }
    unknowns = {
        'b_pts': np.zeros((np.sum(aero_ind[:, 3]), 3), dtype="complex"),
        'mid_b': np.zeros((np.sum(aero_ind[:, 4]), 3), dtype="complex"),
        'c_pts': np.zeros((np.sum(aero_ind[:, 4]), 3), dtype="complex"),
        'widths': np.zeros((np.sum(aero_ind[:, 4]))),
        'normals': np.zeros((np.sum(aero_ind[:, 4]), 3)),
        'S_ref': np.zeros((aero_ind.shape[0]))
    }
    resids = None
    _Component.solve_nonlinear(params, unknowns, resids)
    b_pts = unknowns.get('b_pts')
    mid_b = unknowns.get('mid_b')
    c_pts = unknowns.get('c_pts')
    widths = unknowns.get('widths')
    normals = unknowns.get('normals')
    S_ref = unknowns.get('S_ref')
    return b_pts, mid_b, c_pts, widths, normals, S_ref


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
    _Component = VLMCirculations(aero_ind)
    params = {
        'def_mesh': def_mesh,
        'b_pts': b_pts,
        'c_pts': c_pts,
        'normals': normals,
        'v': v,
        'alpha': alpha
    }
    unknowns = {
        'circulations': np.zeros((np.sum(aero_ind[:, 4])), dtype="complex")
    }
    resids = None
    _Component.solve_nonlinear(params, unknowns, resids)
    circulations = unknowns.get('circulations')
    return circulations


def vlm_forces(def_mesh, aero_ind, b_pts, mid_b, circulations, alpha=3, v=10, rho=3):
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
    _Component = VLMForces(aero_ind)
    params = {
        'def_mesh': def_mesh,
        'b_pts': b_pts,
        'mid_b': mid_b,
        'circulations': circulations,
        'alpha': alpha,
        'v': v,
        'rho': rho
    }
    unknowns = {
        'sec_forces': np.zeros((np.sum(aero_ind[:, 4]), 3))
    }
    resids = None
    _Component.solve_nonlinear(params, unknowns, resids)
    sec_forces = unknowns.get('sec_forces')
    return sec_forces


def transfer_loads(def_mesh, sec_forces, aero_ind, fem_ind, fem_origin=0.35):
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
    _Component = TransferLoads(aero_ind, fem_ind, fem_origin)
    params = {
        'def_mesh': def_mesh,
        'sec_forces': sec_forces
    }
    unknowns = {
        'loads': np.zeros((np.sum(fem_ind[:, 0]), 6))
    }
    resids = None
    _Component.solve_nonlinear(params, unknowns, resids)
    loads = unknowns.get('loads')
    return loads
    

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


def radii(mesh, t_c=0.15):
    """ Obtain the radii of the FEM component based on chord. """
    vectors = mesh[-1, :, :] - mesh[0, :, :]
    # print('sss mesh.shape',mesh.shape)
    #print('vectors.shape',vectors.shape)
    chords = np.sqrt(np.sum(vectors**2, axis=1))
    chords = 0.5 * chords[: -1] + 0.5 * chords[1:]
    return t_c * chords


def spatial_beam_FEM(A, Iy, Iz, J, nodes, loads, aero_ind, fem_ind, E, G, cg_x=5):
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
    _Component = SpatialBeamFEM(aero_ind, fem_ind, E, G, cg_x)
    params = {
        'A': A,
        'Iy': Iy,
        'Iz': Iz,
        'J': J,
        'nodes': nodes,
        'loads': loads
    }
    unknowns = {
        'disp_aug': np.zeros((_Component.size), dtype="complex")
    }
    resids = None
    _Component.solve_nonlinear(params, unknowns, resids)
    disp_aug = unknowns.get('disp_aug')
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
    _Component = SpatialBeamDisp(fem_ind)
    params = {
        'disp_aug': disp_aug
    }
    unknowns = {
        'disp': np.zeros((_Component.tot_n_fem, 6))
    }
    resids = None
    _Component.solve_nonlinear(params, unknowns, resids)
    disp = unknowns.get('disp')
    return disp


def compute_nodes(mesh, fem_ind, aero_ind, fem_origin=0.35):
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
    ComputeNodes_comp = ComputeNodes(fem_ind, aero_ind, fem_origin)
    params = {
        'mesh': mesh
    }
    unknowns = {
        'nodes': np.zeros((ComputeNodes_comp.tot_n_fem, 3))
    }    
    resids = None
    ComputeNodes_comp.solve_nonlinear(params, unknowns, resids)
    nodes = unknowns.get('nodes')
    return nodes


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
    MaterialsTube_comp = MaterialsTube(fem_ind)
    params = {
        'r': r,
        'thickness': thickness
    }
    unknowns = {
        'A': np.zeros((np.sum(fem_ind[:,0]-fem_ind.shape[0]))),
        'Iy': np.zeros((np.sum(fem_ind[:,0]-fem_ind.shape[0]))),
        'Iz': np.zeros((np.sum(fem_ind[:,0]-fem_ind.shape[0]))),
        'J': np.zeros((np.sum(fem_ind[:,0]-fem_ind.shape[0])))
    }
    resids = None
    MaterialsTube_comp.solve_nonlinear(params, unknowns, resids)
    A = unknowns.get('A', None)
    Iy = unknowns.get('Iy', None)
    Iz = unknowns.get('Iz', None)
    J = unknowns.get('J', None)
    return A, Iy, Iz, J
