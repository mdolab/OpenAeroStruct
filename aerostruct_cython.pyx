# Univ. Michigan Aerostructural model.
# Based on OpenAeroStruct by John Hwang, and John Jasa (github.com/johnjasa/OpenAeroStruct)
#
# Cython implementation by Sam Friedman

# Main file for coupled system components

# make compatible Python 2.x to 3.x
from __future__ import print_function, division
# __all__ = ['setup','aerodynamics','structures']
from future.builtins import range  # make compatible Python 2.x to 3.x
import warnings
import sys
import os
import scipy.sparse
from scipy.linalg import lu_factor, lu_solve

try:
    import lib
    fortran_flag = True
except:
    fortran_flag = False
sparse_flag = False  # don't use sparse functions

# import numpy
# DTYPE = numpy.float64  # double precision float
import cython
cimport cython
import numpy
cimport numpy
DTYPE = numpy.float64
ctypedef numpy.float64_t DTYPE_t


"""
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------


                                GEOMETRY / SETUP


--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
From geometry.py: Manipulate geometry mesh based on high-level design parameters """

# @cython.boundscheck(False)
# @cython.wraparound(False)
# @cython.nonecheck(False)

def setup(int num_inboard=3, int num_outboard=4):
    ''' Setup the aerostruct mesh '''
    # Define the aircraft properties, from CRM.py

    cdef double span, M, alpha, rho, a, v, E, G, stress, mrho
    cdef double dihedral, sweep, taper, fem_origin, num_twist
    cdef int num_x, num_y, tot_n_fem

    span = 58.7630524  # [m] baseline CRM
    M = 0.84  # at cruise
    alpha = 3.  # [deg.]
    rho = 0.38  # [kg/m^3] at 35,000 ft
    a = 295.4  # [m/s] at 35,000 ft
    v = a * M
    # Define spatialbeam properties, from aluminum.py
    E = 200.e9  # [Pa]
    G = 30.e9  # [Pa]
    stress = 20.e6  # [Pa]
    mrho = 3.e3  # [kg/m^3]
    # Create the mesh with 3 inboard points and 4 outboard points.
    # This will be mirrored to produce a mesh with ... spanwise points,
    # or ...-1 spanwise panels

    cdef numpy.ndarray[DTYPE_t, ndim=3] mesh1
    mesh1 = gen_crm_mesh(num_inboard, num_outboard, num_x=2)
    num_x, num_y = mesh1.shape[: 2]
    num_twist = numpy.max([((num_y - 1) / 5), 5])
    r = radii(mesh1)
    # Set the number of thickness control points and the initial thicknesses
    num_thickness = num_twist
    t = r / 10

    cdef numpy.ndarray[DTYPE_t, ndim=2] mesh2
    mesh2 = mesh1.reshape(-1, mesh1.shape[2])
    aero_ind1 = numpy.atleast_2d(numpy.array([num_x, num_y]))
    fem_ind1 = [num_y]
    aero_ind, fem_ind = get_inds(aero_ind1, fem_ind1)
    # Set additional mesh parameters
    dihedral = 0.  # dihedral angle in degrees
    sweep = 0.  # shearing sweep angle in degrees
    taper = 1.  # taper ratio
    fem_origin = 0.35
    # Initial displacements of zero
    tot_n_fem = numpy.sum(fem_ind[:, 0])

    cdef numpy.ndarray[DTYPE_t, ndim=2] disp

    disp = numpy.zeros((tot_n_fem, 6))
    # # Define Jacobians for b-spline controls
    tot_n_fem = numpy.sum(fem_ind[:, 0])
    num_surf = fem_ind.shape[0]

    # cdef numpy.ndarray[DTYPE_t, ndim=2] jac_twist, jac_thickness

    jac_twist = get_bspline_mtx(num_twist, num_y)
    jac_thickness = get_bspline_mtx(num_thickness, tot_n_fem - num_surf)
    # # Define ...

    cdef numpy.ndarray[DTYPE_t, ndim=1] twist_cp, thickness_cp, twist, thickness
    twist_cp = numpy.zeros(num_twist)
    thickness_cp = numpy.ones(num_thickness) * numpy.max(t)
    twist = cp2pt(twist_cp, jac_twist)
    thickness = cp2pt(thickness_cp, jac_thickness)
    mesh = geometry_mesh(mesh2, aero_ind, twist, 0, 0, 1, span=58.7630524)
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


def aerodynamics(numpy.ndarray[DTYPE_t, ndim=2] def_mesh, params=None):
    # Unpack variables
    cdef double alpha, rho, v, fem_origin
    aero_ind = params.get('aero_ind')
    alpha = params.get('alpha')
    v = params.get('v')
    rho = params.get('rho')
    fem_ind = params.get('fem_ind')
    fem_origin = params.get('fem_origin', 0.35)

    # # Define Jacobians for b-spline controls
    # tot_n_fem = numpy.sum(fem_ind[:, 0])
    # num_surf = fem_ind.shape[0]
    # jac_twist = get_bspline_mtx(num_twist, num_y)
    # jac_thickness = get_bspline_mtx(num_thickness, tot_n_fem - num_surf)

    cdef numpy.ndarray[DTYPE_t, ndim=2] b_pts, mid_b, c_pts, normals
    cdef numpy.ndarray[DTYPE_t, ndim=1] widths, S_ref
    b_pts, mid_b, c_pts, widths, normals, S_ref = vlm_geometry(aero_ind, def_mesh)

    cdef numpy.ndarray[DTYPE_t, ndim=1] circulations
    circulations = vlm_circulations(aero_ind, def_mesh, b_pts, c_pts, normals, v, alpha)

    cdef numpy.ndarray[DTYPE_t, ndim=2] sec_forces
    sec_forces = vlm_forces(def_mesh, aero_ind, b_pts, mid_b, circulations, alpha, v, rho)

    cdef numpy.ndarray[DTYPE_t, ndim=2] loads
    loads = transfer_loads(def_mesh, sec_forces, aero_ind, fem_ind, fem_origin)

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
    pt = numpy.zeros(jac.shape[0])
    pt = jac.dot(cp)
    return pt


def get_bspline_mtx(int num_cp, int num_pt, int order=4):
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
    cdef int ipt
    knots = numpy.zeros(num_cp + order)
    knots[order - 1:num_cp + 1] = numpy.linspace(0, 1, num_cp - order + 2)
    knots[num_cp + 1:] = 1.0
    t_vec = numpy.linspace(0, 1, num_pt)
    basis = numpy.zeros(order)
    arange = numpy.arange(order)
    data = numpy.zeros((num_pt, order))
    rows = numpy.zeros((num_pt, order), int)
    cols = numpy.zeros((num_pt, order), int)
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
    new_aero_ind = numpy.zeros((aero_ind.shape[0], 8), dtype=int)
    new_aero_ind[:, 0:2] = aero_ind
    for i, row in enumerate(aero_ind):
        nx, ny = aero_ind[i, :]
        new_aero_ind[i, 2] = nx * ny
        new_aero_ind[i, 3] = (nx - 1) * ny
        new_aero_ind[i, 4] = (nx - 1) * (ny - 1)
        new_aero_ind[i, 5] = numpy.sum(numpy.product(aero_ind[:i], axis=1))
        new_aero_ind[i, 6] = numpy.sum((aero_ind[:i, 0] - 1) * aero_ind[:i, 1])
        new_aero_ind[i, 7] = numpy.sum(numpy.product(aero_ind[:i] - 1, axis=1))
    new_fem_ind = numpy.zeros((len(fem_ind), 2), dtype=int)
    new_fem_ind[:, 0] = fem_ind
    for i, row in enumerate(fem_ind):
        new_fem_ind[i, 1] = numpy.sum(fem_ind[:i])
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
    rad_thetas = thetas * numpy.pi / 180.
    mats = numpy.zeros((ny, 3, 3), dtype=DTYPE)
    mats[:, 0, 0] = numpy.cos(rad_thetas)
    mats[:, 0, 2] = numpy.sin(rad_thetas)
    mats[:, 1, 1] = 1
    mats[:, 2, 0] = -numpy.sin(rad_thetas)
    mats[:, 2, 2] = numpy.cos(rad_thetas)
    for ix in range(nx):
        row = mesh[ix]
        row[:] = numpy.einsum("ikj, ij -> ik", mats, row - quarter_chord)
        row += quarter_chord
    return mesh


def sweep(mesh, angle):
    """ Apply shearing sweep. Positive sweeps back. """
    num_x, num_y, _ = mesh.shape
    ny2 = (num_y - 1) / 2
    le = mesh[0]
    y0 = le[ny2, 1]
    p180 = numpy.pi / 180
    tan_theta = numpy.tan(p180 * angle)
    dx_right = (le[ny2:, 1] - y0) * tan_theta
    dx_left = -(le[:ny2, 1] - y0) * tan_theta
    dx = numpy.hstack((dx_left, dx_right))
    for i in range(num_x):
        mesh[i, :, 0] += dx
    return mesh


def dihedral(mesh, angle):
    """ Apply dihedral angle. Positive bends up. """
    num_x, num_y, _ = mesh.shape
    ny2 = (num_y - 1) / 2
    le = mesh[0]
    y0 = le[ny2, 1]
    p180 = numpy.pi / 180
    tan_theta = numpy.tan(p180 * angle)
    dx_right = (le[ny2:, 1] - y0) * tan_theta
    dx_left = -(le[:ny2, 1] - y0) * tan_theta
    dx = numpy.hstack((dx_left, dx_right))
    for i in range(num_x):
        mesh[i, :, 2] += dx
    return mesh


def stretch(mesh, length):
    """ Stretch mesh in spanwise direction to reach specified length. """
    le = mesh[0]
    num_x, num_y, _ = mesh.shape
    span = le[-1, 1] - le[0, 1]
    dy = (length - span) / (num_y - 1) * numpy.arange(1, num_y)
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
    taper = numpy.linspace(1, taper_ratio, ny2)[::-1]
    jac = get_bspline_mtx(ny2, ny2, order=2)
    taper = jac.dot(taper)
    dx = numpy.hstack((taper, taper[::-1][1:]))
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
    new_mesh = numpy.empty((num_x, 2 * num_y - 1, 3))
    mirror_y = numpy.ones(mesh.shape)
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
    raw_crm_points = numpy.array([
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
    # le = numpy.vstack((raw_crm_points[:,1],
    #                 raw_crm_points[:,2],
    #                 raw_crm_points[:,3]))
    # te = numpy.vstack((raw_crm_points[:,1]+raw_crm_points[:,5],
    #                 raw_crm_points[:,2],
    #                 raw_crm_points[:,3]))
    # mesh = numpy.empty((2,20,3))
    # mesh[0,:,:] = le.T
    # mesh[1,:,:] = te.T
    # mesh *= 0.0254 # convert to meters
    # pull out the 3 key y-locations to define the two linear regions of the
    # wing
    crm_base_points = raw_crm_points[(0, 6, 19), :]
    le_base = numpy.vstack((crm_base_points[:, 1],
                            crm_base_points[:, 2],
                            crm_base_points[:, 3]))
    te_base = numpy.vstack((crm_base_points[:, 1] + crm_base_points[:, 5],
                            crm_base_points[:, 2],
                            crm_base_points[:, 3]))
    mesh = numpy.empty((2, 3, 3))
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
    half_mesh = numpy.zeros((2, n_points_total, 3))
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
    full_mesh[:, :, 1] -= numpy.mean(full_mesh[:, :, 1])
    return full_mesh


def add_chordwise_panels(mesh, num_x):
    """ Divide the wing into multiple chordwise panels. """
    le = mesh[0, :, :]
    te = mesh[-1, :, :]
    new_mesh = numpy.zeros((num_x, mesh.shape[1], 3))
    new_mesh[0, :, :] = le
    new_mesh[-1, :, :] = te
    for i in range(1, num_x - 1):
        w = float(i) / (num_x - 1)
        new_mesh[i, :, :] = (1 - w) * le + w * te
    return new_mesh


def gen_mesh(num_x, num_y, span, chord, cosine_spacing=0.):
    """ Generate simple rectangular wing mesh. """
    mesh = numpy.zeros((num_x, num_y, 3))
    ny2 = (num_y + 1) / 2
    beta = numpy.linspace(0, numpy.pi / 2, ny2)
    # mixed spacing with w as a weighting factor
    cosine = .5 * numpy.cos(beta)  # cosine spacing
    uniform = numpy.linspace(0, .5, ny2)[::-1]  # uniform spacing
    half_wing = cosine * cosine_spacing + (1 - cosine_spacing) * uniform
    full_wing = numpy.hstack((-half_wing[:-1], half_wing[::-1])) * span
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
    twist = numpy.zeros(ny)  # default option
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
    cdef int tot_n = numpy.sum(aero_ind[:, 2])
    cdef int tot_n_fem = numpy.sum(fem_ind[:, 0])
    out_def_mesh = numpy.zeros((tot_n, 3), dtype=DTYPE)
    cdef int i_surf, nx, ny, n ,n_bpts, n_panels, i, i_bpts, i_panels
    for i_surf, row in enumerate(fem_ind):
        nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = aero_ind[i_surf, :]
        n_fem, i_fem = row
        mesh2 = mesh[i: i + n, :].reshape(nx, ny, 3)
        disp2 = disp[i_fem: i_fem + n_fem]
        w = fem_origin
        ref_curve = (1 - w) * mesh2[0, :, :] + w * mesh2[-1, :, :]
        Smesh = numpy.zeros(mesh2.shape, dtype=DTYPE)
        for ind in range(nx):
            Smesh[ind, :, :] = mesh2[ind, :, :] - ref_curve
        def_mesh = numpy.zeros(mesh2.shape, dtype=DTYPE)
        cos, sin = numpy.cos, numpy.sin
        for ind in range(ny):
            dx, dy, dz, rx, ry, rz = disp2[ind, :]
            # 1 eye from the axis rotation matrices
            # -3 eye from subtracting Smesh three times
            T = -2 * numpy.eye(3, dtype=DTYPE)
            T[1:,  1:] += [[cos(rx), -sin(rx)], [sin(rx), cos(rx)]]
            T[:: 2, :: 2] += [[cos(ry),  sin(ry)], [-sin(ry), cos(ry)]]
            T[: 2, : 2] += [[cos(rz), -sin(rz)], [sin(rz), cos(rz)]]
            def_mesh[:, ind, :] += Smesh[:, ind, :].dot(T)
            def_mesh[:, ind, 0] += dx
            def_mesh[:, ind, 1] += dy
            def_mesh[:, ind, 2] += dz
        out_def_mesh[i: i + n, :] = (def_mesh + mesh2).reshape(n, 3)
    return out_def_mesh


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
    return (r1_mag + r2_mag) * numpy.cross(r1, r2) / (r1_mag * r2_mag * (r1_mag * r2_mag + r1.dot(r2)))


def get_lengths(A, B, axis):
    return numpy.sqrt(numpy.sum((B - A)**2, axis=axis))


def vlm_geometry(numpy.ndarray[long, ndim=2] aero_ind, numpy.ndarray[DTYPE_t, ndim=2] def_mesh):
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
    cdef int num_surf, tot_n, tot_bpts, tot_panels, i_surf
    num_surf = aero_ind.shape[0]
    tot_n = numpy.sum(aero_ind[:, 2])
    tot_bpts = numpy.sum(aero_ind[:, 3])
    tot_panels = numpy.sum(aero_ind[:, 4])

    B_PTS = numpy.zeros((tot_bpts, 3), dtype=DTYPE)
    MID_B = numpy.zeros((tot_panels, 3), dtype=DTYPE)
    C_PTS = numpy.zeros((tot_panels, 3))
    WIDTHS = numpy.zeros((tot_panels))
    NORMALS = numpy.zeros((tot_panels, 3))
    S_REF = numpy.zeros((num_surf))

    for i_surf, row in enumerate(aero_ind):
        nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = row
        mesh = def_mesh[i: i + n, :].reshape(nx, ny, 3)
        b_pts = mesh[: -1, :, :] * .75 + mesh[1:, :, :] * .25
        mid_b = (b_pts[:, 1:, :] + b_pts[:, : -1, :]) / 2
        c_pts = 0.5 * 0.25 * mesh[: -1, : -1, : ] + \
            0.5 * 0.75 * mesh[1: , : -1, : ] + \
            0.5 * 0.25 * mesh[: -1,  1: , : ] + \
            0.5 * 0.75 * mesh[1:,  1:, :]
        widths = get_lengths(b_pts[:, 1:, :], b_pts[:, : -1, :], 2)
        normals = numpy.cross(
            mesh[: -1,  1:, :] - mesh[1:, : -1, :],
            mesh[: -1, : -1, :] - mesh[1:,  1:, :],
            axis=2)
        norms = numpy.sqrt(numpy.sum(normals**2, axis=2))
        for j in range(3):
            normals[:, :, j] /= norms
        B_PTS[i_bpts: i_bpts + n_bpts, :] = b_pts.reshape(-1, b_pts.shape[-1])
        MID_B[i_panels: i_panels + n_panels,
              :] = mid_b.reshape(-1, mid_b.shape[-1])
        C_PTS[i_panels: i_panels + n_panels,
              :] = c_pts.reshape(-1, c_pts.shape[-1])
        WIDTHS[i_panels: i_panels + n_panels] = widths.flatten()
        NORMALS[i_panels: i_panels + n_panels,
                :] = normals.reshape(-1, normals.shape[-1], order='F')
        S_REF[i_surf] = 0.5 * numpy.sum(norms)
    return B_PTS, MID_B, C_PTS, WIDTHS, NORMALS, S_REF


def assemble_AIC_mtx(mtx, flat_mesh, aero_ind, points, b_pts, alpha, skip=False):
    """
    Compute the aerodynamic influence coefficient matrix
    for either solving the linear system or solving for the drag.

    Parameters
    ----------
    mtx[num_y-1, num_y-1, 3] : array_like
        Aerodynamic influence coefficient (AIC) matrix, or the
        derivative of v w.r.t. circulations.
    flat_mesh[num_x*num_y, 3] : array_like
        Flat array containing nodal coordinates.
    aero_ind[num_surf, 7] : array_like
        Array containing index information for the lifting surfaces.
        See geometry.py/get_inds for more details.
    fem_ind[num_surf, 3] : array_like
        Array containing index information for the lifting surfaces.
        See geometry.py/get_inds for more details.
    points[num_y-1, 3] : array_like
        Collocation points used to find influence coefficient strength.
        Found at 3/4 chord for the linear system and at the midpoint of
        the bound vortices (1/4 chord) for the drag computation.
    b_pts[num_x-1, num_y, 3] : array_like
        Bound vortex coordinates from the 1/4 chord line.
    alpha : float
        Angle of attack.
    skip : boolean
        If false, the bound vortex contributions on the collocation point
        corresponding to the same panel are not included. Used for the drag
        computation.

    Returns
    -------
    mtx[num_y-1, num_y-1, 3] : array_like
        Aerodynamic influence coefficient (AIC) matrix, or the
        derivative of v w.r.t. circulations.
    """
    cdef numpy.ndarray[DTYPE_t, ndim=3] mtx
    mtx[:, :, :] = 0.0
    cdef double cosa, sina
    cosa = numpy.cos(alpha * numpy.pi / 180.)
    sina = numpy.sin(alpha * numpy.pi / 180.)
    cdef numpy.ndarray[DTYPE_t, ndim=1] u
    u = numpy.array([cosa, 0, sina])
    cdef int i_surf, nx_, ny_, n_, n_bpts_, n_panels_, i_, i_bpts_, i_panels_, n
    for i_surf, row in enumerate(aero_ind):
        nx_, ny_, n_, n_bpts_, n_panels_, i_, i_bpts_, i_panels_ = row.copy()
        n = nx_ * ny_
        mesh = flat_mesh[i_: i_ + n_, :].reshape(nx_, ny_, 3)
        bpts = b_pts[i_bpts_: i_bpts_ + n_bpts_].reshape(nx_ - 1, ny_, 3)
        for i_points, row in enumerate(aero_ind):
            nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = row
            pts = points[i_panels: i_panels +
                         n_panels].reshape(nx - 1, ny - 1, 3)
            small_mat = numpy.zeros((n_panels, n_panels_, 3)).astype(DTYPE)
            if fortran_flag:
                small_mat[:, :, :] = lib.assembleaeromtx(ny, nx, ny_, nx_,
                                                         alpha, pts, bpts,
                                                         mesh, skip)
            else:
                # Spanwise loop through horseshoe elements
                for el_j in range(ny_ - 1):
                    el_loc_j = el_j * (nx_ - 1)
                    C_te = mesh[-1, el_j + 1, :]
                    D_te = mesh[-1, el_j + 0, :]
                    # Spanwise loop through control points
                    for cp_j in range(ny - 1):
                        cp_loc_j = cp_j * (nx - 1)
                        # Chordwise loop through control points
                        for cp_i in range(nx - 1):
                            cp_loc = cp_i + cp_loc_j
                            P = pts[cp_i, cp_j]
                            r1 = P - D_te
                            r2 = P - C_te
                            r1_mag = norm(r1)
                            r2_mag = norm(r2)
                            t1 = numpy.cross(u, r2) / \
                                (r2_mag * (r2_mag - u.dot(r2)))
                            t3 = numpy.cross(u, r1) / \
                                (r1_mag * (r1_mag - u.dot(r1)))
                            trailing = t1 - t3
                            edges = 0
                            # Chordwise loop through horseshoe elements
                            for el_i in reversed(range(nx_ - 1)):
                                el_loc = el_i + el_loc_j
                                A = bpts[el_i, el_j + 0, :]
                                B = bpts[el_i, el_j + 1, :]
                                if el_i == nx_ - 2:
                                    C = mesh[-1, el_j + 1, :]
                                    D = mesh[-1, el_j + 0, :]
                                else:
                                    C = bpts[el_i + 1, el_j + 1, :]
                                    D = bpts[el_i + 1, el_j + 0, :]
                                edges += calc_vorticity(B, C, P)
                                edges += calc_vorticity(D, A, P)
                                if skip and el_loc == cp_loc:
                                    small_mat[cp_loc, el_loc,
                                              :] = trailing + edges
                                else:
                                    bound = calc_vorticity(A, B, P)
                                    small_mat[cp_loc, el_loc,
                                              :] = trailing + edges + bound
            mtx[i_panels: i_panels + n_panels,
                i_panels_: i_panels_ + n_panels_, :] = small_mat
    mtx /= 4 * numpy.pi


def assemble_AIC_system(AIC_mtx, def_mesh, aero_ind, c_pts, b_pts, alpha, v, tot_panels, normals):
    mtx = numpy.zeros((tot_panels, tot_panels), dtype=DTYPE)
    rhs = numpy.zeros((tot_panels), dtype=DTYPE)
    assemble_AIC_mtx(AIC_mtx, def_mesh, aero_ind, c_pts, b_pts, alpha)
    for ind in range(3):
        mtx[:, :] += (AIC_mtx[:, :, ind].T * normals[:, ind]).T
    cdef double alpha, cosa, sina
    alpha = alpha * numpy.pi / 180.
    cosa = numpy.cos(alpha)
    sina = numpy.sin(alpha)
    v_inf = v * numpy.array([cosa, 0., sina], dtype=DTYPE)
    rhs[:] = -normals.reshape(-1, normals.shape[-1], order='F').dot(v_inf)
    return mtx, rhs


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
    tot_n = numpy.sum(aero_ind[:, 2])
    tot_bpts = numpy.sum(aero_ind[:, 3])
    tot_panels = numpy.sum(aero_ind[:, 4])
    AIC_mtx = numpy.zeros((tot_panels, tot_panels, 3), dtype=DTYPE)
    mtx = numpy.zeros((tot_panels, tot_panels), dtype=DTYPE)
    rhs = numpy.zeros((tot_panels), dtype=DTYPE)
    mtx, rhs = assemble_AIC_system(
        AIC_mtx, def_mesh, aero_ind, c_pts, b_pts, alpha, v, tot_panels, normals)
    circulations = numpy.linalg.solve(mtx, rhs)
    return circulations


def vlm_forces(def_mesh, aero_ind, b_pts, mid_b, circ, alpha=3, v=10, rho=3):
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
    cdef int tot_n, tot_bpts, tot_panels
    tot_n = numpy.sum(aero_ind[:, 2])
    tot_bpts = numpy.sum(aero_ind[:, 3])
    tot_panels = numpy.sum(aero_ind[:, 4])
    cdef numpy.ndarray[DTYPE, ndim=3] mtx
    mtx = numpy.zeros((tot_panels, tot_panels, 3), dtype=DTYPE)
    cdef numpy.ndarray[DTYPE, ndim=3] vel
    vel = numpy.zeros((tot_panels, 3), dtype=DTYPE)
    sec_forces = numpy.zeros((tot_panels, 3))
    cdef double alpha, cosa, sina
    alpha *= numpy.pi / 180.
    cosa = numpy.cos(alpha)
    sina = numpy.sin(alpha)
    assemble_AIC_mtx(mtx, def_mesh, aero_ind, mid_b, b_pts, alpha, skip=True)
    for i_surf, row in enumerate(aero_ind):
        nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = row

        for ind in range(3):
            vel[:, ind] = mtx[:, :, ind].dot(circ)
        vel[:, 0] += cosa * v
        vel[:, 2] += sina * v

        b_pts = b_pts[i_bpts: i_bpts + n_bpts, :].reshape(nx - 1, ny, 3)

        bound = b_pts[:, 1:, :] - b_pts[:, : -1, :]

        cross = numpy.cross(vel[i_panels: i_panels + n_panels],
                            bound.reshape(-1, bound.shape[-1], order='F'))

        for ind in range(3):
            sec_forces[i_panels: i_panels + n_panels,
                       ind] = (rho * circ[i_panels: i_panels + n_panels] * cross[:, ind])
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
    tot_n = numpy.sum(aero_ind[:, 2])
    tot_panels = numpy.sum(aero_ind[:, 4])
    tot_n_fem = numpy.sum(fem_ind[:, 0])
    output_loads = numpy.zeros((tot_n_fem, 6), dtype=DTYPE)
    for i_surf, row in enumerate(fem_ind):
        nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = aero_ind[i_surf, :]
        n_fem, i_fem = row
        mesh = def_mesh[i: i + n, :].reshape(nx, ny, 3)
        sec_forces = sec_forces[i_panels: i_panels + n_panels, : ]. \
            reshape(nx - 1, ny - 1, 3, order='F')
        sec_forces = numpy.sum(sec_forces, axis=0)
        w = 0.25
        a_pts = 0.5 * (1 - w) * mesh[: -1, : -1, : ] + \
            0.5 *   w   * mesh[1: , : -1, : ] + \
            0.5 * (1 - w) * mesh[: -1,  1: , : ] + \
            0.5 * w * mesh[1:,  1:, :]
        w = fem_origin
        s_pts = 0.5 * (1 - w) * mesh[: -1, : -1, : ] + \
            0.5 *   w   * mesh[1: , : -1, : ] + \
            0.5 * (1 - w) * mesh[: -1,  1: , : ] + \
            0.5 * w * mesh[1:,  1:, :]
        moment = numpy.zeros((ny - 1, 3), dtype=DTYPE)
        for ind in range(ny - 1):
            r = a_pts[0, ind, :] - s_pts[0, ind, :]
            F = sec_forces[ind, :]
            moment[ind, :] = numpy.cross(r, F)
        loads = numpy.zeros((ny, 6), dtype=DTYPE)
        loads[: -1, : 3] += 0.5 * sec_forces[:, :]
        loads[1:, : 3] += 0.5 * sec_forces[:, :]
        loads[: -1, 3:] += 0.5 * moment
        loads[1:, 3:] += 0.5 * moment
        output_loads[i_fem: i_fem + n_fem, :] = loads
    return output_loads


"""
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------


                                    STRUCTURES


--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
From spatialbeam.py: Define the structural analysis component using spatial beam theory. """


def norm(vec):
    return numpy.sqrt(numpy.sum(vec**2))


def unit(vec):
    return vec / norm(vec)


def radii(numpy.ndarray[DTYPE_t, ndim=3] mesh, double t_c=0.15):
    """ Obtain the radii of the FEM component based on chord. """
    # cdef numpy.ndarray[DTYPE_t, ndim=2] vectors
    # cdef double chords
    vectors = mesh[-1, :, :] - mesh[0, :, :]
    chords = numpy.sqrt(numpy.sum(vectors**2, axis=1))
    chords = 0.5 * chords[: -1] + 0.5 * chords[1:]
    return t_c * chords


def assemble_FEM_system(aero_ind, fem_ind, nodes, A, J, Iy, Iz, loads,
                        K_a, K_t, K_y, K_z, elem_IDs, cons, E, G, x_gl, T,
                        K_elem, S_a, S_t, S_y, S_z, T_elem,
                        const2, const_y, const_z, n, size, mtx, rhs):
    """
    Assemble the structural stiffness matrix based on 6 degrees of freedom
    per element.

    Can be run in dense Fortran, Sparse Fortran, or dense
    Python code depending on the flags used. Currently, dense Fortran
    seems to be the fastest version across many matrix sizes.

    """
    data_list = []
    rows_list = []
    cols_list = []
    num_surf = fem_ind.shape[0]
    tot_n_fem = numpy.sum(fem_ind[:, 0])
    size = 6 * tot_n_fem + 6 * num_surf
    for i_surf, row in enumerate(fem_ind):
        n_fem, i_fem = row
        # create truncated versions of the input arrays to assemble
        # smaller matrices that we later assemble into a full matrix
        num_cons = 1  # just one constraint per structural component
        size_ = 6 * (n_fem + num_cons)
        mtx_ = numpy.zeros((size_, size_), dtype=DTYPE)
        rhs_ = numpy.zeros((size_), dtype=DTYPE)
        A_ = A[i_fem - i_surf: i_fem - i_surf + n_fem - 1]
        J_ = J[i_fem - i_surf: i_fem - i_surf + n_fem - 1]
        Iy_ = Iy[i_fem - i_surf: i_fem - i_surf + n_fem - 1]
        Iz_ = Iz[i_fem - i_surf: i_fem - i_surf + n_fem - 1]
        elem_IDs_ = elem_IDs[i_fem - i_surf: i_fem -
                             i_surf + n_fem - 1, :] - i_fem
        loads_ = loads[i_fem: i_fem + n_fem]

        num_elems = elem_IDs_.shape[0]
        E_vec = E * numpy.ones(num_elems)
        G_vec = G * numpy.ones(num_elems)

        # dense Fortran
        if fortran_flag and not sparse_flag:
            mtx_ = lib.assemblestructmtx(nodes, A_, J_, Iy_, Iz_,
                                         K_a, K_t, K_y, K_z,
                                         elem_IDs_ + 1, cons[i_surf],
                                         E_vec, G_vec, x_gl, T,
                                         K_elem, S_a, S_t, S_y, S_z, T_elem,
                                         const2, const_y, const_z, n_fem,
                                         tot_n_fem, size_)
            mtx[(i_fem + i_surf) * 6: (i_fem + n_fem + num_cons + i_surf) * 6,
                (i_fem + i_surf) * 6: (i_fem + n_fem + num_cons + i_surf) * 6] = mtx_

            rhs_[:] = 0.0
            rhs_[: 6 * n_fem] = loads_.reshape((6 * n_fem))
            rhs[6 * (i_fem + i_surf): 6 *
                (i_fem + n_fem + i_surf + num_cons)] = rhs_

        # sparse Fortran
        elif fortran_flag and sparse_flag:
            nnz = 144 * num_elems

            data1, rows1, cols1 = lib.assemblesparsemtx(num_elems, tot_n_fem, nnz, x_gl, E_vec,
                                                        G_vec, A_, J_, Iy_, Iz_,
                                                        nodes, elem_IDs_ + 1, const2, const_y,
                                                        const_z, S_a, S_t, S_y, S_z)

            data2 = numpy.ones(6 * num_cons) * 1.e9
            rows2 = numpy.arange(6 * num_cons) + 6 * n_fem
            cols2 = numpy.zeros(6 * num_cons)
            for ind in range(6):
                cols2[ind:: 6] = 6 * cons[i_surf] + ind

            data = numpy.concatenate([data1, data2, data2])
            rows = numpy.concatenate(
                [rows1, rows2, cols2]) + (i_fem + i_surf) * 6
            cols = numpy.concatenate(
                [cols1, cols2, rows2]) + (i_fem + i_surf) * 6
            data_list.append(data)
            rows_list.append(rows)
            cols_list.append(cols)

            rhs_[:] = 0.0
            rhs_[: 6 * n_fem] = loads_.reshape((6 * n_fem))
            rhs[6 * (i_fem + i_surf): 6 *
                (i_fem + n_fem + i_surf + num_cons)] = rhs_

        # sparse Python
        elif not fortran_flag and sparse_flag:
            data = numpy.concatenate(data_list)
            rows = numpy.concatenate(rows_list)
            cols = numpy.concatenate(cols_list)
            mtx = scipy.sparse.csc_matrix((data, (rows, cols)),
                                          shape=(size, size))

        # dense Python
        else:
            num_nodes = num_elems + 1

            mtx_[:] = 0.
            for ielem in range(num_elems):
                P0 = nodes[elem_IDs_[ielem, 0], :]
                P1 = nodes[elem_IDs_[ielem, 1], :]

                x_loc = unit(P1 - P0)
                y_loc = unit(numpy.cross(x_loc, x_gl))
                z_loc = unit(numpy.cross(x_loc, y_loc))

                T[0, :] = x_loc
                T[1, :] = y_loc
                T[2, :] = z_loc

                for ind in range(4):
                    T_elem[3 * ind: 3 * ind + 3, 3 * ind: 3 * ind + 3] = T

                L = norm(P1 - P0)
                EA_L = E_vec[ielem] * A[ielem] / L
                GJ_L = G_vec[ielem] * J[ielem] / L
                EIy_L3 = E_vec[ielem] * Iy[ielem] / L**3
                EIz_L3 = E_vec[ielem] * Iz[ielem] / L**3

                K_a[:, :] = EA_L * const2
                K_t[:, :] = GJ_L * const2

                K_y[:, :] = EIy_L3 * const_y
                K_y[1, :] *= L
                K_y[3, :] *= L
                K_y[:, 1] *= L
                K_y[:, 3] *= L

                K_z[:, :] = EIz_L3 * const_z
                K_z[1, :] *= L
                K_z[3, :] *= L
                K_z[:, 1] *= L
                K_z[:, 3] *= L

                K_elem[:] = 0
                K_elem += S_a.T.dot(K_a).dot(S_a)
                K_elem += S_t.T.dot(K_t).dot(S_t)
                K_elem += S_y.T.dot(K_y).dot(S_y)
                K_elem += S_z.T.dot(K_z).dot(S_z)

                res = T_elem.T.dot(K_elem).dot(T_elem)

                in0, in1 = elem_IDs[ielem, :]

                mtx_[6 * in0: 6 * in0 + 6, 6 *
                     in0: 6 * in0 + 6] += res[: 6, : 6]
                mtx_[6 * in1: 6 * in1 + 6, 6 * in0: 6 * in0 + 6] += res[6:, : 6]
                mtx_[6 * in0: 6 * in0 + 6, 6 * in1: 6 * in1 + 6] += res[: 6, 6:]
                mtx_[6 * in1: 6 * in1 + 6, 6 * in1: 6 * in1 + 6] += res[6:, 6:]

            for ind in range(num_cons):
                for k in range(6):
                    mtx_[6 * num_nodes + 6 * ind +
                         k, 6 * cons[i_surf] + k] = 1.e9
                    mtx_[6 * cons[i_surf] + k, 6 *
                         num_nodes + 6 * ind + k] = 1.e9

            rhs_[:] = 0.0
            rhs_[: 6 * n_fem] = loads_.reshape((6 * n_fem))
            rhs[6 * (i_fem + i_surf): 6 *
                (i_fem + n_fem + i_surf + num_cons)] = rhs_

            mtx[(i_fem + i_surf) * 6: (i_fem + n_fem + num_cons + i_surf) * 6,
                (i_fem + i_surf) * 6: (i_fem + n_fem + num_cons + i_surf) * 6] = mtx_

    rhs[numpy.abs(rhs) < 1e-6] = 0.  # *should this have lower tolerance?
    return mtx, rhs


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
    n_fem, i_fem = fem_ind[0, :]
    num_surf = fem_ind.shape[0]
    tot_n_fem = numpy.sum(fem_ind[:, 0])
    size = 6 * tot_n_fem + 6 * num_surf
    elem_IDs = numpy.zeros((tot_n_fem - num_surf, 2), int)
    for i_surf, row in enumerate(fem_ind):
        nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = aero_ind[i_surf, :]
        n_fem, i_fem = row
        arange = numpy.arange(n_fem - 1) + i_fem
        elem_IDs[i_fem - i_surf: i_fem - i_surf + n_fem - 1, 0] = arange
        elem_IDs[i_fem - i_surf: i_fem - i_surf + n_fem - 1, 1] = arange + 1
        elem_IDs = elem_IDs
    const2 = numpy.array([
        [1, -1],
        [-1, 1],
    ], dtype=DTYPE)
    const_y = numpy.array([
        [12, -6, -12, -6],
        [-6, 4, 6, 2],
        [-12, 6, 12, 6],
        [-6, 2, 6, 4],
    ], dtype=DTYPE)
    const_z = numpy.array([
        [12, 6, -12, 6],
        [6, 4, -6, 2],
        [-12, -6, 12, -6],
        [6, 2, -6, 4],
    ], dtype=DTYPE)
    x_gl = numpy.array([1, 0, 0], dtype=DTYPE)
    K_elem = numpy.zeros((12, 12), dtype=DTYPE)
    T_elem = numpy.zeros((12, 12), dtype=DTYPE)
    T = numpy.zeros((3, 3), dtype=DTYPE)
    num_nodes = tot_n_fem
    num_cons = num_surf
    size = 6 * num_nodes + 6 * num_cons
    mtx = numpy.zeros((size, size), dtype=DTYPE)
    rhs = numpy.zeros(size, dtype=DTYPE)
    K_a = numpy.zeros((2, 2), dtype=DTYPE)
    K_t = numpy.zeros((2, 2), dtype=DTYPE)
    K_y = numpy.zeros((4, 4), dtype=DTYPE)
    K_z = numpy.zeros((4, 4), dtype=DTYPE)
    S_a = numpy.zeros((2, 12), dtype=DTYPE)
    S_a[(0, 1), (0, 6)] = 1.
    S_t = numpy.zeros((2, 12), dtype=DTYPE)
    S_t[(0, 1), (3, 9)] = 1.
    S_y = numpy.zeros((4, 12), dtype=DTYPE)
    S_y[(0, 1, 2, 3), (2, 4, 8, 10)] = 1.
    S_z = numpy.zeros((4, 12), dtype=DTYPE)
    S_z[(0, 1, 2, 3), (1, 5, 7, 11)] = 1.
    cons = numpy.zeros((num_surf))
    # find constrained nodes based on closeness to specified cg point
    for i_surf, row in enumerate(fem_ind):
        n_fem, i_fem = row
        nodes = nodes[i_fem:i_fem + n_fem]
        dist = nodes - numpy.array([cg_x, 0, 0])
        idx = (numpy.linalg.norm(dist, axis=1)).argmin()
        cons[i_surf] = idx
    mtx, rhs = \
        assemble_FEM_system(aero_ind, fem_ind, nodes, A, J, Iy, Iz, loads,
                            K_a, K_t, K_y, K_z, elem_IDs, cons, E, G, x_gl, T,
                            K_elem, S_a, S_t, S_y, S_z, T_elem, const2, const_y,
                            const_z, n_fem, size, mtx, rhs)
    if type(mtx) == numpy.ndarray:
        # print('AAA')
        disp_aug = numpy.linalg.solve(mtx, rhs)
    else:
        # print('BBB')
        splu = scipy.sparse.linalg.splu(mtx)
        disp_aug = splu.solve(rhs)
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
    num_surf = fem_ind.shape[0]
    tot_n_fem = numpy.sum(fem_ind[:, 0])
    size = 6 * tot_n_fem + 6 * num_surf
    disp = numpy.zeros((tot_n_fem, 6))
    arange = numpy.arange(6 * tot_n_fem)
    for i_surf, row in enumerate(fem_ind):
        n_fem, i_fem = row
        disp[i_fem:i_fem + n_fem] = disp_aug[(i_fem + i_surf) * 6:(
            i_fem + n_fem + i_surf) * 6].reshape((n_fem, 6))
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
    num_surf = fem_ind.shape[0]
    tot_n_fem = numpy.sum(fem_ind[:, 0])
    tot_n = numpy.sum(aero_ind[:, 2])
    nodes = numpy.zeros((tot_n_fem, 3))
    for i_surf, row in enumerate(fem_ind):
        nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = aero_ind[i_surf, :]
        n_fem, i_fem = row
        mesh2 = mesh[i:i + n, :].reshape(nx, ny, 3)
        nodes[i_fem:i_fem + n_fem] = (1 - fem_origin) * \
            mesh2[0, :, :] + fem_origin * mesh2[-1, :, :]
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
    n_fem, i_fem = fem_ind[0, :]
    num_surf = fem_ind.shape[0]
    tot_n_fem = numpy.sum(fem_ind[:, 0])
    size = 6 * tot_n_fem + 6 * num_surf
    A = numpy.zeros((tot_n_fem - num_surf))
    Iy = numpy.zeros((tot_n_fem - num_surf))
    Iz = numpy.zeros((tot_n_fem - num_surf))
    J = numpy.zeros((tot_n_fem - num_surf))
    r1 = r - 0.5 * thickness
    r2 = r + 0.5 * thickness
    A = numpy.pi * (r2**2 - r1**2)
    Iy = numpy.pi * (r2**4 - r1**4) / 4.
    Iz = numpy.pi * (r2**4 - r1**4) / 4.
    J = numpy.pi * (r2**4 - r1**4) / 2.
    return A, Iy, Iz, J
