"""
Define the aerodynamic analysis component using a vortex lattice method.

We input a nodal mesh and properties of the airflow to calculate the
circulations of the horseshoe vortices. We then compute the forces, lift,
and drag acting on the lifting surfaces. Currently we can compute the induced
and viscous drag.

"""

from __future__ import division, print_function
import numpy as np

from openmdao.api import Component, Group
from scipy.linalg import lu_factor, lu_solve

try:
    import OAS_API
    fortran_flag = True
    data_type = float
except:
    fortran_flag = False
    data_type = complex


def norm(vec):
    """ Finds the 2-norm of a vector. """
    return np.sqrt(np.sum(vec**2))


def _calc_vorticity(A, B, P):
    """ Calculates the influence coefficient for a vortex filament.

    Parameters
    ----------
    A[3] : numpy array
        Coordinates for the start point of the filament.
    B[3] : numpy array
        Coordinates for the end point of the filament.
    P[3] : numpy array
        Coordinates for the collocation point where the influence coefficient
        is computed.

    Returns
    -------
    out[3] : numpy array
        Influence coefficient contribution for the described filament.

    """

    r1 = P - A
    r2 = P - B

    r1_mag = norm(r1)
    r2_mag = norm(r2)

    return (r1_mag + r2_mag) * np.cross(r1, r2) / \
           (r1_mag * r2_mag * (r1_mag * r2_mag + r1.dot(r2)))


def _assemble_AIC_mtx(mtx, params, surfaces, skip=False):
    """
    Compute the aerodynamic influence coefficient matrix
    for either solving the linear system or solving for the drag.

    We use a nested for loop structure to loop through the lifting surfaces to
    obtain the corresponding mesh, then for each mesh we again loop through
    the lifting surfaces to obtain the collocation points used to compute
    the horseshoe vortex influence coefficients.

    This creates mtx with blocks corresponding to each lifting surface's
    effects on other lifting surfaces. The block diagonal portions
    correspond to each lifting surface's influencen on itself. For a single
    lifting surface, this is the entire mtx.

    Parameters
    ----------
    mtx[(nx-1)*(ny-1), (nx-1)*(ny-1), 3] : numpy array
        Aerodynamic influence coefficient (AIC) matrix, or the
        derivative of v w.r.t. circulations.
    params : dictionary
        OpenMDAO params dictionary for a given aero problem
    surfaces : dictionary
        Dictionary containing all surfaces in an aero problem.
    skip : boolean
        If false, the bound vortex contributions on the collocation point
        corresponding to the same panel are not included. Used for the drag
        computation.

    Returns
    -------
    mtx[tot_panels, tot_panels, 3] : numpy array
        Aerodynamic influence coefficient (AIC) matrix, or the
        derivative of v w.r.t. circulations.
    """

    alpha = params['alpha']
    mtx[:, :, :] = 0.0
    cosa = np.cos(alpha * np.pi / 180.)
    sina = np.sin(alpha * np.pi / 180.)
    u = np.array([cosa, 0, sina])

    i_ = 0
    i_bpts_ = 0
    i_panels_ = 0

    # Loop over the lifting surfaces to compute their influence on the flow
    # velocity at the collocation points
    for surface_ in surfaces:

        # Variable names with a trailing underscore correspond to the lifting
        # surface being examined, not the collocation point
        name_ = surface_['name']
        nx_ = surface_['num_x']
        ny_ = surface_['num_y']
        n_ = nx_ * ny_
        n_bpts_ = (nx_ - 1) * ny_
        n_panels_ = (nx_ - 1) * (ny_ - 1)

        # Obtain the lifting surface mesh in the form expected by the solver,
        # with shape [nx_, ny_, 3]
        mesh = params[name_+'def_mesh']
        bpts = params[name_+'b_pts']

        # Set counters to know where to index the sub-matrix within the full mtx
        i = 0
        i_bpts = 0
        i_panels = 0

        for surface in surfaces:
            # These variables correspond to the collocation points
            name = surface['name']
            nx = surface['num_x']
            ny = surface['num_y']
            n = nx * ny
            n_bpts = (nx - 1) * ny
            n_panels = (nx - 1) * (ny - 1)
            symmetry = surface['symmetry']

            # Obtain the collocation points used to compute the AIC mtx.
            # If setting up the AIC mtx, we use the collocation points (c_pts),
            # but if setting up the matrix to solve for drag, we use the
            # midpoints of the bound vortices.
            if skip:
                # Find the midpoints of the bound points, used in drag computations
                pts = (params[name+'b_pts'][:, 1:, :] + \
                    params[name+'b_pts'][:, :-1, :]) / 2
            else:
                pts = params[name+'c_pts']

            # Initialize sub-matrix to populate within full mtx
            small_mat = np.zeros((n_panels, n_panels_, 3), dtype=data_type)

            # Dense fortran assembly for the AIC matrix
            if fortran_flag:
                small_mat[:, :, :] = OAS_API.oas_api.assembleaeromtx(alpha, pts, bpts,
                                                         mesh, skip, symmetry)
            # Python matrix assembly
            else:
                # Spanwise loop through horseshoe elements
                for el_j in range(ny_ - 1):
                    el_loc_j = el_j * (nx_ - 1)
                    C_te = mesh[-1, el_j + 1, :]
                    D_te = mesh[-1, el_j + 0, :]

                    # Mirror the horseshoe vortex points
                    if symmetry:
                        C_te_sym = C_te.copy()
                        D_te_sym = D_te.copy()
                        C_te_sym[1] = -C_te_sym[1]
                        D_te_sym[1] = -D_te_sym[1]

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

                            t1 = np.cross(u, r2) / \
                                (r2_mag * (r2_mag - u.dot(r2)))
                            t3 = np.cross(u, r1) / \
                                (r1_mag * (r1_mag - u.dot(r1)))

                            # AIC contribution from trailing vortex filaments
                            # coming off the trailing edge
                            trailing = t1 - t3

                            # Calculate the effects across the symmetry plane
                            if symmetry:
                                r1 = P - D_te_sym
                                r2 = P - C_te_sym

                                r1_mag = norm(r1)
                                r2_mag = norm(r2)

                                t1 = np.cross(u, r2) / \
                                    (r2_mag * (r2_mag - u.dot(r2)))
                                t3 = np.cross(u, r1) / \
                                    (r1_mag * (r1_mag - u.dot(r1)))

                                trailing += t3 - t1

                            edges = 0

                            # Chordwise loop through horseshoe elements in
                            # reversed order, starting with the panel closest
                            # to the leading edge. This is done to sum the
                            # AIC contributions from the side vortex filaments
                            # as we loop through the elements
                            for el_i in reversed(range(nx_ - 1)):
                                el_loc = el_i + el_loc_j

                                A = bpts[el_i, el_j + 0, :]
                                B = bpts[el_i, el_j + 1, :]

                                # Check if this is the last panel; if so, use
                                # the trailing edge mesh points for C & D, else
                                # use the directly aft panel's bound points
                                # for C & D
                                if el_i == nx_ - 2:
                                    C = mesh[-1, el_j + 1, :]
                                    D = mesh[-1, el_j + 0, :]
                                else:
                                    C = bpts[el_i + 1, el_j + 1, :]
                                    D = bpts[el_i + 1, el_j + 0, :]

                                # Calculate and store the contributions from
                                # the vortex filaments on the sides of the
                                # panels, adding as we progress through the
                                # panels
                                edges += _calc_vorticity(B, C, P)
                                edges += _calc_vorticity(D, A, P)

                                # Mirror the horseshoe vortex points and
                                # calculate the effects across
                                # the symmetry plane
                                if symmetry:
                                    A_sym = A.copy()
                                    B_sym = B.copy()
                                    C_sym = C.copy()
                                    D_sym = D.copy()
                                    A_sym[1] = -A_sym[1]
                                    B_sym[1] = -B_sym[1]
                                    C_sym[1] = -C_sym[1]
                                    D_sym[1] = -D_sym[1]

                                    edges += _calc_vorticity(C_sym, B_sym, P)
                                    edges += _calc_vorticity(A_sym, D_sym, P)

                                # If skip, do not include the contributions
                                # from the panel's bound vortex filament, as
                                # this causes a singularity when we're taking
                                # the influence of a panel on its own
                                # collocation point. This true for the drag
                                # computation and false for circulation
                                # computation, due to the different collocation
                                # points.
                                if skip and el_loc == cp_loc:
                                    if symmetry:
                                        bound = _calc_vorticity(B_sym, A_sym, P)
                                    else:
                                        bound = np.zeros((3))
                                    small_mat[cp_loc, el_loc, :] = \
                                        trailing + edges + bound
                                else:
                                    bound = _calc_vorticity(A, B, P)

                                    # Account for symmetry by including the
                                    # mirrored bound vortex
                                    if symmetry:
                                        bound += _calc_vorticity(B_sym, A_sym, P)

                                    small_mat[cp_loc, el_loc, :] = \
                                        trailing + edges + bound

            # Populate the full-size matrix with these surface-surface AICs
            mtx[i_panels:i_panels+n_panels,
                i_panels_:i_panels_+n_panels_, :] = small_mat

            i += n
            i_bpts += n_bpts
            i_panels += n_panels

        i_ += n_
        i_bpts_ += n_bpts_
        i_panels_ += n_panels_

    mtx /= 4 * np.pi

def _assemble_AIC_mtx_d(mtxd, params, dparams, dunknowns, dresids, surfaces, skip=False):

    alpha = params['alpha']
    alphad = dparams['alpha']

    i_ = 0
    i_bpts_ = 0
    i_panels_ = 0

    # Loop over the lifting surfaces to compute their influence on the flow
    # velocity at the collocation points
    for surface_ in surfaces:

        # Variable names with a trailing underscore correspond to the lifting
        # surface being examined, not the collocation point
        name_ = surface_['name']
        nx_ = surface_['num_x']
        ny_ = surface_['num_y']
        n_ = nx_ * ny_
        n_bpts_ = (nx_ - 1) * ny_
        n_panels_ = (nx_ - 1) * (ny_ - 1)

        # Obtain the lifting surface mesh in the form expected by the solver,
        # with shape [nx_, ny_, 3]
        mesh = params[name_+'def_mesh']
        bpts = params[name_+'b_pts']

        meshd = dparams[name_+'def_mesh']
        bptsd = dparams[name_+'b_pts']

        # Set counters to know where to index the sub-matrix within the full mtx
        i = 0
        i_bpts = 0
        i_panels = 0

        for surface in surfaces:
            # These variables correspond to the collocation points
            name = surface['name']
            nx = surface['num_x']
            ny = surface['num_y']
            n = nx * ny
            n_bpts = (nx - 1) * ny
            n_panels = (nx - 1) * (ny - 1)
            symmetry = surface['symmetry']

            # Obtain the collocation points used to compute the AIC mtx.
            # If setting up the AIC mtx, we use the collocation points (c_pts),
            # but if setting up the matrix to solve for drag, we use the
            # midpoints of the bound vortices.
            if skip:
                # Find the midpoints of the bound points, used in drag computations
                pts = (params[name+'b_pts'][:, 1:, :] + \
                    params[name+'b_pts'][:, :-1, :]) / 2
                ptsd = (dparams[name+'b_pts'][:, 1:, :] + \
                    dparams[name+'b_pts'][:, :-1, :]) / 2
            else:
                pts = params[name+'c_pts']
                ptsd = dparams[name+'c_pts']

            _, small_mat = OAS_API.oas_api.assembleaeromtx_d(alpha, alphad, pts, ptsd,
                                                          bpts, bptsd, mesh, meshd,
                                                          skip, symmetry)

            # Populate the full-size matrix with these surface-surface AICs
            mtxd[i_panels:i_panels+n_panels,
                 i_panels_:i_panels_+n_panels_, :] = small_mat

            i += n
            i_bpts += n_bpts
            i_panels += n_panels

        i_ += n_
        i_bpts_ += n_bpts_
        i_panels_ += n_panels_

    mtxd /= 4 * np.pi

def _assemble_AIC_mtx_b(mtxb, params, dparams, dunknowns, dresids, surfaces, skip=False):

    alpha = params['alpha']

    mtxb /= 4 * np.pi

    i_ = 0
    i_bpts_ = 0
    i_panels_ = 0

    # Loop over the lifting surfaces to compute their influence on the flow
    # velocity at the collocation points
    for surface_ in surfaces:

        # Variable names with a trailing underscore correspond to the lifting
        # surface being examined, not the collocation point
        name_ = surface_['name']
        nx_ = surface_['num_x']
        ny_ = surface_['num_y']
        n_ = nx_ * ny_
        n_bpts_ = (nx_ - 1) * ny_
        n_panels_ = (nx_ - 1) * (ny_ - 1)

        # Obtain the lifting surface mesh in the form expected by the solver,
        # with shape [nx_, ny_, 3]
        mesh = params[name_+'def_mesh']
        bpts = params[name_+'b_pts']

        # Set counters to know where to index the sub-matrix within the full mtx
        i = 0
        i_bpts = 0
        i_panels = 0

        for surface in surfaces:
            # These variables correspond to the collocation points
            name = surface['name']
            nx = surface['num_x']
            ny = surface['num_y']
            n = nx * ny
            n_bpts = (nx - 1) * ny
            n_panels = (nx - 1) * (ny - 1)
            symmetry = surface['symmetry']

            # Obtain the collocation points used to compute the AIC mtx.
            # If setting up the AIC mtx, we use the collocation points (c_pts),
            # but if setting up the matrix to solve for drag, we use the
            # midpoints of the bound vortices.
            if skip:
                # Find the midpoints of the bound points, used in drag computations
                pts = (params[name+'b_pts'][:, 1:, :] + \
                    params[name+'b_pts'][:, :-1, :]) / 2
            else:
                pts = params[name+'c_pts']

            small_mtxb = mtxb[i_panels:i_panels+n_panels, i_panels_:i_panels_+n_panels_, :]

            alphab, ptsb, bptsb, meshb, mtx = OAS_API.oas_api.assembleaeromtx_b(alpha, pts, bpts,
                                                     mesh, skip, symmetry, small_mtxb)

            dparams[name_+'def_mesh'] += meshb.real
            dparams[name_+'b_pts'] += bptsb.real

            if skip:
                dparams[name+'b_pts'][:, 1:, :] += ptsb.real / 2
                dparams[name+'b_pts'][:, :-1, :] += ptsb.real / 2
            else:
                dparams[name+'c_pts'] += ptsb.real
            dparams['alpha'] += alphab

            i += n
            i_bpts += n_bpts
            i_panels += n_panels

        i_ += n_
        i_bpts_ += n_bpts_
        i_panels_ += n_panels_



class VLMGeometry(Component):
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

    def __init__(self, surface):
        super(VLMGeometry, self).__init__()

        self.surface = surface

        ny = surface['num_y']
        nx = surface['num_x']

        self.fem_origin = surface['fem_origin']

        self.add_param('def_mesh', val=np.zeros((nx, ny, 3),
                       dtype=data_type))
        self.add_output('b_pts', val=np.zeros((nx-1, ny, 3),
                        dtype=data_type))
        self.add_output('c_pts', val=np.zeros((nx-1, ny-1, 3)))
        self.add_output('widths', val=np.zeros((ny-1)))
        self.add_output('cos_sweep', val=np.zeros((ny-1)))
        self.add_output('lengths', val=np.zeros((ny)))
        self.add_output('normals', val=np.zeros((nx-1, ny-1, 3)))
        self.add_output('S_ref', val=0.)

    def solve_nonlinear(self, params, unknowns, resids):
        mesh = params['def_mesh']

        # Compute the bound points at 1/4 chord
        b_pts = mesh[:-1, :, :] * .75 + mesh[1:, :, :] * .25

        # Compute the collocation points at the midpoints of each
        # panel's 3/4 chord line
        c_pts = 0.5 * 0.25 * mesh[:-1, :-1, :] + \
                0.5 * 0.75 * mesh[1:, :-1, :] + \
                0.5 * 0.25 * mesh[:-1,  1:, :] + \
                0.5 * 0.75 * mesh[1:,  1:, :]

        # Compute the widths of each panel
        qc = 0.25 * mesh[-1] + 0.75 * mesh[0]
        widths = np.linalg.norm(qc[1:, :] - qc[:-1, :], axis=1)

        # Compute the numerator of the cosine of the sweep angle of each panel
        # (we need this for the viscous drag dependence on sweep, and we only compute
        # the numerator because the denominator of the cosine fraction is the width,
        # which we have already computed. They are combined in the viscous drag
        # calculation.)
        cos_sweep = np.linalg.norm(qc[1:, [1,2]] - qc[:-1, [1,2]], axis=1)

        # Compute the cambered length of each chordwise set of mesh points
        dx = mesh[1:, :, 0] - mesh[:-1, :, 0]
        dy = mesh[1:, :, 1] - mesh[:-1, :, 1]
        dz = mesh[1:, :, 2] - mesh[:-1, :, 2]
        lengths = np.sum(np.sqrt(dx**2 + dy**2 + dz**2), axis=0)

        # Compute the normal of each panel by taking the cross-product of
        # its diagonals. Note that this could be a nonplanar surface
        normals = np.cross(
            mesh[:-1,  1:, :] - mesh[1:, :-1, :],
            mesh[:-1, :-1, :] - mesh[1:,  1:, :],
            axis=2)

        norms = np.sqrt(np.sum(normals**2, axis=2))
        for j in range(3):
            normals[:, :, j] /= norms

        if self.surface['S_ref_type'] == 'wetted':
            S_ref = 0.5 * np.sum(norms)

        elif self.surface['S_ref_type'] == 'projected':
            proj_mesh = mesh.copy()
            proj_mesh[:,:,2] = 0.
            proj_normals = np.cross(
                proj_mesh[:-1,  1:, :] - proj_mesh[1:, :-1, :],
                proj_mesh[:-1, :-1, :] - proj_mesh[1:,  1:, :],
                axis=2)

            proj_norms = np.sqrt(np.sum(proj_normals**2, axis=2))
            for j in xrange(3):
                proj_normals[:, :, j] /= proj_norms

            S_ref = 0.5 * np.sum(proj_norms)

        # Store each array
        unknowns['b_pts'] = b_pts
        unknowns['c_pts'] = c_pts
        unknowns['widths'] = widths
        unknowns['cos_sweep'] = cos_sweep
        unknowns['lengths'] = lengths
        unknowns['normals'] = normals
        unknowns['S_ref'] = S_ref

    def linearize(self, params, unknowns, resids):
        """ Jacobian for geometry."""

        jac = self.alloc_jacobian()
        name = self.surface['name']
        mesh = params['def_mesh']

        nx = self.surface['num_x']
        ny = self.surface['num_y']

        if fortran_flag:

            normalsb = np.zeros(unknowns['normals'].shape)
            for i in range(nx-1):
                for j in range(ny-1):
                    for ind in range(3):
                        normalsb[:, :, :] = 0.
                        normalsb[i, j, ind] = 1.
                        meshb, _, _ = OAS_API.oas_api.compute_normals_b(mesh, normalsb, 0.)
                        jac['normals', 'def_mesh'][i*(ny-1)*3 + j*3 + ind, :] = meshb.flatten()

            normalsb[:, :, :] = 0.
            if self.surface['S_ref_type'] == 'wetted':
                seed_mesh = mesh
            elif self.surface['S_ref_type'] == 'projected':
                seed_mesh = mesh.copy()
                seed_mesh[:,:,2] = 0.
            meshb, _, _ = OAS_API.oas_api.compute_normals_b(seed_mesh, normalsb, 1.)
            jac['S_ref', 'def_mesh'] = np.atleast_2d(meshb.flatten())

        else:
            cs_jac = self.complex_step_jacobian(params, unknowns, resids,
                                            fd_params=['def_mesh'],
                                            fd_unknowns=['normals', 'S_ref'],
                                            fd_states=[])
            jac.update(cs_jac)

        for iz, v in zip((0, ny*3), (.75, .25)):
            np.fill_diagonal(jac['b_pts', 'def_mesh'][:, iz:], v)

        for iz, v in zip((0, 3, ny*3, (ny+1)*3),
                         (.125, .125, .375, .375)):
            for ix in range(nx-1):
                np.fill_diagonal(jac['c_pts', 'def_mesh']
                    [(ix*(ny-1))*3:((ix+1)*(ny-1))*3, iz+ix*ny*3:], v)

        jac['widths', 'def_mesh'] = np.zeros_like(jac['widths', 'def_mesh'])
        jac['cos_sweep', 'def_mesh'] = np.zeros_like(jac['cos_sweep', 'def_mesh'])
        qc = 0.25 * mesh[-1] + 0.75 * mesh[0]
        gap = [0, (nx-1)*ny*3]
        factor = [0.75, 0.25]
        for i in range(ny-1):
            w = unknowns['widths'][i]
            cos_sweep = unknowns['cos_sweep'][i]
            dx = (qc[i+1, 0] - qc[i, 0])
            dy = (qc[i+1, 1] - qc[i, 1])
            dz = (qc[i+1, 2] - qc[i, 2])
            for j in range(2):
                jac['widths', 'def_mesh'][i, i*3+gap[j]] -= dx * factor[j] / w
                jac['widths', 'def_mesh'][i, (i+1)*3+gap[j]] += dx * factor[j] / w
                jac['widths', 'def_mesh'][i, i*3+1+gap[j]] -= dy * factor[j] / w
                jac['widths', 'def_mesh'][i, (i+1)*3+1+gap[j]] += dy * factor[j] / w
                jac['widths', 'def_mesh'][i, i*3+2+gap[j]] -= dz * factor[j] / w
                jac['widths', 'def_mesh'][i, (i+1)*3+2+gap[j]] += dz * factor[j] / w
                jac['cos_sweep', 'def_mesh'][i, i*3+1+gap[j]] -= dy / cos_sweep * factor[j]
                jac['cos_sweep', 'def_mesh'][i, (i+1)*3+1+gap[j]] += dy / cos_sweep * factor[j]
                jac['cos_sweep', 'def_mesh'][i, i*3+2+gap[j]] -= dz / cos_sweep * factor[j]
                jac['cos_sweep', 'def_mesh'][i, (i+1)*3+2+gap[j]] += dz / cos_sweep * factor[j]

        jac['lengths', 'def_mesh'] = np.zeros_like(jac['lengths', 'def_mesh'])
        for i in range(ny):
            dx = mesh[1:, i, 0] - mesh[:-1, i, 0]
            dy = mesh[1:, i, 1] - mesh[:-1, i, 1]
            dz = mesh[1:, i, 2] - mesh[:-1, i, 2]
            for j in range(nx-1):
                l = np.sqrt(dx[j]**2 + dy[j]**2 + dz[j]**2)
                jac['lengths', 'def_mesh'][i, (j*ny+i)*3] -= dx[j] / l
                jac['lengths', 'def_mesh'][i, ((j+1)*ny+i)*3] += dx[j] / l
                jac['lengths', 'def_mesh'][i, (j*ny+i)*3 + 1] -= dy[j] / l
                jac['lengths', 'def_mesh'][i, ((j+1)*ny+i)*3 + 1] += dy[j] / l
                jac['lengths', 'def_mesh'][i, (j*ny+i)*3 + 2] -= dz[j] / l
                jac['lengths', 'def_mesh'][i, ((j+1)*ny+i)*3 + 2] += dz[j] / l

        return jac


class AssembleAIC(Component):
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

    def __init__(self, surfaces):
        super(AssembleAIC, self).__init__()

        self.surfaces = surfaces

        tot_panels = 0
        for surface in surfaces:
            self.surface = surface
            ny = surface['num_y']
            nx = surface['num_x']
            name = surface['name']

            self.add_param(name+'def_mesh', val=np.zeros((nx, ny, 3),
                           dtype=data_type))
            self.add_param(name+'b_pts', val=np.zeros((nx-1, ny, 3),
                           dtype=data_type))
            self.add_param(name+'c_pts', val=np.zeros((nx-1, ny-1, 3),
                           dtype=data_type))
            self.add_param(name+'normals', val=np.zeros((nx-1, ny-1, 3)))
            tot_panels += (nx - 1) * (ny - 1)

        self.tot_panels = tot_panels

        self.add_param('v', val=1.)
        self.add_param('alpha', val=0.)

        self.add_output('AIC', val=np.zeros((tot_panels, tot_panels), dtype=data_type))
        self.add_output('rhs', val=np.zeros((tot_panels), dtype=data_type))

        self.AIC_mtx = np.zeros((tot_panels, tot_panels, 3),
                                   dtype=data_type)
        self.mtx = np.zeros((tot_panels, tot_panels),
                                   dtype=data_type)

        if not fortran_flag:
            self.deriv_options['type'] = 'cs'
            self.deriv_options['form'] = 'central'

    def solve_nonlinear(self, params, unknowns, resids):
        # Actually assemble the AIC matrix
        _assemble_AIC_mtx(self.AIC_mtx, params, self.surfaces)

        # Construct an flattened array with the normals of each surface in order
        # so we can do the normals with velocities to set up the right-hand-side
        # of the system.
        flattened_normals = np.zeros((self.tot_panels, 3), dtype=data_type)
        i = 0
        for surface in self.surfaces:
            name = surface['name']
            num_panels = (surface['num_x'] - 1) * (surface['num_y'] - 1)
            flattened_normals[i:i+num_panels, :] = params[name+'normals'].reshape(-1, 3, order='F')
            i += num_panels

        # Construct a matrix that is the AIC_mtx dotted by the normals at each
        # collocation point. This is used to compute the circulations
        self.mtx[:, :] = 0.
        for ind in range(3):
            self.mtx[:, :] += (self.AIC_mtx[:, :, ind].T *
                flattened_normals[:, ind]).T

        # Obtain the freestream velocity direction and magnitude by taking
        # alpha into account
        alpha = params['alpha'] * np.pi / 180.
        cosa = np.cos(alpha)
        sina = np.sin(alpha)
        v_inf = params['v'] * np.array([cosa, 0., sina], dtype=data_type)

        # Populate the right-hand side of the linear system with the
        # expected velocities at each collocation point
        unknowns['rhs'] = -flattened_normals.\
            reshape(-1, flattened_normals.shape[-1], order='F').dot(v_inf)

        unknowns['AIC'] = self.mtx

    def apply_linear(self, params, unknowns, dparams, dunknowns, dresids, mode):

        if mode == 'fwd':

            AIC_mtxd = np.zeros(self.AIC_mtx.shape)

            # Actually assemble the AIC matrix
            _assemble_AIC_mtx_d(AIC_mtxd, params, dparams, dunknowns, dresids, self.surfaces)

            # Construct an flattened array with the normals of each surface in order
            # so we can do the normals with velocities to set up the right-hand-side
            # of the system.
            flattened_normals = np.zeros((self.tot_panels, 3))
            flattened_normalsd = np.zeros((self.tot_panels, 3))
            i = 0
            for surface in self.surfaces:
                name = surface['name']
                num_panels = (surface['num_x'] - 1) * (surface['num_y'] - 1)
                flattened_normals[i:i+num_panels, :] = params[name+'normals'].reshape(-1, 3, order='F')
                flattened_normalsd[i:i+num_panels, :] = dparams[name+'normals'].reshape(-1, 3, order='F')
                i += num_panels

            # Construct a matrix that is the AIC_mtx dotted by the normals at each
            # collocation point. This is used to compute the circulations
            self.mtx[:, :] = 0.
            for ind in range(3):
                self.mtx[:, :] += (AIC_mtxd[:, :, ind].T *
                    flattened_normals[:, ind]).T
                self.mtx[:, :] += (self.AIC_mtx[:, :, ind].T *
                    flattened_normalsd[:, ind]).T

            # Obtain the freestream velocity direction and magnitude by taking
            # alpha into account
            alpha = params['alpha'] * np.pi / 180.
            alphad = dparams['alpha'] * np.pi / 180.
            cosa = np.cos(alpha)
            sina = np.sin(alpha)
            cosad = -sina * alphad
            sinad = cosa * alphad

            freestream_direction = np.array([cosa, 0., sina])
            v_inf = params['v'] * freestream_direction
            v_infd = dparams['v'] * freestream_direction
            v_infd += params['v'] * np.array([cosad, 0., sinad])

            # Populate the right-hand side of the linear system with the
            # expected velocities at each collocation point
            dresids['rhs'] = -flattened_normalsd.\
                reshape(-1, 3, order='F').dot(v_inf)
            dresids['rhs'] += -flattened_normals.\
                reshape(-1, 3, order='F').dot(v_infd)

            dresids['AIC'] = self.mtx

        if mode == 'rev':

            # Construct an flattened array with the normals of each surface in order
            # so we can do the normals with velocities to set up the right-hand-side
            # of the system.
            flattened_normals = np.zeros((self.tot_panels, 3))
            i = 0
            for surface in self.surfaces:
                name = surface['name']
                num_panels = (surface['num_x'] - 1) * (surface['num_y'] - 1)
                flattened_normals[i:i+num_panels, :] = params[name+'normals'].reshape(-1, 3, order='F')
                i += num_panels

            AIC_mtxb = np.zeros((self.tot_panels, self.tot_panels, 3))
            flattened_normalsb = np.zeros(flattened_normals.shape)
            for ind in range(3):
                AIC_mtxb[:, :, ind] = (dresids['AIC'].T * flattened_normals[:, ind]).T
                flattened_normalsb[:, ind] += np.sum(self.AIC_mtx[:, :, ind].real * dresids['AIC'], axis=1).T

            # Actually assemble the AIC matrix
            _assemble_AIC_mtx_b(AIC_mtxb, params, dparams, dunknowns, dresids, self.surfaces)

            # Obtain the freestream velocity direction and magnitude by taking
            # alpha into account
            alpha = params['alpha'] * np.pi / 180.
            cosa = np.cos(alpha)
            sina = np.sin(alpha)
            arr = np.array([cosa, 0., sina])
            v_inf = params['v'] * arr

            fn = flattened_normals
            fnb = np.zeros(fn.shape)
            rhsb = dresids['rhs']

            v_infb = 0.
            for ind in reversed(range(self.tot_panels)):
                fnb[ind, :] -= v_inf * rhsb[ind]
                v_infb -= fn[ind, :] * rhsb[ind]

            dparams['v'] += sum(arr * v_infb)
            arrb = params['v'] * v_infb
            alphab = np.cos(alpha) * arrb[2]
            alphab -= np.sin(alpha) * arrb[0]
            alphab *= np.pi / 180.

            dparams['alpha'] += alphab

            i = 0
            for surface in self.surfaces:
                name = surface['name']
                nx = surface['num_x']
                ny = surface['num_y']
                num_panels = (nx - 1) * (ny - 1)
                dparams[name+'normals'] += flattened_normalsb[i:i+num_panels, :].reshape(nx-1, ny-1, 3, order='F')
                dparams[name+'normals'] += fnb[i:i+num_panels, :].reshape(nx-1, ny-1, 3, order='F')
                i += num_panels

class AeroCirculations(Component):
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

    def __init__(self, size):
        super(AeroCirculations, self).__init__()

        self.add_param('AIC', val=np.zeros((size, size), dtype=data_type))
        self.add_param('rhs', val=np.zeros((size), dtype=data_type))
        self.add_state('circulations', val=np.zeros((size), dtype=data_type))

        self.size = size

        # cache
        self.lup = None
        self.rhs_cache = None

    def solve_nonlinear(self, params, unknowns, resids):
        # lu factorization for use with solve_linear
        self.lup = lu_factor(params['AIC'])

        unknowns['circulations'] = lu_solve(self.lup, params['rhs'])

        resids['circulations'] = params['AIC'].dot(unknowns['circulations']) - params['rhs']

    def apply_nonlinear(self, params, unknowns, resids):
        """Evaluating residual for given state."""

        resids['circulations'] = params['AIC'].dot(unknowns['circulations']) - params['rhs']

    def apply_linear(self, params, unknowns, dparams, dunknowns, dresids, mode):
        """ Apply the derivative of state variable with respect to
        everything."""

        if mode == 'fwd':

            if 'circulations' in dunknowns:
                dresids['circulations'] += params['AIC'].dot(dunknowns['circulations'])
            if 'AIC' in dparams:
                dresids['circulations'] += dparams['AIC'].dot(unknowns['circulations'])
            if 'rhs' in dparams:
                dresids['circulations'] -= dparams['rhs']

        elif mode == 'rev':

            if 'circulations' in dunknowns:
                dunknowns['circulations'] += params['AIC'].T.dot(dresids['circulations'])
            if 'AIC' in dparams:
                dparams['AIC'] += np.outer(unknowns['circulations'], dresids['circulations']).T
            if 'rhs' in dparams:
                dparams['rhs'] -= dresids['circulations']

    def solve_linear(self, dumat, drmat, vois, mode=None):
        """ LU backsubstitution to solve the derivatives of the linear system."""

        if mode == 'fwd':
            sol_vec, rhs_vec = self.dumat, self.drmat
            t=0
        else:
            sol_vec, rhs_vec = self.drmat, self.dumat
            t=1

        if self.rhs_cache is None:
            self.rhs_cache = np.zeros((self.size, ))
        rhs = self.rhs_cache

        for voi in vois:
            rhs[:] = rhs_vec[voi]['circulations']

            sol = lu_solve(self.lup, rhs, trans=t)

            sol_vec[voi]['circulations'] = sol[:]


class VLMForces(Component):
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

    def __init__(self, surfaces):
        super(VLMForces, self).__init__()

        tot_panels = 0
        for surface in surfaces:
            name = surface['name']
            ny = surface['num_y']
            nx = surface['num_x']
            tot_panels += (nx - 1) * (ny - 1)

            self.add_param(name+'def_mesh', val=np.zeros((nx, ny, 3), dtype=data_type))
            self.add_param(name+'b_pts', val=np.zeros((nx-1, ny, 3), dtype=data_type))
            self.add_output(name+'sec_forces', val=np.zeros((nx-1, ny-1, 3), dtype=data_type))

        self.tot_panels = tot_panels

        self.add_param('circulations', val=np.zeros((tot_panels)))
        self.add_param('alpha', val=3.)
        self.add_param('v', val=10.)
        self.add_param('rho', val=3.)
        self.surfaces = surfaces

        self.mtx = np.zeros((tot_panels, tot_panels, 3), dtype=data_type)
        self.v = np.zeros((tot_panels, 3), dtype=data_type)

        if not fortran_flag:
            self.deriv_options['type'] = 'cs'
            self.deriv_options['form'] = 'central'

    def solve_nonlinear(self, params, unknowns, resids):
        circ = params['circulations']
        alpha = params['alpha'] * np.pi / 180.
        rho = params['rho']
        cosa = np.cos(alpha)
        sina = np.sin(alpha)

        # Assemble a different matrix here than the AIC_mtx from above; Note
        # that the collocation points used here are the midpoints of each
        # bound vortex filament, not the collocation points from above
        _assemble_AIC_mtx(self.mtx, params, self.surfaces, skip=True)

        # Compute the induced velocities at the midpoints of the
        # bound vortex filaments
        for ind in range(3):
            self.v[:, ind] = self.mtx[:, :, ind].dot(circ)

        # Add the freestream velocity to the induced velocity so that
        # self.v is the total velocity seen at the point
        self.v[:, 0] += cosa * params['v']
        self.v[:, 2] += sina * params['v']

        i = 0
        for surface in self.surfaces:
            name = surface['name']
            nx = surface['num_x']
            ny = surface['num_y']

            num_panels = (nx - 1) * (ny - 1)

            b_pts = params[name+'b_pts']

            if fortran_flag:
                sec_forces = OAS_API.oas_api.forcecalc(self.v[i:i+num_panels, :], circ[i:i+num_panels], rho, b_pts)
            else:

                bound = b_pts[:, 1:, :] - b_pts[:, :-1, :]

                # Cross the obtained velocities with the bound vortex filament
                # vectors
                cross = np.cross(self.v[i:i+num_panels],
                                    bound.reshape(-1, bound.shape[-1], order='F'))

                sec_forces = np.zeros(((nx-1)*(ny-1), 3), dtype=data_type)
                # Compute the sectional forces acting on each panel
                for ind in range(3):
                    sec_forces[:, ind] = \
                        (params['rho'] * circ[i:i+num_panels] * cross[:, ind])

            unknowns[name+'sec_forces'] = sec_forces.reshape((nx-1, ny-1, 3), order='F')

            i += num_panels

    def apply_linear(self, params, unknowns, dparams, dunknowns, dresids, mode):

        if mode == 'fwd':

            circ = params['circulations']
            alpha = params['alpha'] * np.pi / 180.
            alphad = dparams['alpha'] * np.pi / 180.
            cosa = np.cos(alpha)
            sina = np.sin(alpha)
            cosad = -sina * alphad
            sinad = cosa * alphad
            rho = params['rho']

            mtxd = np.zeros(self.mtx.shape)

            # Actually assemble the AIC matrix
            _assemble_AIC_mtx_d(mtxd, params, dparams, dunknowns, dresids, self.surfaces, skip=True)

            vd = np.zeros(self.v.shape)

            # Compute the induced velocities at the midpoints of the
            # bound vortex filaments
            for ind in range(3):
                vd[:, ind] += mtxd[:, :, ind].dot(circ)
                vd[:, ind] += self.mtx[:, :, ind].real.dot(dparams['circulations'])

            # Add the freestream velocity to the induced velocity so that
            # self.v is the total velocity seen at the point
            vd[:, 0] += cosa * dparams['v']
            vd[:, 2] += sina * dparams['v']
            vd[:, 0] += cosad * params['v']
            vd[:, 2] += sinad * params['v']

            i = 0
            rho = params['rho'].real
            for surface in self.surfaces:
                name = surface['name']
                nx = surface['num_x']
                ny = surface['num_y']

                num_panels = (nx - 1) * (ny - 1)

                b_pts = params[name+'b_pts']

                sec_forces = unknowns[name+'sec_forces'].real

                sec_forces, sec_forcesd = OAS_API.oas_api.forcecalc_d(self.v[i:i+num_panels, :], vd[i:i+num_panels], circ[i:i+num_panels], dparams['circulations'][i:i+num_panels], rho, dparams['rho'], b_pts, dparams[name+'b_pts'])

                dresids[name+'sec_forces'] += sec_forcesd.reshape((nx-1, ny-1, 3), order='F')

                i += num_panels

        if mode == 'rev':

            circ = params['circulations']
            alpha = params['alpha'] * np.pi / 180.
            cosa = np.cos(alpha)
            sina = np.sin(alpha)

            i = 0
            rho = params['rho'].real
            vb = np.zeros(self.v.shape)
            for surface in self.surfaces:
                name = surface['name']
                nx = surface['num_x']
                ny = surface['num_y']
                num_panels = (nx - 1) * (ny - 1)

                b_pts = params[name+'b_pts']
                sec_forcesb = dresids[name+'sec_forces'].reshape((num_panels, 3), order='F')

                sec_forces = unknowns[name+'sec_forces'].real

                v_b, circb, rhob, bptsb, _ = OAS_API.oas_api.forcecalc_b(self.v[i:i+num_panels, :], circ[i:i+num_panels], rho, b_pts, sec_forcesb)

                dparams['circulations'][i:i+num_panels] += circb
                vb[i:i+num_panels] = v_b
                dparams['rho'] += rhob
                dparams[name+'b_pts'] += bptsb

                i += num_panels

            sinab = params['v'] * np.sum(vb[:, 2])
            dparams['v'] += cosa * np.sum(vb[:, 0]) + sina * np.sum(vb[:, 2])
            cosab = params['v'] * np.sum(vb[:, 0])
            ab = np.cos(alpha) * sinab - np.sin(alpha) * cosab
            dparams['alpha'] += np.pi * ab / 180.

            mtxb = np.zeros(self.mtx.shape)
            circb = np.zeros(circ.shape)
            for i in range(3):
                for j in range(self.tot_panels):
                    mtxb[j, :, i] += circ * vb[j, i]
                    circb += self.mtx[j, :, i].real * vb[j, i]

            dparams['circulations'] += circb

            _assemble_AIC_mtx_b(mtxb, params, dparams, dunknowns, dresids, self.surfaces, skip=True)


class VLMLiftDrag(Component):
    """
    Calculate total lift and drag in force units based on section forces.

    Parameters
    ----------
    sec_forces[nx-1, ny-1, 3] : numpy array
        Flattened array containing the sectional forces acting on each panel.
        Stored in Fortran order (only relevant with more than one chordwise
        panel).
    alpha : float
        Angle of attack in degrees.

    Returns
    -------
    L : float
        Total induced lift force for the lifting surface.
    D : float
        Total induced drag force for the lifting surface.

    """

    def __init__(self, surface):
        super(VLMLiftDrag, self).__init__()

        self.surface = surface
        ny = surface['num_y']
        nx = surface['num_x']
        self.num_panels = (nx -1) * (ny - 1)

        self.add_param('sec_forces', val=np.zeros((nx - 1, ny - 1, 3)))
        self.add_param('alpha', val=3.)
        self.add_output('L', val=0.)
        self.add_output('D', val=0.)

    def solve_nonlinear(self, params, unknowns, resids):
        alpha = params['alpha'] * np.pi / 180.
        forces = params['sec_forces'].reshape(-1, 3)
        cosa = np.cos(alpha)
        sina = np.sin(alpha)

        # Compute the induced lift force on each lifting surface
        unknowns['L'] = np.sum(-forces[:, 0] * sina + forces[:, 2] * cosa)

        # Compute the induced drag force on each lifting surface
        unknowns['D'] = np.sum( forces[:, 0] * cosa + forces[:, 2] * sina)

        if self.surface['symmetry']:
            unknowns['D'] *= 2
            unknowns['L'] *= 2

    def linearize(self, params, unknowns, resids):
        """ Jacobian for lift and drag."""

        jac = self.alloc_jacobian()

        # Analytic derivatives for sec_forces
        alpha = params['alpha'] * np.pi / 180.
        cosa = np.cos(alpha)
        sina = np.sin(alpha)

        forces = params['sec_forces']

        if self.surface['symmetry']:
            symmetry_factor = 2.
        else:
            symmetry_factor = 1.

        tmp = np.array([-sina, 0, cosa])
        jac['L', 'sec_forces'] = \
            np.atleast_2d(np.tile(tmp, self.num_panels)) * symmetry_factor
        tmp = np.array([cosa, 0, sina])
        jac['D', 'sec_forces'] = \
            np.atleast_2d(np.tile(tmp, self.num_panels)) * symmetry_factor

        p180 = np.pi / 180.
        jac['L', 'alpha'] = p180 * symmetry_factor * \
            np.sum(-forces[:, :, 0] * cosa - forces[:, :, 2] * sina)
        jac['D', 'alpha'] = p180 * symmetry_factor * \
            np.sum(-forces[:, :, 0] * sina + forces[:, :, 2] * cosa)

        return jac

class ViscousDrag(Component):
    """
    Compute the skin friction drag if the with_viscous option is True.

    Parameters
    ----------
    re : float
        Dimensionalized (1/length) Reynolds number. This is used to compute the
        local Reynolds number based on the local chord length.
    M : float
        Mach number.
    S_ref : float
        The reference area of the lifting surface.
    sweep : float
        The angle (in degrees) of the wing sweep. This is used in the form
        factor calculation.
    widths[ny-1] : numpy array
        The spanwise width of each panel.
    lengths[ny] : numpy array
        The sum of the lengths of each line segment along a chord section.

    Returns
    -------
    CDv : float
        Viscous drag coefficient for the lifting surface computed using flat
        plate skin friction coefficient and a form factor to account for wing
        shape.
    """

    def __init__(self, surface, with_viscous):
        super(ViscousDrag, self).__init__()

        # Percentage of chord with laminar flow
        self.k_lam = surface['k_lam']

        # Thickness over chord for the airfoil
        self.t_over_c = surface['t_over_c']
        self.c_max_t = surface['c_max_t']

        self.ny = surface['num_y']

        self.add_param('re', val=5.e6)
        self.add_param('M', val=.84)
        self.add_param('S_ref', val=0.)
        self.add_param('cos_sweep', val=np.zeros((self.ny-1)))
        self.add_param('widths', val=np.zeros((self.ny-1)))
        self.add_param('lengths', val=np.zeros((self.ny)))
        self.add_output('CDv', val=0.)
        self.with_viscous = with_viscous

    def solve_nonlinear(self, params, unknowns, resids):
        if self.with_viscous:
            re = params['re']
            M = params['M']
            S_ref = params['S_ref']
            widths = params['widths']
            lengths = params['lengths']
            cos_sweep = params['cos_sweep'] / widths

            self.d_over_q = np.zeros((self.ny - 1))

            # Take panel chord length to be average of its edge lengths
            chords = (lengths[1:] + lengths[:-1]) / 2.
            Re_c = re * chords

            cdturb_total = 0.455 / (np.log10(Re_c))**2.58 / \
                (1.0 + 0.144*M**2)**0.65
            cdlam_tr = 1.328 / np.sqrt(Re_c * self.k_lam)

            # Use eq. 12.27 of Raymer for turbulent Cf
            if self.k_lam == 0:
                cdlam_tr = 0.
                cd = cdturb_total

            elif self.k_lam < 1.0:
                cdturb_tr = 0.455 / (np.log10(Re_c*self.k_lam))**2.58 / \
                    (1.0 + 0.144*M**2)**0.65

            else:
                cdturb_total = 0.

            cd = (cdlam_tr - cdturb_tr)*self.k_lam + cdturb_total

            # Multiply by section width to get total normalized drag for section
            # d_over_q = d / 0.5 / rho / v**2
            self.d_over_q = 2 * cd * chords

            # Calculate form factor
            self.k_FF = 1.34 * M**0.18 * \
                (1.0 + 0.6*self.t_over_c/self.c_max_t + 100*self.t_over_c**4)
            FF = self.k_FF * cos_sweep**0.28

            # Sum individual panel drags to get total drag
            self.D_over_q = np.sum(self.d_over_q * widths * FF)

            unknowns['CDv'] = self.D_over_q / S_ref
        else:
            unknowns['CDv'] = 0.0

    def linearize(self, params, unknowns, resids):
        """ Jacobian for viscous drag."""

        jac = self.alloc_jacobian()
        jac['CDv', 'lengths'] = np.zeros_like(jac['CDv', 'lengths'])
        re = params['re']

        if self.with_viscous:
            p180 = np.pi / 180.
            M = params['M']
            S_ref = params['S_ref']
            widths = params['widths']
            lengths = params['lengths']
            cos_sweep = params['cos_sweep'] / widths

            B = (1. + 0.144*M**2)**0.65

            FF = self.k_FF * cos_sweep**0.28

            chords = (lengths[1:] + lengths[:-1]) / 2.
            Re_c = re * chords

            cdl_Re = 0.0
            cdt_Re = 0.0
            cdT_Re = 0.0

            if self.k_lam == 0:
                cdT_Re = 0.455/(np.log10(Re_c))**3.58/B * \
                            -2.58 / np.log(10) / Re_c
            elif self.k_lam < 1.0:

                cdl_Re = 1.328 / (Re_c*self.k_lam)**1.5 * -0.5 * self.k_lam
                cdt_Re = 0.455/(np.log10(Re_c*self.k_lam))**3.58/B * \
                            -2.58 / np.log(10) / Re_c
                cdT_Re = 0.455/(np.log10(Re_c))**3.58/B * \
                            -2.58 / np.log(10) / Re_c
            else:
                cdl_Re = 1.328 / (Re_c*self.k_lam)**1.5 * -0.5 * self.k_lam

            cd_Re = (cdl_Re - cdt_Re)*self.k_lam + cdT_Re

            CDv_lengths = 2 * widths * FF / S_ref * \
                (self.d_over_q / 4 / chords + chords * cd_Re * re / 2.)
            jac['CDv', 'lengths'][0, 1:] += CDv_lengths
            jac['CDv', 'lengths'][0, :-1] += CDv_lengths
            jac['CDv', 'widths'][0, :] = self.d_over_q * FF / S_ref * 0.72
            jac['CDv', 'S_ref'] = - self.D_over_q / S_ref**2
            jac['CDv', 'cos_sweep'][0, :] = 0.28 * self.k_FF * self.d_over_q / S_ref / cos_sweep**0.72

        return jac

class VLMCoeffs(Component):
    """ Compute lift and drag coefficients.

    Parameters
    ----------
    S_ref : float
        The reference areas of the lifting surface.
    L : float
        Total lift for the lifting surface.
    D : float
        Total drag for the lifting surface.
    v : float
        Freestream air velocity in m/s.
    rho : float
        Air density in kg/m^3.

    Returns
    -------
    CL1 : float
        Induced coefficient of lift (CL) for the lifting surface.
    CDi : float
        Induced coefficient of drag (CD) for the lifting surface.
    """

    def __init__(self, surface):
        super(VLMCoeffs, self).__init__()

        self.surface = surface

        self.add_param('S_ref', val=0.)
        self.add_param('L', val=0.)
        self.add_param('D', val=0.)
        self.add_param('v', val=0.)
        self.add_param('rho', val=0.)
        self.add_output('CL1', val=0.)
        self.add_output('CDi', val=0.)

    def solve_nonlinear(self, params, unknowns, resids):
        S_ref = params['S_ref']
        rho = params['rho']
        v = params['v']
        L = params['L']
        D = params['D']

        if self.surface['symmetry']:
            S_ref *= 2

        unknowns['CL1'] = L / (0.5 * rho * v**2 * S_ref)
        unknowns['CDi'] = D / (0.5 * rho * v**2 * S_ref)

    def linearize(self, params, unknowns, resids):
        S_ref = params['S_ref']
        rho = params['rho']
        v = params['v']
        L = params['L']
        D = params['D']

        if self.surface['symmetry']:
            S_ref *= 2

        jac = self.alloc_jacobian()

        jac['CL1', 'L'] = 1. / (0.5 * rho * v**2 * S_ref)
        jac['CDi', 'D'] = 1. / (0.5 * rho * v**2 * S_ref)

        jac['CL1', 'v'] = -2. * L / (0.5 * rho * v**3 * S_ref)
        jac['CDi', 'v'] = -2. * D / (0.5 * rho * v**3 * S_ref)

        jac['CL1', 'rho'] = -L / (0.5 * rho**2 * v**2 * S_ref)
        jac['CDi', 'rho'] = -D / (0.5 * rho**2 * v**2 * S_ref)

        if self.surface['symmetry']:
            jac['CL1', 'S_ref'] = -L / (.25 * rho * v**2 * S_ref**2)
            jac['CDi', 'S_ref'] = -D / (.25 * rho * v**2 * S_ref**2)
        else:
            jac['CL1', 'S_ref'] = -L / (0.5 * rho * v**2 * S_ref**2)
            jac['CDi', 'S_ref'] = -D / (0.5 * rho * v**2 * S_ref**2)

        jac['CL1', 'D'] = 0.
        jac['CDi', 'L'] = 0.

        return jac

class TotalLift(Component):
    """ Calculate total lift in force units.

    Parameters
    ----------
    CL1 : float
        Induced coefficient of lift (CL) for the lifting surface.

    Returns
    -------
    CL : float
        Total coefficient of lift (CL) for the lifting surface.
    CL_wing : float
        CL of the main wing, used for CL constrained optimization.

    """

    def __init__(self, surface):
        super(TotalLift, self).__init__()

        self.add_param('CL1', val=0.)
        self.add_output('CL', val=0.)
        self.CL0 = surface['CL0']

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['CL'] = params['CL1'] + self.CL0

    def linearize(self, params, unknowns, resids):
        jac = self.alloc_jacobian()
        jac['CL', 'CL1'] = 1.
        return jac


class TotalDrag(Component):
    """ Calculate total drag in force units.

    Parameters
    ----------
    CDi : float
        Induced coefficient of drag (CD) for the lifting surface.
    CDv : float
        Calculated viscous drag for the lifting surface..

    Returns
    -------
    CD : float
        Total coefficient of drag (CD) for the lifting surface.

    """

    def __init__(self, surface):
        super(TotalDrag, self).__init__()

        self.add_param('CDi', val=0.)
        self.add_param('CDv', val=0.)
        self.add_output('CD', val=0.)
        self.CD0 = surface['CD0']

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['CD'] = params['CDi'] + params['CDv'] + self.CD0

    def linearize(self, params, unknowns, resids):
        jac = self.alloc_jacobian()
        jac['CD', 'CDi'] = 1.
        jac['CD', 'CDv'] = 1.
        return jac


class VLMStates(Group):
    """ Group that contains the aerodynamic states. """

    def __init__(self, surfaces):
        super(VLMStates, self).__init__()

        tot_panels = 0
        for surface in surfaces:
            ny = surface['num_y']
            nx = surface['num_x']
            tot_panels += (nx - 1) * (ny - 1)

        self.add('assembly',
                 AssembleAIC(surfaces),
                 promotes=['*'])
        self.add('circulations',
                 AeroCirculations(tot_panels),
                 promotes=['*'])
        self.add('forces',
                 VLMForces(surfaces),
                 promotes=['*'])


class VLMFunctionals(Group):
    """ Group that contains the aerodynamic functionals used to evaluate
    performance. """

    def __init__(self, surface, prob_dict):
        super(VLMFunctionals, self).__init__()

        with_viscous = prob_dict['with_viscous']

        self.add('liftdrag',
                 VLMLiftDrag(surface),
                 promotes=['*'])
        self.add('coeffs',
                 VLMCoeffs(surface),
                 promotes=['*'])
        self.add('CL',
                 TotalLift(surface),
                 promotes=['*'])
        self.add('CD',
                 TotalDrag(surface),
                 promotes=['*'])
        self.add('viscousdrag',
                 ViscousDrag(surface, with_viscous),
                 promotes=['*'])
