"""
Define the aerodynamic analysis component using a vortex lattice method.

We input a nodal mesh and properties of the airflow to calculate the
circulations of the horseshoe vortices. We then compute the forces, lift,
and drag acting on the lifting surfaces.

.. todo:: fix depiction of wing
Depiction of wing:

    y -->
 x  -----------------  Leading edge
 |  | + | + | + | + |
 v  |   |   |   |   |
    | o | o | o | o |
    -----------------  Trailing edge

The '+' symbols represent the bound points on the 1/4 chord line, which are
used to compute drag.
The 'o' symbols represent the collocation points on the 3/4 chord line, where
the flow tangency condition is inforced.

For this depiction, num_x = 2 and num_y = 5.
"""

from __future__ import division
import numpy

from openmdao.api import Component, Group
from scipy.linalg import lu_factor, lu_solve
try:
    import lib
    fortran_flag = True
except:
    fortran_flag = False


def view_mat(mat):
    """ Helper function used to visually examine matrices. """
    import matplotlib.pyplot as plt
    if len(mat.shape) > 2:
        mat = numpy.sum(mat, axis=2)
    print "Cond #:", numpy.linalg.cond(mat)
    im = plt.imshow(mat.real, interpolation='none')
    plt.colorbar(im, orientation='horizontal')
    plt.show()


def norm(vec):
    """ Finds the 2-norm of a vector. """
    return numpy.sqrt(numpy.sum(vec**2))


def _calc_vorticity(A, B, P):
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

    return (r1_mag + r2_mag) * numpy.cross(r1, r2) / \
           (r1_mag * r2_mag * (r1_mag * r2_mag + r1.dot(r2)))


def _assemble_AIC_mtx(mtx, flat_mesh, aero_ind, points,
                      b_pts, alpha, skip=False):
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

    mtx[:, :, :] = 0.0
    cosa = numpy.cos(alpha * numpy.pi / 180.)
    sina = numpy.sin(alpha * numpy.pi / 180.)
    u = numpy.array([cosa, 0, sina])

    # Loop over the lifting surfaces to compute their influence on the flow
    # velocity at the collocation points
    for i_surf, row in enumerate(aero_ind):
        # Variable names with a trailing underscore correspond to the lifting
        # surface being examined, not the collocation point
        nx_, ny_, n_, n_bpts_, n_panels_, i_, i_bpts_, i_panels_ = row.copy()
        n = nx_ * ny_

        # Obtain the lifting surface mesh in the form expected by the solver,
        # with shape [nx_, ny_, 3]
        mesh = flat_mesh[i_:i_+n_, :].reshape(nx_, ny_, 3)
        bpts = b_pts[i_bpts_:i_bpts_+n_bpts_].reshape(nx_-1, ny_, 3)

        for i_points, row in enumerate(aero_ind):
            # These variables correspond to the collocation points
            nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = row

            # Obtain the collocation points used to compute the AIC
            pts = points[i_panels:i_panels+n_panels].reshape(nx-1, ny-1, 3)

            small_mat = numpy.zeros((n_panels, n_panels_, 3)).astype("complex")

            # Dense fortran assembly for the AIC matrix
            if fortran_flag:
                small_mat[:, :, :] = lib.assembleaeromtx(ny, nx, ny_, nx_,
                                                         alpha, pts, bpts,
                                                         mesh, skip)
            # Python matrix assembly
            else:
                # Spanwise loop through horseshoe elements
                for el_j in xrange(ny_ - 1):
                    el_loc_j = el_j * (nx_ - 1)
                    C_te = mesh[-1, el_j + 1, :]
                    D_te = mesh[-1, el_j + 0, :]

                    # Spanwise loop through control points
                    for cp_j in xrange(ny - 1):
                        cp_loc_j = cp_j * (nx - 1)

                        # Chordwise loop through control points
                        for cp_i in xrange(nx - 1):
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

                            # AIC contribution from trailing vortex filaments
                            # coming off the trailing edge
                            trailing = t1 - t3
                            edges = 0

                            # Chordwise loop through horseshoe elements in
                            # reversed order, starting with the panel closest
                            # to the leading edge. This is done to sum the
                            # AIC contributions from the side vortex filaments
                            # as we loop through the elements
                            for el_i in reversed(xrange(nx_ - 1)):
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

                                # If skip, do not include the contributions
                                # from the panel's bound vortex filament, as
                                # this causes a singularity when we're taking
                                # the influence of a panel on its own
                                # collocation point. This true for the drag
                                # computation and false for circulation
                                # computation, due to the different collocation
                                # points.
                                if skip and el_loc == cp_loc:
                                    small_mat[cp_loc, el_loc, :] = \
                                        trailing + edges
                                else:
                                    bound = _calc_vorticity(A, B, P)
                                    small_mat[cp_loc, el_loc, :] = \
                                        trailing + edges + bound

            # Populate the full-size matrix with these surface-surface AICs
            mtx[i_panels:i_panels+n_panels,
                i_panels_:i_panels_+n_panels_, :] = small_mat

    mtx /= 4 * numpy.pi


class VLMGeometry(Component):
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

    def __init__(self, aero_ind):
        super(VLMGeometry, self).__init__()

        num_surf = aero_ind.shape[0]
        tot_n = numpy.sum(aero_ind[:, 2])
        tot_bpts = numpy.sum(aero_ind[:, 3])
        tot_panels = numpy.sum(aero_ind[:, 4])
        self.aero_ind = aero_ind

        self.add_param('def_mesh', val=numpy.zeros((tot_n, 3),
                       dtype="complex"))
        self.add_output('b_pts', val=numpy.zeros((tot_bpts, 3),
                        dtype="complex"))
        self.add_output('mid_b', val=numpy.zeros((tot_panels, 3),
                        dtype="complex"))
        self.add_output('c_pts', val=numpy.zeros((tot_panels, 3)))
        self.add_output('widths', val=numpy.zeros((tot_panels)))
        self.add_output('normals', val=numpy.zeros((tot_panels, 3)))
        self.add_output('S_ref', val=numpy.zeros((num_surf)))

    def _get_lengths(self, A, B, axis):
        return numpy.sqrt(numpy.sum((B - A)**2, axis=axis))

    def solve_nonlinear(self, params, unknowns, resids):
        # Loop through each lifting surface
        for i_surf, row in enumerate(self.aero_ind):
            nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = row

            # Reshape the mesh into the expected shape for each surface
            mesh = params['def_mesh'][i:i+n, :].reshape(nx, ny, 3)

            # Compute the bound points at 1/4 chord
            b_pts = mesh[:-1, :, :] * .75 + mesh[1:, :, :] * .25

            # Find the midpoints of the bound points, used in drag computations
            mid_b = (b_pts[:, 1:, :] + b_pts[:, :-1, :]) / 2

            # Compute the collocation points at the midpoints of each
            # panel's 3/4 chord line
            c_pts = 0.5 * 0.25 * mesh[:-1, :-1, :] + \
                    0.5 * 0.75 * mesh[1:, :-1, :] + \
                    0.5 * 0.25 * mesh[:-1,  1:, :] + \
                    0.5 * 0.75 * mesh[1:,  1:, :]

            # Compute the widths of each panel
            widths = self._get_lengths(b_pts[:, 1:, :], b_pts[:, :-1, :], 2)

            # Compute the normal of each panel by taking the cross-product of
            # its diagonals. Note that this could be a nonplanar surface
            normals = numpy.cross(
                mesh[:-1,  1:, :] - mesh[1:, :-1, :],
                mesh[:-1, :-1, :] - mesh[1:,  1:, :],
                axis=2)

            norms = numpy.sqrt(numpy.sum(normals**2, axis=2))

            for j in xrange(3):
                normals[:, :, j] /= norms

            # Flatten and store each array
            unknowns['b_pts'][i_bpts:i_bpts+n_bpts, :] = \
                b_pts.reshape(-1, b_pts.shape[-1])
            unknowns['mid_b'][i_panels:i_panels+n_panels, :] = \
                mid_b.reshape(-1, mid_b.shape[-1])
            unknowns['c_pts'][i_panels:i_panels+n_panels, :] = \
                c_pts.reshape(-1, c_pts.shape[-1])
            unknowns['widths'][i_panels:i_panels+n_panels] = \
                widths.flatten()
            unknowns['normals'][i_panels:i_panels+n_panels, :] = \
                normals.reshape(-1, normals.shape[-1], order='F')
            unknowns['S_ref'][i_surf] = 0.5 * numpy.sum(norms)

    def linearize(self, params, unknowns, resids):
        """ Jacobian for geometry."""

        jac = self.alloc_jacobian()

        fd_jac = self.complex_step_jacobian(params, unknowns, resids,
                                            fd_params=['def_mesh'],
                                            fd_unknowns=['widths', 'normals',
                                                         'S_ref'],
                                            fd_states=[])
        jac.update(fd_jac)

        # Analytic derivatives for some of the parameters
        i_ny = 0.
        for i_surf, row in enumerate(self.aero_ind):
            nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = row

            for iz, v in zip(((i_bpts+i_ny)*3, (i_bpts+i_ny+ny)*3),
                             (.75, .25)):
                numpy.fill_diagonal(jac['b_pts', 'def_mesh']
                    [i_bpts*3:(n_bpts+i_bpts)*3, iz:], v)

            for iz, v in zip((i*3, (i+1)*3, (i+ny)*3, (ny+i+1)*3),
                             (.125, .125, .375, .375)):
                for ix in range(nx-1):
                    numpy.fill_diagonal(jac['c_pts', 'def_mesh']
                        [(i_panels+ix*(ny-1))*3:(i_panels+(ix+1)*(ny-1))*3,
                     iz+ix*ny*3:], v)

            for iz, v in zip((i*3, (i+1)*3, (i+ny)*3, (ny+i+1)*3),
                             (.375, .375, .125, .125)):
                for ix in range(nx-1):
                    numpy.fill_diagonal(jac['mid_b', 'def_mesh']
                        [(i_panels+ix*(ny-1))*3:(i_panels+(ix+1)*(ny-1))*3,
                     iz+ix*ny*3:], v)

            i_ny += ny

        return jac


class VLMCirculations(Component):
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

    def __init__(self, aero_ind):
        super(VLMCirculations, self).__init__()

        tot_n = numpy.sum(aero_ind[:, 2])
        tot_bpts = numpy.sum(aero_ind[:, 3])
        tot_panels = numpy.sum(aero_ind[:, 4])
        self.aero_ind = aero_ind

        self.add_param('def_mesh', val=numpy.zeros((tot_n, 3),
                       dtype="complex"))
        self.add_param('b_pts', val=numpy.zeros((tot_bpts, 3),
                       dtype="complex"))
        self.add_param('c_pts', val=numpy.zeros((tot_panels, 3),
                       dtype="complex"))
        self.add_param('normals', val=numpy.zeros((tot_panels, 3)))
        self.add_param('v', val=10.)
        self.add_param('alpha', val=3.)
        self.add_state('circulations', val=numpy.zeros((tot_panels),
                       dtype="complex"))

        self.deriv_options['linearize'] = True  # only for circulations

        self.AIC_mtx = numpy.zeros((tot_panels, tot_panels, 3),
                                   dtype="complex")
        self.mtx = numpy.zeros((tot_panels, tot_panels), dtype="complex")
        self.rhs = numpy.zeros((tot_panels), dtype="complex")

    def _assemble_system(self, params):
        _assemble_AIC_mtx(self.AIC_mtx, params['def_mesh'], self.aero_ind,
                        params['c_pts'], params['b_pts'], params['alpha'])

        # Construct a matrix that is the AIC_mtx dotted by the normals at each
        # collocation point. This is used to compute the circulations
        self.mtx[:, :] = 0.
        for ind in xrange(3):
            self.mtx[:, :] += (self.AIC_mtx[:, :, ind].T *
                params['normals'][:, ind]).T

        alpha = params['alpha'] * numpy.pi / 180.
        cosa = numpy.cos(alpha)
        sina = numpy.sin(alpha)
        v_inf = params['v'] * numpy.array([cosa, 0., sina], dtype="complex")

        # Populate the right-hand side of the linear system with the
        # expected velocities at each collocation point
        self.rhs[:] = -params['normals'].\
            reshape(-1, params['normals'].shape[-1], order='F').dot(v_inf)

    def solve_nonlinear(self, params, unknowns, resids):
        self._assemble_system(params)
        unknowns['circulations'] = numpy.linalg.solve(self.mtx, self.rhs)

    def apply_nonlinear(self, params, unknowns, resids):
        self._assemble_system(params)

        circ = unknowns['circulations']
        resids['circulations'] = self.mtx.dot(circ) - self.rhs

    def linearize(self, params, unknowns, resids):
        """ Jacobian for circulations."""
        self.lup = lu_factor(self.mtx.real)
        jac = self.alloc_jacobian()

        fd_jac = self.complex_step_jacobian(params, unknowns, resids,
                                         fd_params=['normals', 'alpha',
                                                    'def_mesh', 'b_pts',
                                                    'c_pts'],
                                         fd_states=[])
        jac.update(fd_jac)

        jac['circulations', 'circulations'] = self.mtx.real

        jac['circulations', 'v'][:, 0] = -self.rhs.real / params['v'].real

        return jac

    def solve_linear(self, dumat, drmat, vois, mode=None):

        if mode == 'fwd':
            sol_vec, rhs_vec = self.dumat, self.drmat
            t = 0
        else:
            sol_vec, rhs_vec = self.drmat, self.dumat
            t = 1

        for voi in vois:
            sol_vec[voi].vec[:] = lu_solve(self.lup, rhs_vec[voi].vec, trans=t)


class VLMForces(Component):
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
    circulations : array_like
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

    def __init__(self, aero_ind):
        super(VLMForces, self).__init__()

        tot_n = numpy.sum(aero_ind[:, 2])
        tot_bpts = numpy.sum(aero_ind[:, 3])
        tot_panels = numpy.sum(aero_ind[:, 4])
        self.aero_ind = aero_ind

        self.add_param('def_mesh', val=numpy.zeros((tot_n, 3)))
        self.add_param('b_pts', val=numpy.zeros((tot_bpts, 3)))
        self.add_param('mid_b', val=numpy.zeros((tot_panels, 3)))
        self.add_param('circulations', val=numpy.zeros((tot_panels)))
        self.add_param('alpha', val=3.)
        self.add_param('v', val=10.)
        self.add_param('rho', val=3.)
        self.add_output('sec_forces', val=numpy.zeros((tot_panels, 3)))

        self.mtx = numpy.zeros((tot_panels, tot_panels, 3), dtype="complex")

        self.v = numpy.zeros((tot_panels, 3), dtype="complex")

    def solve_nonlinear(self, params, unknowns, resids):
        circ = params['circulations']
        alpha = params['alpha'] * numpy.pi / 180.
        cosa = numpy.cos(alpha)
        sina = numpy.sin(alpha)

        # Assemble a different matrix here than the AIC_mtx from above; Notes
        # that the collocation points used here are the midpoints of each
        # bound vortex filament, not the collocation points from above
        _assemble_AIC_mtx(self.mtx, params['def_mesh'], self.aero_ind,
                          params['mid_b'], params['b_pts'], params['alpha'],
                          skip=True)

        for i_surf, row in enumerate(self.aero_ind):
            nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = row

            # Compute the velocities at the midpoints of the bound vortex
            # filaments
            for ind in xrange(3):
                self.v[:, ind] = self.mtx[:, :, ind].dot(circ)
            self.v[:, 0] += cosa * params['v']
            self.v[:, 2] += sina * params['v']

            b_pts = params['b_pts'][i_bpts:i_bpts+n_bpts, :].\
                reshape(nx-1, ny, 3)

            bound = b_pts[:, 1:, :] - b_pts[:, :-1, :]

            # Cross the obtained velocities with the bound vortex filament
            # vectors
            cross = numpy.cross(self.v[i_panels:i_panels+n_panels],
                                bound.reshape(-1, bound.shape[-1], order='F'))

            # Compute the sectional forces acting on each panel
            for ind in xrange(3):
                unknowns['sec_forces'][i_panels:i_panels+n_panels, ind] = \
                    (params['rho'] * circ[i_panels:i_panels+n_panels] *
                    cross[:, ind])

    def linearize(self, params, unknowns, resids):
        """ Jacobian for forces."""

        jac = self.alloc_jacobian()

        fd_jac = self.complex_step_jacobian(params, unknowns, resids,
                                         fd_params=['b_pts', 'alpha',
                                                    'circulations', 'v',
                                                    'mid_b', 'def_mesh'],
                                         fd_states=[])
        jac.update(fd_jac)

        rho = params['rho'].real
        sec_forces = unknowns['sec_forces'].real

        jac['sec_forces', 'rho'] = sec_forces.flatten() / rho

        return jac


class VLMLiftDrag(Component):
    """
    Calculate total lift and drag in force units based on section forces.

    If the simulation is given a non-zero Reynolds number, it is used to
    compute the skin friction drag. If Reynolds number == 0, then there is
    no skin friction drag dependence. Currently, the Reynolds number
    must be set by the user in the run script and used as in IndepVarComp.

    Parameters
    ----------
    sec_forces : array_like
        Flattened array containing the sectional forces acting on each panel.
        Stored in Fortran order (only relevant when more than one chordwise
        panel).
    alpha : float
        Angle of attack in degrees.
    Re : float
        Reynolds number.
    M : float
        Mach number.
    v : float
        Freestream air velocity in m/s.
    rho : float
        Air density in kg/m^3.
    S_ref : array_like
        The reference areas of each lifting surface.

    Returns
    -------
    L : array_like
        Total lift force for each lifting surface.
    D : array_like
        Total drag force for each lifting surface.

    """

    def __init__(self, aero_ind):
        super(VLMLiftDrag, self).__init__()

        num_surf = aero_ind.shape[0]
        tot_panels = numpy.sum(aero_ind[:, 4])
        self.aero_ind = aero_ind

        self.add_param('sec_forces', val=numpy.zeros((tot_panels, 3)))
        self.add_param('alpha', val=3.)
        self.add_param('Re', val=5.e6)
        self.add_param('M', val=.84)
        self.add_param('v', val=10.)
        self.add_param('rho', val=3.)
        self.add_param('S_ref', val=numpy.zeros((num_surf)))
        self.add_output('L', val=numpy.zeros((num_surf)))
        self.add_output('D', val=numpy.zeros((num_surf)))

    def solve_nonlinear(self, params, unknowns, resids):
        alpha = params['alpha'] * numpy.pi / 180.
        forces = params['sec_forces']
        cosa = numpy.cos(alpha)
        sina = numpy.sin(alpha)

        Re = params['Re']
        M = params['M']
        v = params['v']
        rho = params['rho']
        S_ref = params['S_ref']

        # Compute the skin friction coefficient
        # Use eq. 12.27 of Raymer for turbulent Cf
        Cf = 0.455 / (numpy.log10(Re)**2.58 * (1 + .144 * M**2)**.65)

        for i_surf, row in enumerate(self.aero_ind):
            nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = row

            # Compute the induced lift force on each lifting surface
            unknowns['L'][i_surf] = \
                numpy.sum(-forces[i_panels:i_panels+n_panels, 0] * sina +
                           forces[i_panels:i_panels+n_panels, 2] * cosa)


            # Compute the induced drag force on each lifting surface
            unknowns['D'][i_surf] = \
                numpy.sum( forces[i_panels:i_panels+n_panels, 0] * cosa +
                           forces[i_panels:i_panels+n_panels, 2] * sina)

            # Compute the drag contribution from skin friction
            D_f = Cf * rho * v**2 / 2. * S_ref[i_surf]
            unknowns['D'][i_surf] += D_f


    def linearize(self, params, unknowns, resids):
        """ Jacobian for lift and drag."""

        jac = self.alloc_jacobian()

        # Analytic derivatives for sec_forces
        alpha = params['alpha'] * numpy.pi / 180.
        cosa = numpy.cos(alpha)
        sina = numpy.sin(alpha)

        for i_surf, row in enumerate(self.aero_ind):
            nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = row

            forces = params['sec_forces'][i_panels:n_panels+i_panels].\
                reshape(nx-1, ny-1, 3)

            tmp = numpy.array([-sina, 0, cosa])
            panel_start = i_panels*3
            panel_end = i_panels*3+n_panels*3
            jac['L', 'sec_forces'][i_surf, panel_start:panel_end] = \
                numpy.atleast_2d(numpy.tile(tmp, n_panels))
            tmp = numpy.array([cosa, 0, sina])
            jac['D', 'sec_forces'][i_surf, panel_start:panel_end] = \
                numpy.atleast_2d(numpy.tile(tmp, n_panels))

            p180 = numpy.pi / 180.
            jac['L', 'alpha'][i_surf] = p180 * \
                numpy.sum(-forces[:, :, 0] * cosa - forces[:, :, 2] * sina)
            jac['D', 'alpha'][i_surf] = p180 * \
                numpy.sum(-forces[:, :, 0] * sina + forces[:, :, 2] * cosa)

        return jac


class VLMCoeffs(Component):
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


    def __init__(self, aero_ind):
        super(VLMCoeffs, self).__init__()

        num_surf = aero_ind.shape[0]
        self.aero_ind = aero_ind
        self.add_param('S_ref', val=numpy.zeros((num_surf)))
        self.add_param('L', val=numpy.zeros((num_surf)))
        self.add_param('D', val=numpy.zeros((num_surf)))
        self.add_param('v', val=0.)
        self.add_param('rho', val=0.)
        self.add_output('CL1', val=numpy.zeros((num_surf)))
        self.add_output('CDi', val=numpy.zeros((num_surf)))

    def solve_nonlinear(self, params, unknowns, resids):
        S_ref = params['S_ref']
        rho = params['rho']
        v = params['v']
        L = params['L']
        D = params['D']
        for i_surf, row in enumerate(self.aero_ind):
            nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = row

            unknowns['CL1'][i_surf] = L[i_surf] / \
                (0.5 * rho * v**2 * S_ref[i_surf])
            unknowns['CDi'][i_surf] = D[i_surf] / \
                (0.5 * rho * v**2 * S_ref[i_surf])

    def linearize(self, params, unknowns, resids):
        """ Jacobian for forces."""

        S_ref = params['S_ref']
        rho = params['rho']
        v = params['v']
        L = params['L']
        D = params['D']

        jac = self.alloc_jacobian()

        for i_surf, row in enumerate(self.aero_ind):
            nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = row

            jac['CL1', 'D'][i_surf] = 0.
            jac['CDi', 'L'][i_surf] = 0.

            tmp = 0.5*rho*v**2*S_ref
            jac['CL1', 'L'][i_surf, i_surf] = 1. / tmp[i_surf]
            jac['CDi', 'D'][i_surf, i_surf] = 1. / tmp[i_surf]

            tmp = -0.5*rho**2*v**2*S_ref
            jac['CL1', 'rho'][i_surf] = L[i_surf] / tmp[i_surf]
            jac['CDi', 'rho'][i_surf] = D[i_surf] / tmp[i_surf]

            tmp = -0.5*rho*v**2*S_ref**2
            jac['CL1', 'S_ref'][i_surf, i_surf] = L[i_surf] / tmp[i_surf]
            jac['CDi', 'S_ref'][i_surf, i_surf] = D[i_surf] / tmp[i_surf]

            tmp = -0.25*rho*v**3*S_ref
            jac['CL1', 'v'][i_surf] = L[i_surf] / tmp[i_surf]
            jac['CDi', 'v'][i_surf] = D[i_surf] / tmp[i_surf]

        return jac


class TotalLift(Component):
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

    def __init__(self, CL0, aero_ind):
        super(TotalLift, self).__init__()

        self.num_surf = aero_ind.shape[0]
        self.add_param('CL1', val=numpy.zeros((self.num_surf)))
        self.add_output('CL', val=numpy.zeros((self.num_surf)))
        self.add_output('CL_wing', val=0.)

        self.CL0 = CL0

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['CL'] = params['CL1'] + self.CL0
        unknowns['CL_wing'] = params['CL1'][0]

    def linearize(self, params, unknowns, resids):
        jac = self.alloc_jacobian()
        jac['CL', 'CL1'][:] = numpy.eye(self.num_surf)
        jac['CL_wing', 'CL1'][0][0] = 1
        return jac


class TotalDrag(Component):
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

    def __init__(self, CD0, aero_ind):
        super(TotalDrag, self).__init__()

        self.num_surf = aero_ind.shape[0]
        self.add_param('CDi', val=numpy.zeros((self.num_surf)))
        self.add_output('CD', val=numpy.zeros((self.num_surf)))
        self.add_output('CD_wing', val=0.)

        self.CD0 = CD0

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['CD'] = params['CDi'] + self.CD0
        unknowns['CD_wing'] = unknowns['CD'][0]

    def linearize(self, params, unknowns, resids):
        jac = self.alloc_jacobian()
        jac['CD', 'CDi'][:] = numpy.eye(self.num_surf)
        jac['CD_wing', 'CDi'][0][0] = 1
        return jac


class VLMStates(Group):
    """ Group that contains the aerodynamic states. """

    def __init__(self, aero_ind):
        super(VLMStates, self).__init__()

        self.add('wgeom',
                 VLMGeometry(aero_ind),
                 promotes=['*'])
        self.add('circ',
                 VLMCirculations(aero_ind),
                 promotes=['*'])
        self.add('forces',
                 VLMForces(aero_ind),
                 promotes=['*'])


class VLMFunctionals(Group):
    """ Group that contains the aerodynamic functionals used to evaluate
    performance. """

    def __init__(self, aero_ind, CL0, CD0):
        super(VLMFunctionals, self).__init__()

        self.add('liftdrag',
                 VLMLiftDrag(aero_ind),
                 promotes=['*'])
        self.add('coeffs',
                 VLMCoeffs(aero_ind),
                 promotes=['*'])
        self.add('CL',
                 TotalLift(CL0, aero_ind),
                 promotes=['*'])
        self.add('CD',
                 TotalDrag(CD0, aero_ind),
                 promotes=['*'])
