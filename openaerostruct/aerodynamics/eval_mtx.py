import numpy as np

import openmdao.api as om

from openaerostruct.utils.vector_algebra import add_ones_axis
from openaerostruct.utils.vector_algebra import compute_dot, compute_dot_deriv
from openaerostruct.utils.vector_algebra import compute_cross, compute_cross_deriv1, compute_cross_deriv2
from openaerostruct.utils.vector_algebra import compute_norm, compute_norm_deriv


tol = 1e-10

def _compute_finite_vortex(r1, r2):
    r1_norm = compute_norm(r1)
    r2_norm = compute_norm(r2)

    r1_x_r2 = compute_cross(r1, r2)
    r1_d_r2 = compute_dot(r1, r2)

    num = (1. / r1_norm + 1. / r2_norm) * r1_x_r2
    den = r1_norm * r2_norm + r1_d_r2

    result = np.divide(num, den * 4 * np.pi, out=np.zeros_like(num), where=np.abs(den)>tol)

    return result

def _compute_finite_vortex_deriv1(r1, r2, r1_deriv):
    r1_norm = add_ones_axis(compute_norm(r1))
    r2_norm = add_ones_axis(compute_norm(r2))
    r1_norm_deriv = compute_norm_deriv(r1, r1_deriv)

    r1_x_r2 = add_ones_axis(compute_cross(r1, r2))
    r1_d_r2 = add_ones_axis(compute_dot(r1, r2))
    r1_x_r2_deriv = compute_cross_deriv1(r1_deriv, r2)
    r1_d_r2_deriv = compute_dot_deriv(r2, r1_deriv)

    num = (1. / r1_norm + 1. / r2_norm) * r1_x_r2
    num_deriv = (-r1_norm_deriv / r1_norm ** 2) * r1_x_r2 \
        + (1. / r1_norm + 1. / r2_norm) * r1_x_r2_deriv

    den = r1_norm * r2_norm + r1_d_r2
    den_deriv = r1_norm_deriv * r2_norm + r1_d_r2_deriv

    result = np.divide(
        num_deriv * den - num * den_deriv,
        den ** 2 * 4 * np.pi,
        out=np.zeros_like(num),
        where=np.abs(den)>tol
        )

    return result

def _compute_finite_vortex_deriv2(r1, r2, r2_deriv):
    r1_norm = add_ones_axis(compute_norm(r1))
    r2_norm = add_ones_axis(compute_norm(r2))
    r2_norm_deriv = compute_norm_deriv(r2, r2_deriv)

    r1_x_r2 = add_ones_axis(compute_cross(r1, r2))
    r1_d_r2 = add_ones_axis(compute_dot(r1, r2))
    r1_x_r2_deriv = compute_cross_deriv2(r1, r2_deriv)
    r1_d_r2_deriv = compute_dot_deriv(r1, r2_deriv)

    num = (1. / r1_norm + 1. / r2_norm) * r1_x_r2
    num_deriv = (-r2_norm_deriv / r2_norm ** 2) * r1_x_r2 \
        + (1. / r1_norm + 1. / r2_norm) * r1_x_r2_deriv

    den = r1_norm * r2_norm + r1_d_r2
    den_deriv = r1_norm * r2_norm_deriv + r1_d_r2_deriv

    result = np.divide(
        num_deriv * den - num * den_deriv,
        den ** 2 * 4 * np.pi,
        out=np.zeros_like(num),
        where=np.abs(den)>tol
        )

    return result

def _compute_semi_infinite_vortex(u, r):
    r_norm = compute_norm(r)
    u_x_r = compute_cross(u, r)
    u_d_r = compute_dot(u, r)

    num = u_x_r
    den = r_norm * (r_norm - u_d_r)
    return num / den / 4 / np.pi

def _compute_semi_infinite_vortex_deriv(u, r, r_deriv):
    r_norm = add_ones_axis(compute_norm(r))
    r_norm_deriv = compute_norm_deriv(r, r_deriv)

    u_x_r = add_ones_axis(compute_cross(u, r))
    u_x_r_deriv = compute_cross_deriv2(u, r_deriv)

    u_d_r = add_ones_axis(compute_dot(u, r))
    u_d_r_deriv = compute_dot_deriv(u, r_deriv)

    num = u_x_r
    num_deriv = u_x_r_deriv

    den = r_norm * (r_norm - u_d_r)
    den_deriv = r_norm_deriv * (r_norm - u_d_r) + r_norm * (r_norm_deriv - u_d_r_deriv)

    return (num_deriv * den - num * den_deriv) / den ** 2 / 4 / np.pi


class EvalVelMtx(om.ExplicitComponent):
    """
    Computes the aerodynamic influence coefficient (AIC) matrix for the VLM
    analysis.

    This component is used in two places a given model, first to
    construct the AIC matrix using the collocation points as evaluation points,
    then to construct the AIC matrix where the force points are the evaluation
    points. The first matrix is used to solve for the circulations, while
    the second matrix is used to compute the forces acting on each panel.

    These calculations are rather complicated for a few reasons.
    Each surface interacts with every other surface, including itself.
    Also, in the general case, we have panel in both the spanwise and chordwise
    directions for all surfaces.
    Because of that, we need to compute the influence of each panel on every
    other panel, which results in rather large arrays for the
    intermediate calculations. Accordingly, the derivatives are complicated.

    The actual calcuations done here vary a fair bit in the case of symmetry.
    Not because the physics change, but because we need to account for a
    "ghost" version of the lifting surface, where we want to add the effects
    from the panels across the symmetry plane, but we don't want to actually
    use any of the evaluation points since we're not interested in the
    performance of this "ghost" version, since it's exactly symmetrical.
    This basically results in us looping through more calculations as if the
    panels were actually there.

    Parameters
    ----------
    alpha : float
        The angle of attack for the aircraft (all lifting surfaces) in degrees.
    vectors[num_eval_points, nx, ny, 3] : numpy array
        The vectors from the aerodynamic meshes to the evaluation points for
        every surface to every surface. For the symmetric case, the third
        dimension is length (2 * ny - 1). There is one of these arrays
        for each lifting surface in the problem.

    Returns
    -------
    vel_mtx[num_eval_points, nx - 1, ny - 1, 3] : numpy array
        The AIC matrix for the all lifting surfaces representing the aircraft.
        This has some sparsity pattern, but it is more dense than the FEM matrix
        and the entries have a wide range of magnitudes. One exists for each
        combination of surface name and evaluation points name.
    """

    def initialize(self):
        self.options.declare('surfaces', types=list)
        self.options.declare('eval_name', types=str)
        self.options.declare('num_eval_points', types=int)

    def setup(self):
        surfaces = self.options['surfaces']
        eval_name = self.options['eval_name']
        num_eval_points = self.options['num_eval_points']

        self.add_input('alpha', val=1., units='deg')

        for surface in surfaces:
            mesh=surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']

            # Get the names for the vectors and vel_mtx. We have the lifting
            # surface name coming in here, as well as the eval_name.
            vectors_name = '{}_{}_vectors'.format(name, eval_name)
            vel_mtx_name = '{}_{}_vel_mtx'.format(name, eval_name)

            # Here we set up the rows and cols for the sparse Jacobians.

            # The logic differs if the surface is symmetric or not, due to the
            # existence of the "ghost" surface; the reflection of the actual.
            if surface['symmetry']:
                self.add_input(vectors_name, shape=(num_eval_points, nx, 2*ny-1, 3), units='m')

                # Get an array of indices representing the number of entries
                # in the vectors array.
                vectors_indices = np.arange(num_eval_points * nx * (2*ny-1) * 3).reshape(
                    (num_eval_points, nx, (2*ny-1), 3))

                # Set up blocks to mannipulate into rows for the sparse indices
                base = np.tile(np.repeat(np.arange(3), 3), ny-1)
                block1 = base + np.repeat(3*np.arange(ny-1), 9)
                block2 = base + np.flip(np.repeat(3*np.arange(ny-1), 9), axis=0)
                block3 = np.concatenate([block1, block2])
                block4 = np.tile(block3, nx-1)
                block5 = block4 + np.repeat(3*(ny-1)*np.arange(nx-1), len(block3))
                block6 = np.tile(block5, num_eval_points)
                row = block6 + np.repeat(3*(ny-1)*(nx-1)*np.arange(num_eval_points), len(block5))

                rows = np.tile(row, 4)

                # Create the columns for each of the tiled out rows based on the
                # previously-assembled vectors_indices.
                cols = np.concatenate([
                    np.einsum('ijkm,l->ijklm', vectors_indices[:, 0:-1, 0:-1, :], np.ones(3, int)).flatten(),
                    np.einsum('ijkm,l->ijklm', vectors_indices[:, 1:  , 0:-1, :], np.ones(3, int)).flatten(),
                    np.einsum('ijkm,l->ijklm', vectors_indices[:, 0:-1, 1:  , :], np.ones(3, int)).flatten(),
                    np.einsum('ijkm,l->ijklm', vectors_indices[:, 1:  , 1:  , :], np.ones(3, int)).flatten(),
                ])

                # Layout logic includes some duplicate entries due to symmetry. Find and remove them.
                nn = len(rows) // 2

                # Determine the repeated indices and store them in an array
                inds = np.arange(nn).reshape((-1, 9))
                to_remove = inds[(ny-1)::2*(ny-1)].flatten()

                # Actually remove the duplicate entries
                rows = np.delete(rows, to_remove)
                cols = np.delete(cols, to_remove)

            # In the nonsymmetric case, the derivative sparsity patterns are
            # much more straightforward.
            else:
                self.add_input(vectors_name, shape=(num_eval_points, nx, ny, 3), units='m')

                vectors_indices = np.arange(num_eval_points * nx * ny * 3).reshape(
                    (num_eval_points, nx, ny, 3))
                vel_mtx_indices = np.arange(num_eval_points * (nx - 1) * (ny - 1) * 3).reshape(
                    (num_eval_points, nx - 1, ny - 1, 3))

                rows = np.tile(np.einsum('ijkl,m->ijklm', vel_mtx_indices, np.ones(3, int)).flatten(), 4)

                cols = np.concatenate([
                    np.einsum('ijkm,l->ijklm', vectors_indices[:, 0:-1, 0:-1, :], np.ones(3, int)).flatten(),
                    np.einsum('ijkm,l->ijklm', vectors_indices[:, 1:  , 0:-1, :], np.ones(3, int)).flatten(),
                    np.einsum('ijkm,l->ijklm', vectors_indices[:, 0:-1, 1:  , :], np.ones(3, int)).flatten(),
                    np.einsum('ijkm,l->ijklm', vectors_indices[:, 1:  , 1:  , :], np.ones(3, int)).flatten(),
                ])

            self.add_output(vel_mtx_name, shape=(num_eval_points, nx - 1, ny - 1, 3), units='1/m')

            self.declare_partials(vel_mtx_name, vectors_name, rows=rows, cols=cols)

            # It's worth the cs cost here because alpha is just a scalar
            self.declare_partials(vel_mtx_name, 'alpha', method='cs')

            self.set_check_partial_options(wrt='*', method='cs')

    def compute(self, inputs, outputs):
        surfaces = self.options['surfaces']
        eval_name = self.options['eval_name']
        num_eval_points = self.options['num_eval_points']

        for surface in surfaces:
            ny = surface['mesh'].shape[1]
            name = surface['name']

            alpha = inputs['alpha'][0]
            cosa = np.cos(alpha * np.pi / 180.)
            sina = np.sin(alpha * np.pi / 180.)

            if surface['symmetry']:
                u = np.einsum('ijk,l->ijkl',
                    np.ones((num_eval_points, 1, 2*(ny - 1))),
                    np.array([cosa, 0, sina]))
            else:
                u = np.einsum('ijk,l->ijkl',
                    np.ones((num_eval_points, 1, ny - 1)),
                    np.array([cosa, 0, sina]))

            vectors_name = '{}_{}_vectors'.format(name, eval_name)
            vel_mtx_name = '{}_{}_vel_mtx'.format(name, eval_name)

            outputs[vel_mtx_name] = 0.

            # Here, we loop through each of the vectors and compute the AIC
            # terms from the four filaments that make up a ring around a single
            # panel. Thus, we are using vortex rings to construct the AIC
            # matrix. Later, we will convert these to horseshoe vortices
            # to compute the panel forces.

            # front vortex
            r1 = inputs[vectors_name][:, 0:-1, 1:  , :]
            r2 = inputs[vectors_name][:, 0:-1, 0:-1, :]
            result1 = _compute_finite_vortex(r1, r2)

            # right vortex
            r1 = inputs[vectors_name][:, 0:-1, 0:-1, :]
            r2 = inputs[vectors_name][:, 1:  , 0:-1, :]
            result2 = _compute_finite_vortex(r1, r2)

            # rear vortex
            r1 = inputs[vectors_name][:, 1:  , 0:-1, :]
            r2 = inputs[vectors_name][:, 1:  , 1:  , :]
            result3 = _compute_finite_vortex(r1, r2)

            # left vortex
            r1 = inputs[vectors_name][:, 1:  , 1:  , :]
            r2 = inputs[vectors_name][:, 0:-1, 1:  , :]
            result4 = _compute_finite_vortex(r1, r2)

            # If the surface is symmetric, mirror the results and add them
            # to the vel_mtx.
            if surface['symmetry']:
                res1 = result1[:, :, :ny-1, :]
                res1 += result1[:, :, ny-1:, :][:, :, ::-1, :]
                res2 = result2[:, :, :ny-1, :]
                res2 += result2[:, :, ny-1:, :][:, :, ::-1, :]
                res3 = result3[:, :, :ny-1, :]
                res3 += result3[:, :, ny-1:, :][:, :, ::-1, :]
                res4 = result4[:, :, :ny-1, :]
                res4 += result4[:, :, ny-1:, :][:, :, ::-1, :]
                outputs[vel_mtx_name] += res1 + res2 + res3 + res4
            else:
                outputs[vel_mtx_name] += result1 + result2 + result3 + result4

            # ----------------- last row -----------------

            r1 = inputs[vectors_name][:, -1:, 1:  , :]
            r2 = inputs[vectors_name][:, -1:, 0:-1, :]
            result1 = _compute_finite_vortex(r1, r2)
            result2 = _compute_semi_infinite_vortex(u, r1)
            result3 = _compute_semi_infinite_vortex(u, r2)

            if surface['symmetry']:
                res1 = result1[:, :, :ny-1, :]
                res1 += result1[:, :, ny-1:, :][:, :, ::-1, :]
                res2 = result2[:, :, :ny-1, :]
                res2 += result2[:, :, ny-1:, :][:, :, ::-1, :]
                res3 = result3[:, :, :ny-1, :]
                res3 += result3[:, :, ny-1:, :][:, :, ::-1, :]
                outputs[vel_mtx_name][:, -1:, :, :] += res1 - res2 + res3
            else:
                outputs[vel_mtx_name][:, -1:, :, :] += result1
                outputs[vel_mtx_name][:, -1:, :, :] -= result2
                outputs[vel_mtx_name][:, -1:, :, :] += result3

    def compute_partials(self, inputs, partials):
        surfaces = self.options['surfaces']
        eval_name = self.options['eval_name']
        num_eval_points = self.options['num_eval_points']

        for surface in surfaces:
            nx = surface['mesh'].shape[0]
            ny = surface['mesh'].shape[1]
            name = surface['name']

            vectors_name = '{}_{}_vectors'.format(name, eval_name)
            vel_mtx_name = '{}_{}_vel_mtx'.format(name, eval_name)

            alpha = inputs['alpha'][0]
            cosa = np.cos(alpha * np.pi / 180.)
            sina = np.sin(alpha * np.pi / 180.)

            if surface['symmetry']:

                u = np.einsum('ijk,l->ijkl',
                    np.ones((num_eval_points, 1, 2*(ny - 1))),
                    np.array([cosa, 0, sina]))

                deriv_array = np.einsum('...,ij->...ij',
                    np.ones((num_eval_points, nx - 1, 2*(ny - 1))),
                    np.eye(3))
                trailing_array = np.einsum('...,ij->...ij',
                    np.ones((num_eval_points, 1, 2*(ny - 1))),
                    np.eye(3))

                derivs0 = np.zeros((num_eval_points, nx - 1, 2*(ny - 1) - 1, 3, 3))
                derivs1 = np.zeros((num_eval_points, nx - 1, 2*(ny - 1) - 1, 3, 3))
                derivs2 = np.zeros((num_eval_points, nx - 1, 2*(ny - 1), 3, 3))
                derivs3 = np.zeros((num_eval_points, nx - 1, 2*(ny - 1), 3, 3))

                # front vortex
                r1 = inputs[vectors_name][:, 0:-1, 1:  , :]
                r2 = inputs[vectors_name][:, 0:-1, 0:-1, :]
                d1 = _compute_finite_vortex_deriv1(r1, r2, deriv_array)
                d2 = _compute_finite_vortex_deriv2(r1, r2, deriv_array)
                derivs2[:, :, :ny-1, :, :] += d1[:, :, :ny-1, :, :]
                derivs0[:, :, :ny-1, :, :] += d2[:, :, :ny-1, :, :]
                derivs2[:, :, ny-1:, :, :] += d1[:, :, ny-1:, :, :]
                derivs0[:, :, ny-1:, :, :] += d2[:, :, ny:, :, :]

                # Formerly duplicated location
                derivs2[:, :, ny-2, :, :] += d2[:, :, ny-1, :, :]

                # right vortex
                r1 = inputs[vectors_name][:, 0:-1, 0:-1, :]
                r2 = inputs[vectors_name][:, 1:  , 0:-1, :]
                d1 = _compute_finite_vortex_deriv1(r1, r2, deriv_array)
                d2 = _compute_finite_vortex_deriv2(r1, r2, deriv_array)
                derivs0[:, :, :ny-1, :, :] += d1[:, :, :ny-1, :, :]
                derivs1[:, :, :ny-1, :, :] += d2[:, :, :ny-1, :, :]
                derivs0[:, :, ny-1:, :, :] += d1[:, :, ny:, :, :]
                derivs1[:, :, ny-1:, :] += d2[:, :, ny:, :, :]

                # Formerly duplicated location
                derivs2[:, :, ny-2, :, :] += d1[:, :, ny-1, :, :]
                derivs3[:, :, ny-2, :, :] += d2[:, :, ny-1, :, :]

                # rear vortex
                r1 = inputs[vectors_name][:, 1:  , 0:-1, :]
                r2 = inputs[vectors_name][:, 1:  , 1:  , :]
                d1 = _compute_finite_vortex_deriv1(r1, r2, deriv_array)
                d2 = _compute_finite_vortex_deriv2(r1, r2, deriv_array)
                derivs1[:, :, :ny-1, :, :] += d1[:, :, :ny-1, :, :]
                derivs3[:, :, :ny-1, :, :] += d2[:, :, :ny-1, :, :]
                derivs1[:, :, ny-1:, :] += d1[:, :, ny:, :, :]
                derivs3[:, :, ny-1:, :] += d2[:, :, ny-1:, :, :]

                # Formerly duplicated location
                derivs3[:, :, ny-2, :, :] += d1[:, :, ny-1, :, :]

                # left vortex
                r1 = inputs[vectors_name][:, 1:  , 1:  , :]
                r2 = inputs[vectors_name][:, 0:-1, 1:  , :]
                d1 = _compute_finite_vortex_deriv1(r1, r2, deriv_array)
                d2 = _compute_finite_vortex_deriv2(r1, r2, deriv_array)
                derivs3[:, :, :ny-1, :, :] += d1[:, :, :ny-1, :, :]
                derivs2[:, :, :ny-1, :, :] += d2[:, :, :ny-1, :, :]
                derivs3[:, :, ny-1:, :] += d1[:, :, ny-1:, :, :]
                derivs2[:, :, ny-1:, :] += d2[:, :, ny-1:, :, :]

                #----------------- last row -----------------

                r1 = inputs[vectors_name][:, -1:, 1:  , :]
                r2 = inputs[vectors_name][:, -1:, 0:-1, :]
                d1 = _compute_finite_vortex_deriv1(r1, r2, trailing_array)
                d2 = _compute_finite_vortex_deriv2(r1, r2, trailing_array)
                d3 = _compute_semi_infinite_vortex_deriv(u, r1, trailing_array)
                d4 = _compute_semi_infinite_vortex_deriv(u, r2, trailing_array)
                derivs3[:, -1:, :ny-1, :] += d1[:, :, :ny-1, :, :]
                derivs1[:, -1:, :ny-1, :] += d2[:, :, :ny-1, :, :]
                derivs3[:, -1:, :ny-1, :] -= d3[:, :, :ny-1, :, :]
                derivs1[:, -1:, :ny-1, :] += d4[:, :, :ny-1, :, :]
                derivs3[:, -1:, ny-1:, :] += d1[:, :, ny-1:, :, :]
                derivs1[:, -1:, ny-1:, :] += d2[:, :, ny:, :, :]
                derivs3[:, -1:, ny-1:, :] -= d3[:, :, ny-1:, :, :]
                derivs1[:, -1:, ny-1:, :] += d4[:, :, ny:, :, :]

                # Formerly duplicated location
                derivs3[:, -1:, ny-2, :, :] += d2[:, :, ny-1, :, :]
                derivs3[:, -1:, ny-2, :, :] += d4[:, :, ny-1, :, :]

                partials[vel_mtx_name, vectors_name] = np.concatenate([
                    derivs0.flatten(),
                    derivs1.flatten(),
                    derivs2.flatten(),
                    derivs3.flatten(),
                ])

            else:
                u = np.einsum('ijk,l->ijkl',
                    np.ones((num_eval_points, 1, ny - 1)),
                    np.array([cosa, 0, sina]))

                deriv_array = np.einsum('...,ij->...ij',
                    np.ones((num_eval_points, nx - 1, ny - 1)),
                    np.eye(3))
                trailing_array = np.einsum('...,ij->...ij',
                    np.ones((num_eval_points, 1, ny - 1)),
                    np.eye(3))

                derivs = np.zeros((4, num_eval_points, nx - 1, ny - 1, 3, 3))

                # front vortex
                r1 = inputs[vectors_name][:, 0:-1, 1:  , :]
                r2 = inputs[vectors_name][:, 0:-1, 0:-1, :]
                derivs[2, :, :, :, :] += _compute_finite_vortex_deriv1(r1, r2, deriv_array)
                derivs[0, :, :, :, :] += _compute_finite_vortex_deriv2(r1, r2, deriv_array)

                # right vortex
                r1 = inputs[vectors_name][:, 0:-1, 0:-1, :]
                r2 = inputs[vectors_name][:, 1:  , 0:-1, :]
                derivs[0, :, :, :, :] += _compute_finite_vortex_deriv1(r1, r2, deriv_array)
                derivs[1, :, :, :, :] += _compute_finite_vortex_deriv2(r1, r2, deriv_array)

                # rear vortex
                r1 = inputs[vectors_name][:, 1:  , 0:-1, :]
                r2 = inputs[vectors_name][:, 1:  , 1:  , :]
                derivs[1, :, :, :, :] += _compute_finite_vortex_deriv1(r1, r2, deriv_array)
                derivs[3, :, :, :, :] += _compute_finite_vortex_deriv2(r1, r2, deriv_array)

                # left vortex
                r1 = inputs[vectors_name][:, 1:  , 1:  , :]
                r2 = inputs[vectors_name][:, 0:-1, 1:  , :]
                derivs[3, :, :, :, :] += _compute_finite_vortex_deriv1(r1, r2, deriv_array)
                derivs[2, :, :, :, :] += _compute_finite_vortex_deriv2(r1, r2, deriv_array)

                # ----------------- last row -----------------

                r1 = inputs[vectors_name][:, -1:, 1:  , :]
                r2 = inputs[vectors_name][:, -1:, 0:-1, :]
                derivs[3, :, -1:, :, :] += _compute_finite_vortex_deriv1(r1, r2, trailing_array)
                derivs[1, :, -1:, :, :] += _compute_finite_vortex_deriv2(r1, r2, trailing_array)
                derivs[3, :, -1:, :, :] -= _compute_semi_infinite_vortex_deriv(u, r1, trailing_array)
                derivs[1, :, -1:, :, :] += _compute_semi_infinite_vortex_deriv(u, r2, trailing_array)

                partials[vel_mtx_name, vectors_name] = derivs.flatten()
