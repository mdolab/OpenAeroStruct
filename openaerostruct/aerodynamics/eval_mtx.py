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

    The calculations also vary when we consider ground effect.
    This is accomplished by mirroring a second copy of the mesh across
    the ground plane. The documentation has more detailed explanations.
    The ground effect is only implemented for symmetric wings.

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

        self.surface_indices_repeated = dict()

        for surface in surfaces:
            mesh=surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']


            ground_effect = surface.get('groundplane', False)

            # Get the names for the vectors and vel_mtx. We have the lifting
            # surface name coming in here, as well as the eval_name.
            vectors_name = '{}_{}_vectors'.format(name, eval_name)
            vel_mtx_name = '{}_{}_vel_mtx'.format(name, eval_name)

            # Here we set up the rows and cols for the sparse Jacobians.

            # The logic differs if the surface is symmetric or not, due to the
            # existence of the "ghost" surface; the reflection of the actual.
            if ground_effect:
                nx_actual = 2*nx
            else:
                nx_actual = nx
            if surface['symmetry']:
                ny_actual = 2*ny - 1
                duplicate_jac_entry_idx_set_1 = np.array([], int)
                duplicate_jac_entry_idx_set_2 = np.array([], int)
                jac_start_ind_running_total = 0
            else:
                ny_actual = ny

            self.add_input(vectors_name, shape=(num_eval_points, nx_actual, ny_actual, 3), units='m')

            # Get an array of indices representing the number of entries
            # in the vectors array.
            vectors_indices = np.arange(num_eval_points * nx_actual * ny_actual * 3).reshape(
                (num_eval_points, nx_actual, ny_actual, 3))
            vel_mtx_indices = np.arange(num_eval_points * (nx - 1) * (ny - 1) * 3).reshape(
                (num_eval_points, nx - 1, ny - 1, 3))
            vel_mtx_idx_expanded = np.arange(num_eval_points * (nx - 1) * (ny - 1) * 3 * 3).reshape(
                (num_eval_points, nx - 1, ny - 1, 3, 3))
            aic_base = np.einsum('ijkl,m->ijklm', vel_mtx_indices, np.ones(3, int))
            aic_len = np.sum(np.product(aic_base.shape))

            if ground_effect:
                # mirrored surface along the x mesh direction
                surfaces_to_compute = [vectors_indices[:,:nx,:], vectors_indices[:,nx:,:]]
            else:
                surfaces_to_compute = [vectors_indices[:,:,:]]

            rows = np.array([], int)
            cols = np.array([], int)

            for surface_to_compute in surfaces_to_compute:
                inds_A = surface_to_compute[:, 0:-1, 1:  , :]
                inds_B = surface_to_compute[:, 0:-1, 0:-1, :]
                inds_C = surface_to_compute[:, 1:  , 0:-1, :]
                inds_D = surface_to_compute[:, 1:  , 1:  , :]
                vertices_to_compute = [inds_A, inds_B, inds_C, inds_D]
                # symmetric meshes end up with duplicated jacobian entries that need to be deleted later
                # vertices A and D duplicate their last entries y-wise
                # vertices B and C duplicate their first entries y-wise
                jac_dup_sets = [1, 2, 2, 1]
                for ivert, vertex_to_compute in enumerate(vertices_to_compute):
                    jac_dup_set = jac_dup_sets[ivert]
                    if surface['symmetry']:
                        rows = np.concatenate([rows, aic_base.flatten()])
                        cols = np.concatenate([cols, np.einsum('ijkm,l->ijklm', vertex_to_compute[:,:,:ny-1,:], np.ones(3, int)).flatten()])
                        if jac_dup_set == 1:
                            duplicate_jac_entry_idx_set_1 = np.concatenate([duplicate_jac_entry_idx_set_1,
                                                                            jac_start_ind_running_total + vel_mtx_idx_expanded[:,:,-1,:,:].flatten()])
                        jac_start_ind_running_total += aic_len

                        rows = np.concatenate([rows, aic_base[:,:,::-1,:].flatten()])
                        cols = np.concatenate([cols, np.einsum('ijkm,l->ijklm', vertex_to_compute[:,:,ny-1:,:], np.ones(3, int)).flatten()])
                        if jac_dup_set == 2:
                            duplicate_jac_entry_idx_set_2 = np.concatenate([duplicate_jac_entry_idx_set_2,
                                                                       jac_start_ind_running_total + vel_mtx_idx_expanded[:,:,0,:,:].flatten()])
                        jac_start_ind_running_total += aic_len

                    else:
                        rows = np.concatenate([rows, aic_base.flatten()])
                        cols = np.concatenate([cols, np.einsum('ijkm,l->ijklm', vertex_to_compute[:,:,:,:], np.ones(3, int)).flatten()])

            if surface['symmetry']:
                # need to determine the location of duplicate indices, knock them out, and save the locations for compute_partials
                self.surface_indices_repeated[name] = [duplicate_jac_entry_idx_set_1.copy(), duplicate_jac_entry_idx_set_2.copy()]

                cols = np.delete(cols, duplicate_jac_entry_idx_set_2)
                rows = np.delete(rows, duplicate_jac_entry_idx_set_2)

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
            nx = surface['mesh'].shape[0]
            ny = surface['mesh'].shape[1]
            name = surface['name']
            ground_effect = surface.get('groundplane', False)

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

            if ground_effect:
                # mirrored surface along the x mesh direction
                surfaces_to_compute = [inputs[vectors_name][:,:nx,:,:], inputs[vectors_name][:,nx:,:,:]]
                vortex_mults = [1.0, -1.0]
            else:
                surfaces_to_compute = [inputs[vectors_name]]
                vortex_mults = [1.0]

            for i_surf, surface_to_compute in enumerate(surfaces_to_compute):
                # vortex vertices:
                #         A ----- B
                #         |       |
                #         |       |
                #         D-------C
                #
                vortex_mult = vortex_mults[i_surf]
                vert_A = surface_to_compute[:, 0:-1, 1:  , :]
                vert_B = surface_to_compute[:, 0:-1, 0:-1, :]
                vert_C = surface_to_compute[:, 1:  , 0:-1, :]
                vert_D = surface_to_compute[:, 1:  , 1:  , :]
                # front vortex
                result1 = _compute_finite_vortex(vert_A, vert_B)
                # right vortex
                result2 = _compute_finite_vortex(vert_B, vert_C)
                # rear vortex
                result3 = _compute_finite_vortex(vert_C, vert_D)
                # left vortex
                result4 = _compute_finite_vortex(vert_D, vert_A)

                # If the surface is symmetric, mirror the results and add them
                # to the vel_mtx.
                if surface['symmetry']:
                    result = vortex_mult*(result1 + result2 + result3 + result4)
                    outputs[vel_mtx_name] += result[:, :, :ny-1, :]
                    outputs[vel_mtx_name] += result[:, :, ny-1:, :][:,:,::-1,:]
                else:
                    outputs[vel_mtx_name] += vortex_mult*(result1 + result2 + result3 + result4)

                # ----------------- last row -----------------

                vert_D_last = vert_D[:,-1:,:,:]
                vert_C_last = vert_C[:, -1:, :,:]
                result1 = _compute_finite_vortex(vert_D_last, vert_C_last)
                result2 = _compute_semi_infinite_vortex(u, vert_D_last)
                result3 = _compute_semi_infinite_vortex(u, vert_C_last)

                if surface['symmetry']:
                    res1 = result1[:, :, :ny-1, :]
                    res1 += result1[:, :, ny-1:, :][:, :, ::-1, :]
                    res2 = result2[:, :, :ny-1, :]
                    res2 += result2[:, :, ny-1:, :][:, :, ::-1, :]
                    res3 = result3[:, :, :ny-1, :]
                    res3 += result3[:, :, ny-1:, :][:, :, ::-1, :]
                    outputs[vel_mtx_name][:, -1:, :, :] += vortex_mult*(res1 - res2 + res3)
                else:
                    outputs[vel_mtx_name][:, -1:, :, :] += vortex_mult*result1
                    outputs[vel_mtx_name][:, -1:, :, :] -= vortex_mult*result2
                    outputs[vel_mtx_name][:, -1:, :, :] += vortex_mult*result3



    def compute_partials(self, inputs, partials):
        surfaces = self.options['surfaces']
        eval_name = self.options['eval_name']
        num_eval_points = self.options['num_eval_points']
        for surface in surfaces:
            nx = surface['mesh'].shape[0]
            ny = surface['mesh'].shape[1]
            name = surface['name']
            ground_effect = surface.get('groundplane', False)


            vectors_name = '{}_{}_vectors'.format(name, eval_name)
            vel_mtx_name = '{}_{}_vel_mtx'.format(name, eval_name)

            alpha = inputs['alpha'][0]
            cosa = np.cos(alpha * np.pi / 180.)
            sina = np.sin(alpha * np.pi / 180.)

            if ground_effect:
                # mirrored surface along the x mesh direction
                surfaces_to_compute = [inputs[vectors_name][:,:nx,:,:], inputs[vectors_name][:,nx:,:,:]]
                vortex_mults = [1.0, -1.0]
            else:
                surfaces_to_compute = [inputs[vectors_name]]
                vortex_mults = [1.0]

            if surface['symmetry']:
                ny_actual = 2*ny-1
            else:
                ny_actual = ny

            assembled_derivs = np.array([], int)

            for i_surf, surface_to_compute in enumerate(surfaces_to_compute):
                    # vortex vertices:
                    #         A ----- B
                    #         |       |
                    #         |       |
                    #         D-------C
                    #
                vortex_mult = vortex_mults[i_surf]


                u = np.einsum('ijk,l->ijkl',
                    np.ones((num_eval_points, 1, ny_actual - 1)),
                    np.array([cosa, 0, sina]))

                deriv_array = np.einsum('...,ij->...ij',
                    np.ones((num_eval_points, nx - 1, ny_actual - 1)),
                    np.eye(3))
                trailing_array = np.einsum('...,ij->...ij',
                    np.ones((num_eval_points, 1, ny_actual - 1)),
                    np.eye(3))

                derivs = np.zeros((4, num_eval_points, nx - 1, ny_actual - 1, 3, 3))

                vert_A = surface_to_compute[:, 0:-1, 1:  , :]
                vert_B = surface_to_compute[:, 0:-1, 0:-1, :]
                vert_C = surface_to_compute[:, 1:  , 0:-1, :]
                vert_D = surface_to_compute[:, 1:  , 1:  , :]

                # front vortex
                derivs[0, :, :, :, :] += _compute_finite_vortex_deriv1(vert_A, vert_B, deriv_array)
                derivs[1, :, :, :, :] += _compute_finite_vortex_deriv2(vert_A, vert_B, deriv_array)

                # right vortex
                derivs[1, :, :, :, :] += _compute_finite_vortex_deriv1(vert_B, vert_C, deriv_array)
                derivs[2, :, :, :, :] += _compute_finite_vortex_deriv2(vert_B, vert_C, deriv_array)

                # rear vortex
                derivs[2, :, :, :, :] += _compute_finite_vortex_deriv1(vert_C, vert_D, deriv_array)
                derivs[3, :, :, :, :] += _compute_finite_vortex_deriv2(vert_C, vert_D, deriv_array)

                # left vortex
                derivs[3, :, :, :, :] += _compute_finite_vortex_deriv1(vert_D, vert_A, deriv_array)
                derivs[0, :, :, :, :] += _compute_finite_vortex_deriv2(vert_D, vert_A, deriv_array)

                # ----------------- last row -----------------
                vert_D_last = vert_D[:,-1:,:,:]
                vert_C_last = vert_C[:,-1:, :,:]

                derivs[3, :, -1:, :, :] += _compute_finite_vortex_deriv1(vert_D_last, vert_C_last, trailing_array)
                derivs[2, :, -1:, :, :] += _compute_finite_vortex_deriv2(vert_D_last, vert_C_last, trailing_array)
                derivs[3, :, -1:, :, :] -= _compute_semi_infinite_vortex_deriv(u, vert_D_last, trailing_array)
                derivs[2, :, -1:, :, :] += _compute_semi_infinite_vortex_deriv(u, vert_C_last, trailing_array)

                derivs = derivs * vortex_mult
                for i in range(4):
                    if surface['symmetry']:
                        assembled_derivs = np.concatenate([assembled_derivs, derivs[i,:,:,:ny-1,:].flatten()])
                        assembled_derivs = np.concatenate([assembled_derivs, derivs[i,:,:,ny-1:,:].flatten()])
                    else:
                        assembled_derivs = np.concatenate([assembled_derivs, derivs[i,:,:,:].flatten()])

            if surface['symmetry']:
                    # now, need to check for duplicate entries and combine / delete
                    first_repeated_index, second_repeated_index = self.surface_indices_repeated[name]
                    assembled_derivs[first_repeated_index] += assembled_derivs[second_repeated_index].copy()
                    assembled_derivs = np.delete(assembled_derivs, second_repeated_index)
            partials[vel_mtx_name, vectors_name] = assembled_derivs
