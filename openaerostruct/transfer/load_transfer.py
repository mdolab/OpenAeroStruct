from __future__ import division, print_function

import numpy as np

import openmdao.api as om


class LoadTransfer(om.ExplicitComponent):
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
    LoadsA[ny,3] : numpy array  = loads[:,:3]
    LoadsB[ny,3]: nympy array   = loads[:,3:]
        Flattened array containing the loads applied on the FEM component,
        computed from the sectional forces.

    Loads[ny, 6] : numpy array = [LoadsA,LoadsB]

    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        self.surface = surface = self.options['surface']

        self.nx = nx = surface['mesh'].shape[0]
        self.ny = ny = surface['mesh'].shape[1]


        if surface['fem_model_type'] == 'tube':
            self.fem_origin = surface['fem_origin']
        else:
            y_upper = surface['data_y_upper']
            x_upper = surface['data_x_upper']
            y_lower = surface['data_y_lower']

            fem_origin = (x_upper[0]  * (y_upper[0]  - y_lower[0]) +
                          x_upper[-1] * (y_upper[-1] - y_lower[-1])) / \
                        ((y_upper[0]  -  y_lower[0]) + (y_upper[-1] - y_lower[-1]))

            # For some reason, surface data is complex in some tests.
            self.fem_origin = np.float(fem_origin)

        self.w1 = 0.25
        self.w2 = self.fem_origin

        self.add_input('def_mesh', val=np.zeros((nx, ny, 3)), units='m')
        self.add_input('sec_forces', val=np.zeros((nx-1, ny-1, 3)), units='N')

        self.add_output('loads', val=np.zeros((self.ny, 6)), units='N') ## WARNING!!! UNITS ARE A MIXTURE OF N & N*m
        # Well, technically the units of this load array are mixed.
        # The first 3 indices are N and the last 3 are N*m.

        # Derivatives

        # First, the direct loads wrt sec_forces terms.
        base_row = np.array([0, 1, 2, 6, 7, 8])
        base_col = np.array([0, 1, 2, 0, 1, 2])
        row = np.tile(base_row, ny-1) + np.repeat(6*np.arange(ny-1), 6)
        col = np.tile(base_col, ny-1) + np.repeat(3*np.arange(ny-1), 6)
        rows1 = np.tile(row, nx-1)
        cols1 = np.tile(col, nx-1) + np.repeat(3*(ny-1)*np.arange(nx-1), 6*(ny-1))

        # Then, the term from the cross product.
        base_row = np.array([3, 3, 4, 4, 5, 5])
        base_col = np.array([1, 2, 0, 2, 0, 1])
        row = np.tile(base_row, ny-1) + np.repeat(6*np.arange(ny-1), 6)
        col = np.tile(base_col, ny-1) + np.repeat(3*np.arange(ny-1), 6)
        row1 = np.tile(row, nx-1)
        col1 = np.tile(col, nx-1) + np.repeat(3*(ny-1)*np.arange(nx-1), 6*(ny-1))
        rows2 = np.tile(row1, 2) + np.repeat(np.array([0, 6]), 6*(nx-1)*(ny-1))
        cols2 = np.tile(col1, 2)

        rows = np.concatenate([rows1, rows2])
        cols = np.concatenate([cols1, cols2])

        self.declare_partials(of='loads', wrt='sec_forces', rows=rows, cols=cols)

        # Top diagonal is forward-most mesh point.
        base_row = np.array([3, 3, 4, 4, 5, 5])
        base_col = np.array([4, 5, 3, 5, 3, 4])
        row = np.tile(base_row, ny-1) + np.repeat(6*np.arange(ny-1), 6)
        col = np.tile(base_col, ny-1) + np.repeat(3*np.arange(ny-1), 6)
        rows1 = np.tile(row, nx)
        cols1 = np.tile(col, nx) + np.repeat(3*ny*np.arange(nx), 6*(ny-1))

        # Bottom diagonal is backward-most mesh point.
        base_row = np.array([9, 9, 10, 10, 11, 11])
        base_col = np.array([1, 2, 0, 2, 0, 1])
        row = np.tile(base_row, ny-1) + np.repeat(6*np.arange(ny-1), 6)
        col = np.tile(base_col, ny-1) + np.repeat(3*np.arange(ny-1), 6)
        rows2 = np.tile(row, nx)
        cols2 = np.tile(col, nx) + np.repeat(3*ny*np.arange(nx), 6*(ny-1))

        # Central Diagonal blocks
        base_row = np.array([3, 3, 4, 4, 5, 5])
        base_col = np.array([1, 2, 0, 2, 0, 1])
        row = np.tile(base_row, ny) + np.repeat(6*np.arange(ny), 6)
        col = np.tile(base_col, ny) + np.repeat(3*np.arange(ny), 6)
        rows3 = np.tile(row, nx)
        cols3 = np.tile(col, nx) + np.repeat(3*ny*np.arange(nx), 6*ny)


        rows = np.concatenate([rows1, rows2, rows3])
        cols = np.concatenate([cols1, cols2, cols3])

        self.declare_partials(of='loads', wrt='def_mesh', rows=rows, cols=cols)

        # -------------------------------- Check Partial Options-------------------------------------
        self.set_check_partial_options('*', method='cs', step=1e-40)

    def compute(self, inputs, outputs):
        mesh = inputs['def_mesh'] #[nx, ny, 3]
        sec_forces = inputs['sec_forces']

        # Compute the aerodynamic centers at the quarter-chord point of each panel
        # a_pts [nx-1, ny-1, 3]
        a_pts = 0.5 * (1-self.w1) * mesh[:-1, :-1, :] + \
                0.5 *   self.w1   * mesh[1:, :-1, :] + \
                0.5 * (1-self.w1) * mesh[:-1,  1:, :] + \
                0.5 *   self.w1   * mesh[1:,  1:, :]

        # Compute the structural midpoints based on the fem_origin location
        # s_pts [ny-1, 3]
        s_pts = 0.5 * (1-self.w2) * mesh[0, :-1, :] + \
                0.5 *   self.w2   * mesh[-1, :-1, :] + \
                0.5 * (1-self.w2) * mesh[0,  1:, :] + \
                0.5 *   self.w2   * mesh[-1,  1:, :]

        # Find the moment arm between the aerodynamic centers of each panel
        # and the FEM elements
        # diff [nx-1, ny-1, 3]
        moment = 0.5 * np.sum(np.cross(a_pts - s_pts, sec_forces), axis=0)

        # Only need to zero out the part that is assigned via +=
        outputs['loads'][-1, :] = 0.

        # Compute the loads based on the xyz forces and the computed moments
        sec_forces_sum = 0.5 * np.sum(sec_forces, axis=0)
        outputs['loads'][:-1, :3] = sec_forces_sum
        outputs['loads'][1:, :3] += sec_forces_sum

        outputs['loads'][:-1, 3:] = moment
        outputs['loads'][1:, 3:] += moment

    def compute_partials(self, inputs, partials):
        mesh = inputs['def_mesh']
        sec_forces = inputs['sec_forces']
        ny = self.ny
        nx = self.nx
        w1 = self.w1
        w2 = self.w2

        # Compute the aerodynamic centers at the quarter-chord point of each panel
        a_pts = 0.5 * (1-w1) * mesh[:-1, :-1, :] + \
                0.5 *   w1   * mesh[1:, :-1, :] + \
                0.5 * (1-w1) * mesh[:-1,  1:, :] + \
                0.5 *   w1   * mesh[1:,  1:, :]

        # Compute the structural midpoints based on the fem_origin location
        s_pts = 0.5 * (1-w2) * mesh[0, :-1, :] + \
                0.5 *   w2   * mesh[-1, :-1, :] + \
                0.5 * (1-w2) * mesh[0,  1:, :] + \
                0.5 *   w2   * mesh[-1,  1:, :]

        diff = 0.5 * (a_pts - s_pts)

        # dmoment__dsec_forces

        dmom_dsec = np.empty((nx-1, ny-1, 6))
        dmom_dsec[:, :, 0] = -diff[:, :, 2]
        dmom_dsec[:, :, 1] = diff[:, :, 1]
        dmom_dsec[:, :, 2] = diff[:, :, 2]
        dmom_dsec[:, :, 3] = -diff[:, :, 0]
        dmom_dsec[:, :, 4] = -diff[:, :, 1]
        dmom_dsec[:, :, 5] = diff[:, :, 0]

        id1 = 6*(ny-1)*(nx-1)
        partials['loads', 'sec_forces'][:id1] = 0.5

        id2 = id1 * 2
        dmom_dsec = dmom_dsec.flatten()
        partials['loads', 'sec_forces'][id1:id2] = dmom_dsec
        partials['loads', 'sec_forces'][id2:] = dmom_dsec

        # dmoment__dmesh

        dmom_ddiff = np.zeros((nx-1, ny-1, 6))
        dmom_ddiff[:, :, 0] = sec_forces[:, :, 2]
        dmom_ddiff[:, :, 1] = -sec_forces[:, :, 1]
        dmom_ddiff[:, :, 2] = -sec_forces[:, :, 2]
        dmom_ddiff[:, :, 3] = sec_forces[:, :, 0]
        dmom_ddiff[:, :, 4] = sec_forces[:, :, 1]
        dmom_ddiff[:, :, 5] = -sec_forces[:, :, 0]

        dmom_ddiff_sum = np.sum(dmom_ddiff, axis=0)

        dmon_ddiff_diag = np.zeros((nx-1, ny, 6))
        dmon_ddiff_diag[:, 1:, :] = dmom_ddiff
        dmon_ddiff_diag[:, :-1, :] += dmom_ddiff
        dmon_ddiff_diag_sum = np.zeros((1, ny, 6))
        dmon_ddiff_diag_sum[:, :-1, :] = dmom_ddiff_sum
        dmon_ddiff_diag_sum[:, 1:, :] += dmom_ddiff_sum

        dmom_ddiff = dmom_ddiff.flatten()
        dmom_ddiff_sum = dmom_ddiff_sum.flatten()
        dmon_ddiff_diag = dmon_ddiff_diag.flatten()
        dmon_ddiff_diag_sum = dmon_ddiff_diag_sum.flatten()

        idy = 6*(ny-1)
        idx = idy*nx
        idw = idy*(nx-1)

        # Need to zero out what's there because our assignments overlap.
        partials['loads','def_mesh'][:] = 0.0

        # Upper diagonal blocks
        partials['loads','def_mesh'][:idw] = dmom_ddiff * ((1-w1) * 0.25)
        partials['loads','def_mesh'][idy:idx] += dmom_ddiff * (w1 * 0.25)
        partials['loads','def_mesh'][:idy] -= dmom_ddiff_sum * ((1-w2) * 0.25)
        partials['loads','def_mesh'][idx-idy:idx] -= dmom_ddiff_sum * (w2 * 0.25)

        # Lower Diagonal blocks
        id2 = idx * 2
        partials['loads','def_mesh'][idx:idx+idw] = dmom_ddiff * ((1-w1) * 0.25)
        partials['loads','def_mesh'][idx+idy:id2] += dmom_ddiff * (w1 * 0.25)
        partials['loads','def_mesh'][idx:idx+idy] -= dmom_ddiff_sum * ((1-w2) * 0.25)
        partials['loads','def_mesh'][id2-idy:id2] -= dmom_ddiff_sum * (w2 * 0.25)

        # Central Diagonal blocks
        idy = 6*ny
        idz = 6*(nx-1)
        id3 = id2 + idw + idz
        partials['loads','def_mesh'][id2:id3] = dmon_ddiff_diag * ((1-w1) * 0.25)
        partials['loads','def_mesh'][id2:id2+idy] -= dmon_ddiff_diag_sum * ((1-w2) * 0.25)

        id2 += idy
        id3 += idy
        partials['loads','def_mesh'][id2:id3] += dmon_ddiff_diag * (w1 * 0.25)
        partials['loads','def_mesh'][id3-idy:id3] -= dmon_ddiff_diag_sum * (w2 * 0.25)
