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
        Array defining the lifting surfaces after deformation.
        Arrays will be flattened in Fortran order (only relevant when more than one chordwise panel).
    sec_forces[nx-1, ny-1, 3] : numpy array
        Array containing the sectional forces acting on each panel.

    Returns
    -------
    loads[ny, 6] : numpy array 
        Array containing the loads applied on the FEM component at each node,
        computed from the sectional forces. The first 3 columns are N, and the last 3 are N*m.
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

            fem_origin = (x_upper[0] * (y_upper[0] - y_lower[0]) + x_upper[-1] * (y_upper[-1] - y_lower[-1])) / (
                (y_upper[0] - y_lower[0]) + (y_upper[-1] - y_lower[-1])
            )

            # For some reason, surface data is complex in some tests.
            self.fem_origin = np.float(fem_origin)

        self.w1 = 0.25
        self.w2 = self.fem_origin

        self.add_input('def_mesh', val=np.zeros((nx, ny, 3)), units='m')
        self.add_input('sec_forces', val=np.zeros((nx - 1, ny - 1, 3)), units='N')

        self.add_output('loads', val=np.zeros((self.ny, 6)), units='N')  ## WARNING!!! UNITS ARE A MIXTURE OF N & N*m
        # Well, technically the units of this load array are mixed.
        # The first 3 indices are N and the last 3 are N*m.

        # Derivatives

        # First, the direct loads wrt sec_forces terms.
        base_row = np.array([0, 1, 2, 6, 7, 8])
        base_col = np.array([0, 1, 2, 0, 1, 2])
        row = np.tile(base_row, ny - 1) + np.repeat(6 * np.arange(ny - 1), 6)
        col = np.tile(base_col, ny - 1) + np.repeat(3 * np.arange(ny - 1), 6)
        rows1 = np.tile(row, nx - 1)
        cols1 = np.tile(col, nx - 1) + np.repeat(3 * (ny - 1) * np.arange(nx - 1), 6 * (ny - 1))

        # Then, the term from the cross product.
        base_row = np.array([3, 3, 4, 4, 5, 5])
        base_col = np.array([1, 2, 0, 2, 0, 1])
        row = np.tile(base_row, ny - 1) + np.repeat(6 * np.arange(ny - 1), 6)
        col = np.tile(base_col, ny - 1) + np.repeat(3 * np.arange(ny - 1), 6)
        row1 = np.tile(row, nx - 1)
        col1 = np.tile(col, nx - 1) + np.repeat(3 * (ny - 1) * np.arange(nx - 1), 6 * (ny - 1))
        rows2 = np.tile(row1, 2) + np.repeat(np.array([0, 6]), 6 * (nx - 1) * (ny - 1))
        cols2 = np.tile(col1, 2)

        rows = np.concatenate([rows1, rows2])
        cols = np.concatenate([cols1, cols2])

        self.declare_partials(of='loads', wrt='sec_forces', rows=rows, cols=cols)

        # Top diagonal is forward-most mesh point.
        base_row = np.array([3, 3, 4, 4, 5, 5])
        base_col = np.array([4, 5, 3, 5, 3, 4])
        row = np.tile(base_row, ny - 1) + np.repeat(6 * np.arange(ny - 1), 6)
        col = np.tile(base_col, ny - 1) + np.repeat(3 * np.arange(ny - 1), 6)
        rows1 = np.tile(row, nx)
        cols1 = np.tile(col, nx) + np.repeat(3 * ny * np.arange(nx), 6 * (ny - 1))

        # Bottom diagonal is backward-most mesh point.
        base_row = np.array([9, 9, 10, 10, 11, 11])
        base_col = np.array([1, 2, 0, 2, 0, 1])
        row = np.tile(base_row, ny - 1) + np.repeat(6 * np.arange(ny - 1), 6)
        col = np.tile(base_col, ny - 1) + np.repeat(3 * np.arange(ny - 1), 6)
        rows2 = np.tile(row, nx)
        cols2 = np.tile(col, nx) + np.repeat(3 * ny * np.arange(nx), 6 * (ny - 1))

        # Central Diagonal blocks
        base_row = np.array([3, 3, 4, 4, 5, 5])
        base_col = np.array([1, 2, 0, 2, 0, 1])
        row = np.tile(base_row, ny) + np.repeat(6 * np.arange(ny), 6)
        col = np.tile(base_col, ny) + np.repeat(3 * np.arange(ny), 6)
        rows3 = np.tile(row, nx)
        cols3 = np.tile(col, nx) + np.repeat(3 * ny * np.arange(nx), 6 * ny)

        rows = np.concatenate([rows1, rows2, rows3])
        cols = np.concatenate([cols1, cols2, cols3])

        self.declare_partials(of='loads', wrt='def_mesh', rows=rows, cols=cols)

        # -------------------------------- Check Partial Options-------------------------------------
        self.set_check_partial_options('*', method='cs', step=1e-40)

    def compute(self, inputs, outputs):
        mesh = inputs['def_mesh']  # [nx, ny, 3]
        sec_forces = inputs['sec_forces']

        # ----- 1. Forces transfer -----
        # Only need to zero out the part that is assigned via +=
        outputs['loads'][-1, :] = 0.0

        # The aero force acting on each panel is evenly transferred to the adjacent FEM nodes.
        sec_forces_sum = 0.5 * np.sum(sec_forces, axis=0)
        outputs['loads'][:-1, :3] = sec_forces_sum
        outputs['loads'][1:, :3] += sec_forces_sum

        # ----- 2. Moments transfer -----
        # Compute the aerodynamic centers at the quarter-chord point of each panel
        # a_pts [nx-1, ny-1, 3]
        a_pts = (
            0.5 * (1 - self.w1) * mesh[:-1, :-1, :]
            + 0.5 * self.w1 * mesh[1:, :-1, :]
            + 0.5 * (1 - self.w1) * mesh[:-1, 1:, :]
            + 0.5 * self.w1 * mesh[1:, 1:, :]
        )

        # Compute the structural nodes based on the fem_origin location (weighted sum of the LE and TE mesh vertices)
        # s_pts [ny, 3]
        s_pts = (1 - self.w2) * mesh[0, :, :] + self.w2 * mesh[-1, :, :]

        # The moment arm is between the aerodynamic centers of each panel and the FEM nodes.
        # Moment contribution of sec_forces (acting on aero center) to the inner/outer adjacent node
        moment_in = np.sum(np.cross(a_pts - s_pts[:-1, :], 0.5 * sec_forces), axis=0)  # [ny-1, 3]
        moment_out = np.sum(np.cross(a_pts - s_pts[1:, :], 0.5 * sec_forces), axis=0)

        # Total moment at each node = sum of moment_in and moment_out, except the edge nodes.s
        outputs['loads'][:-1, 3:] = moment_in
        outputs['loads'][1:, 3:] += moment_out

    def compute_partials(self, inputs, partials):
        mesh = inputs['def_mesh']
        sec_forces = inputs['sec_forces']
        ny = self.ny
        nx = self.nx
        w1 = self.w1
        w2 = self.w2

        # Compute the aerodynamic centers at the quarter-chord point of each panel
        a_pts = (
            0.5 * (1 - w1) * mesh[:-1, :-1, :]
            + 0.5 * w1 * mesh[1:, :-1, :]
            + 0.5 * (1 - w1) * mesh[:-1, 1:, :]
            + 0.5 * w1 * mesh[1:, 1:, :]
        )

        # Compute the structural nodes
        s_pts = (1 - self.w2) * mesh[0, :, :] + self.w2 * mesh[-1, :, :]

        # ----- 1. dmoment__dsec_forces -----
        # Sensitivity of loads (moments) at inner node wrt sec_force
        diff_in = 0.5 * (a_pts - s_pts[:-1, :])  # moment arm from inner node to aero center.

        dmom_dsec_in = np.empty((nx - 1, ny - 1, 6))
        dmom_dsec_in[:, :, 0] = -diff_in[:, :, 2]
        dmom_dsec_in[:, :, 1] = diff_in[:, :, 1]
        dmom_dsec_in[:, :, 2] = diff_in[:, :, 2]
        dmom_dsec_in[:, :, 3] = -diff_in[:, :, 0]
        dmom_dsec_in[:, :, 4] = -diff_in[:, :, 1]
        dmom_dsec_in[:, :, 5] = diff_in[:, :, 0]

        # Repeat for moments at outer node wrt sec_force
        diff_out = 0.5 * (a_pts - s_pts[1:, :])  # moment arm from outer node to aero center.
        dmom_dsec_out = np.empty((nx - 1, ny - 1, 6))
        dmom_dsec_out[:, :, 0] = -diff_out[:, :, 2]
        dmom_dsec_out[:, :, 1] = diff_out[:, :, 1]
        dmom_dsec_out[:, :, 2] = diff_out[:, :, 2]
        dmom_dsec_out[:, :, 3] = -diff_out[:, :, 0]
        dmom_dsec_out[:, :, 4] = -diff_out[:, :, 1]
        dmom_dsec_out[:, :, 5] = diff_out[:, :, 0]

        id1 = 6 * (ny - 1) * (nx - 1)
        partials['loads', 'sec_forces'][:id1] = 0.5

        id2 = id1 * 2
        dmom_dsec_in = dmom_dsec_in.flatten()
        dmom_dsec_out = dmom_dsec_out.flatten()
        partials['loads', 'sec_forces'][id1:id2] = dmom_dsec_in
        partials['loads', 'sec_forces'][id2:] = dmom_dsec_out

        # ----- 2. dmoment__dmesh -----
        # Sensitivity of moments at inner nodes wrt diff_in (upper diagonal)
        dmom_ddiff_in = np.zeros((nx - 1, ny - 1, 6))
        dmom_ddiff_in[:, :, 0] = sec_forces[:, :, 2]
        dmom_ddiff_in[:, :, 1] = -sec_forces[:, :, 1]
        dmom_ddiff_in[:, :, 2] = -sec_forces[:, :, 2]
        dmom_ddiff_in[:, :, 3] = sec_forces[:, :, 0]
        dmom_ddiff_in[:, :, 4] = sec_forces[:, :, 1]
        dmom_ddiff_in[:, :, 5] = -sec_forces[:, :, 0]
        dmom_ddiff_in_sum = np.sum(dmom_ddiff_in, axis=0)

        # Sensitivity of moments at outer nodes wrt diff_out (lower diagonal)
        dmom_ddiff_out = np.zeros((nx - 1, ny - 1, 6))
        dmom_ddiff_out[:, :, 0] = sec_forces[:, :, 2]
        dmom_ddiff_out[:, :, 1] = -sec_forces[:, :, 1]
        dmom_ddiff_out[:, :, 2] = -sec_forces[:, :, 2]
        dmom_ddiff_out[:, :, 3] = sec_forces[:, :, 0]
        dmom_ddiff_out[:, :, 4] = sec_forces[:, :, 1]
        dmom_ddiff_out[:, :, 5] = -sec_forces[:, :, 0]
        dmom_ddiff_out_sum = np.sum(dmom_ddiff_out, axis=0)

        dmon_ddiff_diag = np.zeros((nx - 1, ny, 6))
        dmon_ddiff_diag[:, 1:, :] = dmom_ddiff_out
        dmon_ddiff_diag[:, :-1, :] += dmom_ddiff_in
        dmon_ddiff_diag_sum = np.zeros((1, ny, 6))
        dmon_ddiff_diag_sum[:, :-1, :] = dmom_ddiff_in_sum
        dmon_ddiff_diag_sum[:, 1:, :] += dmom_ddiff_out_sum

        dmom_ddiff_in = dmom_ddiff_in.flatten()
        dmom_ddiff_out = dmom_ddiff_out.flatten()
        dmon_ddiff_diag = dmon_ddiff_diag.flatten()
        dmon_ddiff_diag_sum = dmon_ddiff_diag_sum.flatten()

        idy = 6 * (ny - 1)
        idx = idy * nx
        idw = idy * (nx - 1)

        # Need to zero out what's there because our assignments overlap.
        partials['loads', 'def_mesh'][:] = 0.0

        # Upper diagonal blocks
        partials['loads', 'def_mesh'][:idw] = dmom_ddiff_in * ((1 - w1) * 0.25)
        partials['loads', 'def_mesh'][idy:idx] += dmom_ddiff_in * (w1 * 0.25)

        # Lower Diagonal blocks
        id2 = idx * 2
        partials['loads', 'def_mesh'][idx : idx + idw] = dmom_ddiff_out * ((1 - w1) * 0.25)
        partials['loads', 'def_mesh'][idx + idy : id2] += dmom_ddiff_out * (w1 * 0.25)

        # Central Diagonal blocks
        idy = 6 * ny
        idz = 6 * (nx - 1)
        id3 = id2 + idw + idz
        partials['loads', 'def_mesh'][id2:id3] = dmon_ddiff_diag * ((1 - w1) * 0.25)
        partials['loads', 'def_mesh'][id2 : id2 + idy] -= dmon_ddiff_diag_sum * ((1 - w2) * 0.5)

        id2 += idy
        id3 += idy
        partials['loads', 'def_mesh'][id2:id3] += dmon_ddiff_diag * (w1 * 0.25)
        partials['loads', 'def_mesh'][id3 - idy : id3] -= dmon_ddiff_diag_sum * (w2 * 0.5)
