from __future__ import division, print_function
import numpy as np

from openmdao.api import ExplicitComponent

try:
    from openaerostruct.fortran import OAS_API
    fortran_flag = True
    data_type = float
except:
    fortran_flag = False
    data_type = complex


class LoadTransfer(ExplicitComponent):
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

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        self.surface = surface = self.options['surface']

        self.ny = surface['num_y']
        self.nx = surface['num_x']

        if surface['fem_model_type'] == 'tube':
            self.fem_origin = surface['fem_origin']
        else:
            y_upper = surface['data_y_upper']
            x_upper = surface['data_x_upper']
            y_lower = surface['data_y_lower']

            self.fem_origin = (x_upper[0]  * (y_upper[0]  - y_lower[0]) +
                               x_upper[-1] * (y_upper[-1] - y_lower[-1])) / \
                             ((y_upper[0]  -  y_lower[0]) + (y_upper[-1] - y_lower[-1]))

        self.add_input('def_mesh', val=np.ones((self.nx, self.ny, 3)), units='m')
        self.add_input('sec_forces', val=np.ones((self.nx-1, self.ny-1, 3)), units='N')

        # Well, technically the units of this load array are mixed.
        # The first 3 indices are N and the last 3 are N*m.
        self.add_output('loads', val=np.zeros((self.ny, 6)), units='N')

        self.declare_partials('*', '*')

        if not fortran_flag:
            self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        mesh = inputs['def_mesh'].copy()
        sec_forces = inputs['sec_forces'].copy()

        if fortran_flag:
            loads = OAS_API.oas_api.transferloads(mesh, sec_forces, self.fem_origin)
        else:
            # Compute the aerodynamic centers at the quarter-chord point of each panel
            w = 0.25
            a_pts = 0.5 * (1-w) * mesh[:-1, :-1, :] + \
                    0.5 *   w   * mesh[1:, :-1, :] + \
                    0.5 * (1-w) * mesh[:-1,  1:, :] + \
                    0.5 *   w   * mesh[1:,  1:, :]

            # Compute the structural midpoints based on the fem_origin location
            w = self.fem_origin
            s_pts = 0.5 * (1-w) * mesh[0, :-1, :] + \
                    0.5 *   w   * mesh[-1, :-1, :] + \
                    0.5 * (1-w) * mesh[0,  1:, :] + \
                    0.5 *   w   * mesh[-1,  1:, :]

            # Find the moment arm between the aerodynamic centers of each panel
            # and the FEM elements
            diff = a_pts - s_pts
            moment = np.zeros((self.ny - 1, 3))
            for ind in range(self.nx-1):
                moment += np.cross(diff[ind, :, :], sec_forces[ind, :, :], axis=1)

            # Compute the loads based on the xyz forces and the computed moments
            loads = np.zeros((self.ny, 6))
            sec_forces_sum = np.sum(sec_forces, axis=0)
            loads[:-1, :3] += 0.5 * sec_forces_sum[:, :]
            loads[ 1:, :3] += 0.5 * sec_forces_sum[:, :]
            loads[:-1, 3:] += 0.5 * moment
            loads[ 1:, 3:] += 0.5 * moment

        outputs['loads'] = loads

    if fortran_flag:
        def compute_partials(self, inputs, partials):

            self.compute(inputs, {})

            ny = self.surface['num_y']
            d_loads = np.zeros((ny, 6))

            def_mesh = inputs['def_mesh']
            sec_forces = inputs['sec_forces']

            for j, val in enumerate(np.array(d_loads).flatten()):
                d_out_b = np.array(d_loads).flatten()
                d_out_b[:] = 0.
                d_out_b[j] = 1.
                d_loads = d_out_b.reshape(d_loads.shape)

                d_def_mesh, d_sec_forces, _ = OAS_API.oas_api.transferloads_b(def_mesh, sec_forces, self.fem_origin, d_loads)

                partials['loads', 'def_mesh'][j, :] = d_def_mesh.flatten()
                partials['loads', 'sec_forces'][j, :] = d_sec_forces.flatten()
