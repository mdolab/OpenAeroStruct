import numpy as np

import openmdao.api as om


class EvalVelocities(om.ExplicitComponent):
    """
    Compute the total velocities at each of the evaluation points for every
    panel in the entire system. This is the sum of the freestream and induced
    velocities caused by the circulations.

    Parameters
    ----------
    freestream_velocities[system_size, 3] : numpy array
        The rotated freestream velocities at each evaluation point for all
        lifting surfaces. system_size is the sum of the count of all panels
        for all lifting surfaces.
    circulations[system_size] : numpy array
        The vortex ring circulations obtained from solving the AIC linear
        system.
    vel_mtx[num_eval_points, nx - 1, ny - 1, 3] : numpy array
        The AIC matrix for the all lifting surfaces representing the aircraft.
        This has some sparsity pattern, but it is more dense than the FEM matrix
        and the entries have a wide range of magnitudes. One exists for each
        combination of surface name and evaluation points name.

    Returns
    -------
    velocities[num_eval_points, 3] : numpy array
        The actual velocities experienced at the evaluation points for each
        lifting surface in the system. This is the summation of the freestream
        velocities and the induced velocities caused by the circulations.

    """

    def initialize(self):
        self.options.declare('surfaces', types=list)
        self.options.declare('eval_name', types=str)
        self.options.declare('num_eval_points', types=int)

    def setup(self):
        surfaces = self.options['surfaces']
        eval_name = self.options['eval_name']
        num_eval_points = self.options['num_eval_points']

        system_size = 0

        # Determine system_size by looping through each surface and summing
        # the number of panels.
        for surface in surfaces:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            system_size += (nx - 1) * (ny - 1)

        self.system_size = system_size

        self.add_input('freestream_velocities', shape=(system_size, 3), units='m/s')
        self.add_input('circulations', shape=system_size, units='m**2/s')

        # Get the correct output name; the velocities output depends on which
        # set of evaluation points we use, either collocation or force.
        velocities_name = '{}_velocities'.format(eval_name)
        self.add_output(velocities_name, shape=(num_eval_points, 3), units='m/s')

        # Set up indices to create the sparsity pattern for the derivatives.
        circulations_indices = np.arange(system_size)
        velocities_indices = np.arange(num_eval_points * 3).reshape((num_eval_points, 3))

        self.declare_partials(velocities_name, 'circulations',
            rows=np.einsum('ik,j->ijk',
                velocities_indices, np.ones(system_size, int)).flatten(),
            cols=np.einsum('ik,j->ijk',
                np.ones((num_eval_points, 3), int), circulations_indices).flatten(),
        )

        # These derivatives are linear and don't change so we set the val here
        self.declare_partials(velocities_name, 'freestream_velocities', val=1.,
            rows=np.arange(3 * num_eval_points),
            cols=np.arange(3 * num_eval_points),
        )

        # For each surface we need to correctly set up the sparsity pattern
        # based on the vel_mtx. This is pretty hairy due to the highly
        # dimensional nature of the vel_mtx.
        ind_1 = 0
        ind_2 = 0
        for surface in surfaces:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']
            num = (nx - 1) * (ny - 1)

            ind_2 += num

            vel_mtx_name = '{}_{}_vel_mtx'.format(name, eval_name)

            self.add_input(vel_mtx_name,
                shape=(num_eval_points, nx - 1, ny - 1, 3), units='1/m')

            vel_mtx_indices = np.arange(num_eval_points * num * 3).reshape(
                (num_eval_points, num, 3))

            self.declare_partials(velocities_name, vel_mtx_name,
                rows=np.einsum('ik,j->ijk', velocities_indices, np.ones(num, int)).flatten(),
                cols=vel_mtx_indices.flatten(),
            )

            ind_1 += num

    def compute(self, inputs, outputs):
        surfaces = self.options['surfaces']
        eval_name = self.options['eval_name']
        num_eval_points = self.options['num_eval_points']

        velocities_name = '{}_velocities'.format(eval_name)

        # Start with just the freestream velocities as the base for the output
        # velocities.
        outputs[velocities_name] = inputs['freestream_velocities']

        ind_1 = 0
        ind_2 = 0
        for surface in surfaces:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']
            num = (nx - 1) * (ny - 1)

            ind_2 += num

            vel_mtx_name = '{}_{}_vel_mtx'.format(name, eval_name)

            # Add the induced velocities; the contributions from the
            # circulations.
            outputs[velocities_name] += np.einsum('ijk,j->ik',
                inputs[vel_mtx_name].reshape((num_eval_points, num, 3)),
                inputs['circulations'][ind_1:ind_2],
            )

            ind_1 += num

    def compute_partials(self, inputs, partials):
        surfaces = self.options['surfaces']
        eval_name = self.options['eval_name']
        num_eval_points = self.options['num_eval_points']

        system_size = self.system_size

        velocities_name = '{}_velocities'.format(eval_name)

        dv_dcirc = np.zeros((num_eval_points, system_size, 3))

        ind_1 = 0
        ind_2 = 0
        for surface in surfaces:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']
            num = (nx - 1) * (ny - 1)

            ind_2 += num

            vel_mtx_name = '{}_{}_vel_mtx'.format(name, eval_name)

            partials[velocities_name, vel_mtx_name] = np.einsum('ijk,j->ijk',
                np.ones((num_eval_points, num, 3)),
                inputs['circulations'][ind_1:ind_2],
            ).flatten()

            dv_dcirc[:, ind_1:ind_2, :] = inputs[vel_mtx_name].reshape((num_eval_points, num, 3))

            ind_1 += num

        partials[velocities_name, 'circulations'] = dv_dcirc.flatten()
