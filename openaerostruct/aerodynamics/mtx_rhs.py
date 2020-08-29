import numpy as np

import openmdao.api as om


class VLMMtxRHSComp(om.ExplicitComponent):
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
    vel_mtx[num_eval_points, nx - 1, ny - 1, 3] : numpy array
        The AIC matrix for the all lifting surfaces representing the aircraft.
        This has some sparsity pattern, but it is more dense than the FEM matrix
        and the entries have a wide range of magnitudes. One exists for each
        combination of surface name and evaluation points name.
    normals[nx-1, ny-1, 3] : numpy array
        The normal vector for each panel, computed as the cross of the two
        diagonals from the mesh points.

    Returns
    -------
    mtx[system_size, system_size] : numpy array
        Final fully assembled AIC matrix that is used to solve for the
        circulations.
    rhs[system_size] : numpy array
        Right-hand side of the AIC linear system, constructed from the
        freestream velocities and panel normals.
    """

    def initialize(self):
        self.options.declare('surfaces', types=list)

    def setup(self):
        surfaces = self.options['surfaces']

        system_size = 0

        # Loop through the surfaces to compute the total number of panels;
        # the system_size
        for surface in surfaces:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            system_size += (nx - 1) * (ny - 1)

        self.system_size = system_size

        self.add_input('freestream_velocities', shape=(system_size, 3), units='m/s')
        self.add_output('mtx', shape=(system_size, system_size), units='1/m')
        self.add_output('rhs', shape=system_size, units='m/s')

        # Set up indicies arrays for sparse Jacobians
        vel_indices = np.arange(system_size * 3).reshape((system_size, 3))
        mtx_indices = np.arange(system_size * system_size).reshape((system_size, system_size))
        rhs_indices = np.arange(system_size)

        self.declare_partials('rhs', 'freestream_velocities',
            rows=np.einsum('i,j->ij', rhs_indices, np.ones(3, int)).flatten(),
            cols=vel_indices.flatten()
        )

        ind_1 = 0
        ind_2 = 0

        # Loop through each surface to add inputs and set up derivatives.
        # We keep track of the surface's indices within the total system's
        # indices to access the matrix in the correct locations for the derivs.
        # This is because the AIC linear system has information for all surfaces
        # together.
        for surface in surfaces:
            mesh=surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']
            num = (nx - 1) * (ny - 1)

            ind_2 += num

            # Get the correct names for each vel_mtx and normals, then
            # add them to the component
            vel_mtx_name = '{}_{}_vel_mtx'.format(name, 'coll_pts')
            normals_name = '{}_normals'.format(name)

            self.add_input(vel_mtx_name,
                shape=(system_size, nx - 1, ny - 1, 3), units='1/m')
            self.add_input(normals_name, shape=(nx - 1, ny - 1, 3))

            velocities_indices = np.arange(system_size * num * 3).reshape(
                (system_size, nx - 1, ny - 1, 3)
            )
            normals_indices = np.arange(num * 3).reshape((num, 3))

            # Declare each set of partials based on the indices, ind_1 and ind_2
            self.declare_partials('mtx', vel_mtx_name,
                rows=np.einsum('ij,k->ijk', mtx_indices[:, ind_1:ind_2], np.ones(3, int)).flatten(),
                cols=velocities_indices.flatten(),
            )
            self.declare_partials('mtx', normals_name,
                rows=np.einsum('ij,k->ijk', mtx_indices[ind_1:ind_2, :], np.ones(3, int)).flatten(),
                cols=np.einsum('ik,j->ijk', normals_indices, np.ones(system_size, int)).flatten(),
            )
            self.declare_partials('rhs', normals_name,
                rows=np.outer(rhs_indices[ind_1:ind_2], np.ones(3, int)).flatten(),
                cols=normals_indices.flatten(),
            )

            ind_1 += num

        self.mtx_n_n_3 = np.zeros((system_size, system_size, 3))
        self.normals_n_3 = np.zeros((system_size, 3))
        self.set_check_partial_options(wrt='*', method='fd', step=1e-5)

    def compute(self, inputs, outputs):
        surfaces = self.options['surfaces']

        system_size = self.system_size

        ind_1 = 0
        ind_2 = 0
        for surface in surfaces:
            nx = surface['mesh'].shape[0]
            ny = surface['mesh'].shape[1]
            name = surface['name']
            num = (nx - 1) * (ny - 1)

            ind_2 += num

            vel_mtx_name = '{}_{}_vel_mtx'.format(name, 'coll_pts')
            normals_name = '{}_normals'.format(name)

            # Construct the full matrix and all of the lifting surfaces
            # together
            self.mtx_n_n_3[:, ind_1:ind_2, :] = inputs[vel_mtx_name].reshape((system_size, num, 3))
            self.normals_n_3[ind_1:ind_2, :] = inputs[normals_name].reshape((num, 3))

            ind_1 += num

        # Actually obtain the final matrix by multiplying through with the
        # normals. Also create the rhs based on v dot n.
        outputs['mtx'] = np.einsum('ijk,ik->ij', self.mtx_n_n_3, self.normals_n_3)
        outputs['rhs'] = -np.einsum('ij,ij->i', inputs['freestream_velocities'], self.normals_n_3)

    def compute_partials(self, inputs, partials):
        surfaces = self.options['surfaces']

        system_size = self.system_size

        ind_1 = 0
        ind_2 = 0
        for surface in surfaces:
            nx = surface['mesh'].shape[0]
            ny = surface['mesh'].shape[1]
            name = surface['name']
            num = (nx - 1) * (ny - 1)

            ind_2 += num

            vel_mtx_name = '{}_{}_vel_mtx'.format(name, 'coll_pts')
            normals_name = '{}_normals'.format(name)

            partials['mtx', vel_mtx_name] = np.einsum('ijk,ik->ijk',
                np.ones((system_size, num, 3)),
                self.normals_n_3,
            ).flatten()

            partials['mtx', normals_name] = self.mtx_n_n_3[ind_1:ind_2, :, :].flatten()

            partials['rhs', normals_name] = -inputs['freestream_velocities'][ind_1:ind_2, :].flatten()

            ind_1 += num

        partials['rhs', 'freestream_velocities'] = -self.normals_n_3.flatten()
