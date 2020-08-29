import numpy as np

import openmdao.api as om

from openaerostruct.utils.vector_algebra import compute_cross, compute_cross_deriv1, compute_cross_deriv2


class PanelForces(om.ExplicitComponent):
    """
    Compute the panel forces acting on all surfaces in the system.

    Parameters
    ----------
    rho : float
        Air density at the flight condition.
    horseshoe_circulations[system_size] : numpy array
        The equivalent horseshoe circulations obtained by intelligently summing
        the vortex ring circulations, accounting for overlaps between rings.
    bound_vecs[system_size, 3] : numpy array
        The vectors representing the bound vortices for each panel in the
        problem.
        This array contains points for all lifting surfaces in the problem.
    force_pts_velocities[system_size, 3] : numpy array
        The actual velocities experienced at the evaluation points for each
        lifting surface in the system. This is the summation of the freestream
        velocities and the induced velocities caused by the circulations.

    Returns
    -------
    panel_forces[system_size, 3] : numpy array
        All of the forces acting on all panels in the total system.
    """

    def initialize(self):
        self.options.declare('surfaces', types=list)

    def setup(self):
        surfaces = self.options['surfaces']

        system_size = 0

        for surface in surfaces:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]

            system_size += (nx - 1) * (ny - 1)

        self.system_size = system_size

        self.add_input('rho', units='kg/m**3')
        self.add_input('horseshoe_circulations', shape=system_size, units='m**2/s')
        self.add_input('force_pts_velocities', shape=(system_size, 3), units='m/s')
        self.add_input('bound_vecs', shape=(system_size, 3), units='m')

        self.add_output('panel_forces', shape=(system_size, 3), units='N')

        # Set up all the sparse Jacobians
        self.declare_partials('panel_forces', 'rho',
            rows=np.arange(3 * system_size),
            cols=np.zeros(3 * system_size, int),
        )
        self.declare_partials('panel_forces', 'horseshoe_circulations',
            rows=np.arange(3 * system_size),
            cols=np.outer(np.arange(system_size), np.ones(3, int)).flatten(),
        )
        self.declare_partials('panel_forces', 'force_pts_velocities',
            rows=np.einsum('ij,k->ijk',
                np.arange(3 * system_size).reshape((system_size, 3)),
                np.ones(3, int),
            ).flatten(),
            cols=np.einsum('ik,j->ijk',
                np.arange(3 * system_size).reshape((system_size, 3)),
                np.ones(3, int),
            ).flatten(),
        )
        self.declare_partials('panel_forces', 'bound_vecs',
            rows=np.einsum('ij,k->ijk',
                np.arange(3 * system_size).reshape((system_size, 3)),
                np.ones(3, int),
            ).flatten(),
            cols=np.einsum('ik,j->ijk',
                np.arange(3 * system_size).reshape((system_size, 3)),
                np.ones(3, int),
            ).flatten(),
        )

    def compute(self, inputs, outputs):
        rho = inputs['rho'][0]
        horseshoe_circulations = np.outer(inputs['horseshoe_circulations'], np.ones(3))
        velocities = inputs['force_pts_velocities']
        bound_vecs = inputs['bound_vecs']

        # Actually compute the forces by taking the cross of velocities acting
        # at the force points with the bound vortex filament vector.
        outputs['panel_forces'] = \
            rho * horseshoe_circulations * compute_cross(velocities, bound_vecs)

    def compute_partials(self, inputs, partials):
        rho = inputs['rho'][0]
        horseshoe_circulations = np.outer(inputs['horseshoe_circulations'], np.ones(3))
        velocities = inputs['force_pts_velocities']
        bound_vecs = inputs['bound_vecs']

        horseshoe_circulations_ones = np.einsum('i,jk->ijk', inputs['horseshoe_circulations'], np.ones((3, 3)))

        deriv_array = np.einsum('i,jk->ijk',
            np.ones(self.system_size),
            np.eye(3))

        partials['panel_forces', 'rho'] = \
            (horseshoe_circulations * compute_cross(velocities, bound_vecs)).flatten()
        partials['panel_forces', 'horseshoe_circulations'] = \
            (rho * compute_cross(velocities, bound_vecs)).flatten()
        partials['panel_forces', 'force_pts_velocities'] = \
            (rho * horseshoe_circulations_ones * compute_cross_deriv1(deriv_array, bound_vecs)).flatten()
        partials['panel_forces', 'bound_vecs'] = \
            (rho * horseshoe_circulations_ones * compute_cross_deriv2(velocities, deriv_array)).flatten()
