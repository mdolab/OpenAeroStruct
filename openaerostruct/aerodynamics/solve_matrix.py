import numpy as np
from scipy.linalg import lu_factor, lu_solve

import openmdao.api as om


class SolveMatrix(om.ImplicitComponent):
    """
    Solve the AIC linear system to obtain the vortex ring circulations.

    Parameters
    ----------
    mtx[system_size, system_size] : numpy array
        Final fully assembled AIC matrix that is used to solve for the
        circulations.
    rhs[system_size] : numpy array
        Right-hand side of the AIC linear system, constructed from the
        freestream velocities and panel normals.

    Returns
    -------
    circulations[system_size] : numpy array
        The vortex ring circulations obtained by solving the AIC linear system.

    """

    def initialize(self):
        self.options.declare('surfaces', types=list)

    def setup(self):
        system_size = 0

        for surface in self.options['surfaces']:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]

            system_size += (nx - 1) * (ny - 1)

        self.system_size = system_size

        self.add_input('mtx', shape=(system_size, system_size), units='1/m')
        self.add_input('rhs', shape=system_size, units='m/s')
        self.add_output('circulations', shape=system_size, units='m**2/s')

        self.declare_partials('circulations', 'circulations',
            rows=np.outer(np.arange(system_size), np.ones(system_size, int)).flatten(),
            cols=np.outer(np.ones(system_size, int), np.arange(system_size)).flatten(),
        )
        self.declare_partials('circulations', 'mtx',
            rows=np.outer(np.arange(system_size), np.ones(system_size, int)).flatten(),
            cols=np.arange(system_size ** 2),
        )
        self.declare_partials('circulations', 'rhs', val=-1.,
            rows=np.arange(system_size),
            cols=np.arange(system_size),
        )

    def apply_nonlinear(self, inputs, outputs, residuals):
        residuals['circulations'] = inputs['mtx'].dot(outputs['circulations']) - inputs['rhs']

    def solve_nonlinear(self, inputs, outputs):
        self.lu = lu_factor(inputs['mtx'])

        outputs['circulations'] = lu_solve(self.lu, inputs['rhs'])

    def linearize(self, inputs, outputs, partials):
        system_size = self.system_size
        self.lu = lu_factor(inputs['mtx'])

        partials['circulations', 'circulations'] = inputs['mtx'].flatten()
        partials['circulations', 'mtx'] = \
            np.outer(np.ones(system_size), outputs['circulations']).flatten()

    def solve_linear(self, d_outputs, d_residuals, mode):
        if mode == 'fwd':
            d_outputs['circulations'] = lu_solve(self.lu, d_residuals['circulations'], trans=0)
        else:
            d_residuals['circulations'] = lu_solve(self.lu, d_outputs['circulations'], trans=1)
