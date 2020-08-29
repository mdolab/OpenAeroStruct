import numpy as np

import openmdao.api as om


class LiftDrag(om.ExplicitComponent):
    """
    Calculate total lift and drag in force units based on section forces.
    This is for one given lifting surface.

    parameters
    ----------
    sec_forces[nx-1, ny-1, 3] : numpy array
        Contains the sectional forces acting on each panel.
    alpha : float
        Angle of attack in degrees.
    beta : float
        Sideslip angle in degrees.

    Returns
    -------
    L : float
        Total induced lift force for the lifting surface.
    D : float
        Total induced drag force for the lifting surface.

    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        self.surface = surface = self.options['surface']

        self.nx = nx = surface['mesh'].shape[0]
        self.ny = ny = surface['mesh'].shape[1]
        self.num_panels = (nx - 1) * (ny - 1)

        self.add_input('sec_forces', val=np.ones((nx-1, ny-1, 3)), units='N')
        self.add_input('alpha', val=3., units='deg')
        self.add_input('beta', val=0., units='deg')
        self.add_output('L', val=0., units='N')
        self.add_output('D', val=0., units='N')

        self.declare_partials('*', ['sec_forces', 'alpha', 'L', 'D'])
        self.declare_partials('D', 'beta')

    def compute(self, inputs, outputs):
        alpha = inputs['alpha'] * np.pi / 180.
        beta = inputs['beta'] * np.pi / 180.
        forces = inputs['sec_forces'].reshape(-1, 3)

        cosa = np.cos(alpha)
        sina = np.sin(alpha)
        cosb = np.cos(beta)
        sinb = np.sin(beta)

        # Compute the induced lift force on each lifting surface
        outputs['L'] = np.sum(-forces[:, 0] * sina + forces[:, 2] * cosa)

        # Compute the induced drag force on each lifting surface
        outputs['D'] = np.sum( forces[:, 0] * cosa * cosb \
                               - forces[:, 1] * sinb \
                               + forces[:, 2] * sina * cosb)

        if self.surface['symmetry']:
            outputs['D'] *= 2.
            outputs['L'] *= 2.

    def compute_partials(self, inputs, partials):
        """ Jacobian for lift and drag."""

        # Analytic derivatives for sec_forces
        p180 = np.pi / 180.
        alpha = inputs['alpha'][0] * p180
        beta = inputs['beta'][0] * p180

        cosa = np.cos(alpha)
        sina = np.sin(alpha)
        cosb = np.cos(beta)
        sinb = np.sin(beta)

        forces = inputs['sec_forces']

        if self.surface['symmetry']:
            symmetry_factor = 2.
        else:
            symmetry_factor = 1.

        tmp = np.array([-sina, 0, cosa])
        partials['L', 'sec_forces'] = \
            np.atleast_2d(np.tile(tmp, self.num_panels)) * symmetry_factor

        tmp = np.array([cosa * cosb, -sinb, sina * cosb])
        partials['D', 'sec_forces'] = \
            np.atleast_2d(np.tile(tmp, self.num_panels)) * symmetry_factor

        partials['L', 'alpha'] = p180 * symmetry_factor * \
            np.sum(-forces[:, :, 0] * cosa - forces[:, :, 2] * sina)
        partials['D', 'alpha'] = p180 * symmetry_factor * \
            np.sum(-forces[:, :, 0] * sina * cosb + \
                   forces[:, :, 2] * cosa * cosb)

        partials['D', 'beta'] = p180 * symmetry_factor * \
            np.sum(-forces[:, :, 0] * cosa * sinb + \
                   -forces[:, :, 1] * cosb + \
                   -forces[:, :, 2] * sina * sinb)