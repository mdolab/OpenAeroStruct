import numpy as np

import openmdao.api as om


class ConvertVelocity(om.ExplicitComponent):
    """
    Convert the freestream velocity magnitude into a velocity vector at each
    evaluation point. In this case, each of the panels sees the same velocity.
    This really just helps us set up the velocities for use in the VLM analysis.

    Parameters
    ----------
    alpha : float
        The angle of attack for the aircraft (all lifting surfaces) in degrees.
    beta : float
        The sideslip angle for the aircraft (all lifting surfaces) in degrees.
    v : float
        The freestream velocity magnitude.
    rotational_velocities[system_size, 3] : numpy array
        The rotated freestream velocities at each evaluation point for all
        lifting surfaces. system_size is the sum of the count of all panels
        for all lifting surfaces.

    Returns
    -------
    freestream_velocities[system_size, 3] : numpy array
        The rotated freestream velocities at each evaluation point for all
        lifting surfaces. system_size is the sum of the count of all panels
        for all lifting surfaces.
    """

    def initialize(self):
        self.options.declare('surfaces', types=list)
        self.options.declare('rotational', False, types=bool,
                             desc="Set to True to turn on support for computing angular velocities")

    def setup(self):
        surfaces = self.options['surfaces']
        rotational = self.options['rotational']

        system_size = 0
        sizes = []

        # Loop through each surface and cumulatively add the number of panels
        # to obtain system_size.
        for surface in surfaces:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            size = (nx - 1) * (ny - 1)
            system_size += size
            sizes.append(size)

        self.system_size = system_size

        self.add_input('alpha', val=0., units='deg')
        self.add_input('beta', val=0., units='deg')
        self.add_input('v', val=1., units='m/s')

        if rotational:
            self.add_input('rotational_velocities', shape=(system_size, 3), units='m/s')

        self.add_output('freestream_velocities', shape=(system_size, 3), units='m/s')

        self.declare_partials('freestream_velocities', 'alpha')
        self.declare_partials('freestream_velocities', 'beta')
        self.declare_partials('freestream_velocities', 'v')

        if rotational:
            nn = 3 * system_size
            row_col = np.arange(nn)
            val = np.ones((nn, ))
            self.declare_partials('freestream_velocities', 'rotational_velocities',
                                  rows=row_col, cols=row_col, val=val)

    def compute(self, inputs, outputs):
        # Rotate the freestream velocities based on the angle of attack and the sideslip angle.
        alpha = inputs['alpha'][0] * np.pi / 180.
        beta = inputs['beta'][0] * np.pi / 180.

        cosa = np.cos(alpha)
        sina = np.sin(alpha)
        cosb = np.cos(beta)
        sinb = np.sin(beta)

        v_inf = inputs['v'][0] * np.array([cosa*cosb, -sinb, sina*cosb])
        outputs['freestream_velocities'][:, :] = v_inf

        if self.options['rotational']:
            outputs['freestream_velocities'][:, :] += inputs['rotational_velocities']

    def compute_partials(self, inputs, J):
        alpha = inputs['alpha'][0] * np.pi / 180.
        beta = inputs['beta'][0] * np.pi / 180.

        cosa = np.cos(alpha)
        sina = np.sin(alpha)
        cosb = np.cos(beta)
        sinb = np.sin(beta)

        J['freestream_velocities','v'] = np.tile(np.array([cosa*cosb, -sinb, sina*cosb]), self.system_size)
        J['freestream_velocities','alpha'] = np.tile(inputs['v'][0] * np.array([-sina*cosb, 0., cosa*cosb]) * np.pi/180.,
                                                     self.system_size)
        J['freestream_velocities','beta'] = np.tile(inputs['v'][0] * np.array([-cosa*sinb, -cosb, -sina*sinb]) * np.pi/180.,
                                                    self.system_size)

