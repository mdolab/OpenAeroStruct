import openmdao.api as om


class Coeffs(om.ExplicitComponent):
    """ Compute lift and drag coefficients for each individual lifting surface.

    Parameters
    ----------
    S_ref : float
        The reference areas of the lifting surface.
    L : float
        Total lift for the lifting surface.
    D : float
        Total drag for the lifting surface.
    v : float
        Freestream air velocity in m/s.
    rho : float
        Air density in kg/m^3.

    Returns
    -------
    CL1 : float
        Induced coefficient of lift (CL) for the lifting surface.
    CDi : float
        Induced coefficient of drag (CD) for the lifting surface.
    """

    def setup(self):
        self.add_input('S_ref', val=1., units='m**2')
        self.add_input('L', val=1., units='N')
        self.add_input('D', val=1., units='N')
        self.add_input('v', val=1., units='m/s')
        self.add_input('rho', val=1., units='kg/m**3')

        self.add_output('CL1', val=0.)
        self.add_output('CDi', val=0.)

        self.declare_partials('CL1', 'L')
        self.declare_partials('CDi', 'D')
        self.declare_partials('CL1', 'v')
        self.declare_partials('CDi', 'v')
        self.declare_partials('CL1', 'rho')
        self.declare_partials('CDi', 'rho')
        self.declare_partials('CL1', 'S_ref')
        self.declare_partials('CDi', 'S_ref')

    def compute(self, inputs, outputs):
        S_ref = inputs['S_ref']
        rho = inputs['rho']
        v = inputs['v']
        L = inputs['L']
        D = inputs['D']

        outputs['CL1'] = L / (0.5 * rho * v**2 * S_ref)
        outputs['CDi'] = D / (0.5 * rho * v**2 * S_ref)

    def compute_partials(self, inputs, partials):
        S_ref = inputs['S_ref']
        rho = inputs['rho']
        v = inputs['v']
        L = inputs['L']
        D = inputs['D']

        partials['CL1', 'L'] = 1. / (0.5 * rho * v**2 * S_ref)
        partials['CDi', 'D'] = 1. / (0.5 * rho * v**2 * S_ref)

        partials['CL1', 'v'] = -2. * L / (0.5 * rho * v**3 * S_ref)
        partials['CDi', 'v'] = -2. * D / (0.5 * rho * v**3 * S_ref)

        partials['CL1', 'rho'] = -L / (0.5 * rho**2 * v**2 * S_ref)
        partials['CDi', 'rho'] = -D / (0.5 * rho**2 * v**2 * S_ref)

        partials['CL1', 'S_ref'] = -L / (0.5 * rho * v**2 * S_ref**2)
        partials['CDi', 'S_ref'] = -D / (0.5 * rho * v**2 * S_ref**2)
