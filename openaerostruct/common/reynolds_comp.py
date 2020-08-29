import openmdao.api as om


class ReynoldsComp(om.ExplicitComponent):

    def setup(self):
        self.add_input('rho', val=1., units='slug/ft**3')
        self.add_input('mu', val=1., units='lbf*s/ft**2')
        self.add_input('v', val=1., units='ft/s')

        self.add_output('re', val=1., units='1/ft')

        self.declare_partials('re', ['rho', 'mu', 'v'])

    def compute(self, inputs, outputs):
        outputs['re'] = inputs['rho'] * inputs['v'] / inputs['mu']

    def compute_partials(self, inputs, partials):
        partials['re', 'rho'] = inputs['v'] / inputs['mu']
        partials['re', 'v'] = inputs['rho'] / inputs['mu']
        partials['re', 'mu'] = - inputs['rho'] * inputs['v'] / inputs['mu']**2
