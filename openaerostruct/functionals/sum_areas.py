import openmdao.api as om

class SumAreas(om.ExplicitComponent):
    """
    Compute the total surface area of the entire aircraft as a sum of its
    individual surfaces' surface areas.

    Parameters
    ----------
    S_ref : float
        Surface area for one lifting surface.

    Returns
    -------
    S_ref_total : float
        Total surface area of the aircraft based on the sum of individual
        surface areas.

    """

    def initialize(self):
        self.options.declare('surfaces', types=list)

    def setup(self):
        for surface in self.options['surfaces']:
            name = surface['name']
            self.add_input(name + '_S_ref', val=1., units='m**2')

        self.add_output('S_ref_total', val=0., units='m**2')

        self.declare_partials('*', '*', val=1.)

    def compute(self, inputs, outputs):
        outputs['S_ref_total'] = 0.
        for surface in self.options['surfaces']:
            name = surface['name']
            S_ref = inputs[name + '_S_ref']
            outputs['S_ref_total'] += S_ref
