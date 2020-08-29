import openmdao.api as om

class TotalDrag(om.ExplicitComponent):
    """ Calculate total drag in force units.

    parameters
    ----------
    CDi : float
        Induced coefficient of drag (CD) for the lifting surface.
    CDv : float
        Calculated coefficient of viscous drag for the lifting surface.

    Returns
    -------
    CD : float
        Total coefficient of drag (CD) for the lifting surface.
    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        surface = self.options['surface']

        self.add_input('CDi', val=1.)
        self.add_input('CDv', val=1.)
        self.add_input('CDw', val=1.)

        self.add_output('CD', val=1.)

        self.CD0 = surface['CD0']

        self.declare_partials('CD', 'CDi', val=1.)
        self.declare_partials('CD', 'CDv', val=1.)
        self.declare_partials('CD', 'CDw', val=1.)

    def compute(self, inputs, outputs):
        outputs['CD'] = inputs['CDi'] + inputs['CDv'] + inputs['CDw'] + self.CD0
