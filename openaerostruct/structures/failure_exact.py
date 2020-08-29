import numpy as np

import openmdao.api as om


class FailureExact(om.ExplicitComponent):
    """
    Output individual failure constraints on each FEM element.

    Parameters
    ----------
    vonmises[ny-1, 2] : numpy array
        von Mises stress magnitudes for each FEM element.

    Returns
    -------
    failure[ny-1, 2] : numpy array
        Array of failure conditions. Positive if element has failed.

    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        surface = self.options['surface']

        if surface['fem_model_type'] == 'tube':
            num_failure_criteria = 2
        elif surface['fem_model_type'] == 'wingbox':
            num_failure_criteria = 4

        self.ny = surface['mesh'].shape[1]
        self.sigma = surface['yield']

        self.add_input('vonmises', val=np.zeros((self.ny-1, num_failure_criteria)), units='N/m**2')
        self.add_output('failure', val=np.zeros((self.ny-1, num_failure_criteria)))

        self.declare_partials('failure', 'vonmises', val=np.eye(((self.ny-1)*num_failure_criteria)) / self.sigma)

    def compute(self, inputs, outputs):
        outputs['failure'] = inputs['vonmises'] / self.sigma - 1
