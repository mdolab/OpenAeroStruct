import numpy as np

import openmdao.api as om


class Energy(om.ExplicitComponent):
    """ Compute strain energy.

    Parameters
    ----------
    disp[ny, 6] : numpy array
        Actual displacement array formed by truncating disp_aug.
    loads[ny, 6] : numpy array
        Array containing the loads applied on the FEM component,
        computed from the sectional forces.

    Returns
    -------
    energy : float
        Total strain energy of the structural component.

    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        surface = self.options['surface']

        ny = surface['mesh'].shape[1]

        self.add_input('disp', val=np.zeros((ny, 6)), units='m')
        self.add_input('loads', val=np.zeros((ny, 6)), units='N')
        self.add_output('energy', val=0., units='N*m')

        self.declare_partials('*', '*')

    def compute(self, inputs, outputs):
        outputs['energy'] = np.sum(inputs['disp'] * inputs['loads'])

    def compute_partials(self, inputs, partials):
        partials['energy', 'disp'][0, :] = inputs['loads'].real.flatten()
        partials['energy', 'loads'][0, :] = inputs['disp'].real.flatten()
