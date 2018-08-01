from __future__ import print_function, division
import numpy as np

from openmdao.api import ExplicitComponent

try:
    from openaerostruct.fortran import OAS_API
    fortran_flag = True
    data_type = float
except:
    fortran_flag = False
    data_type = complex

class Disp(ExplicitComponent):
    """
    Reshape the flattened displacements from the linear system solution into
    a 2D array so we can more easily use the results.

    The solution to the linear system has meaingless entires due to the
    constraints on the FEM model. The displacements from this portion of
    the linear system are not needed, so we select only the relevant
    portion of the displacements for further calculations.

    Parameters
    ----------
    disp_aug[6*(ny+1)] : numpy array
        Augmented displacement array. Obtained by solving the system
        K * disp_aug = forces, where forces is a flattened version of loads.

    Returns
    -------
    disp[6*ny] : numpy array
        Actual displacement array formed by truncating disp_aug.

    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        surface = self.options['surface']

        self.ny = surface['num_y']

        self.add_input('disp_aug', val=np.zeros(((self.ny+1)*6)), units='m')
        self.add_output('disp', val=np.zeros((self.ny, 6)), units='m')

        n = self.ny * 6
        disp_disp_aug = np.zeros((n, n+6))
        disp_disp_aug[:n, :n] = np.eye((n))
        self.declare_partials('disp', 'disp_aug', val=disp_disp_aug)

    def compute(self, inputs, outputs):
        # Obtain the relevant portions of disp_aug and store the reshaped
        # displacements in disp
        outputs['disp'] = inputs['disp_aug'][:-6].reshape((-1, 6))
