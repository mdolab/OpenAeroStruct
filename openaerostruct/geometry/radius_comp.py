""" Manipulate geometry mesh based on high-level design parameters. """

from __future__ import division, print_function
import numpy as np

from openmdao.api import ExplicitComponent
from openaerostruct.structures.utils import radii

try:
    from openaerostruct.fortran import OAS_API
    fortran_flag = True
    data_type = float
except:
    fortran_flag = False
    data_type = complex

def view_mat(mat):
    """ Helper function used to visually examine matrices. """
    import matplotlib.pyplot as plt
    if len(mat.shape) > 2:
        mat = np.sum(mat, axis=2)
    im = plt.imshow(mat.real, interpolation='none')
    plt.colorbar(im, orientation='horizontal')
    plt.show()

class RadiusComp(ExplicitComponent):
    """

    Parameters
    ----------

    Returns
    -------
    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        surface = self.options['surface']

        self.nx, self.ny = surface['num_x'], surface['num_y']
        self.add_input('mesh', val=np.zeros((self.nx, self.ny, 3)), units='m')
        self.add_output('radius', val=np.zeros((self.ny - 1)), units='m')

        self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        outputs['radius'] = radii(inputs['mesh'], self.options['surface']['t_over_c'])

    def compute_partials(self, inputs, partials):
        """
        Obtain the radii of the FEM element based on local chord.
        """
        mesh = inputs['mesh']
        t_c = self.options['surface']['t_over_c']
        vectors = mesh[-1, :, :] - mesh[0, :, :]
        chords = np.sqrt(np.sum(vectors**2, axis=1))
        mean_chords = 0.5 * chords[:-1] + 0.5 * chords[1:]
        radii_output = t_c * mean_chords / 2.

        for iy in range(self.ny-1):
            partials['radius', 'mesh'][iy, iy*3:(iy+1)*3] = -vectors[iy, :] / chords[iy] * t_c / 4
            partials['radius', 'mesh'][iy, (iy+1)*3:(iy+2)*3] = -vectors[iy+1, :] / chords[iy+1] * t_c / 4
