import numpy as np

import unittest

import openmdao.api as om

from openaerostruct.geometry.geometry_mesh import GeometryMesh

from openaerostruct.utils.testing import run_test


class Test(unittest.TestCase):

    def test(self):
        surface = {}
        surface['symmetry'] = False
        surface['type'] = 'aero'
        span = 1
        ny = 5
        mesh = np.zeros((2, ny, 3))
        mesh[0, :, 0] = -1.
        mesh[:, :, 1] = np.linspace(0, span, ny)
        surface['mesh'] = mesh
        surface['name'] = 'wing'

        # The way this is currently set up, we don't actually use the values here
        surface['twist_cp'] = np.zeros((5))
        surface['chord_cp'] = np.zeros((5))
        surface['xshear_cp'] = np.zeros((5))
        surface['yshear_cp'] = np.zeros((5))
        surface['zshear_cp'] = np.zeros((5))

        surface['span'] = span
        surface['t_over_c_cp'] = np.array([.12])

        comp = GeometryMesh(surface=surface)

        run_test(self, comp, complex_flag=True, method='cs')


if __name__ == '__main__':
    unittest.main()
