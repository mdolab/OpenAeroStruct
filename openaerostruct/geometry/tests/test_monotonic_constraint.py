import unittest

import numpy as np

from openaerostruct.geometry.monotonic_constraint import MonotonicConstraint

from openaerostruct.utils.testing import run_test


class Test(unittest.TestCase):

    def test_sym1(self):
        surface = {'symmetry' : True,
                   'mesh' : np.zeros((1,5,3))}

        comp = MonotonicConstraint(var_name='x', surface=surface)
        run_test(self, comp)

    def test_sym2(self):
        surface = {'symmetry' : True,
                   'mesh' : np.zeros((1,16,3))}
        comp = MonotonicConstraint(var_name='x', surface=surface)

        run_test(self, comp)

    def test_assym1(self):
        surface = {'symmetry' : False,
                   'mesh' : np.zeros((1,5,3))}
        comp = MonotonicConstraint(var_name='x', surface=surface)
        run_test(self, comp)

    def test_assym2(self):
        surface = {'symmetry' : False,
                   'mesh' : np.zeros((1,16,3))}
        comp = MonotonicConstraint(var_name='x', surface=surface)
        run_test(self, comp)

if __name__ == '__main__':
    unittest.main()
