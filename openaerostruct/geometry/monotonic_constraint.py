from __future__ import print_function, division
import numpy as np

from openmdao.api import ExplicitComponent


class MonotonicConstraint(ExplicitComponent):
    """
    Produce a constraint that is violated if a user-chosen measure on the
    wing does not decrease monotonically from the root to the tip.

    Parameters
    ----------
    var_name : string
        The variable to which the user would like to apply the monotonic constraint.

    Returns
    -------
    monotonic[ny-1] : numpy array
        Values are greater than 0 if the constraint is violated.

    """

    def initialize(self):
        self.options.declare('var_name', types=str)
        self.options.declare('surface', types=dict)

    def setup(self):

        self.surface = surface = self.options['surface']
        self.var_name = self.options['var_name']
        self.con_name = 'monotonic_' + self.var_name

        self.symmetry = surface['symmetry']
        self.ny = surface['num_y']

        self.add_input(self.var_name, val=np.zeros(self.ny))
        self.add_output(self.con_name, val=np.zeros(self.ny-1))

        self.declare_partials('*', '*')

    def compute(self, inputs, outputs):
        # Compute the difference between adjacent variable values
        diff = inputs[self.var_name][:-1] - inputs[self.var_name][1:]
        if self.symmetry:
            outputs[self.con_name] = diff
        else:
            ny2 = (self.ny - 1) // 2
            outputs[self.con_name][:ny2] = diff[:ny2]
            outputs[self.con_name][ny2:] = -diff[ny2:]

    def compute_partials(self, inputs, partials):
        np.fill_diagonal(partials[self.con_name, self.var_name][:, :], 1)
        np.fill_diagonal(partials[self.con_name, self.var_name][:, 1:], -1)

        if not self.symmetry:
            ny2 = (self.ny - 1) // 2
            partials[self.con_name, self.var_name][ny2:, :] *= -1
