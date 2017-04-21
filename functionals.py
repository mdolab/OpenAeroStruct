from __future__ import division, print_function
import numpy as np

from openmdao.api import Component

class FunctionalBreguetRange(Component):
    """ Computes the fuel burn using the Breguet range equation using
        the computed CL, CD, weight, and provided specific fuel consumption, speed of sound,
        Mach number, initial weight, and range.

        Note that we add information from each lifting surface.

    Parameters
    ----------
    CL : float
        Total coefficient of lift (CL) for the lifting surface.
    CD : float
        Total coefficient of drag (CD) for the lifting surface.
    weight : float
        Total weight of the structural spar.
    Returns
    -------
    fuelburn : float
        Computed fuel burn in kg based on the Breguet range equation.
    """

    def __init__(self, surfaces, prob_dict):
        super(FunctionalBreguetRange, self).__init__()

        self.surfaces = surfaces
        self.prob_dict = prob_dict

        for surface in surfaces:
            name = surface['name']

            self.add_param(name+'CL', val=0.)
            self.add_param(name+'CD', val=0.)
            self.add_param(name+'structural_weight', val=0.)

        self.add_output('fuelburn', val=0.)

        self.deriv_options['type'] = 'cs'
        self.deriv_options['form'] = 'central'

    def solve_nonlinear(self, params, unknowns, resids):
        CT = self.prob_dict['CT']
        a = self.prob_dict['a']
        R = self.prob_dict['R']
        M = self.prob_dict['M']
        W0 = self.prob_dict['W0'] * self.prob_dict['g']
        fuelburn = 0.

        for surface in self.surfaces:
            name = surface['name']

            CL = params[name+'CL']
            CD = params[name+'CD']
            Ws = params[name+'structural_weight']

            fuelburn += np.sum((W0 + Ws) * (np.exp(R * CT / a / M * CD / CL) - 1))

        # Convert fuelburn from N to kg
        unknowns['fuelburn'] = fuelburn / self.prob_dict['g']
        if unknowns['fuelburn'].imag == 0.:
            print(unknowns['fuelburn'])

class FunctionalEquilibrium(Component):
    """ Lift = weight constraint.

    Note that we add information from each lifting surface.

    Parameters
    ----------
    L : float
        Total lift for the lifting surface.
    weight : float
        Total weight of the structural spar.
    fuelburn : float
        Computed fuel burn in kg based on the Breguet range equation.

    Returns
    -------
    eq_con : float
        Equality constraint for L=W. eq_con = 0 for the constraint to be satisfied.
    """

    def __init__(self, surfaces, prob_dict):
        super(FunctionalEquilibrium, self).__init__()

        self.surfaces = surfaces
        self.prob_dict = prob_dict

        for surface in surfaces:
            name = surface['name']

            self.add_param(name+'L', val=0.)
            self.add_param(name+'structural_weight', val=0.)

        self.add_param('fuelburn', val=0.)
        self.add_output('eq_con', val=0.)

        self.deriv_options['type'] = 'cs'
        self.deriv_options['form'] = 'central'

    def solve_nonlinear(self, params, unknowns, resids):
        structural_weight = 0.
        L = 0.
        W0 = 0.
        for surface in self.surfaces:
            name = surface['name']
            structural_weight += params[name+'structural_weight']
            L += params[name+'L']
            W0 += (self.prob_dict['W0'] * self.prob_dict['g'])

        unknowns['eq_con'] = (structural_weight + params['fuelburn'] * self.prob_dict['g'] + W0 - L) / W0
