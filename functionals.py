from __future__ import division
import numpy

from openmdao.api import Component



class FunctionalBreguetRange(Component):
    """ Computes the fuel burn using the Breguet range equation """

    def __init__(self, W0, CT, a, R, M, aero_ind):
        super(FunctionalBreguetRange, self).__init__()

        n_surf = aero_ind.shape[0]

        self.add_param('CL', val=numpy.zeros((n_surf)))
        self.add_param('CD', val=numpy.zeros((n_surf)))
        self.add_param('weight', val=0.)

        self.add_output('fuelburn', val=0.)

        self.W0 = W0
        self.CT = CT
        self.a = a
        self.R = R
        self.M = M

        self.deriv_options['type'] = 'cs'
        self.deriv_options['form'] = 'central'
        #self.deriv_options['extra_check_partials_form'] = "central"

    def solve_nonlinear(self, params, unknowns, resids):
        W0 = self.W0
        CT = self.CT
        a = self.a
        R = self.R
        M = self.M

        CL = params['CL']
        CD = params['CD']
        Ws = params['weight']

        unknowns['fuelburn'] = numpy.sum((W0 + Ws) * (numpy.exp(R * CT / a / M * CD / CL) - 1))


class FunctionalEquilibrium(Component):
    """ L = W constraint """

    def __init__(self, W0, aero_ind):
        super(FunctionalEquilibrium, self).__init__()

        n_surf = aero_ind.shape[0]

        self.add_param('L', val=numpy.zeros((n_surf)))
        self.add_param('weight', val=1.)
        self.add_param('fuelburn', val=1.)

        self.add_output('eq_con', val=0.)

        self.deriv_options['type'] = 'cs'
        self.deriv_options['form'] = 'central'
        #self.deriv_options['extra_check_partials_form'] = "central"

        self.W0 = W0

    def solve_nonlinear(self, params, unknowns, resids):
        W0 = self.W0

        unknowns['eq_con'] = (params['weight'] + params['fuelburn'] + W0 - numpy.sum(params['L'])) / W0
