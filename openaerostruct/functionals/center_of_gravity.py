import numpy as np

import openmdao.api as om
from openaerostruct.utils.constants import grav_constant


class CenterOfGravity(om.ExplicitComponent):
    """
    Compute the center of gravity of the entire aircraft based on the inputted W0
    and its corresponding cg and the weighted sum of each surface's structural
    weight and location. We assume the fuel mass is acting at the cg point.

    Note that we add information from each lifting surface.

    Parameters
    ----------
    structural_mass : float
        Total weight of the structural spar for a given surface.
    cg_location[3] : numpy array
        Location of the structural spar's cg for a given surface.

    total_weight : float
        Total weight of the entire aircraft, including W0, all structural weights,
        and fuel.
    fuelburn : float
        Computed fuel burn in kg based on the Breguet range equation.
    W0 : float
        The operating empty weight of the aircraft, without fuel or structural
        mass. Supplied in kg despite being a 'weight' due to convention.
    load_factor : float
        Multiplicative factor on gravity. 1.0 is normal flight; 2.5 would be
        for a 2.5g manuever.
    empty_cg[3] : numpy array
        The location of the cg of the empty aircraft, without considering
        the structural spar or fuel mass's contribution to the cg location.

    Returns
    -------
    cg[3] : numpy array
        The x, y, z coordinates of the center of gravity for the entire aircraft.

    """

    def initialize(self):
        self.options.declare('surfaces', types=list)

    def setup(self):

        arange = np.arange(3)
        for surface in self.options['surfaces']:
            name = surface['name']
            self.add_input(name + '_structural_mass', val=1., units='kg')
            self.add_input(name + '_cg_location', val=np.ones(3), units='m')

            self.declare_partials('cg', name + '_cg_location', rows=arange, cols=arange)
            self.declare_partials('cg', name + '_structural_mass')

        self.add_input('total_weight', val=1000., units='N')
        self.add_input('fuelburn', val=1.5, units='kg')
        self.add_input('W0', val=123., units='kg')
        self.add_input('load_factor', val=1.05)
        self.add_input('empty_cg', val=np.ones((3)), units='m')

        self.add_output('cg', val=np.ones(3), units='m')

        self.declare_partials('cg', 'total_weight')
        self.declare_partials('cg', 'W0')
        self.declare_partials('cg', 'fuelburn')
        self.declare_partials('cg', 'load_factor')
        self.declare_partials('cg', 'empty_cg', rows=arange, cols=arange)

    def compute(self, inputs, outputs):

        g = grav_constant * inputs['load_factor']
        W0_cg = inputs['W0'] * inputs['empty_cg']

        spar_cg = np.zeros(3)

        # Loop through the surfaces and compute the weighted cg location
        # of all structural spars
        for surface in self.options['surfaces']:
            name = surface['name']
            spar_cg = spar_cg + inputs[name + '_cg_location'] * inputs[name + '_structural_mass']

        # Compute the total cg of the aircraft based on the empty weight cg and
        # the structures cg. Here we assume the fuel weight is at the cg.
        outputs['cg'] = (W0_cg + spar_cg) / (inputs['total_weight'] / g - inputs['fuelburn'])

    def compute_partials(self, inputs, partials):

        g = grav_constant * inputs['load_factor']
        W0 = inputs['W0']
        cg = inputs['empty_cg']
        fb = inputs['fuelburn']
        tw = inputs['total_weight']
        W0_cg = W0 * cg

        spar_cg = np.zeros(3)

        # Loop through the surfaces and compute the weighted cg location
        # of all structural spars
        for surface in self.options['surfaces']:
            name = surface['name']
            spar_cg = spar_cg + inputs[name + '_cg_location'] * inputs[name + '_structural_mass']

        partials['cg', 'total_weight'] = - g * (W0_cg + spar_cg) / (tw - fb * g) ** 2
        partials['cg', 'fuelburn'] = g**2 * (W0_cg + spar_cg) / (tw - fb * g) ** 2
        partials['cg', 'load_factor'] =  grav_constant * tw * (W0_cg + spar_cg) / (tw - fb * g)**2
        partials['cg', 'empty_cg'] = W0 / (tw / g - fb)

        partials['cg', 'W0'] = cg / (tw / g - fb)

        for surface in self.options['surfaces']:
            name = surface['name']
            partials['cg', name + '_cg_location'] = -g * inputs[name + '_structural_mass'] / \
                (fb * g - tw)
            partials['cg', name + '_structural_mass'] = - inputs[name + '_cg_location'] * \
                g / (fb * g - tw)
