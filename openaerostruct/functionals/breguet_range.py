import numpy as np

import openmdao.api as om


class BreguetRange(om.ExplicitComponent):
    """
    Computes the fuel burn using the Breguet range equation using
    the computed CL, CD, weight, and provided specific fuel consumption, speed of sound,
    Mach number, initial weight, and range.

    Note that we add information from each lifting surface.

    Parameters
    ----------
    CL : float
        Total coefficient of lift (CL) for the lifting surface.
    CD : float
        Total coefficient of drag (CD) for the lifting surface.
    CT : float
        Specific fuel consumption for the entire aircraft.
    speed_of_sound : float
        The Mach speed, speed of sound, at the specified flight condition.
    R : float
        The total range of the aircraft, used to backcalculate the fuel mass.
    Mach_number : float
        The Mach number of the aircraft at the specified flight condition.
    W0 : float
        The operating empty weight of the aircraft, without fuel or structural
        mass. Supplied in kg despite being a 'weight' due to convention.
    load_factor : float
        Multiplicative factor on gravity. 1.0 is normal flight; 2.5 would be
        for a 2.5g manuever.
    _structural_mass : float
        Weight of a single lifting surface's structural spar.

    Returns
    -------
    fuelburn : float
        Computed fuel burn in kg based on the Breguet range equation.

    """

    def initialize(self):
        self.options.declare('surfaces', types=list)

    def setup(self):
        for surface in self.options['surfaces']:
            name = surface['name']
            self.add_input(name + '_structural_mass', val=1., units='kg')

        self.add_input('CT', val=0.25, units='1/s')
        self.add_input('CL', val=0.7)
        self.add_input('CD', val=0.02)
        self.add_input('speed_of_sound', val=100., units='m/s')
        self.add_input('R', val=3000., units='m')
        self.add_input('Mach_number', val=1.2)
        self.add_input('W0', val=200., units='kg')

        self.add_output('fuelburn', val=1., units='kg')

        self.declare_partials('*', '*')
        self.set_check_partial_options(wrt='*', method='cs', step=1e-30)

    def compute(self, inputs, outputs):

        CT = inputs['CT']
        a = inputs['speed_of_sound']
        R = inputs['R']
        M = inputs['Mach_number']
        W0 = inputs['W0']

        # Loop through the surfaces and add up the structural weights
        # to get the total structural weight.
        Ws = 0.
        for surface in self.options['surfaces']:
            name = surface['name']
            Ws += inputs[name + '_structural_mass']

        CL = inputs['CL']
        CD = inputs['CD']

        outputs['fuelburn'] = (W0 + Ws) * (np.exp(R * CT / a / M * CD / CL) - 1)

    def compute_partials(self, inputs, partials):

        CT = inputs['CT']
        a = inputs['speed_of_sound']
        R = inputs['R']
        M = inputs['Mach_number']
        W0 = inputs['W0']

        Ws = 0.
        for surface in self.options['surfaces']:
            name = surface['name']
            Ws += inputs[name + '_structural_mass']

        CL = inputs['CL']
        CD = inputs['CD']

        dfb_dCL = -(W0 + Ws) * np.exp(R * CT / a / M * CD / CL) \
            * R * CT / a / M * CD / CL ** 2
        dfb_dCD = (W0 + Ws) * np.exp(R * CT / a / M * CD / CL) \
            * R * CT / a / M / CL
        dfb_dCT = (W0 + Ws) * np.exp(R * CT / a / M * CD / CL) \
            * R / a / M / CL * CD
        dfb_dR = (W0 + Ws) * np.exp(R * CT / a / M * CD / CL) \
            / a / M / CL * CD * CT
        dfb_da = -(W0 + Ws) * np.exp(R * CT / a / M * CD / CL) \
            * R * CT / a**2 / M * CD / CL
        dfb_dM = -(W0 + Ws) * np.exp(R * CT / a / M * CD / CL) \
            * R * CT / a / M**2 * CD / CL

        dfb_dW = np.exp(R * CT / a / M * CD / CL) - 1

        partials['fuelburn', 'CL'] = dfb_dCL
        partials['fuelburn', 'CD'] = dfb_dCD
        partials['fuelburn', 'CT'] = dfb_dCT
        partials['fuelburn', 'speed_of_sound'] = dfb_da
        partials['fuelburn', 'R'] = dfb_dR
        partials['fuelburn', 'Mach_number'] = dfb_dM
        partials['fuelburn', 'W0'] = dfb_dW

        for surface in self.options['surfaces']:
            name = surface['name']
            inp_name = name + '_structural_mass'
            partials['fuelburn', inp_name] = dfb_dW
