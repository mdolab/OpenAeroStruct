from __future__ import division, print_function
import numpy as np

from openmdao.api import Component, Group

try:
    import OAS_API
    fortran_flag = True
    data_type = float
except:
    fortran_flag = False
    data_type = complex

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
        self.add_output('total_weight', val=0.)

        self.deriv_options['type'] = 'cs'
        self.deriv_options['form'] = 'central'

    def solve_nonlinear(self, params, unknowns, resids):
        structural_weight = 0.
        L = 0.
        W0 = self.prob_dict['W0'] * self.prob_dict['g']
        for surface in self.surfaces:
            name = surface['name']
            structural_weight += params[name+'structural_weight']
            L += params[name+'L']

        tot_weight = structural_weight + params['fuelburn'] * self.prob_dict['g'] + W0

        unknowns['total_weight'] = tot_weight
        unknowns['eq_con'] = (tot_weight - L) / W0

class ComputeCM(Component):
    """

    Parameters
    ----------

    Returns
    -------
    """

    def __init__(self, surfaces, prob_dict):
        super(ComputeCM, self).__init__()

        tot_panels = 0
        for surface in surfaces:
            name = surface['name']
            ny = surface['num_y']
            nx = surface['num_x']
            tot_panels += (nx - 1) * (ny - 1)

            self.add_param(name+'b_pts', val=np.zeros((nx-1, ny, 3), dtype=data_type))
            self.add_param(name+'widths', val=np.zeros((ny-1), dtype=data_type))
            self.add_param(name+'lengths', val=np.zeros((ny), dtype=data_type))
            self.add_param(name+'S_ref', val=0.)
            self.add_param(name+'sec_forces', val=np.zeros((nx-1, ny-1, 3), dtype=data_type))

        self.tot_panels = tot_panels

        self.add_param('cg', val=np.zeros((3), dtype=data_type))
        self.add_param('v', val=10.)
        self.add_param('rho', val=3.)

        self.add_output('CM', val=np.zeros((3), dtype=data_type))

        self.surfaces = surfaces

        self.mtx = np.zeros((tot_panels, tot_panels, 3), dtype=data_type)
        self.v = np.zeros((tot_panels, 3), dtype=data_type)

        if not fortran_flag:
            self.deriv_options['type'] = 'cs'
            self.deriv_options['form'] = 'central'

    def solve_nonlinear(self, params, unknowns, resids):
        rho = params['rho']
        cg = params['cg']

        i = 0
        S_ref_tot = 0.
        M = np.zeros((3), dtype=data_type)
        for surface in self.surfaces:
            name = surface['name']
            nx = surface['num_x']
            ny = surface['num_y']

            num_panels = (nx - 1) * (ny - 1)

            b_pts = params[name+'b_pts']
            widths = params[name+'widths']
            lengths = params[name+'lengths']
            S_ref = params[name+'S_ref']
            sec_forces = params[name+'sec_forces']

            if fortran_flag:
                M_tmp = OAS_API.oas_api.momentcalc(b_pts, cg, lengths, widths, S_ref, sec_forces, surface['symmetry'])
                M += M_tmp

            else:
                panel_chords = (lengths[1:] + lengths[:-1]) / 2.
                MAC = 1. / S_ref * np.sum(panel_chords**2 * widths)

                if surface['symmetry']:
                    MAC *= 2

                pts = (params[name+'b_pts'][:, 1:, :] + \
                    params[name+'b_pts'][:, :-1, :]) / 2
                diff = (pts - cg) / MAC
                moment = np.zeros((ny - 1, 3), dtype=data_type)
                for ind in range(nx-1):
                    moment += np.cross(diff[ind, :, :], sec_forces[ind, :, :], axis=1)

                if surface['symmetry']:
                    moment[:, 0] = 0.
                    moment[:, 2] = 0.
                M += np.sum(moment, axis=0)

            S_ref_tot += S_ref
            i += num_panels

        self.M = M
        self.S_ref_tot = S_ref_tot

        unknowns['CM'] = M / (0.5 * rho * params['v']**2 * S_ref_tot)

    def apply_linear(self, params, unknowns, dparams, dunknowns, dresids, mode):

        if mode == 'fwd':

            rho = params['rho']
            v = params['v']
            i = 0
            Md = 0.
            S_ref_tot = 0.
            S_ref_totd = 0.
            for surface in self.surfaces:
                name = surface['name']
                nx = surface['num_x']
                ny = surface['num_y']

                num_panels = (nx - 1) * (ny - 1)

                b_pts = params[name+'b_pts']
                widths = params[name+'widths']
                lengths = params[name+'lengths']
                S_ref = params[name+'S_ref']
                sec_forces = params[name+'sec_forces'].real

                M, Md_tmp = OAS_API.oas_api.momentcalc_d(
                                            b_pts, dparams[name+'b_pts'],
                                            params['cg'], dparams['cg'],
                                            lengths, dparams[name+'lengths'],
                                            widths, dparams[name+'widths'],
                                            S_ref, dparams[name+'S_ref'],
                                            sec_forces, dparams[name+'sec_forces'],
                                            surface['symmetry'])

                Md += Md_tmp
                S_ref_tot += S_ref
                S_ref_totd += dparams[name+'S_ref']
                i += num_panels

            dresids['CM'] = (Md*0.5*rho*v**2*S_ref_tot - M*0.5*((dparams['rho']*S_ref_tot + rho*S_ref_totd)*v**2 + rho*S_ref_tot*2*v*dparams['v'])) / (.5*rho*v**2*S_ref_tot)**2

        if mode == 'rev':
            i = 0
            cg = params['cg']
            rho = params['rho']
            v = params['v']
            vb = np.zeros(self.v.shape)

            temp0 = 0.5*rho*self.S_ref_tot
            temp = temp0*v**2
            tempb = np.sum(-(self.M * dresids['CM'] / temp))/temp
            tempb0 = v**2*tempb
            Mb = dresids['CM']/temp
            dparams['rho'] += self.S_ref_tot*0.5*tempb0
            sb = 0.5*rho*tempb0
            dparams['v'] += temp0*2*v*tempb

            for surface in self.surfaces:
                name = surface['name']
                nx = surface['num_x']
                ny = surface['num_y']
                num_panels = (nx - 1) * (ny - 1)

                b_pts = params[name+'b_pts']
                widths = params[name+'widths']
                lengths = params[name+'lengths']
                S_ref = params[name+'S_ref']
                sec_forces = params[name+'sec_forces']

                bptsb, cgb, lengthsb, widthsb, S_refb, sec_forcesb, _ = OAS_API.oas_api.momentcalc_b(b_pts, cg, lengths, widths, S_ref, sec_forces, surface['symmetry'], Mb)

                dparams['cg'] += cgb
                dparams[name+'b_pts'] += bptsb
                dparams[name+'lengths'] += lengthsb
                dparams[name+'widths'] += widthsb
                dparams[name+'sec_forces'] += sec_forcesb
                dparams[name+'S_ref'] += S_refb + sb

                i += num_panels

class ComputeCG(Component):
    """ Lift = weight constraint.

    Note that we add information from each lifting surface.

    Parameters
    ----------

    Returns
    -------
    """

    def __init__(self, surfaces, prob_dict):
        super(ComputeCG, self).__init__()

        self.prob_dict = prob_dict
        self.surfaces = surfaces

        for surface in surfaces:
            name = surface['name']

            self.add_param(name+'nodes', val=0.)
            self.add_param(name+'structural_weight', val=0.)
            self.add_param(name+'cg_location', val=np.zeros((3), dtype=data_type))

        self.add_param('total_weight', val=0.)
        self.add_param('fuelburn', val=0.)
        self.add_output('cg', val=np.zeros((3), dtype=complex))

        self.deriv_options['type'] = 'cs'
        self.deriv_options['form'] = 'central'

    def solve_nonlinear(self, params, unknowns, resids):
        g = self.prob_dict['g']
        W0_cg = self.prob_dict['W0'] * self.prob_dict['cg'] * g

        spar_cg = 0.
        for surface in self.surfaces:
            name = surface['name']
            spar_cg = params[name + 'cg_location'] * params[name + 'structural_weight']

        unknowns['cg'] = (W0_cg + spar_cg) / (params['total_weight'] - params['fuelburn'] * g)

class TotalPerformance(Group):

    def __init__(self, surfaces, prob_dict):
        super(TotalPerformance, self).__init__()

        self.add('fuelburn',
                 FunctionalBreguetRange(surfaces, prob_dict),
                 promotes=['*'])
        self.add('eq_con',
                 FunctionalEquilibrium(surfaces, prob_dict),
                 promotes=['*'])
        self.add('moment',
                 ComputeCM(surfaces, prob_dict),
                 promotes=['*'])
        self.add('CG',
                 ComputeCG(surfaces, prob_dict),
                 promotes=['*'])

class TotalAeroPerformance(Group):
    "Placeholder"

    def __init__(self, surfaces, prob_dict):
        super(TotalAeroPerformance, self).__init__()

        self.add('moment',
                 ComputeCM(surfaces, prob_dict),
                 promotes=['*'])
