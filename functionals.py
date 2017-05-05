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
            self.add_param(name+'structural_weight', val=0.)

        self.add_param('CL', val=0.)
        self.add_param('CD', val=0.)

        self.add_output('fuelburn', val=0.)
        self.add_output('weighted_obj', val=0.)

        self.deriv_options['type'] = 'cs'
        self.deriv_options['form'] = 'central'

    def solve_nonlinear(self, params, unknowns, resids):
        CT = self.prob_dict['CT']
        a = self.prob_dict['a']
        R = self.prob_dict['R']
        M = self.prob_dict['M']
        W0 = self.prob_dict['W0'] * self.prob_dict['g']

        beta = self.prob_dict['beta']

        # Loop through the surfaces and add up the structural weights
        # to get the total structural weight.
        Ws = 0.
        for surface in self.surfaces:
            name = surface['name']
            Ws += params[name+'structural_weight']

        CL = params['CL']
        CD = params['CD']
        fuelburn = np.sum((W0 + Ws) * (np.exp(R * CT / a / M * CD / CL) - 1))

        # Convert fuelburn from N to kg
        unknowns['fuelburn'] = fuelburn / self.prob_dict['g']

        # This lines makes the 'weight' the total aircraft weight
        unknowns['weighted_obj'] = (beta * fuelburn + (1 - beta) * (W0 + Ws + fuelburn)) / self.prob_dict['g']

        # Whereas this line only considers the structural weight
        # unknowns['weighted_obj'] = (beta * fuelburn + (1 - beta) * Ws) / self.prob_dict['g']

class FunctionalEquilibrium(Component):
    """
    Lift = weight constraint.

    Note that we add information from each lifting surface.

    Parameters
    ----------
    L : float
        Total lift for the lifting surface.
    structural_weight : float
        Total weight of the structural spar.

    fuelburn : float
        Computed fuel burn in kg based on the Breguet range equation.

    Returns
    -------
    L_equals_W : float
        Equality constraint for lift = total weight. L_equals_W = 0 for the constraint to be satisfied.
    total_weight : float
        Total weight of the entire aircraft, including W0, all structural weights,
        and fuel.

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
        self.add_output('L_equals_W', val=0.)
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
        unknowns['L_equals_W'] = (tot_weight - L) / tot_weight

class ComputeCG(Component):
    """
    Compute the center of gravity of the entire aircraft based on the inputted W0
    and its corresponding cg and the weighted sum of each surface's structural
    weight and location.

    Note that we add information from each lifting surface.

    Parameters
    ----------
    nodes[ny, 3] : numpy array
        Flattened array with coordinates for each FEM node.
    structural_weight : float
        Total weight of the structural spar for a given surface.
    cg_location[3] : numpy array
        Location of the structural spar's cg for a given surface.

    total_weight : float
        Total weight of the entire aircraft, including W0, all structural weights,
        and fuel.
    fuelburn : float
        Computed fuel burn in kg based on the Breguet range equation.

    Returns
    -------
    cg[3] : numpy array
        The x, y, z coordinates of the center of gravity for the entire aircraft.
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
        W0 = self.prob_dict['W0']
        W0_cg = W0 * self.prob_dict['cg'] * g

        spar_cg = 0.
        structural_weight = 0.
        for surface in self.surfaces:
            name = surface['name']
            spar_cg = params[name + 'cg_location'] * params[name + 'structural_weight']
            structural_weight += params[name + 'structural_weight']
        tot_weight = structural_weight + params['fuelburn'] * g + W0

        unknowns['cg'] = (W0_cg + spar_cg) / (params['total_weight'] - params['fuelburn'] * g)

class ComputeCM(Component):
    """
    Compute the coefficient of moment (CM) for the entire aircraft.

    Parameters
    ----------
    b_pts[nx-1, ny, 3] : numpy array
        Bound points for the horseshoe vortices, found along the 1/4 chord.
    widths[ny-1] : numpy array
        The spanwise widths of each individual panel.
    chords[ny] : numpy array
        The chordwise length of the entire airfoil following the camber line.
    S_ref : float
        The reference area of the lifting surface.
    sec_forces[nx-1, ny-1, 3] : numpy array
        Contains the sectional forces acting on each panel.
        Stored in Fortran order (only relevant with more than one chordwise
        panel).

    cg[3] : numpy array
        The x, y, z coordinates of the center of gravity for the entire aircraft.
    v : float
        Freestream air velocity in m/s.
    rho : float
        Air density in kg/m^3.

    Returns
    -------
    CM[3] : numpy array
        The coefficient of moment around the x-, y-, and z-axes at the cg point.
    """

    def __init__(self, surfaces, prob_dict):
        super(ComputeCM, self).__init__()

        tot_panels = 0
        for surface in surfaces:
            name = surface['name']
            ny = surface['num_y']
            nx = surface['num_x']

            self.add_param(name+'b_pts', val=np.zeros((nx-1, ny, 3), dtype=data_type))
            self.add_param(name+'widths', val=np.zeros((ny-1), dtype=data_type))
            self.add_param(name+'chords', val=np.zeros((ny), dtype=data_type))
            self.add_param(name+'S_ref', val=0.)
            self.add_param(name+'sec_forces', val=np.zeros((nx-1, ny-1, 3), dtype=data_type))

        self.add_param('cg', val=np.zeros((3), dtype=data_type))
        self.add_param('v', val=10.)
        self.add_param('rho', val=3.)

        self.add_output('CM', val=np.zeros((3), dtype=data_type))

        self.surfaces = surfaces
        self.prob_dict = prob_dict

        if not fortran_flag:
            self.deriv_options['type'] = 'cs'
            self.deriv_options['form'] = 'central'

    def solve_nonlinear(self, params, unknowns, resids):
        rho = params['rho']
        cg = params['cg']

        S_ref_tot = 0.
        M = np.zeros((3), dtype=data_type)

        for j, surface in enumerate(self.surfaces):
            name = surface['name']
            nx = surface['num_x']
            ny = surface['num_y']

            b_pts = params[name+'b_pts']
            widths = params[name+'widths']
            chords = params[name+'chords']
            S_ref = params[name+'S_ref']
            sec_forces = params[name+'sec_forces']

            panel_chords = (chords[1:] + chords[:-1]) / 2.
            MAC = 1. / S_ref * np.sum(panel_chords**2 * widths)

            if fortran_flag:
                M_tmp = OAS_API.oas_api.momentcalc(b_pts, cg, chords, widths, S_ref, sec_forces, surface['symmetry'])
                M += M_tmp

            else:
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
                    moment[:, 1] *= 2
                    moment[:, 2] = 0.
                M += np.sum(moment, axis=0)

            # For the first (main) lifting surface, we save the MAC
            if j == 0:
                MAC_wing = MAC
            S_ref_tot += S_ref

        self.M = M

        # Use the user-provided reference area; otherwise compute the total
        # area of all lifting surfaces.
        if self.prob_dict['S_ref_total'] is None:
            self.S_ref_tot = S_ref_tot
        else:
            self.S_ref_tot = self.prob_dict['S_ref_total']

        unknowns['CM'] = M / (0.5 * rho * params['v']**2 * self.S_ref_tot * MAC_wing)

    def linearize(self, params, unknowns, resids):

        jac = self.alloc_jacobian()

        cg = params['cg']
        rho = params['rho']
        v = params['v']

        for j in range(3):
            CMb = np.zeros((3))
            CMb[j] = 1.

            for i, surface in enumerate(self.surfaces):
                name = surface['name']
                ny = surface['num_y']

                b_pts = params[name+'b_pts']
                widths = params[name+'widths']
                chords = params[name+'chords']
                S_ref = params[name+'S_ref']
                sec_forces = params[name+'sec_forces']

                if i == 0:
                    panel_chords = (chords[1:] + chords[:-1]) / 2.
                    MAC = 1. / S_ref * np.sum(panel_chords**2 * widths)
                    temp1 = self.S_ref_tot * MAC
                    temp0 = 0.5 * rho * v**2
                    temp = temp0 * temp1
                    tempb = np.sum(-(self.M * CMb / temp)) / temp
                    Mb = CMb / temp
                    Mb_master = Mb.copy()
                    jac['CM', 'rho'][j] += v**2*temp1*0.5*tempb
                    jac['CM', 'v'][j] += 0.5*rho*temp1*2*v*tempb
                    s_totb = temp0 * MAC * tempb
                    macb = temp0 * self.S_ref_tot * tempb
                    chordsb = np.zeros((ny))
                    tempb0 = macb / S_ref
                    panel_chordsb = 2*panel_chords*widths*tempb0
                    widthsb = panel_chords**2*tempb0
                    sb = -np.sum(panel_chords**2*widths)*tempb0/S_ref
                    chordsb[1:] += panel_chordsb/2.
                    chordsb[:-1] += panel_chordsb/2.
                    cb = chordsb
                    wb = widthsb

                else:
                    cb = 0.
                    wb = 0.
                    sb = 0.
                    Mb = Mb_master.copy()

                bptsb, cgb, chordsb, widthsb, S_refb, sec_forcesb, _ = OAS_API.oas_api.momentcalc_b(b_pts, cg, chords, widths, S_ref, sec_forces, surface['symmetry'], Mb)

                jac['CM', 'cg'][j, :] += cgb
                jac['CM', name+'b_pts'][j, :] += bptsb.flatten()
                jac['CM', name+'chords'][j, :] += chordsb + cb
                jac['CM', name+'widths'][j, :] += widthsb + wb
                jac['CM', name+'sec_forces'][j, :] += sec_forcesb.flatten()
                jac['CM', name+'S_ref'][j, :] += S_refb + sb + s_totb

        return jac

class ComputeTotalCLCD(Component):
    """
    Compute the coefficients of lift (CL) and drag (CD) for the entire aircraft.

    Parameters
    ----------
    CL : float
        Coefficient of lift (CL) for one lifting surface.
    CD : float
        Coefficient of drag (CD) for one lifting surface.
    S_ref : float
        Surface area for one lifting surface.
    v : float
        Fresstream air velocity.
    rho : float
        Air density in kg/m^3.

    Returns
    -------
    CL : float
        Total coefficient of lift (CL) for the entire aircraft.
    CD : float
        Total coefficient of drag (CD) for the entire aircraft.

    """

    def __init__(self, surfaces, prob_dict):
        super(ComputeTotalCLCD, self).__init__()

        for surface in surfaces:
            name = surface['name']
            self.add_param(name+'CL', val=0.)
            self.add_param(name+'CD', val=0.)
            self.add_param(name+'S_ref', val=0.)

        self.add_param('v', val=10.)
        self.add_param('rho', val=3.)

        self.add_output('CL', val=0.)
        self.add_output('CD', val=0.)

        self.surfaces = surfaces
        self.prob_dict = prob_dict

    def solve_nonlinear(self, params, unknowns, resids):
        rho = params['rho']
        v = params['v']

        CL = 0.
        CD = 0.
        computed_total_S_ref = 0.
        for surface in self.surfaces:
            name = surface['name']
            S_ref = params[name+'S_ref']
            CL += params[name+'CL'] * S_ref
            CD += params[name+'CD'] * S_ref
            computed_total_S_ref += S_ref

        if self.prob_dict['S_ref_total'] is not None:
            S_ref_total = self.prob_dict['S_ref_total']
        else:
            S_ref_total = computed_total_S_ref

        unknowns['CL'] = CL / S_ref_total
        unknowns['CD'] = CD / S_ref_total
        self.S_ref_total = S_ref_total

    def linearize(self, params, unknowns, resids):
        jac = self.alloc_jacobian()

        for surface in self.surfaces:
            name = surface['name']
            S_ref = params[name+'S_ref']
            jac['CL', name+'CL'] = S_ref / self.S_ref_total
            jac['CD', name+'CD'] = S_ref / self.S_ref_total
            jac['CL', name+'S_ref'] = params[name+'CL'] / self.S_ref_total
            jac['CD', name+'S_ref'] = params[name+'CD'] / self.S_ref_total

        return jac


class TotalPerformance(Group):
    """
    Group to contain the total aerostructural performance components.
    """

    def __init__(self, surfaces, prob_dict):
        super(TotalPerformance, self).__init__()

        self.add('fuelburn',
                 FunctionalBreguetRange(surfaces, prob_dict),
                 promotes=['*'])
        self.add('L_equals_W',
                 FunctionalEquilibrium(surfaces, prob_dict),
                 promotes=['*'])
        self.add('CG',
                 ComputeCG(surfaces, prob_dict),
                 promotes=['*'])
        self.add('moment',
                 ComputeCM(surfaces, prob_dict),
                 promotes=['*'])
        self.add('CL_CD',
                 ComputeTotalCLCD(surfaces, prob_dict),
                 promotes=['*'])

class TotalAeroPerformance(Group):
    """
    Group to contain the total aerodynamic performance components.
    """

    def __init__(self, surfaces, prob_dict):
        super(TotalAeroPerformance, self).__init__()

        self.add('moment',
                 ComputeCM(surfaces, prob_dict),
                 promotes=['*'])
        self.add('CL_CD',
                 ComputeTotalCLCD(surfaces, prob_dict),
                 promotes=['*'])
