import numpy as np

import openmdao.api as om


class MomentCoefficient(om.ExplicitComponent):
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
    S_ref_total : float
        Total surface area of the aircraft based on the sum of individual
        surface areas.

    Returns
    -------
    CM[3] : numpy array
        The coefficient of moment around the x-, y-, and z-axes at the cg point.
    """

    def initialize(self):
        self.options.declare('surfaces', types=list)

    def setup(self):
        for surface in self.options['surfaces']:
            name = surface['name']
            nx = surface['mesh'].shape[0]
            ny = surface['mesh'].shape[1]

            self.add_input(name + '_b_pts', val=np.ones((nx-1, ny, 3)), units='m')
            self.add_input(name + '_widths', val=np.ones((ny-1)), units='m')
            self.add_input(name + '_chords', val=np.ones((ny)), units='m')
            self.add_input(name + '_S_ref', val=1., units='m**2')
            self.add_input(name + '_sec_forces', val=np.ones((nx-1, ny-1, 3)), units='N')

        self.add_input('cg', val=np.ones((3)), units='m')
        self.add_input('v', val=10., units='m/s')
        self.add_input('rho', val=3., units='kg/m**3')
        self.add_input('S_ref_total', val=1., units='m**2')

        self.add_output('CM', val=np.ones((3)))

        self.declare_partials(of='*', wrt='*')

    def compute(self, inputs, outputs):
        cg = inputs['cg']

        M = np.zeros((3))

        # Loop through each surface and find its contributions to the moment
        # of the aircraft based on the section forces and their location
        for j, surface in enumerate(self.options['surfaces']):
            name = surface['name']

            b_pts = inputs[name + '_b_pts']
            widths = inputs[name + '_widths']
            chords = inputs[name + '_chords']
            S_ref = inputs[name + '_S_ref']
            sec_forces = inputs[name + '_sec_forces']

            # Compute the average chord for each panel and then the
            # mean aerodynamic chord (MAC) based on these chords and the
            # computed area
            panel_chords = (chords[1:] + chords[:-1]) * 0.5
            MAC = 1. / S_ref * np.sum(panel_chords**2 * widths)

            # If the surface is symmetric, then the previously computed MAC
            # is half what it should be
            if surface['symmetry']:
                MAC *= 2.0

            # Get the moment arm acting on each panel, relative to the cg
            pts = (b_pts[:, 1:, :] + b_pts[:, :-1, :]) * 0.5
            diff = (pts - cg)

            # Compute the moment based on the previously computed moment
            # arm and the section forces
            moment = np.sum(np.cross(diff, sec_forces, axis=2), axis=0)

            # If the surface is symmetric, set the x- and z-direction moments
            # to 0 and double the y-direction moment
            if surface['symmetry']:
                moment[:, 0] = 0.
                moment[:, 1] *= 2.0
                moment[:, 2] = 0.

            # Note: a scalar can be factored from a cross product, so I moved the division by MAC
            # down here for efficiency of calc and derivs.
            M = M + np.sum(moment, axis=0)

            # For the first (main) lifting surface, we save the MAC to correctly
            # normalize CM
            if j == 0:
                self.MAC_wing = MAC

        self.M = M

        # Compute the normalized CM
        rho = inputs['rho']
        outputs['CM'] = M / (0.5 * rho * inputs['v']**2 * inputs['S_ref_total'] * self.MAC_wing)

    def compute_partials(self, inputs, partials):
        cg = inputs['cg']
        rho = inputs['rho']
        S_ref_total = inputs['S_ref_total']
        v = inputs['v']

        # Cached values
        M = self.M
        MAC_wing = self.MAC_wing

        fact = 1.0 / (0.5 * rho * v**2 * S_ref_total * MAC_wing)

        partials['CM', 'rho'] = -M * fact**2 * 0.5 * v**2 * S_ref_total * MAC_wing
        partials['CM', 'v'] = -M * fact**2 * rho * v * S_ref_total * MAC_wing
        partials['CM', 'S_ref_total'] = -M * fact**2 * 0.5 * rho * v**2 * MAC_wing

        partials['CM', 'cg'][:] = 0.0

        # Loop through each surface.
        for j, surface in enumerate(self.options['surfaces']):
            name = surface['name']
            nx = surface['mesh'].shape[0]
            ny = surface['mesh'].shape[1]

            partials['CM', name + '_sec_forces'][:] = 0.0
            partials['CM', name + '_b_pts'][:] = 0.0

            b_pts = inputs[name + '_b_pts']
            widths = inputs[name + '_widths']
            chords = inputs[name + '_chords']
            S_ref = inputs[name + '_S_ref']
            sec_forces = inputs[name + '_sec_forces']

            # MAC derivs
            panel_chords = (chords[1:] + chords[:-1]) * 0.5
            MAC = 1. / S_ref * np.sum(panel_chords**2 * widths)

            # This transformation is used for multiple derivatives
            dpc_dc = np.zeros((ny-1, ny))
            idx = np.arange(ny-1)
            dpc_dc[idx, idx] = 0.5
            dpc_dc[idx, idx+1] = 0.5

            dMAC_dc = (2.0 / S_ref) * np.einsum('i,ij', panel_chords * widths, dpc_dc)
            dMAC_dw = (1.0 / S_ref) * panel_chords**2
            dMAC_dS = -MAC / S_ref

            # If the surface is symmetric, then the previously computed MAC
            # is half what it should be
            if surface['symmetry']:
                MAC *= 2.0
                dMAC_dc *= 2.0
                dMAC_dw *= 2.0
                dMAC_dS *= 2.0

            # diff derivs
            pts = (b_pts[:, 1:, :] + b_pts[:, :-1, :]) * 0.5
            diff = (pts - cg)

            c = np.cross(diff, sec_forces, axis=2)
            moment = np.sum(c, axis=0)

            dcda = np.zeros((3, nx-1, ny-1, 3))
            dcda[0, :, :, 1] = sec_forces[:, :, 2]
            dcda[0, :, :, 2] = -sec_forces[:, :, 1]
            dcda[1, :, :, 0] = -sec_forces[:, :, 2]
            dcda[1, :, :, 2] = sec_forces[:, :, 0]
            dcda[2, :, :, 0] = sec_forces[:, :, 1]
            dcda[2, :, :, 1] = -sec_forces[:, :, 0]

            dcdb = np.zeros((3, nx-1, ny-1, 3))
            dcdb[0, :, :, 1] = -diff[:, :, 2]
            dcdb[0, :, :, 2] = diff[:, :, 1]
            dcdb[1, :, :, 0] = diff[:, :, 2]
            dcdb[1, :, :, 2] = -diff[:, :, 0]
            dcdb[2, :, :, 0] = -diff[:, :, 1]
            dcdb[2, :, :, 1] = diff[:, :, 0]

            partials['CM', name + '_sec_forces'] += dcdb.reshape((3, 3*(nx-1)*(ny-1))) * fact

            dc_dchord = np.einsum('ijkl,km->ijml', dcda, dpc_dc)
            partials['CM', name + '_b_pts'] += dc_dchord.reshape((3, 3*(nx-1)*ny)) * fact

            dcda = np.einsum('ijkl->il', dcda)

            # If the surface is symmetric, set the x- and z-direction moments
            # to 0 and double the y-direction moment
            if surface['symmetry']:
                moment[:, 0] = 0.
                moment[:, 1] *= 2.0
                moment[:, 2] = 0.
                partials['CM', name + '_sec_forces'][0, :] = 0.
                partials['CM', name + '_sec_forces'][1, :] *= 2.0
                partials['CM', name + '_sec_forces'][2, :] = 0.
                partials['CM', name + '_b_pts'][0, :] = 0.
                partials['CM', name + '_b_pts'][1, :] *= 2.0
                partials['CM', name + '_b_pts'][2, :] = 0.
                dcda[0, :] = 0.
                dcda[1, :] *= 2.0
                dcda[2, :] = 0.

            partials['CM', 'cg'] -= dcda * fact

            M_j = np.sum(moment, axis=0)
            term = fact / MAC
            partials['CM', name + '_chords'] = -np.outer(M_j * term, dMAC_dc)
            partials['CM', name + '_widths'] = -np.outer(M_j * term, dMAC_dw)
            partials['CM', name + '_S_ref'] = -np.outer(M_j, dMAC_dS * term)

            # For first surface, we need to save the deriv results
            if j == 0:
                base_name = name
                base_dMAC_dc = dMAC_dc
                base_dMAC_dw = dMAC_dw
                base_dMAC_dS = dMAC_dS

            else:
                # Apply this surface's portion of the moment to the MAC_wing term.
                term = fact / (MAC_wing * MAC)
                partials['CM', base_name + '_chords'] -= np.outer(M_j * term, base_dMAC_dc)
                partials['CM', base_name + '_widths'] -= np.outer(M_j * term, base_dMAC_dw)
                partials['CM', base_name + '_S_ref'] -= np.outer(M_j, base_dMAC_dS * term)
