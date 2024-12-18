import numpy as np

import openmdao.api as om


class MomentCoefficient(om.ExplicitComponent):
    """
    Compute the coefficient of moment (CM) and moment (M) for the entire aircraft.

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
    M[3] : numpy array
        The moment around the x-, y-, and z-axes at the cg point.
    CM[3] : numpy array
        The coefficient of moment around the x-, y-, and z-axes at the cg point.
    """

    def initialize(self):
        self.options.declare("surfaces", types=list)

    def setup(self):
        for surface in self.options["surfaces"]:
            name = surface["name"]
            nx = surface["mesh"].shape[0]
            ny = surface["mesh"].shape[1]

            self.add_input(name + "_b_pts", val=np.ones((nx - 1, ny, 3)), units="m", tags=["mphys_coupling"])
            self.add_input(name + "_widths", val=np.ones((ny - 1)), units="m", tags=["mphys_coupling"])
            self.add_input(name + "_chords", val=np.ones((ny)), units="m", tags=["mphys_coupling"])
            self.add_input(name + "_S_ref", val=1.0, units="m**2", tags=["mphys_coupling"])
            self.add_input(name + "_sec_forces", val=np.ones((nx - 1, ny - 1, 3)), units="N", tags=["mphys_coupling"])

        self.add_input("cg", val=np.ones((3)), units="m", tags=["mphys_input"])
        self.add_input("v", val=10.0, units="m/s", tags=["mphys_input"])
        self.add_input("rho", val=3.0, units="kg/m**3", tags=["mphys_input"])
        self.add_input("S_ref_total", val=1.0, units="m**2", tags=["mphys_input"])

        self.add_output("CM", val=np.ones((3)), tags=["mphys_result"])
        self.add_output("M", val=np.ones((3)), units="N*m", tags=["mphys_result"])

        self.declare_partials(of="*", wrt="*")

    def compute(self, inputs, outputs):
        cg = inputs["cg"]

        M = np.zeros((3))

        # Loop through each surface and find its contributions to the moment
        # of the aircraft based on the section forces and their location
        for j, surface in enumerate(self.options["surfaces"]):
            name = surface["name"]

            b_pts = inputs[name + "_b_pts"]
            widths = inputs[name + "_widths"]
            chords = inputs[name + "_chords"]
            S_ref = inputs[name + "_S_ref"]
            sec_forces = inputs[name + "_sec_forces"]

            # Compute the average chord for each panel and then the
            # mean aerodynamic chord (MAC) based on these chords and the
            # computed area
            panel_chords = (chords[1:] + chords[:-1]) * 0.5
            MAC = 1.0 / S_ref * np.sum(panel_chords**2 * widths)

            # If the surface is symmetric, then the previously computed MAC
            # is half what it should be
            if surface["symmetry"]:
                MAC *= 2.0

            # Get the moment arm acting on each panel, relative to the cg
            pts = (b_pts[:, 1:, :] + b_pts[:, :-1, :]) * 0.5
            diff = pts - cg

            # Compute the moment based on the previously computed moment
            # arm and the section forces
            moment = np.sum(np.cross(diff, sec_forces, axis=2), axis=0)

            # If the surface is symmetric, set the x- and z-direction moments
            # to 0 and double the y-direction moment
            if surface["symmetry"]:
                moment[:, 0] = 0.0
                moment[:, 1] *= 2.0
                moment[:, 2] = 0.0

            # Note: a scalar can be factored from a cross product, so I moved the division by MAC
            # down here for efficiency of calc and derivs.
            M = M + np.sum(moment, axis=0)

            # For the first (main) lifting surface, we save the MAC to correctly
            # normalize CM
            if j == 0:
                self.MAC_wing = MAC
                self.S_ref_wing = S_ref

        self.M = M

        # Output the moment vector
        outputs["M"] = M

        # Compute the normalized CM
        outputs["CM"] = M / (0.5 * inputs["rho"] * inputs["v"] ** 2 * inputs["S_ref_total"] * self.MAC_wing)

    def compute_partials(self, inputs, partials):
        cg = inputs["cg"]
        rho = inputs["rho"]
        S_ref_total = inputs["S_ref_total"]
        v = inputs["v"]

        # Cached values
        M = self.M
        MAC_wing = self.MAC_wing
        S_ref_wing = self.S_ref_wing

        # Scaling factor of one over the dynamic pressure times sum of reference areas times the wing MAC
        fact = 1.0 / (0.5 * rho * v**2 * S_ref_total * MAC_wing)

        partials["CM", "rho"] = -M * fact / rho
        partials["CM", "v"] = -2 * M * fact / v
        partials["CM", "S_ref_total"] = -M * fact / S_ref_total

        partials["CM", "cg"][:] = 0.0

        # Loop through each surface.
        for j, surface in enumerate(self.options["surfaces"]):
            name = surface["name"]
            nx = surface["mesh"].shape[0]
            ny = surface["mesh"].shape[1]

            partials["CM", name + "_sec_forces"][:] = 0.0
            partials["CM", name + "_b_pts"][:] = 0.0

            b_pts = inputs[name + "_b_pts"]
            widths = inputs[name + "_widths"]
            chords = inputs[name + "_chords"]
            S_ref = inputs[name + "_S_ref"]
            sec_forces = inputs[name + "_sec_forces"]

            # MAC derivs
            panel_chords = (chords[1:] + chords[:-1]) * 0.5
            MAC = 1.0 / S_ref * np.sum(panel_chords**2 * widths)

            # This produces a bi-diagonal matrix for the derivative of panel_chords with respect to chords
            # This transformation matrix is further used for multiple derivatives later
            dpc_dc = np.zeros((ny - 1, ny))
            idx = np.arange(ny - 1)
            dpc_dc[idx, idx] = 0.5
            dpc_dc[idx, idx + 1] = 0.5

            dMAC_dc = (2.0 / S_ref) * np.einsum("i,ij", panel_chords * widths, dpc_dc)
            dMAC_dw = (1.0 / S_ref) * panel_chords**2
            dMAC_dS = -MAC / S_ref

            # If the surface is symmetric, then the previously computed MAC
            # is half what it should be
            if surface["symmetry"]:
                MAC *= 2.0
                dMAC_dc *= 2.0
                dMAC_dw *= 2.0
                dMAC_dS *= 2.0

            # Compute the bound vortex(quarter chord) points at mid-panel
            pts = (b_pts[:, 1:, :] + b_pts[:, :-1, :]) * 0.5

            # Compute the vectors between the cg and the mid-panel bound vortex points
            diff = pts - cg

            # Compute the cross product of the panel bound vortex vectors from cg and the panel forces
            c = np.cross(diff, sec_forces, axis=2)

            # Compute the spanwise moment vector distribution by summing over each resulting column
            moment = np.sum(c, axis=0)

            # Compute the derviative of the moment vectors(c) with respect to the diff vectors(a) multiplied by -1
            dcda = np.zeros((3, nx - 1, ny - 1, 3))

            # Compute the derivative wrt to the first element of diff
            dcda[0, :, :, 1] = sec_forces[:, :, 2]
            dcda[0, :, :, 2] = -sec_forces[:, :, 1]

            # Compute the derivative wrt to the second element of diff
            dcda[1, :, :, 0] = -sec_forces[:, :, 2]
            dcda[1, :, :, 2] = sec_forces[:, :, 0]

            # Compute the derivative wrt to the third element of diff
            dcda[2, :, :, 0] = sec_forces[:, :, 1]
            dcda[2, :, :, 1] = -sec_forces[:, :, 0]

            # Compute the derviative of the moment vectors(c) with respect to the sec_forces vectors(b) multiplied by -1
            dcdb = np.zeros((3, nx - 1, ny - 1, 3))

            # Compute the derivative wrt to the first element of sec_forces
            dcdb[0, :, :, 1] = -diff[:, :, 2]
            dcdb[0, :, :, 2] = diff[:, :, 1]

            # Compute the derivative wrt to the second element of sec_forces
            dcdb[1, :, :, 0] = diff[:, :, 2]
            dcdb[1, :, :, 2] = -diff[:, :, 0]

            # Compute the derivative wrt to the third element of sec_forces
            dcdb[2, :, :, 0] = -diff[:, :, 1]
            dcdb[2, :, :, 1] = diff[:, :, 0]

            # Compute derivative of CM wrt to the sec_forces of the section by reshaping to 3 rows and multiplying by fact.
            partials["CM", name + "_sec_forces"] += dcdb.reshape((3, 3 * (nx - 1) * (ny - 1))) * fact

            # Compute derivative of M wrt to the sec_forces of the section by reshaping to 3 rows
            partials["M", name + "_sec_forces"] += dcdb.reshape((3, 3 * (nx - 1) * (ny - 1)))

            # Project the derviative of the moment vectors(c) with respect to the diff vectors(a)
            # onto the derivative of mid-panel chord distribution wrt to the chord distribution giving the derivative of
            # the moment vectors(c) with respect to the chord distribution(dc_dchord). This works because the diff
            # vectors are difference between the mid-panel bound vortex(quarter chord) points and cg which is static in this derivative.
            # The spanwise component of the mid-panel bound vortex(quarter chord) points have the same derivatrive wrt to the chord
            # distribution as the mid-panel chord distribution does wrt the the chord distribution.
            dc_dchord = np.einsum("ijkl,km->ijml", dcda, dpc_dc)

            # Compute the derivative of CM wrt to the bound vortex points(b_pts) by reshaping dc_dchord to three rows
            # and multiplying by fact.
            partials["CM", name + "_b_pts"] += dc_dchord.reshape((3, 3 * (nx - 1) * ny)) * fact

            # Compute the derivative of M wrt to the bound vortex points(b_pts) by reshaping dc_dchord to three rows
            partials["M", name + "_b_pts"] += dc_dchord.reshape((3, 3 * (nx - 1) * ny))

            # Reduce the derivative of the moment vectors(c) with respect to the diff vectors(a) by summing over all
            # chordwise and spanwise panels(j and k). Reduces to a 3x3 matrix for the whole surface by summing over all
            # panels.
            dcda = np.einsum("ijkl->il", dcda)

            # If the surface is symmetric, set the x- and z-direction moments
            # to 0 and double the y-direction moment
            if surface["symmetry"]:
                moment[:, 0] = 0.0
                moment[:, 1] *= 2.0
                moment[:, 2] = 0.0
                partials["CM", name + "_sec_forces"][0, :] = 0.0
                partials["CM", name + "_sec_forces"][1, :] *= 2.0
                partials["CM", name + "_sec_forces"][2, :] = 0.0
                partials["CM", name + "_b_pts"][0, :] = 0.0
                partials["CM", name + "_b_pts"][1, :] *= 2.0
                partials["CM", name + "_b_pts"][2, :] = 0.0
                partials["M", name + "_sec_forces"][0, :] = 0.0
                partials["M", name + "_sec_forces"][1, :] *= 2.0
                partials["M", name + "_sec_forces"][2, :] = 0.0
                partials["M", name + "_b_pts"][0, :] = 0.0
                partials["M", name + "_b_pts"][1, :] *= 2.0
                partials["M", name + "_b_pts"][2, :] = 0.0
                dcda[0, :] = 0.0
                dcda[1, :] *= 2.0
                dcda[2, :] = 0.0

            # Compute the derivative of CM wrt to the cg position which is negative dcda since diff = pts - cg times fact
            # Accumlate the derivative over each surface as the total moment vector is sum over all surfaces.
            partials["CM", "cg"] -= dcda * fact

            # Compute the derivative of M wrt to the cg position which is negative dcda since diff = pts - cg
            # Accumlate the derivative over each surface as the total moment vector is sum over all surfaces.
            partials["M", "cg"] -= dcda

            # Compute the total surface moment vector by summing spanwise
            M_j = np.sum(moment, axis=0)

            # For first surface, we need to save the deriv results
            if j == 0:
                # Compute a term by dividing fact by MAC. Note that MAC is the mean aerodynamic chord for the surface and
                # the MAC_wing terms already factored into fact is of the main wing surface
                term = fact / MAC

                # Compute the derivative of CM wrt to the chord distribution by taking the negative outer product of the
                # moment vector(M_j) time the term with the derivative of MAC wrt to the chord distribution. We only do
                # this for the main wing since CM only depends on the MAC of the main wing and the chord distribution of
                # the main wing is the only chord distribution of all the surfaces that can impact the MAC of the main wing.
                partials["CM", name + "_chords"] = -np.outer(M_j * term, dMAC_dc)

                # Compute the derivative of CM wrt to the width distribution by taking the negative outer product of the
                # moment vector(M_j) time the term with the derivative of MAC wrt to the width distribution. We only do
                # this for the main wing since CM only depends on the MAC of the main wing and the panel width distribution of
                # the main wing is the only panel width distribution of all the surfaces that can impact the MAC of the main wing.
                partials["CM", name + "_widths"] = -np.outer(M_j * term, dMAC_dw)

                # Compute the derivative of CM wrt to the surface S_ref by taking the negative outer product of the
                # moment vector(M_j) time the term with the derivative of MAC wrt to the surface S_ref. The CM depends on
                # the total references area of all surfaces including the main wing and the MAC of them main wing itself
                # As result, this derivative has two parts only for the main wing.
                # partials["CM", name + "_S_ref"] = -np.outer(M_j, dMAC_dS * term)
                partials["CM", name + "_S_ref"] = np.outer(M_j * fact, (1 / S_ref))

                # Cache the main wing's MAC derivatives
                base_name = name
                base_dMAC_dc = dMAC_dc
                base_dMAC_dw = dMAC_dw
                # base_dMAC_dS = dMAC_dS
            else:
                # Apply this surface's portion of the moment to the MAC_wing term.
                # We need to do this because CM is normalized by the MAC of the main wing
                term = fact / MAC_wing
                partials["CM", base_name + "_chords"] -= np.outer(M_j * term, base_dMAC_dc)
                partials["CM", base_name + "_widths"] -= np.outer(M_j * term, base_dMAC_dw)
                # partials["CM", base_name + "_S_ref"] -= np.outer(M_j, base_dMAC_dS * term)
                partials["CM", base_name + "_S_ref"] += np.outer(M_j * fact, (1 / S_ref_wing))
