import numpy as np

import openmdao.api as om

from openaerostruct.structures.utils import norm, unit


class TsaiWuWingbox(om.ExplicitComponent):
    """Compute the von Mises stresses for each element.
    See Chauhan et al. (https://doi.org/10.1007/978-3-319-97773-7_38) for more.

    NOTE: Compute the Tsai-Wu failure criteria for each element, and return an array of the Strength Ratios (SR) for each ply (4) for each element (4).

    Parameters
    ----------
    nodes[ny, 3] : numpy array
        Flattened array with coordinates for each FEM node.
    disp[ny, 6] : numpy array
        Displacements of each FEM node.
    Qz[ny-1] : numpy array
        First moment of area above the neutral axis parallel to the local
        z-axis (for each wingbox segment).
    J[ny-1] : numpy array
        Torsion constants for each wingbox segment.
    A_enc[ny-1] : numpy array
        Cross-sectional enclosed area (measured using the material midlines) of
        each wingbox segment.
    spar_thickness[ny-1] : numpy array
        Material thicknesses of the front and rear spars for each wingbox segment.
    htop[ny-1] : numpy array
        Distance to the point on the top skin that is the farthest away from
        the local-z neutral axis (for each wingbox segment).
    hbottom[ny-1] : numpy array
        Distance to the point on the bottom skin that is the farthest away from
        the local-z neutral axis (for each wingbox segment).
    hfront[ny-1] : numpy array
        Distance to the point on the front spar that is the farthest away from
        the local-y neutral axis (for each wingbox segment).
    hrear[ny-1] : numpy array
        Distance to the point on the rear spar that is the farthest away
        from the local-y neutral axis (for each wingbox segment).

    Returns
    -------
    vonmises[ny-1, 4] : numpy array
        von Mises stresses for 4 stress combinations for each FEM element.

    NOTE: tsaiwu_sr[ny-1, 16] : numpy array
        Tsai-Wu Strength Ratios for 16 plys for each FEM element.

    """

    def initialize(self):
        self.options.declare("surface", types=dict)

    def setup(self):
        self.surface = surface = self.options["surface"]

        self.ny = surface["mesh"].shape[1]

        self.add_input("nodes", val=np.zeros((self.ny, 3)), units="m")
        self.add_input("disp", val=np.zeros((self.ny, 6)), units="m")
        self.add_input("Qz", val=np.zeros((self.ny - 1)), units="m**3")
        self.add_input("J", val=np.zeros((self.ny - 1)), units="m**4")
        self.add_input("A_enc", val=np.zeros((self.ny - 1)), units="m**2")
        self.add_input("spar_thickness", val=np.zeros((self.ny - 1)), units="m")
        self.add_input("htop", val=np.zeros((self.ny - 1)), units="m")
        self.add_input("hbottom", val=np.zeros((self.ny - 1)), units="m")
        self.add_input("hfront", val=np.zeros((self.ny - 1)), units="m")
        self.add_input("hrear", val=np.zeros((self.ny - 1)), units="m")

        # self.add_output("vonmises", val=np.zeros((self.ny - 1, 4)), units="N/m**2")
        self.add_output("tsaiwu_sr", val=np.zeros((self.ny - 1, 16))) # NOTE: To be verified again

        self.E = surface["E"]
        self.G = surface["G"]

        self.tssf = surface["strength_factor_for_upper_skin"]

        

        self.declare_partials("*", "*", method="cs")

    def compute(self, inputs, outputs):
        disp = inputs["disp"]
        nodes = inputs["nodes"]
        A_enc = inputs["A_enc"]
        Qy = inputs["Qz"]
        J = inputs["J"]
        htop = inputs["htop"]
        hbottom = inputs["hbottom"]
        hfront = inputs["hfront"]
        hrear = inputs["hrear"]
        spar_thickness = inputs["spar_thickness"]
        # vonmises = outputs["vonmises"]
        tsaiwu_sr = outputs["tsaiwu_sr"]

        # Only use complex type for these arrays if we're using cs to check derivs
        dtype = type(disp[0, 0])
        T = np.zeros((3, 3), dtype=dtype)
        x_gl = np.array([1, 0, 0], dtype=dtype)

        E = self.E
        G = self.G
        

        num_elems = self.ny - 1
        for ielem in range(num_elems):
            P0 = nodes[ielem, :]
            P1 = nodes[ielem + 1, :]
            L = norm(P1 - P0)

            x_loc = unit(P1 - P0)
            y_loc = unit(np.cross(x_loc, x_gl))
            z_loc = unit(np.cross(x_loc, y_loc))

            T[0, :] = x_loc
            T[1, :] = y_loc
            T[2, :] = z_loc

            u0x, u0y, u0z = T.dot(disp[ielem, :3])
            r0x, r0y, r0z = T.dot(disp[ielem, 3:])
            u1x, u1y, u1z = T.dot(disp[ielem + 1, :3])
            r1x, r1y, r1z = T.dot(disp[ielem + 1, 3:])

            # # this is stress = modulus * strain; positive is tensile
            # axial_stress = E * (u1x - u0x) / L

            # # this is Torque / (2 * thickness_min * Area_enclosed)
            # torsion_stress = G * J[ielem] / L * (r1x - r0x) / 2 / spar_thickness[ielem] / A_enc[ielem]

            # # this is moment * h / I
            # top_bending_stress = E / (L**2) * (6 * u0y + 2 * r0z * L - 6 * u1y + 4 * r1z * L) * htop[ielem]

            # # this is moment * h / I
            # bottom_bending_stress = -E / (L**2) * (6 * u0y + 2 * r0z * L - 6 * u1y + 4 * r1z * L) * hbottom[ielem]

            # # this is moment * h / I
            # front_bending_stress = -E / (L**2) * (-6 * u0z + 2 * r0y * L + 6 * u1z + 4 * r1y * L) * hfront[ielem]

            # # this is moment * h / I
            # rear_bending_stress = E / (L**2) * (-6 * u0z + 2 * r0y * L + 6 * u1z + 4 * r1y * L) * hrear[ielem]

            # # shear due to bending (VQ/It) note: the I used to get V cancels the other I
            # vertical_shear = (
            #     E
            #     / (L**3)
            #     * (-12 * u0y - 6 * r0z * L + 12 * u1y - 6 * r1z * L)
            #     * Qy[ielem]
            #     / (2 * spar_thickness[ielem])
            # )

            # ==============================================================================
            # strain equations
            # ==============================================================================
            # this is stress = modulus * strain; positive is tensile
            axial_strain = (u1x - u0x) / L

            # this is Torque / (2 * thickness_min * Area_enclosed)
            torsion_shear_strain = J[ielem] / L * (r1x - r0x) / 2 / spar_thickness[ielem] / A_enc[ielem]

            # this is moment * h / I
            top_bending_strain = 1.0 / (L**2) * (6 * u0y + 2 * r0z * L - 6 * u1y + 4 * r1z * L) * htop[ielem]

            # this is moment * h / I
            bottom_bending_strain = -1.0 / (L**2) * (6 * u0y + 2 * r0z * L - 6 * u1y + 4 * r1z * L) * hbottom[ielem]

            # this is moment * h / I
            front_bending_strain = -1.0 / (L**2) * (-6 * u0z + 2 * r0y * L + 6 * u1z + 4 * r1y * L) * hfront[ielem]

            # this is moment * h / I
            rear_bending_strain = 1.0 / (L**2) * (-6 * u0z + 2 * r0y * L + 6 * u1z + 4 * r1y * L) * hrear[ielem]

            # shear due to bending (VQ/It) note: the I used to get V cancels the other I
            vertical_shear_strain = (
                1.0
                / (L**3)
                * (-12 * u0y - 6 * r0z * L + 12 * u1y - 6 * r1z * L)
                * Qy[ielem]
                / (2 * spar_thickness[ielem])
            )

            # print("==========",ielem,"================")
            # print("vertical_shear", vertical_shear)
            # print("top",top_bending_stress)
            # print("bottom",bottom_bending_stress)
            # print("front",front_bending_stress)
            # print("rear",rear_bending_stress)
            # print("axial", axial_stress)
            # print("torsion", torsion_stress)

            # # The 4 stress combinations:
            # vonmises[ielem, 0] = (
            #     np.sqrt((top_bending_stress + rear_bending_stress + axial_stress) ** 2 + 3 * torsion_stress**2)
            #     / self.tssf
            # )
            # vonmises[ielem, 1] = np.sqrt(
            #     (bottom_bending_stress + front_bending_stress + axial_stress) ** 2 + 3 * torsion_stress**2
            # )
            # vonmises[ielem, 2] = np.sqrt(
            #     (front_bending_stress + axial_stress) ** 2 + 3 * (torsion_stress - vertical_shear) ** 2
            # )
            # vonmises[ielem, 3] = (
            #     np.sqrt((rear_bending_stress + axial_stress) ** 2 + 3 * (torsion_stress + vertical_shear) ** 2)
            #     / self.tssf
            # )

            # The strain combinations for the 4 elements under consideration: (split into epsilonx, epsilony and gammatau)
            # Defining the epsilon_elem array for the epsionx, epsiony and gammatau for each element
            epsilon_elem = np.zeros((4, 3), dtype=dtype)

            # Element 1:
            epsilon_elem[0, 0] = top_bending_strain + rear_bending_strain + axial_strain
            epsilon_elem[0, 1] = 0
            epsilon_elem[0, 2] = torsion_shear_strain

            # Element 2:
            epsilon_elem[1, 0] = bottom_bending_strain + front_bending_strain + axial_strain
            epsilon_elem[1, 1] = 0
            epsilon_elem[1, 2] = torsion_shear_strain

            # Element 3:
            epsilon_elem[2, 0] = front_bending_strain + axial_strain
            epsilon_elem[2, 1] = 0
            epsilon_elem[2, 2] = torsion_shear_strain + vertical_shear_strain # NOTE: To be verified again

            # Element 4:
            epsilon_elem[3, 0] = rear_bending_strain + axial_strain
            epsilon_elem[3, 1] = 0
            epsilon_elem[3, 2] = - torsion_shear_strain + vertical_shear_strain # NOTE: To be verified again

            # defining the array for ply-orientation angles:
            theta = np.array([0, 90, 45, -45])

            # defining the epsilon_elem_ply array for the epsilon1, epsilon2 and gamma12 for each ply
            epsilon_elem_ply = np.zeros((4, 4, 3), dtype=dtype)
            sigma_elem_ply = np.zeros((4, 4, 3), dtype=dtype)

            # running a loop over the 4 elements and 4 plies to calculate the epsilon_elem_ply array
            for elem_num in range(4):
                for ply_num in range(4):
                    epsilon_elem_ply[elem_num, ply_num, 0] = epsilon_elem[elem_num, 0] * np.cos(np.radians(theta[ply_num]))**2 + epsilon_elem[elem_num, 1] * np.sin(np.radians(theta[ply_num]))**2 + 2 * epsilon_elem[elem_num, 2] * np.sin(np.radians(theta[ply_num])) * np.cos(np.radians(theta[ply_num]))
                    epsilon_elem_ply[elem_num, ply_num, 1] = epsilon_elem[elem_num, 0] * np.sin(np.radians(theta[ply_num]))**2 + epsilon_elem[elem_num, 1] * np.cos(np.radians(theta[ply_num]))**2 - 2 * epsilon_elem[elem_num, 2] * np.sin(np.radians(theta[ply_num])) * np.cos(np.radians(theta[ply_num]))
                    epsilon_elem_ply[elem_num, ply_num, 2] = - epsilon_elem[elem_num, 0] * np.sin(np.radians(theta[ply_num])) * np.cos(np.radians(theta[ply_num])) + epsilon_elem[elem_num, 1] * np.sin(np.radians(theta[ply_num])) * np.cos(np.radians(theta[ply_num])) + epsilon_elem[elem_num, 2] * (np.cos(np.radians(theta[ply_num]))**2 - np.sin(np.radians(theta[ply_num]))**2)

            # defining the strain-stress relations for the material:
            E1 = 117.7e9
            E2 = 9.7e9
            G12 = 4.8e9
            nu_12 = 0.35
            nu_21 = nu_12 * E2 / E1

            Q11 = E1 / (1 - nu_12 * nu_21)
            Q22 = E2 / (1 - nu_12 * nu_21)
            Q12 = nu_12 * E2 / (1 - nu_12 * nu_21)
            Q21 = Q12
            Q66 = G12

            # converting the strains to stresses using strain-stress relations
            for elem_num in range(4):
                for ply_num in range(4):
                    epsilon1 = epsilon_elem_ply[elem_num, ply_num, 0]
                    epsilon2 = epsilon_elem_ply[elem_num, ply_num, 1]
                    gamma12 = epsilon_elem_ply[elem_num, ply_num, 2]

                    sigma1 = Q11 * epsilon1 + Q12 * epsilon2
                    sigma2 = Q21 * epsilon1 + Q22 * epsilon2
                    sigma12 = Q66 * gamma12

                    # assigning the stresses to the sigma_elem_ply array
                    sigma_elem_ply[elem_num, ply_num, 0] = sigma1
                    sigma_elem_ply[elem_num, ply_num, 1] = sigma2
                    sigma_elem_ply[elem_num, ply_num, 2] = sigma12


            # defining the Tsai-Wu constants for the material:
            sigma_c1 = 1034.0e6
            sigma_t1 = 1648.0e6
            sigma_c2 = 228.0e6
            sigma_t2 = 64.0e6
            sigma_12_max = 71.0e6

            F1 =  1 / sigma_t1 - 1 / sigma_c1 # NOTE: To be verified again
            F11 = 1 / (sigma_t1 * sigma_c1)
            F2 =  1 / sigma_t2 - 1 / sigma_c2 # NOTE: To be verified again
            F22 = 1 / (sigma_t2 * sigma_c2)
            F66 = 1 / (sigma_12_max**2)

            # Find the Tsai-Wu Strength Ratios for each ply in each element and store them in the tsaiwu_sr array
            for elem_num in range(4):
                for ply_num in range(4):
                    a = F1 * sigma_elem_ply[elem_num, ply_num, 0] + F2 * sigma_elem_ply[elem_num, ply_num, 1]
                    b = F11 * sigma_elem_ply[elem_num, ply_num, 0]**2 + F22 * sigma_elem_ply[elem_num, ply_num, 1]**2 + F66 * sigma_elem_ply[elem_num, ply_num, 2]**2
                    tsaiwu_sr[ielem, elem_num * 4 + ply_num] = 0.5 * (a + np.sqrt(a**2 + 4 * b))

                    # a = F1 * epsilon_elem_ply[elem_num, ply_num, 0] + F2 * epsilon_elem_ply[elem_num, ply_num, 1]
                    # b = F11 * epsilon_elem_ply[elem_num, ply_num, 0]**2 + F22 * epsilon_elem_ply[elem_num, ply_num, 1]**2 + F66 * epsilon_elem_ply[elem_num, ply_num, 0] * epsilon_elem_ply[elem_num, ply_num, 1]
                    # tsaiwu_sr[ielem, elem_num * 4 + ply_num] = 0.5 * (a + np.sqrt(a**2 + 4 * b))
