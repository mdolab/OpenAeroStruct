import numpy as np

import openmdao.api as om


class SectionPropertiesWingbox(om.ExplicitComponent):
    """
    Compute geometric cross-section properties for the wingbox elements.
    See Chauhan et al. (https://doi.org/10.1007/978-3-319-97773-7_38) for more.

    Parameters
    ----------
    streamwise_chords[ny-1] : numpy array
        Average streamwise chord lengths for each streamwise VLM panel.
    fem_chords[ny-1] : numpy array
        Effective chord lengths normal to the FEM elements.
    fem_twists[ny-1] : numpy array
        Twist angles in planes normal to the FEM elements.
    spar_thickness[ny-1] : numpy array
        Material thicknesses of the front and rear spars for each wingbox segment.
    skin_thickness[ny-1] : numpy array
        Material thicknesses of the top and bottom skins for each wingbox segment.
    t_over_c[ny-1] : numpy array
        Streamwise thickness-to-chord ratios for each wingbox segment.


    Returns
    -------
    A[ny-1] : numpy array
        Cross-sectional area of each wingbox segment.
    A_enc[ny-1] : numpy array
        Cross-sectional enclosed area (measured using the material midlines) of
        each wingbox segment.
    A_int[ny-1] : numpy array
        Cross-sectional internal area of each wingbox segment (used for fuel
        volume).
    Iy[ny-1] : numpy array
        Second moment of area about the neutral axis parallel to the local
        y-axis (for each wingbox segment).
    Qz[ny-1] : numpy array
        First moment of area above the neutral axis parallel to the local
        z-axis (for each wingbox segment).
    Iz[ny-1] : numpy array
        Second moment of area about the neutral axis parallel to the local
        z-axis (for each wingbox segment).
    J[ny-1] : numpy array
        Torsion constants for each wingbox segment.
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
    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        self.surface = surface = self.options['surface']

        self.mesh = surface['mesh']
        self.ny = self.mesh.shape[1]

        # original thickness-to-chord ratio of the airfoil provided by the user
        self.orig_wb_af_t_over_c = surface['original_wingbox_airfoil_t_over_c']

        # airfoil coordinates provided by the user
        self.data_x_upper = surface['data_x_upper']
        self.data_x_lower = surface['data_x_lower']
        self.data_y_upper = surface['data_y_upper']
        self.data_y_lower = surface['data_y_lower']

        self.add_input('streamwise_chords', val=np.ones((self.ny - 1)),units='m')
        self.add_input('fem_chords', val=np.ones((self.ny - 1)),units='m')
        self.add_input('fem_twists', val=np.ones((self.ny - 1)),units='deg')

        self.add_input('spar_thickness', val=np.ones((self.ny - 1)),units='m')
        self.add_input('skin_thickness', val=np.ones((self.ny - 1)),units='m')
        self.add_input('t_over_c', val=np.ones((self.ny - 1)))

        self.add_output('A', val=np.ones((self.ny - 1)),units='m**2')
        self.add_output('A_enc', val=np.ones((self.ny - 1)),units='m**2')
        self.add_output('A_int', val=np.ones((self.ny - 1)),units='m**2')
        self.add_output('Iy', val=np.ones((self.ny - 1)),units='m**4')
        self.add_output('Qz', val=np.ones((self.ny - 1)),units='m**3')
        self.add_output('Iz', val=np.ones((self.ny - 1)),units='m**4')
        self.add_output('J', val=np.ones((self.ny - 1)),units='m**4')
        self.add_output('htop', val=np.ones((self.ny - 1)),units='m')
        self.add_output('hbottom', val=np.ones((self.ny - 1)),units='m')
        self.add_output('hfront', val=np.ones((self.ny - 1)),units='m')
        self.add_output('hrear', val=np.ones((self.ny - 1)),units='m')

        self.declare_partials('*', '*', method='cs')

    def compute(self, inputs, outputs):

        # NOTE: In the code below, the x- and y-axes correspond to the element
        # local z- and y-axes, respectively.

        chord = inputs['fem_chords']
        spar_thickness = inputs['spar_thickness']
        skin_thickness = inputs['skin_thickness']
        t_over_c_original = self.orig_wb_af_t_over_c
        t_over_c = inputs['t_over_c']
        streamwise_chord = inputs['streamwise_chords']
        theta = inputs['fem_twists']

        # Scale data points with chord
        data_x_upper = np.outer(self.data_x_upper, chord)
        data_y_upper = np.outer(self.data_y_upper, chord)
        data_x_lower = np.outer(self.data_x_lower, chord)
        data_y_lower = np.outer(self.data_y_lower, chord)

        # Scale y-coordinates by t/c design variable which is streamwise t/c
        data_y_upper *= t_over_c / t_over_c_original * streamwise_chord / chord
        data_y_lower *= t_over_c / t_over_c_original * streamwise_chord / chord

        # Compute enclosed area for torsion constant
        # This currently does not change with twist
        # Also compute internal area for internal volume calculation for fuel
        x_up_diff = data_x_upper[1:] - data_x_upper[:-1]
        x_low_diff = data_x_lower[1:] - data_x_lower[:-1]
        y_up_diff = data_y_upper[1:] - data_y_upper[:-1]
        y_low_diff = data_y_lower[1:] - data_y_lower[:-1]

        y_up_add = data_y_upper[1:] + data_y_upper[:-1]
        y_low_add = data_y_lower[1:] + data_y_lower[:-1]

        A_enc = (x_up_diff) * (y_up_add - skin_thickness ) / 2 # area above 0 line
        A_enc += (x_low_diff) * (-y_low_add - skin_thickness ) / 2 # area below 0 line
        A_int = (x_up_diff) * (y_up_add - 2*skin_thickness ) / 2 # area above 0 line
        A_int += (x_low_diff) * (-y_low_add - 2*skin_thickness ) / 2 # area below 0 line

        A_enc = np.sum(A_enc, axis=0)
        A_int = np.sum(A_int, axis=0)

        A_enc -= (data_y_upper[0] - data_y_lower[0]) * spar_thickness / 2 # area of spars
        A_enc -= (data_y_upper[-1] - data_y_lower[-1]) * spar_thickness / 2 # area of spars
        A_int -= (data_y_upper[0] - data_y_lower[0]) * spar_thickness # area of spars
        A_int -= (data_y_upper[-1] - data_y_lower[-1]) * spar_thickness # area of spars

        outputs['A_enc'] = A_enc
        outputs['A_int'] = A_int

        # Compute perimeter to thickness ratio for torsion constant
        # This currently does not change with twist
        p_by_t_1 = (x_up_diff**2 + y_up_diff**2)**0.5 / skin_thickness # length / thickness of caps
        p_by_t_2 = (x_low_diff**2 + y_low_diff**2)**0.5 / skin_thickness # length / thickness of caps

        p_by_t = np.sum(p_by_t_1 + p_by_t_2, axis=0)

        p_by_t += (data_y_upper[0] - data_y_lower[0] - skin_thickness) / spar_thickness # length / thickness of spars
        p_by_t += (data_y_upper[-1] - data_y_lower[-1] - skin_thickness) / spar_thickness # length / thickness of spars

        # Torsion constant
        J = 4 * A_enc**2 / p_by_t

        outputs['J'] = J

        # Rotate the wingbox
        rot_mat = np.array([[np.cos(theta), np.sin(theta)],
                            [-np.sin(theta), np.cos(theta)]])

        data_x_upper_2 = rot_mat[0, 0] * data_x_upper + rot_mat[0, 1] * data_y_upper
        data_y_upper_2 = rot_mat[1, 0] * data_x_upper + rot_mat[1, 1] * data_y_upper
        data_x_lower_2 = rot_mat[0, 0] * data_x_lower + rot_mat[0, 1] * data_y_lower
        data_y_lower_2 = rot_mat[1, 0] * data_x_lower + rot_mat[1, 1] * data_y_lower

        data_x_upper = data_x_upper_2.copy()
        data_y_upper = data_y_upper_2.copy()
        data_x_lower = data_x_lower_2.copy()
        data_y_lower = data_y_lower_2.copy()

        x_up_diff = data_x_upper[1:] - data_x_upper[:-1]
        x_low_diff = data_x_lower[1:] - data_x_lower[:-1]
        y_up_diff = data_y_upper[1:] - data_y_upper[:-1]
        y_low_diff = data_y_lower[1:] - data_y_lower[:-1]

        y_up_add = data_y_upper[1:] + data_y_upper[:-1]
        y_low_add = data_y_lower[1:] + data_y_lower[:-1]

        # Compute area moment of inertia about x axis
        # First compute centroid and area
        first_moment_area_upper = (y_up_add / 2 - (skin_thickness/2) ) * skin_thickness * x_up_diff
        upper_area = skin_thickness * x_up_diff

        first_moment_area_lower = (y_low_add / 2 + (skin_thickness/2) ) * skin_thickness * x_low_diff
        lower_area = skin_thickness * x_low_diff

        first_moment_area_front_spar = (data_y_upper[0] - data_y_lower[0] - 2 * skin_thickness) * spar_thickness * (data_y_upper[0] + data_y_lower[0]) / 2
        first_moment_area_rear_spar = (data_y_upper[-1] - data_y_lower[-1] - 2 * skin_thickness) * spar_thickness * (data_y_upper[-1] + data_y_lower[-1]) / 2
        area_spars = ((data_y_upper[0] - data_y_lower[0] - 2 * skin_thickness) + (data_y_upper[-1] - data_y_lower[-1] - 2 * skin_thickness)) * spar_thickness

        area = np.sum(upper_area, axis=0) + np.sum(lower_area, axis=0) + area_spars

        outputs['A'] = area

        centroid = (np.sum(first_moment_area_upper, axis=0) + np.sum(first_moment_area_lower, axis=0) + first_moment_area_front_spar + first_moment_area_rear_spar) / area

        # Then compute area moment of inertia for upward bending
        # This is calculated using derived analytical expression assuming linear interpolation between airfoil data points
        a = y_up_diff / x_up_diff
        b = (y_up_diff + skin_thickness) / 2
        x2 = x_up_diff

        I_horiz_1 = 2 * ((1./12. * a**3 * x2**4 + 1./3. * a**2 * x2**3 * b + 1./2. * a * x2**2 * b**2 + 1./3. * b**3 * x2))
        I_horiz_2 = x2 * skin_thickness * (y_up_add / 2 - skin_thickness / 2 - centroid)**2

        I_horiz = np.sum(I_horiz_1 + I_horiz_2, axis=0)

        # Compute area moment of inertia about y axis
        a = -y_low_diff / x_low_diff
        b = (-y_low_diff + skin_thickness) / 2
        x2 = x_low_diff

        I_horiz += np.sum(2 * ((1./12. * a**3 * x2**4 + 1./3. * a**2 * x2**3 * b + 1./2. * a * x2**2 * b**2 + 1./3. * b**3 * x2)), axis=0)
        I_horiz += np.sum(x2 * skin_thickness * ((-y_low_add)/2 - skin_thickness/2 + centroid)**2, axis=0)

        # Contribution from the forward spar
        I_horiz += 1./12. * spar_thickness * (data_y_upper[0] - data_y_lower[0] - 2 * skin_thickness)**3 + spar_thickness * (data_y_upper[0] - data_y_lower[0] - 2 * skin_thickness) * ((data_y_upper[0] + data_y_lower[0]) / 2 - centroid)**2
        # Contribution from the rear spar
        I_horiz += 1./12. * spar_thickness * (data_y_upper[-1] - data_y_lower[-1] - 2 * skin_thickness)**3 + spar_thickness * (data_y_upper[-1] - data_y_lower[-1] - 2 * skin_thickness) * ((data_y_upper[-1] + data_y_lower[-1]) / 2 - centroid)**2

        outputs['Iz'] = I_horiz

        # Compute the Q required for transverse shear stress due to upward bending
        Q_upper = np.sum(((y_up_add / 2 - (skin_thickness/2) ) - centroid) * skin_thickness * x_up_diff, axis=0)

        Q_upper += (data_y_upper[0] - skin_thickness - centroid)**2 / 2 * (spar_thickness)
        Q_upper += (data_y_upper[-1] - skin_thickness - centroid)**2 / 2 * (spar_thickness)

        outputs['Qz'] = Q_upper

        # Compute area moment of inertia for backward bending
        I_vert = 0
        first_moment_area_front = (data_y_upper[0] - data_y_lower[0]) * spar_thickness * (data_x_upper[0] + spar_thickness / 2)
        first_moment_area_rear = (data_y_upper[-1] - data_y_lower[-1]) * spar_thickness * (data_x_upper[-1] - spar_thickness / 2)
        centroid_Ivert = (first_moment_area_front + first_moment_area_rear) / \
                        ( ((data_y_upper[0] - data_y_lower[0]) + (data_y_upper[-1] - data_y_lower[-1])) * spar_thickness)

        I_vert += 1./12. * (data_y_upper[0] - data_y_lower[0]) * spar_thickness**3 + (data_y_upper[0] - data_y_lower[0]) * spar_thickness * (centroid_Ivert - (data_x_upper[0] + spar_thickness/2))**2
        I_vert += 1./12. * (data_y_upper[-1] - data_y_lower[-1]) * spar_thickness**3 + (data_y_upper[-1] - data_y_lower[-1]) * spar_thickness * (data_x_upper[-1] - spar_thickness/2 - centroid_Ivert)**2

        # Add contribution of skins
        I_vert += 2 * ( 1./12. * skin_thickness * (data_x_upper[-1] - data_x_upper[0] - 2 * spar_thickness)**3 + skin_thickness * (data_x_upper[-1] - data_x_upper[0] - 2 * spar_thickness) * (centroid_Ivert - (data_x_upper[-1] + data_x_upper[0]) / 2)**2 )

        outputs['Iy'] = I_vert

        # Distances for calculating max bending stresses (KS function used)
        ks_rho = 500. # Hard coded, see Martins and Poon 2005 for more
        fmax_upper = np.max(data_y_upper, axis=0)
        htop = fmax_upper + 1 / ks_rho * np.log(np.sum(np.exp(ks_rho * (data_y_upper - fmax_upper)), axis=0)) - centroid

        fmax_lower = np.max(-data_y_lower, axis=0)
        hbottom = fmax_lower + 1 / ks_rho * np.log(np.sum(np.exp(ks_rho * (-data_y_lower - fmax_lower)), axis=0)) + centroid

        hfront =  centroid_Ivert - data_x_upper[0]
        hrear = data_x_upper[-1] - centroid_Ivert

        outputs['htop'] = htop
        outputs['hbottom'] = hbottom
        outputs['hfront'] = hfront
        outputs['hrear'] = hrear
