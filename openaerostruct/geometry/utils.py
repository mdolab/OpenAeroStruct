import numpy as np
from numpy import cos, sin, tan

from openaerostruct.geometry.CRM_definitions import get_crm_points


def rotate(mesh, theta_y, symmetry, rotate_x=True):
    """
    Compute rotation matrices given mesh and rotation angles in degrees.

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the initial aerodynamic surface.
    theta_y[ny] : numpy array
        1-D array of rotation angles about y-axis for each wing slice in degrees.
    symmetry : boolean
        Flag set to True if surface is reflected about y=0 plane.
    rotate_x : boolean
        Flag set to True if the user desires the twist variable to always be
        applied perpendicular to the wing (say, in the case of a winglet).

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the twisted aerodynamic surface.

    """
    te = mesh[-1]
    le = mesh[ 0]
    quarter_chord = 0.25 * te + 0.75 * le

    nx, ny, _ = mesh.shape

    if rotate_x:
        # Compute spanwise z displacements along quarter chord
        if symmetry:
            dz_qc = quarter_chord[:-1,2] - quarter_chord[1:,2]
            dy_qc = quarter_chord[:-1,1] - quarter_chord[1:,1]
            theta_x = np.arctan(dz_qc/dy_qc)

            # Prepend with 0 so that root is not rotated
            rad_theta_x = np.append(theta_x, 0.0)
        else:
            root_index = int((ny - 1) / 2)
            dz_qc_left = quarter_chord[:root_index,2] - quarter_chord[1:root_index+1,2]
            dy_qc_left = quarter_chord[:root_index,1] - quarter_chord[1:root_index+1,1]
            theta_x_left = np.arctan(dz_qc_left/dy_qc_left)
            dz_qc_right = quarter_chord[root_index+1:,2] - quarter_chord[root_index:-1,2]
            dy_qc_right = quarter_chord[root_index+1:,1] - quarter_chord[root_index:-1,1]
            theta_x_right = np.arctan(dz_qc_right/dy_qc_right)

            # Concatenate thetas
            rad_theta_x = np.concatenate((theta_x_left, np.zeros(1), theta_x_right))

    else:
        rad_theta_x = 0.0

    rad_theta_y = theta_y * np.pi / 180.

    mats = np.zeros((ny, 3, 3), dtype=type(rad_theta_y[0]))

    cos_rtx = cos(rad_theta_x)
    cos_rty = cos(rad_theta_y)
    sin_rtx = sin(rad_theta_x)
    sin_rty = sin(rad_theta_y)

    mats[:, 0, 0] = cos_rty
    mats[:, 0, 2] = sin_rty
    mats[:, 1, 0] = sin_rtx * sin_rty
    mats[:, 1, 1] = cos_rtx
    mats[:, 1, 2] = -sin_rtx * cos_rty
    mats[:, 2, 0] = -cos_rtx * sin_rty
    mats[:, 2, 1] = sin_rtx
    mats[:, 2, 2] = cos_rtx*cos_rty

    mesh[:] = np.einsum("ikj, mij -> mik", mats, mesh - quarter_chord) + quarter_chord


def scale_x(mesh, chord_dist):
    """
    Modify the chords along the span of the wing by scaling only the x-coord.

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the initial aerodynamic surface.
    chord_dist[ny] : numpy array
        Chord length for each panel edge.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh with the new chord lengths.
    """
    te = mesh[-1]
    le = mesh[ 0]
    quarter_chord = 0.25 * te + 0.75 * le

    ny = mesh.shape[1]

    for i in range(ny):
        mesh[:, i, 0] = (mesh[:, i, 0] - quarter_chord[i, 0]) * chord_dist[i] + \
            quarter_chord[i, 0]

def shear_x(mesh, xshear):
    """
    Shear the wing in the x direction (distributed sweep).

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the initial aerodynamic surface.
    xshear[ny] : numpy array
        Distance to translate wing in x direction.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh with the new chord lengths.
    """
    mesh[:, :, 0] += xshear

def shear_y(mesh, yshear):
    """ Shear the wing in the y direction (distributed span).

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the initial aerodynamic surface.
    yshear[ny] : numpy array
        Distance to translate wing in y direction.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh with the new span widths.
    """
    mesh[:, :, 1] += yshear

def shear_z(mesh, zshear):
    """
    Shear the wing in the z direction (distributed dihedral).

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the initial aerodynamic surface.
    zshear[ny] : numpy array
        Distance to translate wing in z direction.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh with the new chord lengths.
    """
    mesh[:, :, 2] += zshear

def sweep(mesh, sweep_angle, symmetry):
    """
    Apply shearing sweep. Positive sweeps back.

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the initial aerodynamic surface.
    sweep_angle : float
        Shearing sweep angle in degrees.
    symmetry : boolean
        Flag set to true if surface is reflected about y=0 plane.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the swept aerodynamic surface.

    """

    # Get the mesh parameters and desired sweep angle
    num_x, num_y, _ = mesh.shape
    le = mesh[0]
    p180 = np.pi / 180
    tan_theta = tan(p180*sweep_angle)

    # If symmetric, simply vary the x-coord based on the distance from the
    # center of the wing
    if symmetry:
        y0 = le[-1, 1]
        dx = -(le[:, 1] - y0) * tan_theta

    # Else, vary the x-coord on either side of the wing
    else:
        ny2 = (num_y - 1) // 2
        y0 = le[ny2, 1]

        dx_right = (le[ny2:, 1] - y0) * tan_theta
        dx_left = -(le[:ny2, 1] - y0) * tan_theta
        dx = np.hstack((dx_left, dx_right))

    # dx added spanwise.
    mesh[:, :, 0] += dx

def dihedral(mesh, dihedral_angle, symmetry):
    """
    Apply dihedral angle. Positive angles up.

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the initial aerodynamic surface.
    dihedral_angle : float
        Dihedral angle in degrees.
    symmetry : boolean
        Flag set to true if surface is reflected about y=0 plane.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the aerodynamic surface with dihedral angle.

    """

    # Get the mesh parameters and desired sweep angle
    num_x, num_y, _ = mesh.shape
    le = mesh[0]
    p180 = np.pi / 180
    tan_theta = tan(p180*dihedral_angle)

    # If symmetric, simply vary the z-coord based on the distance from the
    # center of the wing
    if symmetry:
        y0 = le[-1, 1]
        dz = -(le[:, 1] - y0) * tan_theta

    else:
        ny2 = (num_y-1) // 2
        y0 = le[ny2, 1]
        dz_right = (le[ny2:, 1] - y0) * tan_theta
        dz_left = -(le[:ny2, 1] - y0) * tan_theta
        dz = np.hstack((dz_left, dz_right))

    # dz added spanwise.
    mesh[:, :, 2] += dz


def stretch(mesh, span, symmetry):
    """
    Stretch mesh in spanwise direction to reach specified span.

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the initial aerodynamic surface.
    span : float
        Relative stetch ratio in the spanwise direction.
    symmetry : boolean
        Flag set to true if surface is reflected about y=0 plane.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the stretched aerodynamic surface.

    """

    # Set the span along the quarter-chord line
    le = mesh[0]
    te = mesh[-1]
    quarter_chord = 0.25 * te + 0.75 * le

    # The user always deals with the full span, so if they input a specific
    # span value and have symmetry enabled, we divide this value by 2.
    if symmetry:
        span /= 2.

    # Compute the previous span and determine the scalar needed to reach the
    # desired span
    prev_span = quarter_chord[-1, 1] - quarter_chord[0, 1]
    s = quarter_chord[:,1] / prev_span
    mesh[:, :, 1] = s * span

def taper(mesh, taper_ratio, symmetry):
    """
    Alter the spanwise chord linearly to produce a tapered wing. Note that
    we apply taper around the quarter-chord line.

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the initial aerodynamic surface.
    taper_ratio : float
        Taper ratio for the wing; 1 is untapered, 0 goes to a point.
    symmetry : boolean
        Flag set to true if surface is reflected about y=0 plane.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the tapered aerodynamic surface.

    """

    # Get mesh parameters and the quarter-chord
    le = mesh[0]
    te = mesh[-1]
    num_x, num_y, _ = mesh.shape
    quarter_chord = 0.25 * te + 0.75 * le
    x = quarter_chord[:, 1]
    span = x[-1] - x[0]

    # If symmetric, solve for the correct taper ratio, which is a linear
    # interpolation problem
    if symmetry:
        xp = np.array([-span, 0.])
        fp = np.array([taper_ratio, 1.])

    # Otherwise, we set up an interpolation problem for the entire wing, which
    # consists of two linear segments
    else:
        xp = np.array([-span/2, 0., span/2])
        fp = np.array([taper_ratio, 1., taper_ratio])

    taper = np.interp(x.real, xp.real, fp.real)

    # Modify the mesh based on the taper amount computed per spanwise section
    mesh[:] = np.einsum('ijk, j->ijk', mesh - quarter_chord, taper) + quarter_chord


def gen_rect_mesh(num_x, num_y, span, chord, span_cos_spacing=0., chord_cos_spacing=0.):
    """
    Generate simple rectangular wing mesh.

    Parameters
    ----------
    num_x : float
        Desired number of chordwise node points for the final mesh.
    num_y : float
        Desired number of chordwise node points for the final mesh.
    span : float
        Total wingspan.
    chord : float
        Root chord.
    span_cos_spacing : float (optional)
        Blending ratio of uniform and cosine spacing in the spanwise direction.
        A value of 0. corresponds to uniform spacing and a value of 1.
        corresponds to regular cosine spacing. This increases the number of
        spanwise node points near the wingtips.
    chord_cos_spacing : float (optional)
        Blending ratio of uniform and cosine spacing in the chordwise direction.
        A value of 0. corresponds to uniform spacing and a value of 1.
        corresponds to regular cosine spacing. This increases the number of
        chordwise node points near the wingtips.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Rectangular nodal mesh defining the final aerodynamic surface with the
        specified parameters.
    """

    mesh = np.zeros((num_x, num_y, 3))
    ny2 = (num_y + 1) // 2

    # Hotfix a special case for spacing bunched at the root and tips
    if span_cos_spacing == 2.:
        beta = np.linspace(0, np.pi, ny2)

        # mixed spacing with span_cos_spacing as a weighting factor
        # this is for the spanwise spacing
        cosine = .25 * (1 - np.cos(beta)) # cosine spacing
        uniform = np.linspace(0, .5, ny2)[::-1]  # uniform spacing
        half_wing = cosine[::-1] * span_cos_spacing + (1 - span_cos_spacing) * uniform
        full_wing = np.hstack((-half_wing[:-1], half_wing[::-1])) * span

    else:
        beta = np.linspace(0, np.pi/2, ny2)

        # mixed spacing with span_cos_spacing as a weighting factor
        # this is for the spanwise spacing
        cosine = .5 * np.cos(beta)  # cosine spacing
        uniform = np.linspace(0, .5, ny2)[::-1]  # uniform spacing
        half_wing = cosine * span_cos_spacing + (1 - span_cos_spacing) * uniform
        full_wing = np.hstack((-half_wing[:-1], half_wing[::-1])) * span

    nx2 = (num_x + 1) // 2
    beta = np.linspace(0, np.pi/2, nx2)

    # mixed spacing with span_cos_spacing as a weighting factor
    # this is for the chordwise spacing
    cosine = .5 * np.cos(beta)  # cosine spacing
    uniform = np.linspace(0, .5, nx2)[::-1]  # uniform spacing
    half_wing = cosine * chord_cos_spacing + (1 - chord_cos_spacing) * uniform
    full_wing_x = np.hstack((-half_wing[:-1], half_wing[::-1])) * chord

    # Special case if there are only 2 chordwise nodes
    if num_x <= 2:
        full_wing_x = np.array([0., chord])

    for ind_x in range(num_x):
        for ind_y in range(num_y):
            mesh[ind_x, ind_y, :] = [full_wing_x[ind_x], full_wing[ind_y], 0]

    return mesh


def gen_crm_mesh(num_x, num_y, span_cos_spacing=0., chord_cos_spacing=0., wing_type="CRM:jig"):
    """
    Generate Common Research Model wing mesh.

    Parameters
    ----------
    num_x : float
        Desired number of chordwise node points for the final mesh.
    num_y : float
        Desired number of chordwise node points for the final mesh.
    span : float
        Total wingspan.
    chord : float
        Root chord.
    span_cos_spacing : float (optional)
        Blending ratio of uniform and cosine spacing in the spanwise direction.
        A value of 0. corresponds to uniform spacing and a value of 1.
        corresponds to regular cosine spacing. This increases the number of
        spanwise node points near the wingtips.
    chord_cos_spacing : float (optional)
        Blending ratio of uniform and cosine spacing in the chordwise direction.
        A value of 0. corresponds to uniform spacing and a value of 1.
        corresponds to regular cosine spacing. This increases the number of
        chordwise node points near the wingtips.
    wing_type : string (optional)
        Describes the desired CRM shape. Current options are:
        "CRM:jig" (undeformed jig shape),
        "CRM:alpha_2.75" (shape from wind tunnel testing at a=2.75 from DPW6)

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Rectangular nodal mesh defining the final aerodynamic surface with the
        specified parameters.
    eta : numpy array
        Spanwise locations of the airfoil slices. Later used in the
        interpolation function to obtain correct twist values at
        points along the span that are not aligned with these slices.
    twist : numpy array
        Twist along the span at the spanwise eta locations. We use these twists
        as training points for interpolation to obtain twist values at
        arbitrary points along the span.

    """

    # Call an external function to get the data points for the specific CRM
    # type requested. See `CRM_definitions.py` for more information and the
    # raw data.
    raw_crm_points = get_crm_points(wing_type)

    # If this is a jig shape, remove all z-deflection to create a
    # poor person's version of the undeformed CRM.
    if 'jig' in wing_type or 'CRM' == wing_type:
        raw_crm_points[:, 3] = 0.

    # Get the leading edge of the raw crm points
    le = np.vstack((raw_crm_points[:,1],
                    raw_crm_points[:,2],
                    raw_crm_points[:,3]))

    # Get the chord, twist(in correct order), and eta values from the points
    chord = raw_crm_points[:, 5]
    twist = raw_crm_points[:, 4][::-1]
    eta = raw_crm_points[:, 0]

    # Get the trailing edge of the crm points, based on the chord + le distance.
    # Note that we do not account for twist here; instead we set that using
    # the twist design variable later in run_classes.py.
    te = np.vstack((raw_crm_points[:,1] + chord,
                       raw_crm_points[:,2],
                       raw_crm_points[:,3]))

    # Get the number of points that define this CRM shape and create a mesh
    # array based on this size
    n_raw_points = raw_crm_points.shape[0]
    mesh = np.empty((2, n_raw_points, 3))

    # Set the leading and trailing edges of the mesh matrix
    mesh[0, :, :] = le.T
    mesh[1, :, :] = te.T

    # Convert the mesh points to meters from inches.
    raw_mesh = mesh * 0.0254

    # Create the blended spacing using the user input for span_cos_spacing
    ny2 = (num_y + 1) // 2
    beta = np.linspace(0, np.pi/2, ny2)

    # Distribution for cosine spacing
    cosine = np.cos(beta)

    # Distribution for uniform spacing
    uniform = np.linspace(0, 1., ny2)[::-1]

    # Combine the two distrubtions using span_cos_spacing as the weighting factor.
    # span_cos_spacing == 1. is for fully cosine, 0. for uniform
    lins = cosine * span_cos_spacing + (1 - span_cos_spacing) * uniform

    # Populate a mesh object with the desired num_y dimension based on
    # interpolated values from the raw CRM points.
    mesh = np.empty((2, ny2, 3))
    for j in range(2):
        for i in range(3):
            mesh[j, :, i] = np.interp(lins[::-1], eta, raw_mesh[j, :, i].real)

    # That is just one half of the mesh and we later expect the full mesh,
    # even if we're using symmetry == True.
    # So here we mirror and stack the two halves of the wing.
    full_mesh = getFullMesh(right_mesh=mesh)

    # If we need to add chordwise panels, do so
    if num_x > 2:
        full_mesh = add_chordwise_panels(full_mesh, num_x, chord_cos_spacing)

    return full_mesh, eta, twist


def add_chordwise_panels(mesh, num_x, chord_cos_spacing):
    """
    Generate a new mesh with multiple chordwise panels.

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the initial aerodynamic surface with only
        the leading and trailing edges defined.
    num_x : float
        Desired number of chordwise node points for the final mesh.
    chord_cos_spacing : float
        Blending ratio of uniform and cosine spacing in the chordwise direction.
        A value of 0. corresponds to uniform spacing and a value of 1.
        corresponds to regular cosine spacing. This increases the number of
        chordwise node points near the wingtips.

    Returns
    -------
    new_mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the final aerodynamic surface with the
        specified number of chordwise node points.

    """

    # Obtain mesh and num properties
    num_y = mesh.shape[1]
    nx2 = (num_x + 1) // 2

    # Create beta, an array of linear sampling points to pi/2
    beta = np.linspace(0, np.pi/2, nx2)

    # Obtain the two spacings that we will use to blend
    cosine = .5 * np.cos(beta)  # cosine spacing
    uniform = np.linspace(0, .5, nx2)[::-1]  # uniform spacing

    # Create half of the wing in the chordwise direction
    half_wing = cosine * chord_cos_spacing + (1 - chord_cos_spacing) * uniform

    if chord_cos_spacing == 0.:
        full_wing_x = np.linspace(0, 1., num_x)

    else:
        # Mirror this half wing into a full wing; offset by 0.5 so it goes 0 to 1
        full_wing_x = np.hstack((-half_wing[:-1], half_wing[::-1])) + .5

    # Obtain the leading and trailing edges
    le = mesh[ 0, :, :]
    te = mesh[-1, :, :]

    # Create a new mesh with the desired num_x and set the leading and trailing edge values
    new_mesh = np.zeros((num_x, num_y, 3))
    new_mesh[ 0, :, :] = le
    new_mesh[-1, :, :] = te

    for i in range(1, num_x-1):
        w = full_wing_x[i]
        new_mesh[i, :, :] = (1 - w) * le + w * te

    return new_mesh

def get_default_geo_dict():
    """
    Obtain the default settings for the surface descriptions. Note that
    these defaults are overwritten based on user input for each surface.
    Each dictionary describes one surface.

    Returns
    -------
    defaults : dict
        A python dict containing the default surface-level settings.
    """

    defaults = {
                # Wing definition
                'num_x' : 3,            # number of chordwise points
                'num_y' : 5,            # number of spanwise points
                'span_cos_spacing' : 0, # 0 for uniform spanwise panels
                                        # 1 for cosine-spaced panels
                                        # any value between 0 and 1 for
                                        # a mixed spacing
                'chord_cos_spacing' : 0.,   # 0 for uniform chordwise panels
                                        # 1 for cosine-spaced panels
                                        # any value between 0 and 1 for
                                        # a mixed spacing
                'wing_type' : 'rect',   # initial shape of the wing
                                        # either 'CRM' or 'rect'
                                        # 'CRM' can have different options
                                        # after it, such as 'CRM:alpha_2.75'
                                        # for the CRM shape at alpha=2.75
                'symmetry' : True,     # if true, model one half of wing
                                        # reflected across the plane y = 0
                'offset' : np.zeros((3)), # coordinates to offset
                                # the surface from its default location

                # Simple Geometric Variables
                'span' : 10.,           # full wingspan, even for symmetric cases
                'root_chord' : 1.,      # root chord
                'dihedral' : 0.,        # wing dihedral angle in degrees
                                        # positive is upward
                'sweep' : 0.,           # wing sweep angle in degrees
                                        # positive sweeps back
                'taper' : 1.,           # taper ratio; 1. is uniform chord
                }

    return defaults

def generate_mesh(input_dict):

    # Get defaults and update surface with the user-provided input
    surf_dict = get_default_geo_dict()
    surf_dict.update(input_dict)

    num_x = surf_dict['num_x']
    num_y = surf_dict['num_y']
    span = surf_dict['span']
    chord = surf_dict['root_chord']
    span_cos_spacing = surf_dict['span_cos_spacing']
    chord_cos_spacing = surf_dict['chord_cos_spacing']

    # Check to make sure that an odd number of spanwise points (num_y) was provided
    if not num_y % 2:
        raise ValueError('num_y must be an odd number.')

    # Check to make sure that an odd number of chordwise points (num_x) was provided
    if not num_x % 2 and not num_x==2:
        raise ValueError('num_x must be an odd number.')

    # Generate rectangular mesh
    if surf_dict['wing_type'] == 'rect':
        mesh = gen_rect_mesh(num_x, num_y, span, chord,
            span_cos_spacing, chord_cos_spacing)

    # Generate CRM mesh. Note that this outputs twist information
    # based on the data from the CRM definition paper, so we save
    # this twist information to the surf_dict.
    elif 'CRM' in surf_dict['wing_type']:
        mesh, eta, twist = gen_crm_mesh(num_x, num_y,
            span_cos_spacing, chord_cos_spacing, surf_dict['wing_type'])
        surf_dict['crm_twist'] = twist

    else:
        raise NameError('wing_type option not understood. Must be either a type of ' +
              '"CRM" or "rect".')

    # Chop the mesh in half if using symmetry during analysis.
    # Note that this means that the provided mesh should be the full mesh
    if surf_dict['symmetry']:
        num_y = int((num_y+1)/2)
        mesh = mesh[:, :num_y, :]

    # Apply the user-provided coordinate offset to position the mesh
    mesh = mesh + surf_dict['offset']

    # If CRM wing, then compute the jig twist values.
    # Interpolate the twist values from the CRM wing definition to the twist
    # control points.
    if 'CRM' in surf_dict['wing_type']:
        num_twist = surf_dict['num_twist_cp']

        # If the surface is symmetric, simply interpolate the initial
        # twist_cp values based on the mesh data
        if surf_dict['symmetry']:
            twist = np.interp(np.linspace(0, 1, num_twist), eta, surf_dict['crm_twist'])
        else:

            # If num_twist is odd, create the twist vector and mirror it
            # then stack the two together, but remove the duplicated twist
            # value.
            if num_twist % 2:
                twist = np.interp(np.linspace(0, 1, (num_twist+1)//2), eta, surf_dict['crm_twist'])
                twist = np.hstack((twist[:-1], twist[::-1]))

            # If num_twist is even, mirror the twist vector and stack
            # them together
            else:
                twist = np.interp(np.linspace(0, 1, num_twist//2), eta, surf_dict['crm_twist'])
                twist = np.hstack((twist, twist[::-1]))

        return mesh, twist

    else:

        return mesh

def write_FFD_file(surface, mx, my):

    mesh = surface['mesh']
    nx, ny = mesh.shape[:2]

    half_ffd = np.zeros((mx, my, 3))

    LE = mesh[0, :, :]
    TE = mesh[-1, :, :]

    half_ffd[0, :, 0] = np.interp(np.linspace(0, 1, my), np.linspace(0, 1, ny), LE[:, 0])
    half_ffd[0, :, 1] = np.interp(np.linspace(0, 1, my), np.linspace(0, 1, ny), LE[:, 1])
    half_ffd[0, :, 2] = np.interp(np.linspace(0, 1, my), np.linspace(0, 1, ny), LE[:, 2])

    half_ffd[-1, :, 0] = np.interp(np.linspace(0, 1, my), np.linspace(0, 1, ny), TE[:, 0])
    half_ffd[-1, :, 1] = np.interp(np.linspace(0, 1, my), np.linspace(0, 1, ny), TE[:, 1])
    half_ffd[-1, :, 2] = np.interp(np.linspace(0, 1, my), np.linspace(0, 1, ny), TE[:, 2])

    for i in range(my):
        half_ffd[:, i, 0] = np.linspace(half_ffd[0, i, 0], half_ffd[-1, i, 0], mx)
        half_ffd[:, i, 1] = np.linspace(half_ffd[0, i, 1], half_ffd[-1, i, 1], mx)
        half_ffd[:, i, 2] = np.linspace(half_ffd[0, i, 2], half_ffd[-1, i, 2], mx)

    cushion = .5

    half_ffd[0, :, 0] -= cushion
    half_ffd[-1, :, 0] += cushion
    half_ffd[:, 0, 1] -= cushion
    half_ffd[:, -1, 1] += cushion

    bottom_ffd = half_ffd.copy()
    bottom_ffd[:, :, 2] -= cushion

    top_ffd = half_ffd.copy()
    top_ffd[:, :, 2] += cushion

    ffd = np.vstack((bottom_ffd, top_ffd))

    # ### Uncomment this to plot the FFD points
    # import matplotlib.pyplot as plt
    # from mpl_toolkits.mplot3d import Axes3D
    #
    # fig = plt.figure()
    # axes = []
    #
    # axes.append(fig.add_subplot(221, projection='3d'))
    # axes.append(fig.add_subplot(222, projection='3d'))
    # axes.append(fig.add_subplot(223, projection='3d'))
    # axes.append(fig.add_subplot(224, projection='3d'))
    #
    # for i, ax in enumerate(axes):
    #     xs = ffd[:, :, 0].flatten()
    #     ys = ffd[:, :, 1].flatten()
    #     zs = ffd[:, :, 2].flatten()
    #
    #     ax.scatter(xs, ys, zs, c='red', alpha=1., clip_on=False)
    #
    #     xs = ffd[:, :, 0].flatten()
    #     ys = ffd[:, :, 1].flatten()
    #     zs = ffd[:, :, 2].flatten()
    #
    #     ax.scatter(xs, ys, zs, c='blue', alpha=1.)
    #
    #     xs = mesh[:, :, 0]
    #     ys = mesh[:, :, 1]
    #     zs = mesh[:, :, 2]
    #
    #     ax.plot_wireframe(xs, ys, zs, color='k')
    #
    #     ax.set_xlim([-5, 5])
    #     ax.set_ylim([-5, 5])
    #     ax.set_zlim([-5, 5])
    #
    #     ax.set_xlim([20, 40])
    #     ax.set_ylim([-25, -5.])
    #     ax.set_zlim([-10, 10])
    #
    #     ax.set_xlabel('x')
    #     ax.set_ylabel('y')
    #     ax.set_zlabel('z')
    #
    #     ax.set_axis_off()
    #
    #     ax.set_axis_off()
    #
    #     if i == 0:
    #         ax.view_init(elev=0, azim=180)
    #     elif i == 1:
    #         ax.view_init(elev=0, azim=90)
    #     elif i == 2:
    #         ax.view_init(elev=100000, azim=0)
    #     else:
    #         ax.view_init(elev=40, azim=-30)
    #
    # plt.tight_layout()
    # plt.subplots_adjust(wspace=0, hspace=0)
    #
    # plt.show()

    filename = surface['name'] + '_ffd.fmt'

    with open(filename, 'w') as f:
        f.write('1\n')
        f.write('{} {} {}\n'.format(mx, 2, my))
        x = np.array_str(ffd[:, :, 0].flatten(order='F'))[1:-1] + '\n'
        y = np.array_str(ffd[:, :, 1].flatten(order='F'))[1:-1] + '\n'
        z = np.array_str(ffd[:, :, 2].flatten(order='F'))[1:-1] + '\n'

        f.write(x)
        f.write(y)
        f.write(z)

    return filename

def writeMesh(mesh,filename):
    """
    Writes the OAS mesh in Tecplot .dat file format, for visualization and debugging purposes.

    Parameters
    ----------
    mesh[nx,ny,3] : numpy array
        The OAS mesh to be written.
    filename : str
        The file name including the .dat extension.
    """
    num_y = mesh.shape[0]
    num_x = mesh.shape[1]
    f = open(filename, 'w')
    f.write('\t\t1\n')
    f.write('\t\t%d\t\t%d\t\t%d\n' % (num_y, num_x, 1))

    x = mesh[:, :, 0]
    y = mesh[:, :, 1]
    z = mesh[:, :, 2]

    for dim in [x, y, z]:
        for iy in range(num_x):
            row = dim[:, iy]
            for val in row:
                f.write('\t{: 3.6f}'.format(val))
            f.write('\n')
    f.close()


def getFullMesh(left_mesh=None, right_mesh=None):
    """
    For a symmetric wing, OAS only keeps and does computation on the left half.
    This script mirros the OAS mesh and attaches it to the existing mesh to
    obtain the full mesh.

    Parameters
    ----------
    left_mesh[nx,ny,3] or right_mesh : numpy array
        The half mesh to be mirrored.

    Returns
    -------
    full_mesh[nx,2*ny-1,3] : numpy array
        The computed full mesh.
    """
    if left_mesh is None and right_mesh is None:
        raise ValueError("Either the left or right mesh need to be supplied.")
    elif left_mesh is not None and right_mesh is not None:
        raise ValueError("Please only provide either left or right mesh, not both.")
    elif left_mesh is not None:
        right_mesh = np.flip(left_mesh,axis=1).copy()
        right_mesh[:,:,1] *= -1
    else:
        left_mesh = np.flip(right_mesh,axis=1).copy()
        left_mesh[:,:,1] *= -1
    full_mesh = np.concatenate((left_mesh,right_mesh[:,1:,:]),axis=1)
    return full_mesh


def plot3D_meshes(file_name, zero_tol=0):
    """
    Reads in multi-surface meshes from a Plot3D mesh file for VLM analysis.

    Parameters
    ----------
    fileName : str
        Plot3D file name to be read in.
    zero_tol : float
        If a node location read in the file is below this magnitude we will just
        make it zero. This is useful for getting rid of noise in the surface
        that may be due to the meshing tools geometry tolerance.

    Returns
    -------
    mesh_dict : dict
        Dictionary holding the mesh of every surface included in the plot3D
        sorted by surface name.
    """
    file_handle = open(file_name, 'r')
    num_panels = int(file_handle.readline())
    # Get the multi-block dimensions of every included surface
    block_dims = file_handle.readline().split()

    # Now loop through remainder of file and pluck out mesh node locations
    mesh_list = []
    mesh_dict = {}
    for i in range(num_panels):
        [nx, ny, nz] = block_dims[3*i:3*i+3]
        # Use nx and ny to intialize mesh. Since these are surfaces nz always
        # equals 1, so no need to use it
        mesh = np.zeros(int(nx)*int(ny)*3)

        for j in range(mesh.size):
            line = file_handle.readline()
            val = float(line)
            if np.abs(val) < zero_tol:
                val = 0
            mesh[j] = val

        # Restructure mesh as 3D array,
        # Plot3D files are always written using Fortran order
        mesh_list.append(mesh.reshape([int(nx), int(ny), 3], order='f'))

    # Now read in names for each surface mesh
    for i in range(num_panels):
        name = file_handle.readline()[:-1]
        mesh_dict[name] = mesh_list[i]

    return mesh_dict
