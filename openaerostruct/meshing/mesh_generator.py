import numpy as np
import warnings

from openaerostruct.meshing.CRM_definitions import get_crm_points
from openaerostruct.meshing.utils import getFullMesh, regen_chordwise_panels


def gen_rect_mesh(num_x, num_y, span, chord, span_cos_spacing=0.0, chord_cos_spacing=0.0):
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
        spanwise node points near the wingtips. A value between 2 and 3 will
        create cosine spacing at both the root and tips.
    chord_cos_spacing : float (optional)
        Blending ratio of uniform and cosine spacing in the chordwise direction.
        A value of 0. corresponds to uniform spacing and a value of 1.
        corresponds to regular cosine spacing. This increases the number of
        chordwise node points near the leading/trailing edge.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Rectangular nodal mesh defining the final aerodynamic surface with the
        specified parameters.
    """

    # Preallocate the mesh array
    mesh = np.zeros((num_x, num_y, 3))
    # Compute the symmetry index of the mesh
    ny2 = (num_y + 1) // 2

    # --- spanwise discretization ---
    # Spacings >= 2.0 bunch panels at both the root and tips
    if span_cos_spacing >= 2.0:
        beta = np.linspace(0, np.pi, ny2)

        # mixed spacing with span_cos_spacing as a weighting factor
        # this is for the spanwise spacing
        cosine = 0.5 * (1 - np.cos(beta))  # cosine spacing
        uniform = np.linspace(0, 0.5, ny2)[::-1]  # uniform spacing
        half_wing = cosine[::-1] * (span_cos_spacing - 2.0) + (1 - (span_cos_spacing - 2.0)) * uniform
        full_wing = np.hstack((-half_wing[:-1], half_wing[::-1])) * span

    else:
        beta = np.linspace(0, np.pi / 2, ny2)

        # mixed spacing with span_cos_spacing as a weighting factor
        # this is for the spanwise spacing
        cosine = 0.5 * np.cos(beta)  # cosine spacing
        uniform = np.linspace(0, 0.5, ny2)[::-1]  # uniform spacing
        half_wing = cosine * span_cos_spacing + (1 - span_cos_spacing) * uniform
        full_wing = np.hstack((-half_wing[:-1], half_wing[::-1])) * span

    # --- chordwise discretization ---
    cosine = 0.5 * (1 - np.cos(np.linspace(0, np.pi, num_x)))  # cosine spacing from 0 to 1
    uniform = np.linspace(0, 1, num_x)  # uniform spacing
    # mixed spacing with chord_cos_spacing as a weighting factor
    wing_x = cosine * chord_cos_spacing + (1 - chord_cos_spacing) * uniform
    wing_x *= chord  # apply chord length

    # --- form 3D mesh array ---
    for ind_x in range(num_x):
        for ind_y in range(num_y):
            mesh[ind_x, ind_y, :] = [wing_x[ind_x], full_wing[ind_y], 0]

    return mesh


def gen_crm_mesh(num_x, num_y, span_cos_spacing=0.0, chord_cos_spacing=0.0, wing_type="CRM:jig"):
    """
    Generate Common Research Model wing mesh.

    Parameters
    ----------
    num_x : float
        Desired number of chordwise node points for the final mesh.
    num_y : float
        Desired number of spanwise node points for the final mesh.
    span : float
        Total wingspan.
    chord : float
        Root chord.
    span_cos_spacing : float (optional)
        Blending ratio of uniform and cosine spacing in the spanwise direction.
        A value of 0. corresponds to uniform spacing and a value of 1.
        corresponds to regular cosine spacing. This increases the number of
        spanwise node points near the wingtips. A value between 2 and 3 will
        create cosine spacing at both the root and tips.
    chord_cos_spacing : float (optional)
        Blending ratio of uniform and cosine spacing in the chordwise direction.
        A value of 0. corresponds to uniform spacing and a value of 1.
        corresponds to regular cosine spacing. This increases the number of
        chordwise node points near the leading/trailing edge.
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
    if "jig" in wing_type or "CRM" == wing_type:
        raw_crm_points[:, 3] = 0.0

    # Get the leading edge of the raw crm points
    le = np.vstack((raw_crm_points[:, 1], raw_crm_points[:, 2], raw_crm_points[:, 3]))

    # Get the chord, twist(in correct order), and eta values from the points
    chord = raw_crm_points[:, 5]
    twist = raw_crm_points[:, 4][::-1]
    eta = raw_crm_points[:, 0]

    # Get the trailing edge of the crm points, based on the chord + le distance.
    # Note that we do not account for twist here; instead we set that using
    # the twist design variable later in run_classes.py.
    te = np.vstack((raw_crm_points[:, 1] + chord, raw_crm_points[:, 2], raw_crm_points[:, 3]))

    # Get the number of points that define this CRM shape and create a mesh
    # array based on this size
    n_raw_points = raw_crm_points.shape[0]
    mesh = np.empty((2, n_raw_points, 3))

    # Set the leading and trailing edges of the mesh matrix
    mesh[0, :, :] = le.T
    mesh[1, :, :] = te.T

    # Convert the mesh points to meters from inches.
    raw_mesh = mesh * 0.0254

    # Index of symmetry line
    ny2 = (num_y + 1) // 2

    # --- spanwise discretization ---
    # Create the blended spacing using the user input for span_cos_spacing
    # Spacings >= 2.0 bunch panels at both the root and tips
    if span_cos_spacing >= 2.0:
        beta = np.linspace(0, np.pi, ny2)

        # mixed spacing with span_cos_spacing as a weighting factor
        # this is for the spanwise spacing
        cosine = 1 - np.cos(beta)  # cosine spacing
        uniform = np.linspace(0, 1.0, ny2)[::-1]  # uniform spacing
        lins = cosine[::-1] * (span_cos_spacing - 2.0) + (1 - (span_cos_spacing - 2.0)) * uniform
    else:
        beta = np.linspace(0, np.pi / 2, ny2)

        # mixed spacing with span_cos_spacing as a weighting factor
        # this is for the spanwise spacing
        cosine = np.cos(beta)  # cosine spacing
        uniform = np.linspace(0, 1.0, ny2)[::-1]  # uniform spacing
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
        full_mesh = regen_chordwise_panels(full_mesh, num_x, chord_cos_spacing)

    return full_mesh, eta, twist


def get_default_geo_dict():
    """
    Obtain the default settings for the mesh descriptions. Note that
    these defaults are overwritten based on user input for each mesh.
    Each dictionary describes one mesh.

    Returns
    -------
    defaults : dict
        A python dict containing the default surface-level settings.
    """

    defaults = {
        # Wing definition
        "num_x": 3,  # number of chordwise points
        "num_y": 5,  # number of spanwise points
        "span_cos_spacing": 0,  # 0 for uniform spanwise panels
        # 1 for cosine-spaced panels
        # any value between 0 and 1 for
        # a mixed spacing
        "chord_cos_spacing": 0.0,  # 0 for uniform chordwise panels
        # 1 for cosine-spaced panels
        # any value between 0 and 1 for
        # a mixed spacing
        "wing_type": "rect",  # initial shape of the wing
        # either 'CRM' or 'rect'
        # 'CRM' can have different options
        # after it, such as 'CRM:alpha_2.75'
        # for the CRM shape at alpha=2.75
        "symmetry": True,  # if true, model one half of wing
        # reflected across the plane y = 0
        "offset": np.zeros((3)),  # coordinates to offset
        # the surface from its default location
        # Simple Geometric Variables
        "span": 10.0,  # full wingspan, even for symmetric cases
        "root_chord": 1.0,  # root chord
        "num_twist_cp": 2,  # number of twist controling point, only relevant for CRM wings.
    }

    return defaults


def generate_mesh(input_dict):
    """
    Generate a rectangular or CRM OAS mesh an input mesh dictionary.

    Parameters
    ----------
    input_dict : dict
        Dictionary containing user-provided parameters for the surface definition.
        See the following for more information:
        https://mdolab-openaerostruct.readthedocs-hosted.com/en/latest/user_reference/mesh_surface_dict.html#mesh-dict

    Returns
    -------
    mesh : numpy array
        Nodal coordinates defining the mesh.
        shape = (nx, ny, 3),
        where nx is the number of chordwise discretization nodes, ny is the number of spanwise discretization nodes.
        If input_dict["symmetry"] is True, mesh defines left half of wing.
    twist : numpy array, optional
        Only for CRM wing (input_dict["wing_type"] == "CRM").
        Twist values at the spanwise locations.

    """

    # Get defaults and update surface with the user-provided input
    surf_dict = get_default_geo_dict()

    # Warn if a user provided a key that is not implemented
    user_defined_keys = input_dict.keys()
    for key in user_defined_keys:
        if key not in surf_dict:
            warnings.warn(
                "Key `{}` in mesh_dict is not implemented and will be ignored".format(key),
                category=RuntimeWarning,
                stacklevel=2,
            )
    # Warn if a user did not define important keys
    for key in ["num_x", "num_y", "wing_type", "symmetry"]:
        if key not in user_defined_keys:
            warnings.warn(
                "Missing `{}` in mesh_dict. The default value of {} will be used.".format(key, surf_dict[key]),
                category=RuntimeWarning,
                stacklevel=2,
            )

    # Apply user-defined options
    surf_dict.update(input_dict)

    # Warn if a user defined span and root_chord for CRM
    if surf_dict["wing_type"] == "CRM":
        if "span" in user_defined_keys or "root_chord" in user_defined_keys:
            warnings.warn(
                "`span` and `root_chord` in mesh_dict will be ignored for the CRM wing.",
                category=RuntimeWarning,
                stacklevel=2,
            )

    num_x = surf_dict["num_x"]
    num_y = surf_dict["num_y"]
    span_cos_spacing = surf_dict["span_cos_spacing"]
    chord_cos_spacing = surf_dict["chord_cos_spacing"]

    # Check to make sure that an odd number of spanwise points (num_y) was provided
    if not num_y % 2:
        raise ValueError("num_y must be an odd number.")

    # Generate rectangular mesh
    if surf_dict["wing_type"] == "rect":
        span = surf_dict["span"]
        chord = surf_dict["root_chord"]
        mesh = gen_rect_mesh(num_x, num_y, span, chord, span_cos_spacing, chord_cos_spacing)

    # Generate CRM mesh. Note that this outputs twist information
    # based on the data from the CRM definition paper, so we save
    # this twist information to the surf_dict.
    elif "CRM" in surf_dict["wing_type"]:
        mesh, eta, twist = gen_crm_mesh(num_x, num_y, span_cos_spacing, chord_cos_spacing, surf_dict["wing_type"])
        surf_dict["crm_twist"] = twist

    else:
        raise NameError("wing_type option not understood. Must be either a type of " + '"CRM" or "rect".')

    # Chop the mesh in half if using symmetry during analysis.
    # Note that this means that the provided mesh should be the full mesh
    if surf_dict["symmetry"]:
        num_y = int((num_y + 1) / 2)
        mesh = mesh[:, :num_y, :]

    # Apply the user-provided coordinate offset to position the mesh
    mesh = mesh + surf_dict["offset"]

    # If CRM wing, then compute the jig twist values.
    # Interpolate the twist values from the CRM wing definition to the twist
    # control points.
    if "CRM" in surf_dict["wing_type"]:
        num_twist = surf_dict["num_twist_cp"]

        # If the surface is symmetric, simply interpolate the initial
        # twist_cp values based on the mesh data
        if surf_dict["symmetry"]:
            twist = np.interp(np.linspace(0, 1, num_twist), eta, surf_dict["crm_twist"])
        else:
            # If num_twist is odd, create the twist vector and mirror it
            # then stack the two together, but remove the duplicated twist
            # value.
            if num_twist % 2:
                twist = np.interp(np.linspace(0, 1, (num_twist + 1) // 2), eta, surf_dict["crm_twist"])
                twist = np.hstack((twist[:-1], twist[::-1]))

            # If num_twist is even, mirror the twist vector and stack
            # them together
            else:
                twist = np.interp(np.linspace(0, 1, num_twist // 2), eta, surf_dict["crm_twist"])
                twist = np.hstack((twist, twist[::-1]))

        return mesh, twist

    else:
        return mesh
