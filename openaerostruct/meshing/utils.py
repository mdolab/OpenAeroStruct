import numpy as np


def regen_chordwise_panels(mesh, num_x, chord_cos_spacing):
    """
    Generates a new mesh based on an existing mesh with the specified number of
    chordwise points.

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

    # chordwise discretization
    cosine = 0.5 * (1 - np.cos(np.linspace(0, np.pi, num_x)))  # cosine spacing from 0 to 1
    uniform = np.linspace(0, 1, num_x)  # uniform spacing
    # mixed spacing with chord_cos_spacing as a weighting factor
    wing_x = cosine * chord_cos_spacing + (1 - chord_cos_spacing) * uniform

    # Obtain the leading and trailing edges
    le = mesh[0, :, :]
    te = mesh[-1, :, :]

    # Create a new mesh with the desired num_x and set the leading and trailing edge values
    new_mesh = np.zeros((num_x, num_y, 3))
    new_mesh[0, :, :] = le
    new_mesh[-1, :, :] = te

    for i in range(1, num_x - 1):
        w = wing_x[i]
        new_mesh[i, :, :] = (1 - w) * le + w * te

    return new_mesh


def regen_spanwise_panels(mesh, num_y, span_cos_spacing):
    """
    Generates a new mesh based on an existing mesh with the specified number of
    spanwise points.

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the initial aerodynamic surface with only
        the leading and trailing edges defined.
    num_y : float
        Desired number of spanwise node points for the final mesh.
    span_cos_spacing : float
        Blending ratio of uniform and cosine spacing in the chordwise direction.
        A value of 0. corresponds to uniform spacing and a value of 1.
        corresponds to regular cosine spacing. This increases the number of
        chordwise node points near the wingtips. A value between 2 and 3 will
        create cosine spacing at both the root and tips.

    Returns
    -------
    new_mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the final aerodynamic surface with the
        specified number of spanwise node points.

    """
    # Determine if we have a half mesh or a full mesh
    if np.all(mesh[:, :, 1] >= -1e-8) or np.all(mesh[:, :, 1] <= 1e-8):
        symmetry = True
        print("Converting half mesh to full mesh")
        right_wing = abs(mesh[0, 0, 1]) < abs(mesh[0, -1, 1])
        if right_wing:
            mesh = getFullMesh(right_mesh=mesh)
        else:
            mesh = getFullMesh(left_mesh=mesh)
    else:
        symmetry = False

    # Obtain mesh and num properties
    num_x = mesh.shape[0]
    ny2 = (num_y + 1) // 2

    # Get le y coordds
    le = mesh[0, :, 1]
    span = np.abs(le[-1] - le[0])

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

    # Populate a mesh object with the desired num_y dimension based on
    # interpolated values from the raw CRM points.
    new_mesh = np.empty((num_x, num_y, 3))
    for j in range(num_x):
        for i in range(3):
            new_mesh[j, :, i] = np.interp(full_wing, le, mesh[j, :, i].real)

    if symmetry:
        return new_mesh[:, :ny2, :]

    return new_mesh


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
        right_mesh = np.flip(left_mesh, axis=1).copy()
        right_mesh[:, :, 1] *= -1
    else:
        left_mesh = np.flip(right_mesh, axis=1).copy()
        left_mesh[:, :, 1] *= -1
    full_mesh = np.concatenate((left_mesh, right_mesh[:, 1:, :]), axis=1)
    return full_mesh


def write_tecplot(mesh, filename):
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
    f = open(filename, "w")
    f.write("\t\t1\n")
    f.write("\t\t%d\t\t%d\t\t%d\n" % (num_y, num_x, 1))

    x = mesh[:, :, 0]
    y = mesh[:, :, 1]
    z = mesh[:, :, 2]

    for dim in [x, y, z]:
        for iy in range(num_x):
            row = dim[:, iy]
            for val in row:
                f.write("\t{: 3.6f}".format(val))
            f.write("\n")
    f.close()
