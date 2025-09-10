import numpy as np


def norm(vec, axis=None):
    return np.sqrt(np.sum(vec**2, axis=axis))


def unit(vec):
    return vec / norm(vec)


def norm_d(vec):
    vec_d = vec / norm(vec)
    return vec_d


def unit_d(vec):
    n_d = norm_d(vec)
    normvec = norm(vec)
    vec_d = np.outer((-vec / (normvec * normvec)), n_d) + 1 / normvec * np.eye(len(vec))

    return vec_d


# This is a limited cross product definition for 3 vectors
def cross_d(a, b):
    if not isinstance(a, np.ndarray):
        a = np.array(a)
        if a.shape != (3,):
            raise ValueError("a must be a (3,) nd array")
    if not isinstance(b, np.ndarray):
        b = np.array(b)
        if b.shape != (3,):
            raise ValueError("b must be a (3,) nd array")

    dcda = np.zeros([3, 3])
    dcdb = np.zeros([3, 3])

    dcda[0, 1] = b[2]
    dcda[0, 2] = -b[1]
    dcda[1, 0] = -b[2]
    dcda[1, 2] = b[0]
    dcda[2, 0] = b[1]
    dcda[2, 1] = -b[0]

    dcdb[0, 1] = -a[2]
    dcdb[0, 2] = a[1]
    dcdb[1, 0] = a[2]
    dcdb[1, 2] = -a[0]
    dcdb[2, 0] = -a[1]
    dcdb[2, 1] = a[0]

    return dcda, dcdb


def radii(mesh, t_c=0.15):
    """
    Obtain the radii of the FEM element based on local chord.
    """
    vectors = mesh[-1, :, :] - mesh[0, :, :]
    chords = np.sqrt(np.sum(vectors**2, axis=1))
    mean_chords = 0.5 * chords[:-1] + 0.5 * chords[1:]
    return t_c * mean_chords * 0.5


def transformation_matrix(theta):
    """
    Function to find the transformation matrix for a given angle.

    Parameters
    ----------
    theta : float
        Angle in degrees.

    Returns
    -------
    T : numpy array
        Transformation matrix.
    """
    theta = theta * np.pi / 180
    c = np.cos(theta)
    s = np.sin(theta)
    T = np.zeros((3, 3))
    T[0, 0] = c**2
    T[0, 1] = s**2
    T[0, 2] = 2 * s * c
    T[1, 0] = s**2
    T[1, 1] = c**2
    T[1, 2] = -2 * s * c
    T[2, 0] = -s * c
    T[2, 1] = s * c
    T[2, 2] = c**2 - s**2
    return T


def compute_composite_stiffness(surface):
    """
    Function to compute the effective E and G stiffness values for a composite material,
    based on the ply_fractions, ply angles and individual fiber and matrix properties.

    Parameters
    ----------
    surface : dict
        OpenAeroStruct surface dictionary.
    """
    E1 = surface["E1"]
    E2 = surface["E2"]
    v12 = surface["nu12"]
    G12 = surface["G12"]
    v21 = (E2 / E1) * v12
    ply_fractions = surface["ply_fractions"]
    ply_angles = surface["ply_angles"]
    num_plies = len(ply_fractions)

    # check inputs
    if len(ply_fractions) != len(ply_angles):
        raise ValueError("Length of ply_fractions and ply_angles must be equal")
    if sum(ply_fractions) != 1:
        raise ValueError("Sum of ply_fractions must be 1")

    # finding the Q matrix
    Q = np.zeros((3, 3))
    Q[0, 0] = E1 / (1 - v12 * v21)
    Q[0, 1] = v12 * E2 / (1 - v12 * v21)
    Q[0, 2] = 0
    Q[1, 0] = v21 * E1 / (1 - v12 * v21)
    Q[1, 1] = E2 / (1 - v12 * v21)
    Q[1, 2] = 0
    Q[2, 0] = 0
    Q[2, 1] = 0
    Q[2, 2] = G12

    # finding the Q_bar matrix for each ply in the form of a 3D Array
    Q_bar = np.zeros((num_plies, 3, 3))
    Q_bar_eff = np.zeros((3, 3))
    for i in range(num_plies):
        theta = ply_angles[i]
        T = transformation_matrix(theta)
        Q_bar[i] = np.dot(np.dot(np.linalg.inv(T), Q), T)
        Q_bar_eff += ply_fractions[i] * Q_bar[i]

    S_bar_eff = np.linalg.inv(Q_bar_eff)
    E_eff = 1 / S_bar_eff[0, 0]
    G_eff = 1 / S_bar_eff[2, 2]

    # replacing the values in the surface dictionary
    surface["E"] = E_eff
    surface["G"] = G_eff

    # no need to return anything as the values are updated in the surface dictionary (call by reference)