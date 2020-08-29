import numpy as np


def norm(vec, axis=None):
    return np.sqrt(np.sum(vec**2, axis=axis))

def unit(vec):
    return vec / norm(vec)

def norm_d(vec):
    vec_d = vec/norm(vec)
    return vec_d

def unit_d(vec):
    n_d = norm_d(vec)
    normvec = norm(vec)
    vec_d = np.outer((-vec/(normvec*normvec)),n_d) + 1/normvec*np.eye(len(vec))

    return vec_d

# This is a limited cross product definition for 3 vectors
def cross_d(a,b):
    if not isinstance(a, np.ndarray) :
        a = np.array(a)
        if a.shape!=(3,):
            raise ValueError('a must be a (3,) nd array')
    if not isinstance(b, np.ndarray):
        b = np.array(b)
        if b.shape!=(3,):
            raise ValueError('b must be a (3,) nd array')

    dcda = np.zeros([3,3])
    dcdb = np.zeros([3,3])

    dcda[0,1]=b[2]
    dcda[0,2]=-b[1]
    dcda[1,0]=-b[2]
    dcda[1,2]=b[0]
    dcda[2,0]=b[1]
    dcda[2,1]=-b[0]

    dcdb[0,1]=-a[2]
    dcdb[0,2]=a[1]
    dcdb[1,0]=a[2]
    dcdb[1,2]=-a[0]
    dcdb[2,0]=-a[1]
    dcdb[2,1]=a[0]

    return dcda,dcdb


def radii(mesh, t_c=0.15):

    """
    Obtain the radii of the FEM element based on local chord.
    """
    vectors = mesh[-1, :, :] - mesh[0, :, :]
    chords = np.sqrt(np.sum(vectors**2, axis=1))
    mean_chords = 0.5 * chords[:-1] + 0.5 * chords[1:]
    return t_c * mean_chords * 0.5
