import numpy as np


def get_array_indices(*shape):
    return np.arange(np.prod(shape)).reshape(shape)


def add_ones_axis(array):
    return np.einsum("...,l->...l", array, np.ones(3))


def compute_dot(array1, array2):
    """
    Parameters
    ----------
    array1 : numpy array[..., 3]
        First argument in the dot product.
        The dot product axis is the last one.
    array2 : numpy array[..., 3]
        Second argument in the dot product.
        The dot product axis is the last one.
    """
    return np.einsum("...,i->...i", np.einsum("...i,...i->...", array1, array2), np.ones(3))


def compute_dot_deriv(array, deriv_array):
    """
    Parameters
    ----------
    array : numpy array[..., 3]
        The argument in the dot product we are not taking the derivatives for.
        The dot product axis is the last one.
    deriv_array : numpy array[..., 3, 3]
        The derivatives of the argument in the dot product of interest.
        The dot product axis is the last one.
    """
    return np.einsum(
        "...j,i->...ij",
        np.einsum("...ij,...i->...j", deriv_array, array),
        np.ones(3),
    )


def compute_cross(array1, array2):
    """
    Parameters
    ----------
    array1 : numpy array[..., 3]
        First argument in the cross product (order matters).
        The cross product axis is the last one.
    array2 : numpy array[..., 3]
        Second argument in the cross product (order matters).
        The cross product axis is the last one.
    """
    return np.cross(array1, array2, axis=-1)


def compute_cross_deriv1(deriv_array, array):
    """
    Parameters
    ----------
    deriv_array : numpy array[..., 3, 3]
        Derivatives of the first argument in the cross product.
        The cross product axis is the second last one.
    array : numpy array[..., 3]
        Second argument in the cross product (order matters).
        The cross product axis is the last one.
    """
    tmp_0 = np.einsum("...i,j->...ij", compute_cross(deriv_array[..., 0], array), np.array([1.0, 0.0, 0.0]))
    tmp_1 = np.einsum("...i,j->...ij", compute_cross(deriv_array[..., 1], array), np.array([0.0, 1.0, 0.0]))
    tmp_2 = np.einsum("...i,j->...ij", compute_cross(deriv_array[..., 2], array), np.array([0.0, 0.0, 1.0]))
    return tmp_0 + tmp_1 + tmp_2


def compute_cross_deriv2(array, deriv_array):
    """
    Parameters
    ----------
    array : numpy array[..., 3]
        First argument in the cross product (order matters).
        The cross product axis is the last one.
    deriv_array : numpy array[..., 3, 3]
        Derivatives of the second argument in the cross product.
        The cross product axis is the second last one.
    """
    tmp_0 = np.einsum("...i,j->...ij", compute_cross(array, deriv_array[..., 0]), np.array([1.0, 0.0, 0.0]))
    tmp_1 = np.einsum("...i,j->...ij", compute_cross(array, deriv_array[..., 1]), np.array([0.0, 1.0, 0.0]))
    tmp_2 = np.einsum("...i,j->...ij", compute_cross(array, deriv_array[..., 2]), np.array([0.0, 0.0, 1.0]))
    return tmp_0 + tmp_1 + tmp_2


def compute_norm(array):
    """
    Parameters
    ----------
    array : numpy array[..., 3]
        Array we are taking the norm of in the last axis.
    """
    return np.einsum("...,k->...k", np.sum(array**2, axis=-1) ** 0.5, np.ones(3))


def compute_norm_deriv(array, deriv_array):
    """
    Parameters
    ----------
    array : numpy array[..., 3]
        Array we are taking the norm of in the last axis.
    deriv_array : numpy array[..., 3, 3]
        Derivatives of the argument.
    """
    return np.einsum(
        "...j,...i->...ij",
        np.einsum("...i,...ij->...j", array, deriv_array),
        1.0 / compute_norm(array),
    )
