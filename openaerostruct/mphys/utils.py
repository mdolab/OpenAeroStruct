import numpy as np


def get_number_of_nodes(surfaces):
    """
    Get the total number of nodes over all surfaces.

    Parameters
    ----------
    surfaces : list[dict]
        List of all surfaces.

    Returns
    -------
    nnodes : int
        Total number of nodes across all surfaces.

    """
    nnodes = 0
    for surface in surfaces:
        nnodes += surface["mesh"].size // 3
    return nnodes


def get_src_indices(surfaces):
    """
    Get src indices for each surface that will project each mesh into a single flattened array.

    Parameters
    ----------
    surfaces : list[dict]
        List of all surfaces.

    Returns
    -------
    src_indices : dict
        Dictionary holding source indices sorted by surface name.

    """
    src_indices = {}
    nindices = 0
    for surface in surfaces:
        surf_name = surface["name"]
        mesh = surface["mesh"]
        nx, ny, _ = mesh.shape
        surf_indices = np.arange(mesh.size) + nindices
        src_indices[surf_name] = surf_indices.reshape(nx, ny, 3)
        nindices += mesh.size
    return src_indices


def get_node_indices(surfaces):
    """
    Get node indices for each surface that define the ordering of each node for all surfaces.

    Parameters
    ----------
    surfaces : list[dict]
        List of all surfaces.

    Returns
    -------
    node_indices : dict
        Dictionary holding node indices sorted by surface name.

    """
    node_indices = {}
    nnodes = 0
    for surface in surfaces:
        surf_name = surface["name"]
        mesh = surface["mesh"]
        nx, ny, _ = mesh.shape
        surf_indices = np.arange(nx * ny) + nnodes
        node_indices[surf_name] = surf_indices.reshape(nx, ny)
        nnodes += nx * ny
    return node_indices
