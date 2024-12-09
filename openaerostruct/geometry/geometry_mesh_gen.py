""" Utility for quickly generating a multi-section user specificed OAS mesh
"""

import numpy as np
import matplotlib.pyplot as plt


def generate_mesh(surface):
    """
    Automatically generates meshes for single and multi-section wings using the OpenAeroStruct surface dictionary.

    Parameters
    ----------
    surface : dict
        OpenAeroStruct surface or multi-section surface dictionary.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the aerodynamic surface.

    sec_mesh : list
        List of nodal meshes corresponding to each section in a multi-section surface.
    """

    # Get Data
    num_sections = surface["num_sections"]
    symmetry = surface["symmetry"]

    if symmetry or num_sections == 1:
        root_section = num_sections - 1
    else:
        if "root_section" not in surface.keys():
            raise Exception("The root section of an asymmetrical mesh needs to be identified")
        else:
            root_section = surface["root_section"]

    # Geometry data dictionary
    section_data = {
        "taper": surface["taper"],
        "sweep": surface["sweep"],
        "span": surface["span"],
        "root_chord": surface["root_chord"],
    }

    if "bpanels" in surface.keys():
        ny = surface["bpanels"] + 1
    else:
        ny = surface["ny"]

    if "cpanels" in surface.keys():
        nx = surface["cpanels"] + 1
    else:
        nx = surface["nx"]

    # Generate the mesh data
    panel_gx, panel_gy = generate_section_geometry(num_sections, symmetry, section_data, ny, nx, root_section)

    # Sitch the mesh into a unified mesh and out in OAS format
    panel_geom_x, panel_geom_y = stitch_section_geometry(num_sections, panel_gy, panel_gx)
    mesh = output_oas_mesh(panel_geom_x, panel_geom_y)

    # Produce meshes for each section in OAS format
    sec_meshes = []
    for section in range(num_sections):
        sec_mesh = output_oas_mesh(panel_gx[section], panel_gy[section])
        sec_meshes.append(sec_mesh)

    return mesh, sec_meshes


def generate_section_geometry(sections, symmetry, section_data, ny, nx, root_section):
    """
    Constructs the multi-section wing geometry specified by the user and generates a mesh for each section.

    Parameters
    ----------
    sections : int
        Integer for number of wing sections specified
    symmetry : bool
        Bool inidicating if the funciton should only generate the left span of the wing
    section_data : dict
        Dictionary with arrays corresponding the taper, span, sweep, and root chord of each section
    ny : numpy array
        Array with ints correponding to the number of spanwise points per section
    nx : int
        Number of chordwise points
    root_section:
        The section number that should be treated as the root section(y=0 origin)


    Returns
    -------
    panel_gx : List
         List containing the mesh x-coordinates for each section
    panel_gy : List
        List containing the mesh y-coordinates for each section
    """

    # Preallocate the lists
    panel_gy = [None] * sections
    panel_gx = [None] * sections

    # Jump to root section and build left wing
    for sec in np.arange(root_section, -1, -1):
        taper = section_data["taper"][sec]
        b = section_data["span"][sec]
        le_lambda = section_data["sweep"][sec]

        if sec == root_section:
            root_c = section_data["root_chord"]
            root_te = 0
            root_y = 0
        else:
            root_c = np.abs(panel_gx[sec + 1][0, 0] - panel_gx[sec + 1][nx - 1, 0])
            root_te = panel_gx[sec + 1][nx - 1, 0]
            root_y = panel_gy[sec + 1][0]
        tip_c = root_c * taper

        root_le = root_c + root_te
        tip_le = root_le - b * np.tan(le_lambda)
        tip_te = tip_le - tip_c

        root_x = np.linspace(root_le, root_te, nx)

        if tip_le == tip_te:
            tip_x = tip_le * np.ones(len(root_x))
        else:
            tip_x = np.linspace(tip_le, tip_te, len(root_x))

        panel_geom_y = np.zeros(ny[sec])
        panel_geom_x = np.zeros([nx, ny[sec]])

        panel_geom_y = np.linspace(root_y - b, root_y, ny[sec])

        for i in range(len(root_x)):
            panel_geom_x[i, :] = root_x[i] - ((tip_x[i] - root_x[i]) / b) * (panel_geom_y - root_y)

        panel_gy[sec] = panel_geom_y
        panel_gx[sec] = panel_geom_x

    # Build the right wing if asymmetrical
    if not symmetry:
        for sec in np.arange(root_section + 1, sections):
            taper = section_data["taper"][sec]
            b = section_data["span"][sec]
            le_lambda = section_data["sweep"][sec]

            root_c = np.abs(panel_gx[sec - 1][0, -1] - panel_gx[sec - 1][nx - 1, -1])
            tip_c = root_c * taper
            root_y = panel_gy[sec - 1][-1]

            root_te = panel_gx[sec - 1][nx - 1, 0]
            root_le = root_c + root_te
            tip_le = root_le + b * np.tan(le_lambda)
            tip_te = tip_le - tip_c

            root_x = np.linspace(root_le, root_te, nx)

            if tip_le == tip_te:
                tip_x = tip_le * np.ones(len(root_x))
            else:
                tip_x = np.linspace(tip_le, tip_te, len(root_x))

            panel_geom_y = np.zeros(ny[sec])
            panel_geom_x = np.zeros([nx, ny[sec]])

            panel_geom_y[0 : ny[sec]] = np.linspace(root_y, root_y + b, ny[sec])

            for i in range(len(root_x)):
                panel_geom_x[i, :] = root_x[i] + ((tip_x[i] - root_x[i]) / (b / 2)) * (panel_geom_y - root_y)

        panel_gy[sec] = panel_geom_y
        panel_gx[sec] = panel_geom_x

    return panel_gx, panel_gy


def stitch_section_geometry(sections, panel_gy, panel_gx):
    """
    Combines the split section array into singular unified mesh

    Parameters
    ----------
    sections : int
        Integer for number of wing sections specified
    panel_gx : List
         List containing the mesh x-coordinates for each section
    panel_gy : List
        List containing the mesh y-coordinates for each section


    Returns
    -------
    panel_geom_x : numpy array
         Array of the mesh x-coordinates
    panel_geom_y : numpy array
        Array of the mesh y-coordinates
    """
    # Stitch the results into a singular mesh

    if sections > 1:
        panel_geom_y = panel_gy[0][:-1]
        panel_geom_x = panel_gx[0][:, :-1]
        for i in np.arange(1, sections - 1):
            panel_geom_y = np.concatenate((panel_geom_y, panel_gy[i][:-1]))
            panel_geom_x = np.concatenate((panel_geom_x, panel_gx[i][:, :-1]), axis=1)
        panel_geom_y = np.concatenate((panel_geom_y, panel_gy[sections - 1]))
        panel_geom_x = np.concatenate((panel_geom_x, panel_gx[sections - 1]), axis=1)
    else:
        panel_geom_y = panel_gy[0]
        panel_geom_x = panel_gx[0]
    return panel_geom_x, panel_geom_y


def reflect_symmetric(panel_geom_x, panel_geom_y):
    """
    Reflects the mesh over y=0

    Parameters
    ----------
    panel_geom_x : numpy array
         Array of the mesh x-coordinates
    panel_geom_y : numpy array
        Array of the mesh y-coordinates

    Returns
    -------
    panel_geom_x : numpy array
         Array of the mesh x-coordinates
    panel_geom_y : numpy array
        Array of the mesh y-coordinates
    """
    panel_geom_x = np.hstack((panel_geom_x, panel_geom_x))
    panel_geom_y = np.hstack((panel_geom_y, -panel_geom_y))
    return panel_geom_x, panel_geom_y


def output_oas_mesh(panel_geom_x, panel_geom_y):
    """
    Outputs the mesh in OAS format

    Parameters
    ----------
    panel_geom_x : numpy array
        2D array of the mesh x-coordinates
    panel_geom_y : numpy array
        1D array of the mesh y-coordinates

    Returns
    -------
    mesh : numpy array
         3-D array with the OAS format mesh
    """
    panel_geom_y = np.broadcast_to(panel_geom_y, (panel_geom_x.shape[0], len(panel_geom_y)))
    mesh = np.zeros((panel_geom_x.shape[0], panel_geom_y.shape[1], 3))
    mesh[:, :, 0] = panel_geom_x[::-1]
    mesh[:, :, 1] = panel_geom_y
    return mesh


if __name__ == "__main__":
    """Runs the mesh generator to be run independently of OAS. Example 2 section mesh provided."""

    surface = {
        # Wing definition
        # Basic surface parameters
        "name": "surface",
        "num_sections": 2,  # The number of sections in the multi-section surface
        "sec_name": ["sec0", "sec1"],  # names of the individual sections
        "symmetry": True,  # if true, model one half of wing. reflected across the midspan of the root section
        "S_ref_type": "wetted",  # how we compute the wing area,
        # can be 'wetted' or 'projected'
        "root_section": 1,
        # Geometry Parameters
        "taper": np.array([1.0, 1.0]),  # Wing taper for each section
        "span": np.array([1.0, 1.0]),  # Wing span for each section
        "sweep": np.array([0.0, 0.0]),  # Wing sweep for each section
        "chord_cp": [np.array([1, 1]), np.array([1.0, 0.2])],
        # "sec_chord_cp": [np.ones(1),2*np.ones(1),3*np.ones(1)], #Chord B-spline control points for each section
        "root_chord": 1.0,  # Wing root chord for each section
        # Mesh Parameters
        "meshes": "gen-meshes",  # Supply a mesh for each section or "gen-meshes" for automatic mesh generation
        "nx": 2,  # Number of chordwise points. Same for all sections
        "ny": np.array([2, 2]),  # Number of spanwise points for each section
        # Aerodynamic Parameters
        "CL0": 0.0,  # CL of the surface at alpha=0
        "CD0": 0.015,  # CD of the surface at alpha=0
        # Airfoil properties for viscous drag calculation
        "k_lam": 0.05,  # percentage of chord with laminar
        # flow, used for viscous drag
        "sec_t_over_c_cp": [np.array([0.15]), np.array([0.15])],  # thickness over chord ratio (NACA0015)
        "sec_c_max_t": 0.303,  # chordwise location of maximum (NACA0015)
        # thickness
        "with_viscous": False,  # if true, compute viscous drag
        "with_wave": False,  # if true, compute wave drag
        "groundplane": False,
    }

    meshT, sec_meshes = generate_mesh(surface)

    def plot_meshes(meshes):
        """this function plots to plot the mesh"""
        plt.figure(figsize=(8, 4))
        for i, mesh in enumerate(meshes):
            mesh_x = mesh[:, :, 0]
            mesh_y = mesh[:, :, 1]
            color = "k"
            for i in range(mesh_x.shape[0]):
                plt.plot(mesh_y[i, :], mesh_x[i, :], color, lw=1)
                # plt.plot(-mesh_y[i, :], mesh_x[i, :], color, lw=1)   # plots the other side of symmetric wing
            for j in range(mesh_x.shape[1]):
                plt.plot(mesh_y[:, j], mesh_x[:, j], color, lw=1)
                # plt.plot(-mesh_y[:, j], mesh_x[:, j], color, lw=1)   # plots the other side of symmetric wing
        plt.axis("equal")
        plt.xlabel("y (m)")
        plt.ylabel("x (m)")

    plot_meshes(sec_meshes)
    plt.show()
