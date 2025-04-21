import numpy as np
from numpy import cos, sin, tan
import warnings
import copy

# openvsp python interface
try:
    import openvsp as vsp
    import degen_geom as dg
except ImportError:
    vsp = None
    dg = None

import openmdao.api as om
from openaerostruct.meshing.section_mesh_generator import generate_mesh as generate_section_mesh

# import functions for backward compatibility with old scripts
from openaerostruct.meshing.mesh_generator import (
    generate_mesh as _generate_mesh,
    gen_rect_mesh as _gen_rect_mesh,
    gen_crm_mesh as _gen_crm_mesh,
    get_default_geo_dict as _get_default_geo_dict,
)

from openaerostruct.meshing.utils import (
    regen_chordwise_panels as _regen_chordwise_panels,
    getFullMesh as _getFullMesh,
    write_tecplot as _write_tecplot,
)


def rotate(mesh, theta_y, symmetry, rotate_x=True, ref_axis_pos=0.25):
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
        applied perpendicular to the wing. To clarify, this is to ensure that non-planar surfaces
        such as winglets are twisted correctly about their axis. The winglets themselves have
        to be created with zshear or a user created mesh.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the twisted aerodynamic surface.

    """

    # Get trailing edge coordinates (ny, 3)
    te = mesh[-1]
    # Get leading edge coordinates (ny, 3)
    le = mesh[0]
    # Linear interpolation to compute the ref_axis coordinates (ny, 3)
    ref_axis = ref_axis_pos * te + (1 - ref_axis_pos) * le

    # Get number of spanwise stations (ny)
    nx, ny, _ = mesh.shape

    # Option to include mesh rotations about x-axis
    if rotate_x:
        # Compute x-axis rotation angle distribution using spanwise z displacements along quarter chord
        if symmetry:
            # This computes the change in dihedral angle along the references axis
            dz_qc = ref_axis[:-1, 2] - ref_axis[1:, 2]
            dy_qc = ref_axis[:-1, 1] - ref_axis[1:, 1]
            theta_x = np.arctan(dz_qc / dy_qc)

            # Prepend with 0 so that root is not rotated
            rad_theta_x = np.append(theta_x, 0.0)
        else:
            root_index = int((ny - 1) / 2)
            dz_qc_left = ref_axis[:root_index, 2] - ref_axis[1 : root_index + 1, 2]
            dy_qc_left = ref_axis[:root_index, 1] - ref_axis[1 : root_index + 1, 1]
            theta_x_left = np.arctan(dz_qc_left / dy_qc_left)
            dz_qc_right = ref_axis[root_index + 1 :, 2] - ref_axis[root_index:-1, 2]
            dy_qc_right = ref_axis[root_index + 1 :, 1] - ref_axis[root_index:-1, 1]
            theta_x_right = np.arctan(dz_qc_right / dy_qc_right)

            # Concatenate thetas with 0 at the root so it's not rotated
            rad_theta_x = np.concatenate((theta_x_left, np.zeros(1), theta_x_right))

    else:
        # If there is no rotation about x applied then the angle is 0
        rad_theta_x = 0.0

    rad_theta_y = theta_y * np.pi / 180.0

    # Initialize rotation matrix
    # Each spanwise (ny) station needs it's own 3x3 rotation matrix so this is 3D array of size (ny, 3, 3)
    mats = np.zeros((ny, 3, 3), dtype=type(rad_theta_y[0]))

    # Compute sin and cos of angles for the matrix
    cos_rtx = cos(rad_theta_x)
    cos_rty = cos(rad_theta_y)
    sin_rtx = sin(rad_theta_x)
    sin_rty = sin(rad_theta_y)

    # Each rotation matrix is 3x3 and is the product Rx(rad_theta_x)Ry(rad_theta_y)
    # Rx = [[0, 0, 0], [0, cos(rad_theta_x), -sin(rad_theta_x)], [0, sin(rad_theta_x), cos(rad_theta_x)]]
    # Ry = [[cos(rad_theta_y),0,-sin(rad_theta_y)], [0, 0, 0], [-sin(rad_theta_y), 0, cos(rad_theta_y)]]
    # RxRy = [[cos(rad_theta_y), 0, sin(rad_theta_y)],[sin(rad_theta_x)*sin(rad_theta_y), cos(rad_theta_x), -sin(rad_theta_x)*cos(rad_theta_y)], ...
    # [-cos(rad_theta_x)*sin(rad_theta_y), sin(rad_theta_x), cos(rad_theta_x)*cos(rad_theta_y)]]

    mats[:, 0, 0] = cos_rty
    mats[:, 0, 2] = sin_rty
    mats[:, 1, 0] = sin_rtx * sin_rty
    mats[:, 1, 1] = cos_rtx
    mats[:, 1, 2] = -sin_rtx * cos_rty
    mats[:, 2, 0] = -cos_rtx * sin_rty
    mats[:, 2, 1] = sin_rtx
    mats[:, 2, 2] = cos_rtx * cos_rty

    # Multiply each point on the mesh by the rotation matrix associated with its spanwise station
    # i - spanwise station index (ny)
    # m - chordwise station index
    # k - output vector(After rotation)
    # j - inputs vector(Before rotation)
    mesh[:] = np.einsum("ikj, mij -> mik", mats, mesh - ref_axis) + ref_axis


def scale_x(mesh, chord_dist, ref_axis_pos=0.25):
    """
    Modify the chords along the span of the wing by scaling only the x-coord.

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the initial aerodynamic surface.
    chord_dist[ny] : numpy array
        Spanwise distribution of the chord scaler.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh with the new chord lengths.
    """
    # Get trailing edge coordinates (ny, 3)
    te = mesh[-1]
    # Get leading edge coordinates (ny, 3)
    le = mesh[0]
    # Linear interpolation to compute the reference axis coordinates (ny, 3)
    ref_axis = ref_axis_pos * te + (1 - ref_axis_pos) * le

    # Get number of spanwise stations (ny)
    ny = mesh.shape[1]

    # Loop over each spanwise station and scale its x coodinates by chord_dist[i]
    for i in range(ny):
        mesh[:, i, 0] = (mesh[:, i, 0] - ref_axis[i, 0]) * chord_dist[i] + ref_axis[i, 0]


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
    # Add the xshear distribution to all x coordinates
    mesh[:, :, 0] += xshear


def shear_y(mesh, yshear):
    """Shear the wing in the y direction (distributed span).

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
    # Add the yshear distribution to all x coordinates
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
    # Add the zshear distribution to all x coordinates
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
    tan_theta = tan(p180 * sweep_angle)

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

    # dx added to mesh x coordinates spanwise.
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
    tan_theta = tan(p180 * dihedral_angle)

    # If symmetric, simply vary the z-coord based on the distance from the
    # center of the wing
    if symmetry:
        y0 = le[-1, 1]
        dz = -(le[:, 1] - y0) * tan_theta

    else:
        ny2 = (num_y - 1) // 2
        y0 = le[ny2, 1]
        dz_right = (le[ny2:, 1] - y0) * tan_theta
        dz_left = -(le[:ny2, 1] - y0) * tan_theta
        dz = np.hstack((dz_left, dz_right))

    # dz added to z coordinates spanwise.
    mesh[:, :, 2] += dz


def stretch(mesh, span, symmetry, ref_axis_pos=0.25):
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

    # Set the span along the reference axis
    le = mesh[0]
    te = mesh[-1]
    ref_axis = ref_axis_pos * te + (1 - ref_axis_pos) * le

    # The user always deals with the full span, so if they input a specific
    # span value and have symmetry enabled, we divide this value by 2.
    if symmetry:
        span /= 2.0

    # Compute the previous span and determine the scalar needed to reach the
    # desired span
    prev_span = ref_axis[-1, 1] - ref_axis[0, 1]
    s = ref_axis[:, 1] / prev_span
    mesh[:, :, 1] = s * span


def taper(mesh, taper_ratio, symmetry, ref_axis_pos=0.25):
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
    ref_axis = ref_axis_pos * te + (1 - ref_axis_pos) * le
    x = ref_axis[:, 1]

    # If symmetric, solve for the correct taper ratio, which is a linear
    # interpolation problem (assume symmetry axis is not necessarily at y = 0)
    if symmetry:
        xp = np.array([x[0], x[-1]])
        fp = np.array([taper_ratio, 1.0])

    # Otherwise, we set up an interpolation problem for the entire wing, which
    # consists of two linear segments (assume symmetry axis is not necessarily at y = 0)
    else:
        xp = np.array([x[0], x[((len(x) + 1) // 2) - 1], x[-1]])
        fp = np.array([taper_ratio, 1.0, taper_ratio])

    # Interpolate over quarter chord line to compute the taper at each spanwise stations
    taper = np.interp(x, xp, fp)

    # Modify the mesh based on the taper amount computed per spanwise section
    # j - spanwise station index (ny)
    # Broadcast taper array over the mesh along spanwise(j) index multiply it by the x and z coordinates
    mesh[:] = np.einsum("ijk, j->ijk", mesh - ref_axis, taper) + ref_axis


def generate_vsp_surfaces(vsp_file, symmetry=False, include=None, scale=1.0):
    """
    Generate a series of VLM surfaces based on geometries in an OpenVSP model.

    Parameters
    ----------
    vsp_file : str
        OpenVSP file to generate meshes from.
    symmetry : bool
        Flag specifying if the full model should be read in (False) or only half (True).
        Half model only reads in right side surfaces.
        Defaults to full model.
    include : list[str]
        List of body names defined in OpenVSP model that should be included in VLM mesh output.
        Defaults to all bodies found in model.
    scale: float
        A global scale factor from the OpenVSP geometry to incoming VLM mesh
        geometry. For example, if the OpenVSP model is in inches, and the VLM
        in meters, scale=0.0254. Defaults to 1.0.

    Returns
    -------
    surfaces : list[dict]
        List of surfaces dictionaries, one (two if symmetry==False) for each body requested in include.
        This is a relatively empty surface dictionary that contains only basic information about the VLM mesh
        (i.e. name, symmetry, mesh).

    """

    if vsp is None:
        raise ImportError("The OpenVSP Python API is required in order to use generate_vsp_surfaces")

    # Check if VSPVehicle class exits
    if hasattr(vsp, "VSPVehicle"):
        # Create a private vehicle geometry instance
        vsp_model = vsp.VSPVehicle()
    # Otherwise use module level API
    # This is less safe since any python module that loads
    # the OpenVSP module has access to our geometry instance
    else:
        vsp_model = vsp

    # Read in file
    vsp_model.ReadVSPFile(vsp_file)

    # Find all vsp bodies
    all_geoms = vsp_model.FindGeoms()

    # If surfaces to include were not specified, we'll output all of them
    if include is None:
        include = []
        for geom_id in all_geoms:
            geom_name = vsp_model.GetContainerName(geom_id)
            if geom_name not in include:
                include.append(geom_name)

    # Create a VSP set that we'll use to identify surfaces we want to output
    for geom_id in all_geoms:
        geom_name = vsp_model.GetContainerName(geom_id)
        if geom_name in include:
            set_flag = True
        else:
            set_flag = False
        vsp_model.SetSetFlag(geom_id, 3, set_flag)

    # Create a degengeom set that will have our VLM surfaces in it
    vsp_model.SetAnalysisInputDefaults("DegenGeom")
    vsp_model.SetIntAnalysisInput("DegenGeom", "WriteCSVFlag", [0], 0)
    vsp_model.SetIntAnalysisInput("DegenGeom", "WriteMFileFlag", [0], 0)
    vsp_model.SetIntAnalysisInput("DegenGeom", "Set", [3], 0)

    # Export all degengeoms to a list
    degen_results_id = vsp_model.ExecAnalysis("DegenGeom")

    # Get all of the degen geom results managers ids
    degen_ids = vsp_model.GetStringResults(degen_results_id, "Degen_DegenGeoms")

    # Create a list of all degen surfaces
    degens = []
    # loop over all degen objects
    for degen_id in degen_ids:
        res = vsp_model.parse_results_object(degen_id)
        degen_obj = dg.DegenGeom(res)

        # Create a degengeom object for the cambersurface
        plate_ids = vsp_model.GetStringResults(degen_id, "plates")
        for plate_id in plate_ids:
            res = vsp_model.parse_results_object(plate_id)
            degen_obj.plates.append(dg.DegenPlate(res))

        degens.append(degen_obj)

    # Loop through each included body and generate a surface dict
    surfaces = {}
    symm_surfaces = []
    for degen in degens:
        if degen.name in include:
            # We found a right surface or a full model was requested
            if degen.surf_index == 0 or symmetry is False:
                flip_normal = degen.flip_normal
                for plate_idx, plate in enumerate(degen.plates):
                    # Some vsp bodies (fuselages) have two surfaces associated with them
                    if len(degen.plates) > 1:
                        surf_name = f"{degen.name}_{plate_idx}"
                    # If there's only one surface (wings) we don't need to append plate id
                    else:
                        surf_name = degen.name
                    # Remove any spaces from name to be OpenMDAO-compatible
                    surf_name = surf_name.replace(" ", "_")
                    # For now, set symmetry to false, we'll update in next step if user requested a half model
                    surf_dict = {"name": surf_name, "symmetry": False}

                    nx = (plate.num_pnts + 1) // 2
                    ny = plate.num_secs
                    mesh = np.zeros([nx, ny, 3])

                    # Extract camber-surface from plate info
                    x = np.array(plate.x) + np.array(plate.nCamber_x) * np.array(plate.zCamber)
                    y = np.array(plate.y) + np.array(plate.nCamber_y) * np.array(plate.zCamber)
                    z = np.array(plate.z) + np.array(plate.nCamber_z) * np.array(plate.zCamber)

                    # Make sure VLM mesh is ordered in right direction
                    if not flip_normal:
                        x = np.flipud(x)
                        y = np.flipud(y)
                        z = np.flipud(z)

                    mesh[:, :, 0] = np.flipud(x.T)
                    mesh[:, :, 1] = np.flipud(y.T)
                    mesh[:, :, 2] = np.flipud(z.T)
                    mesh *= scale

                    # Check if the surface has already been added (i.e. symmetry == False)
                    if surf_name not in surfaces:
                        surf_dict["mesh"] = mesh
                        surfaces[surf_name] = surf_dict
                    # If so, this surface has a left and right segment that must be concatonated
                    else:
                        if degen.surf_index == 0:
                            right_mesh = mesh
                            left_mesh = surfaces[surf_name]["mesh"]
                        else:
                            right_mesh = surfaces[surf_name]["mesh"]
                            left_mesh = mesh
                        new_mesh = np.hstack((left_mesh[:, :-1, :], right_mesh))
                        surfaces[surf_name]["mesh"] = new_mesh

            # We found a left surface, but a half-model was requested, flag the surface as symmetrical
            elif degen.surf_index == 1 and symmetry is True:
                surf_name = degen.name
                surf_name = surf_name.replace(" ", "_")
                symm_surfaces.append(surf_name)

    # If a half-model was requested, go through and flag each surface as symmetrical
    # if a left and right surface was found.
    # NOTE: We don't necessarily want to mark every surface as symmetrical,
    # even if a half-model is requested, since some surfaces, like vertical tails,
    # might lie perfectly on the symmetry plane.
    if symmetry:
        for surf_name in surfaces:
            if surf_name in symm_surfaces:
                surfaces[surf_name]["symmetry"] = True

    # Make sure vsp model is cleared before exit
    vsp_model.ClearVSPModel()

    # Return surfaces as list
    return list(surfaces.values())


def write_FFD_file(surface, mx, my):
    mesh = surface["mesh"]
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

    cushion = 0.5

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

    filename = surface["name"] + "_ffd.fmt"

    with open(filename, "w") as f:
        f.write("1\n")
        f.write("{} {} {}\n".format(mx, 2, my))
        x = np.array_str(ffd[:, :, 0].flatten(order="F"))[1:-1] + "\n"
        y = np.array_str(ffd[:, :, 1].flatten(order="F"))[1:-1] + "\n"
        z = np.array_str(ffd[:, :, 2].flatten(order="F"))[1:-1] + "\n"

        f.write(x)
        f.write(y)
        f.write(z)

    return filename


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
    file_handle = open(file_name, "r")
    num_panels = int(file_handle.readline())
    # Get the multi-block dimensions of every included surface
    block_dims = file_handle.readline().split()

    # Now loop through remainder of file and pluck out mesh node locations
    mesh_list = []
    mesh_dict = {}
    for i in range(num_panels):
        [nx, ny, nz] = block_dims[3 * i : 3 * i + 3]
        # Use nx and ny to intialize mesh. Since these are surfaces nz always
        # equals 1, so no need to use it
        mesh = np.zeros(int(nx) * int(ny) * 3)

        for j in range(mesh.size):
            line = file_handle.readline()
            val = float(line)
            if np.abs(val) < zero_tol:
                val = 0
            mesh[j] = val

        # Restructure mesh as 3D array,
        # Plot3D files are always written using Fortran order
        mesh_list.append(mesh.reshape([int(nx), int(ny), 3], order="f"))

    # Now read in names for each surface mesh
    for i in range(num_panels):
        name = file_handle.readline()[:-1]
        mesh_dict[name] = mesh_list[i]

    return mesh_dict


def build_section_dicts(surface):
    """This utility function takes a multi-section surface dictionary and outputs a list
    of individual section surface dictionaries so the geometry group for each individual
    section can be initialized.

    Parameters
    ----------
    surface: dict
        OpenAeroStruct multi-section surface dictionary

    Returns
    -------
    section_surfaces : list
        List of OpenAeroStruct surface dictionaries for each individual surface

    """
    # Get number of sections
    num_sections = surface["num_sections"]

    if surface["meshes"] == "gen-meshes":
        # Verify that all required inputs for automatic mesh generation are provided for each section
        if len(surface["ny"]) != num_sections:
            raise ValueError("Number of spanwise points needs to be provided for each section")
        if len(surface["taper"]) != num_sections:
            raise ValueError("Taper needs to be provided for each section")
        if len(surface["span"]) != num_sections:
            raise ValueError("Span needs to be provided for each section")
        if len(surface["sweep"]) != num_sections:
            raise ValueError("Sweep needs to be provided for each section")

        # Generate unified and individual section meshes
        _, sec_meshes = generate_section_mesh(surface)
    else:
        # Allow user to provide mesh for each section
        if len(surface["meshes"]) != num_sections:
            raise ValueError("A mesh needs to be provided for each section.")
        sec_meshes = surface["meshes"]

    if len(surface["sec_name"]) != num_sections:
        raise ValueError("A name needs to be provided for each section.")

    # List of support keys for multi-section wings
    # NOTE: make sure this is consistent to the documentation's surface dict page
    target_keys = [
        # Essential Info
        "num_section",
        "symmetry",
        "S_ref_type",
        "ref_axis_pos",
        # wing definition
        "span",
        "taper",
        "sweep",
        "dihedral",
        "twist_cp",
        "chord_cp",
        "xshear_cp",
        "yshear_cp",
        "zshear_cp",
        # aerodynamics
        "CL0",
        "CD0",
        "with_viscous",
        "with_wave",
        "groundplane",
        "k_lam",
        "t_over_c_cp",
        "c_max_t",
    ]

    # Constructs a list of section dictionaries and adds the specified supported keys and values from the mult-section surface dictionary.
    surface_sections = []
    num_sections = surface["num_sections"]

    for i in range(num_sections):
        section = {}
        for k in set(surface).intersection(target_keys):
            if type(surface[k]) is list:
                # Reset taper, sweep, and span so that OAS doesn't apply the the transformations again
                if k == "taper":
                    section[k] = 1.0
                elif k == "sweep":
                    section[k] = 0.0
                elif k == "span":
                    if surface["symmetry"]:
                        section[k] = 2.0 * surface[k][i]
                    else:
                        section[k] = surface[k][i]
                else:
                    section[k] = surface[k][i]
            else:
                section[k] = surface[k]
        section["mesh"] = sec_meshes[i]
        section["name"] = surface["sec_name"][i]
        surface_sections.append(section)
    return surface_sections


def unify_mesh(sections, shift_uni_mesh=True):
    """
    Function that produces a unified mesh from all the individual wing section meshes.

    Parameters
    ----------
    sections : list
        List of section OpenAeroStruct surface dictionaries

    shift_uni_mesh : bool
        Flag that shifts sections so that their leading edges are coincident. Intended to keep sections from seperating
        or intersecting during scalar span or sweep operations without the use of the constraint component.

    Returns
    -------
    uni_mesh : numpy array
        Unfied surface mesh in OAS format
    """
    for i_sec in np.arange(0, len(sections) - 1):
        mesh = sections[i_sec]["mesh"]

        if i_sec == 0:
            uni_mesh = copy.deepcopy(mesh[:, :-1, :])
        else:
            if shift_uni_mesh:
                # translate or shift uni_mesh (outer sections) to align leading edge at unification boundary
                last_mesh = sections[i_sec - 1]["mesh"]
                uni_mesh = uni_mesh - last_mesh[0, -1, :] + mesh[0, 0, :]

            uni_mesh = np.concatenate([uni_mesh, mesh[:, :-1, :]], axis=1)

    # Stitch the results into a singular mesh
    mesh = sections[len(sections) - 1]["mesh"]
    if len(sections) == 1:
        uni_mesh = copy.deepcopy(mesh)
    else:
        uni_mesh = np.concatenate([uni_mesh, mesh], axis=1)

    return uni_mesh


def build_multi_spline(out_name, num_sections, control_points):
    """This function returns an OpenMDAO Independent Variable Component with an output vector appropriately
    named and sized to function as an unified set of B-spline control poitns that join multiple sections by construction.

    Parameters
    ----------
    out_name: string
        Name of the output to assign to the B-spline
    num_sections : int
        Number of sections
    control_points: list
        List of B-spline control point arrays corresponding to each section

    Returns
    -------
    spline_control : OpenMDAO component object
        The unified B-spline control point indpendent variable component

    """
    if len(control_points) != num_sections:
        raise Exception("Target sections need to match with control points!")

    single_sections = len([cp for cp in control_points if len(cp) == 1])

    control_poin_vec = np.ones(len(np.concatenate(control_points)) - (num_sections - 1 - single_sections))

    spline_control = om.IndepVarComp()
    spline_control.add_output("{}_spline".format(out_name), val=control_poin_vec)

    return spline_control


def connect_multi_spline(prob, section_surfaces, sec_cp, out_name, comp_name, geom_name, return_bind_inds=False):
    """This function connects the the unified B-spline component with the individual B-splines
    of each section. There is a point of overlap at each section so that each edge control point control the edge
    controls points of each section's B-spline. This is how section joining by consturction is acheived.
    An issue occurs however when a B-spline in a particular section only has one control point. In this case the one
    section control point is bound to the left edge B-spline component control point. As result, there is nothing to
    maintain C0 continuity with the next section. As result a constraint will need to be manually set. To facilitate this,
    the array bind_inds will contain a list of the B-spline control point indicies that will need to be manually constrained to
    their previous sections.


    Parameters
    ----------
    prob : OpenMDAO problem object
        The OpenAeroStruct problem object with the unified B-spline component added.
    section_surfaces : list
        List of the surface dictionaries for each section.
    sec_cp : list
        List of B-spline control point arrays for each section.
    out_name: string
        Name of the unified B-spline component output to connect from
    comp_name: string
        Name of the unified B-spline component added to the problem object
    geom_name : string
        Name of the multi-section geometry group
    return_bind_inds: bool
        Return list of unjoined unified B-spline inidices. Default is False.

    Returns
    -------
    bind_inds : list
        List of unified B-spline control point indicies not connected due to the presence of a single control point section.(Only if return bind_inds specified)

    """
    acc = 0
    bind_inds = []
    for i, section in enumerate(section_surfaces):
        point_count = len(sec_cp[i])
        src_inds = np.arange(acc, acc + point_count)
        acc += point_count - 1
        if point_count == 1:
            acc += 1
            bind_inds.append(acc)
        prob.model.connect(
            "{}.{}".format(comp_name, out_name) + "_spline",
            geom_name + "." + section["name"] + ".{}".format(out_name),
            src_indices=src_inds,
        )

    if return_bind_inds:
        return bind_inds


def generate_mesh(input_dict):
    warnings.warn(
        "generate_mesh has been moved to mesh_generator.py. Importing from utils.py is deprecated and will be removed in a future release.",
        DeprecationWarning,
        stacklevel=2,
    )
    return _generate_mesh(input_dict)


def gen_rect_mesh(num_x, num_y, span, chord, span_cos_spacing=0.0, chord_cos_spacing=0.0):
    warnings.warn(
        "gen_rect_mesh has been moved to mesh_generator.py. Importing from utils.py is deprecated and will be removed in a future release.",
        DeprecationWarning,
        stacklevel=2,
    )
    return _gen_rect_mesh(num_x, num_y, span, chord, span_cos_spacing=0.0, chord_cos_spacing=0.0)


def gen_crm_mesh(num_x, num_y, span_cos_spacing=0.0, chord_cos_spacing=0.0, wing_type="CRM:jig"):
    warnings.warn(
        "gen_crm_mesh has been moved to mesh_generator.py. Importing from utils.py is deprecated and will be removed in a future release.",
        DeprecationWarning,
        stacklevel=2,
    )
    return _gen_crm_mesh(num_x, num_y, span_cos_spacing=0.0, chord_cos_spacing=0.0, wing_type="CRM:jig")


def add_chordwise_panels(num_x, num_y, span_cos_spacing=0.0, chord_cos_spacing=0.0, wing_type="CRM:jig"):
    warnings.warn(
        "add_chordwise_panels has been moved to mesh_generator.py and renamed to regen_chordwise_panels. Importing from utils.py is deprecated and will be removed in a future release.",
        DeprecationWarning,
        stacklevel=2,
    )
    return _regen_chordwise_panels(num_x, num_y, span_cos_spacing=0.0, chord_cos_spacing=0.0, wing_type="CRM:jig")


def add_chordwise_panels(num_x, num_y, span_cos_spacing=0.0, chord_cos_spacing=0.0, wing_type="CRM:jig"):
    warnings.warn(
        "add_chordwise_panels has been moved to mesh_generator.py and renamed to regen_chordwise_panels. Importing from utils.py is deprecated and will be removed in a future release.",
        DeprecationWarning,
        stacklevel=2,
    )
    return _regen_chordwise_panels(num_x, num_y, span_cos_spacing=0.0, chord_cos_spacing=0.0, wing_type="CRM:jig")


def get_default_geo_dict():
    warnings.warn(
        "get_default_geo_dict has been moved to mesh_generator.py and renamed to regen_chordwise_panels. Importing from utils.py is deprecated and will be removed in a future release.",
        DeprecationWarning,
        stacklevel=2,
    )
    return _get_default_geo_dict()


def writeMesh():
    warnings.warn(
        "writeMesh has been moved to mesh_generator.py and renamed to write_tecplot. Importing from utils.py is deprecated and will be removed in a future release.",
        DeprecationWarning,
        stacklevel=2,
    )
    return _write_tecplot()


def getFullMesh():
    warnings.warn(
        "getFullMesh has been moved to mesh_generator.py. Importing from utils.py is deprecated and will be removed in a future release.",
        DeprecationWarning,
        stacklevel=2,
    )
    return _getFullMesh()
