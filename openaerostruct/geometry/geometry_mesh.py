""" Group that manipulates geometry mesh based on high-level design parameters. """

import numpy as np

import openmdao.api as om

from openaerostruct.geometry.geometry_mesh_transformations import \
     Taper, ScaleX, Sweep, ShearX, Stretch, ShearY, Dihedral, \
     ShearZ, Rotate


class GeometryMesh(om.Group):
    """
    OpenMDAO group that performs mesh manipulation functions. It reads in
    the initial mesh from the surface dictionary and outputs the altered
    mesh based on the geometric design variables.

    Depending on the design variables selected or the supplied geometry information,
    only some of the follow parameters will actually be given to this component.
    If parameters are not active (they do not deform the mesh), then
    they will not be given to this component.

    Parameters
    ----------
    sweep : float
        Shearing sweep angle in degrees.
    dihedral : float
        Dihedral angle in degrees.
    twist[ny] : numpy array
        1-D array of rotation angles for each wing slice in degrees.
    chord_dist[ny] : numpy array
        Chord length for each panel edge.
    taper : float
        Taper ratio for the wing; 1 is untapered, 0 goes to a point at the tip.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Modified mesh based on the initial mesh in the surface dictionary and
        the geometric design variables.
    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        surface = self.options['surface']

        mesh = surface['mesh']
        ny = mesh.shape[1]
        mesh_shape = mesh.shape
        symmetry = surface['symmetry']

        # This flag determines whether or not changes in z (dihedral) add an
        # additional rotation matrix to modify the twist direction
        self.rotate_x = True

        # 1. Taper

        if 'taper' in surface:
            val = surface['taper']
            promotes = ['taper']
        else:
            val = 1.
            promotes = []

        self.add_subsystem('taper', Taper(val=val, mesh=mesh, symmetry=symmetry),
                           promotes_inputs=promotes)

        # 2. Scale X

        val = np.ones(ny)
        if 'chord_cp' in surface:
            promotes = ['chord']
        else:
            promotes = []

        self.add_subsystem('scale_x', ScaleX(val=val, mesh_shape=mesh_shape),
                           promotes_inputs=promotes)

        # 3. Sweep

        if 'sweep' in surface:
            val = surface['sweep']
            promotes = ['sweep']
        else:
            val = 0.
            promotes = []

        self.add_subsystem('sweep', Sweep(val=val, mesh_shape=mesh_shape, symmetry=symmetry),
                           promotes_inputs=promotes)

        # 4. Shear X

        val = np.zeros(ny)
        if 'xshear_cp' in surface:
            promotes = ['xshear']
        else:
            promotes = []

        self.add_subsystem('shear_x', ShearX(val=val, mesh_shape=mesh_shape),
                           promotes_inputs=promotes)

        # 5. Stretch

        if 'span' in surface:
            promotes = ['span']
            val = surface['span']
        else:
            # Compute span. We need .real to make span to avoid OpenMDAO warnings.
            quarter_chord = 0.25 * mesh[-1, :, :] + 0.75 * mesh[0, :, :]
            span = max(quarter_chord[:, 1]).real - min(quarter_chord[:, 1]).real
            if symmetry:
                span *= 2.
            val = span
            promotes = []

        self.add_subsystem('stretch', Stretch(val=val, mesh_shape=mesh_shape, symmetry=symmetry),
                           promotes_inputs=promotes)

        # 6. Shear Y

        val = np.zeros(ny)
        if 'yshear_cp' in surface:
            promotes = ['yshear']
        else:
            promotes = []

        self.add_subsystem('shear_y', ShearY(val=val, mesh_shape=mesh_shape),
                           promotes_inputs=promotes)

        # 7. Dihedral

        if 'dihedral' in surface:
            val = surface['dihedral']
            promotes = ['dihedral']
        else:
            val = 0.
            promotes = []

        self.add_subsystem('dihedral', Dihedral(val=val, mesh_shape=mesh_shape, symmetry=symmetry),
                           promotes_inputs=promotes)

        # 8. Shear Z

        val = np.zeros(ny)
        if 'zshear_cp' in surface:
            promotes = ['zshear']
        else:
            promotes = []

        self.add_subsystem('shear_z', ShearZ(val=val, mesh_shape=mesh_shape),
                           promotes_inputs=promotes)

        # 9. Rotate

        val = np.zeros(ny)
        if 'twist_cp' in surface:
            promotes = ['twist']
        else:
            val = np.zeros(ny)
            promotes = []

        self.add_subsystem('rotate', Rotate(val=val, mesh_shape=mesh_shape, symmetry=symmetry),
                           promotes_inputs=promotes, promotes_outputs=['mesh'])


        names = ['taper', 'scale_x', 'sweep', 'shear_x', 'stretch', 'shear_y', 'dihedral',
                 'shear_z', 'rotate']

        for j in np.arange(len(names) - 1):
            self.connect(names[j] + '.mesh', names[j+1] + '.in_mesh')
