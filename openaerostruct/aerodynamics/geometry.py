import numpy as np

import openmdao.api as om


np.random.seed(314)

class VLMGeometry(om.ExplicitComponent):
    """
    Compute various geometric properties for VLM analysis.
    These are used primarily to help compute postprocessing quantities,
    such as wave CD, viscous CD, etc.
    Some of the quantities, like `normals`, are used to compute the RHS
    of the AIC linear system.

    parameters
    ----------
    def_mesh[nx, ny, 3] : numpy array
        Array defining the nodal coordinates of the lifting surface.

    Returns
    -------
    b_pts[nx-1, ny, 3] : numpy array
        Bound points for the horseshoe vortices, found along the 1/4 chord.
    widths[ny-1] : numpy array
        The spanwise widths of each individual panel.
    lengths[ny] : numpy array
        The chordwise length of the entire airfoil following the camber line.
    chords[ny] : numpy array
        The chordwise distance between the leading and trailing edges.
    normals[nx-1, ny-1, 3] : numpy array
        The normal vector for each panel, computed as the cross of the two
        diagonals from the mesh points.
    S_ref : float
        The reference area of the lifting surface.
    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        self.surface = surface = self.options['surface']

        mesh=surface['mesh']
        nx = self.nx = mesh.shape[0]
        ny = self.ny = mesh.shape[1]

        # All of these computations only need the deformed mesh
        self.add_input('def_mesh', val=np.zeros((nx, ny, 3)), units='m')

        self.add_output('b_pts', val=np.random.random((nx-1, ny, 3)), units='m')
        self.add_output('widths', val=np.ones((ny-1)), units='m')
        self.add_output('cos_sweep', val=np.zeros((ny-1)), units='m')
        self.add_output('lengths', val=np.zeros((ny)), units='m')
        self.add_output('chords', val=np.zeros((ny)), units='m')
        self.add_output('normals', val=np.zeros((nx-1, ny-1, 3)))
        self.add_output('S_ref', val=1., units='m**2')

        # Next up we have a lot of rows and cols settings for the sparse
        # Jacobians. Each set of partials needs a different rows/cols setup

        # b_pts
        size = (nx-1) * ny * 3
        base = np.arange(size)
        rows = np.tile(base, 2)
        cols = rows + np.repeat([0, ny*3], len(base))
        val = np.empty((2*size, ))
        val[:size] = 0.75
        val[size:] = 0.25
        self.declare_partials('b_pts', 'def_mesh', rows=rows, cols=cols, val=val)

        # widths
        size = ny - 1
        base = np.arange(size)
        rows = np.tile(base, 12)
        col = np.tile(3*base, 6) + np.repeat(np.arange(6), len(base))
        cols = np.tile(col, 2) + np.repeat([0, (nx-1)*ny*3], len(col))
        self.declare_partials('widths', 'def_mesh', rows=rows, cols=cols)

        # cos_sweep
        rows = np.tile(base, 8)
        col = np.tile(3*base, 4) + np.repeat([1, 2, 4, 5], len(base))
        cols = np.tile(col, 2) + np.repeat([0, (nx-1)*ny*3], len(col))
        self.declare_partials('cos_sweep', 'def_mesh', rows=rows, cols=cols)

        # lengths
        size = ny
        base = np.arange(size)
        rows = np.tile(base, nx * 3)
        col = np.tile(3*base, 3) + np.repeat(np.arange(3), len(base))
        cols = np.tile(col, nx) + np.repeat(3*ny*np.arange(nx), len(col))
        self.declare_partials('lengths', 'def_mesh', rows=rows, cols=cols)

        # chords
        rows = np.tile(base, 6)
        col = np.tile(3*base, 3) + np.repeat(np.arange(3), len(base))
        cols = np.tile(col, 2) + np.repeat([0, (nx-1)*ny*3], len(col))
        self.declare_partials('chords', 'def_mesh', rows=rows, cols=cols)

        # normals
        size = (ny-1)*(nx-1)*3
        row = np.tile(np.arange(size).reshape((size, 1)), 3).flatten()
        rows = np.tile(row, 4)
        base = np.tile(np.arange(3), size) + np.repeat(3*np.arange(size//3), 9)
        base += np.repeat(3*np.arange(nx-1), 9*(ny-1))
        cols = np.concatenate([
            base + 3,
            base + ny*3,
            base,
            base + (ny+1)*3
        ])
        self.declare_partials('normals', 'def_mesh', rows=rows, cols=cols)

        # And here actually all parts of the mesh influence the area, so it's
        # fully dense
        self.declare_partials('S_ref', 'def_mesh')
        self.set_check_partial_options(wrt='def_mesh', method='fd', step=1e-6)

    def compute(self, inputs, outputs):
        mesh = inputs['def_mesh']

        # Compute the bound points at quarter-chord
        b_pts = mesh[:-1, :, :] * .75 + mesh[1:, :, :] * .25

        # Compute the widths of each panel at the quarter-chord line
        quarter_chord = 0.25 * mesh[-1] + 0.75 * mesh[0]
        widths = np.linalg.norm(quarter_chord[1:, :] - quarter_chord[:-1, :], axis=1)

        # Compute the numerator of the cosine of the sweep angle of each panel
        # (we need this for the viscous drag dependence on sweep, and we only compute
        # the numerator because the denominator of the cosine fraction is the width,
        # which we have already computed. They are combined in the viscous drag
        # calculation.)
        cos_sweep = np.linalg.norm(quarter_chord[1:, [1,2]] - quarter_chord[:-1, [1,2]], axis=1)

        # Compute the length of each chordwise set of mesh points through the camber line.
        dx = mesh[1:, :, 0] - mesh[:-1, :, 0]
        dy = mesh[1:, :, 1] - mesh[:-1, :, 1]
        dz = mesh[1:, :, 2] - mesh[:-1, :, 2]
        lengths = np.sum(np.sqrt(dx**2 + dy**2 + dz**2), axis=0)

        # Compute the normal of each panel by taking the cross-product of
        # its diagonals. Note that this could be a nonplanar surface
        normals = np.cross(
            mesh[:-1,  1:, :] - mesh[1:, :-1, :],
            mesh[:-1, :-1, :] - mesh[1:,  1:, :],
            axis=2)

        # Normalize the normal vectors
        norms = np.sqrt(np.sum(normals**2, axis=2))
        for j in range(3):
            normals[:, :, j] /= norms

        # Compute the wetted surface area
        if self.surface['S_ref_type'] == 'wetted':
            S_ref = 0.5 * np.sum(norms)

        # Compute the projected surface area
        elif self.surface['S_ref_type'] == 'projected':
            proj_mesh = mesh.copy()
            proj_mesh[: , :, 2] = 0.
            proj_normals = np.cross(
                proj_mesh[:-1,  1:, :] - proj_mesh[1:, :-1, :],
                proj_mesh[:-1, :-1, :] - proj_mesh[1:,  1:, :],
                axis=2)

            proj_norms = np.sqrt(np.sum(proj_normals**2, axis=2))
            for j in range(3):
                proj_normals[:, :, j] /= proj_norms

            S_ref = 0.5 * np.sum(proj_norms)

        # Multiply the surface area by 2 if symmetric to get consistent area measures
        if self.surface['symmetry']:
            S_ref *= 2

        # Compute the chord for each spanwise portion.
        # This is the distance from the leading to trailing edge.
        chords = np.linalg.norm(mesh[0, :, :] - mesh[-1, :, :], axis=1)

        # Store each array in the outputs dict
        outputs['b_pts'] = b_pts
        outputs['widths'] = widths
        outputs['cos_sweep'] = cos_sweep
        outputs['lengths'] = lengths
        outputs['normals'] = normals
        outputs['S_ref'] = S_ref
        outputs['chords'] = chords

    def compute_partials(self, inputs, partials):
        """ Jacobian for VLM geometry."""

        nx = self.nx
        ny = self.ny
        mesh = inputs['def_mesh']

        # Compute the widths of each panel at the quarter-chord line
        quarter_chord = 0.25 * mesh[-1] + 0.75 * mesh[0]
        widths = np.linalg.norm(quarter_chord[1:, :] - quarter_chord[:-1, :], axis=1)

        # Compute the cosine of the sweep angle of each panel
        cos_sweep_array = np.linalg.norm(quarter_chord[1:, [1, 2]] - \
                                         quarter_chord[:-1, [1, 2]], axis=1)

        delta = np.diff(quarter_chord, axis=0).T
        d1 = delta / widths
        partials['widths', 'def_mesh'] = np.outer([-0.75, 0.75, -0.25, 0.25],
                                                  d1.flatten()).flatten()
        d1 = delta[1:, :] / cos_sweep_array
        partials['cos_sweep', 'def_mesh'] = np.outer([-0.75, 0.75, -0.25, 0.25],
                                                     d1.flatten()).flatten()

        partials['lengths', 'def_mesh'][:] = 0.0
        dmesh = np.diff(mesh, axis=0)
        l = np.sqrt(np.sum(dmesh**2, axis=2))
        dmesh = dmesh / l[:, :, np.newaxis]
        derivs = np.transpose(dmesh, axes=[0, 2, 1]).flatten()
        nn = len(derivs)
        partials['lengths', 'def_mesh'][:nn] -= derivs
        partials['lengths', 'def_mesh'][-nn:] += derivs

        dfullmesh = mesh[0, :] - mesh[-1, :]
        l = np.sqrt(np.sum(dfullmesh**2, axis=1))
        derivs = (dfullmesh.T/l).flatten()
        partials['chords', 'def_mesh'] = np.concatenate([derivs, -derivs])

        # f = c / n
        a = mesh[:-1, 1:, :] - mesh[1:, :-1, :]
        b = mesh[:-1, :-1, :] - mesh[1:, 1:, :]
        c = np.cross(a, b, axis=2)
        n = np.sqrt(np.sum(c**2, axis=2))

        # Now let's work backwards to get derivative
        # dfdc = (dcdc * n - c * dndc) / n**2
        dndc = c / n[:, :, np.newaxis]
        dcdc = np.zeros((nx-1, ny-1, 3, 3))
        dcdc[:, :, 0, 0] = 1.0
        dcdc[:, :, 1, 1] = 1.0
        dcdc[:, :, 2, 2] = 1.0

        # dfdc is now a 3x3 jacobian with f along the rows and c along the columns
        dfdc = (dcdc*n[:, :, np.newaxis, np.newaxis] - np.einsum('ijk,ijl->ijkl', c, dndc)) / (n**2)[:, :, np.newaxis, np.newaxis]

        # The next step is to get dcda and dcdb, both of which will be
        # 3x3 jacobians with c along the rows

        dcda = np.zeros((nx-1, ny-1, 3, 3))
        dcda[:, :, 0, 1] = b[:, :, 2]
        dcda[:, :, 0, 2] = -b[:, :, 1]
        dcda[:, :, 1, 0] = -b[:, :, 2]
        dcda[:, :, 1, 2] = b[:, :, 0]
        dcda[:, :, 2, 0] = b[:, :, 1]
        dcda[:, :, 2, 1] = -b[:, :, 0]

        dcdb = np.zeros((nx-1, ny-1, 3, 3))
        dcdb[:, :, 0, 1] = -a[:, :, 2]
        dcdb[:, :, 0, 2] = a[:, :, 1]
        dcdb[:, :, 1, 0] = a[:, :, 2]
        dcdb[:, :, 1, 2] = -a[:, :, 0]
        dcdb[:, :, 2, 0] = -a[:, :, 1]
        dcdb[:, :, 2, 1] = a[:, :, 0]

        # Now let's do some matrix multiplication to get dfda and dfdb
        dfda = np.einsum('ijkl,ijlm->ijkm', dfdc, dcda)
        dfdb = np.einsum('ijkl,ijlm->ijkm', dfdc, dcdb)

        # Aside: preparation for surface area deriv computation under 'projected' option.
        if self.surface['S_ref_type'] == 'projected':
            # Compute the projected surface area by zeroing out z dimension.
            proj_mesh = mesh.copy()
            proj_mesh[: , :, 2] = 0.

            a = proj_mesh[:-1, 1:, :] - proj_mesh[1:, :-1, :]
            b = proj_mesh[:-1, :-1, :] - proj_mesh[1:, 1:, :]
            c = np.cross(a, b, axis=2)
            n = np.sqrt(np.sum(c**2, axis=2))

            dcda[:, :, 0, 1] = b[:, :, 2]
            dcda[:, :, 0, 2] = -b[:, :, 1]
            dcda[:, :, 1, 0] = -b[:, :, 2]
            dcda[:, :, 1, 2] = b[:, :, 0]
            dcda[:, :, 2, 0] = b[:, :, 1]
            dcda[:, :, 2, 1] = -b[:, :, 0]

            dcdb[:, :, 0, 1] = -a[:, :, 2]
            dcdb[:, :, 0, 2] = a[:, :, 1]
            dcdb[:, :, 1, 0] = a[:, :, 2]
            dcdb[:, :, 1, 2] = -a[:, :, 0]
            dcdb[:, :, 2, 0] = -a[:, :, 1]
            dcdb[:, :, 2, 1] = a[:, :, 0]

        # Need these for wetted surface area derivs.
        dsda = np.einsum("ijk,ijkl->ijl", c, dcda) / n[:, :, np.newaxis]
        dsdb = np.einsum("ijk,ijkl->ijl", c, dcdb) / n[:, :, np.newaxis]

        # Note: this is faster than np.concatenate for large meshes.
        nn = (nx-1)*(ny-1)*9
        dfda_flat = dfda.flatten()
        dfdb_flat = dfdb.flatten()
        partials['normals', 'def_mesh'][:nn] = dfda_flat
        partials['normals', 'def_mesh'][nn:2*nn] = -dfda_flat
        partials['normals', 'def_mesh'][2*nn:3*nn] = dfdb_flat
        partials['normals', 'def_mesh'][3*nn:4*nn] = -dfdb_flat

        # At this point, same calculation for wetted and projected surface.
        dsda_flat = 0.5*dsda.flatten()
        dsdb_flat = 0.5*dsdb.flatten()
        idx = np.arange((nx-1)*(ny-1)*3) + np.repeat(3*np.arange(nx-1), 3*(ny-1))
        partials['S_ref', 'def_mesh'][:] = 0.0
        partials['S_ref', 'def_mesh'][:, idx+3] += dsda_flat
        partials['S_ref', 'def_mesh'][:, idx+ny*3] -= dsda_flat
        partials['S_ref', 'def_mesh'][:, idx] += dsdb_flat
        partials['S_ref', 'def_mesh'][:, idx+(ny+1)*3] -= dsdb_flat

        # Multiply the surface area by 2 if symmetric to get consistent area measures
        if self.surface['symmetry']:
            partials['S_ref', 'def_mesh'] *= 2.0
