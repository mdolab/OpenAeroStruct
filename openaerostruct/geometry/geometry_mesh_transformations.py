"""A set of components that manipulate geometry mesh
based on high-level design parameters."""

import numpy as np

import openmdao.api as om


class Taper(om.ExplicitComponent):
    """
    OpenMDAO component that manipulates the mesh by altering the spanwise chord linearly to produce
    a tapered wing. Note that we apply taper around the reference axis line which is the quarter-chord by default.

    Parameters
    ----------
    taper : float
        Taper ratio for the wing; 1 is untapered, 0 goes to a point at the tip.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the tapered aerodynamic surface.
    """

    def initialize(self):
        """
        Declare options.
        """
        self.options.declare("val", desc="Initial value for the taper ratio.")
        self.options.declare("mesh", desc="Nodal mesh defining the initial aerodynamic surface.")
        self.options.declare(
            "symmetry", default=False, desc="Flag set to true if surface is reflected about y=0 plane."
        )
        self.options.declare(
            "ref_axis_pos",
            default=0.25,
            desc="Fraction of the chord to use as the reference axis",
        )

    def setup(self):
        mesh = self.options["mesh"]
        val = self.options["val"]
        self.ref_axis_pos = self.options["ref_axis_pos"]

        self.add_input("taper", val=val)

        self.add_output("mesh", val=mesh, units="m")

        self.declare_partials("*", "*")

    def compute(self, inputs, outputs):
        mesh = self.options["mesh"]
        symmetry = self.options["symmetry"]
        taper_ratio = inputs["taper"][0]

        # Get mesh parameters and the quarter-chord
        le = mesh[0]
        te = mesh[-1]
        ref_axis = self.ref_axis_pos * te + (1 - self.ref_axis_pos) * le
        x = ref_axis[:, 1]

        # Spanwise(j) index of wing centerline
        n_sym = (len(x) + 1) // 2 - 1

        # If symmetric, solve for the correct taper ratio, which is a linear
        # interpolation problem (assume symmetry axis is not necessarily at y = 0)
        if symmetry:
            xp = np.array([x[0], x[-1]])
            fp = np.array([taper_ratio, 1.0])

        # Otherwise, we set up an interpolation problem for the entire wing, which
        # consists of two linear segments (assume symmetry axis is not necessarily at y = 0)
        else:
            xp = np.array([x[0], x[n_sym], x[-1]])
            fp = np.array([taper_ratio, 1.0, taper_ratio])

        # Interpolate over quarter chord line to compute the taper at each spanwise stations
        taper = np.interp(x, xp, fp)

        # Modify the mesh based on the taper amount computed per spanwise section
        # j - spanwise station index (ny)
        # Broadcast taper array over the mesh along spanwise(j) index multiply it by the x and z coordinates
        outputs["mesh"] = np.einsum("ijk,j->ijk", mesh - ref_axis, taper) + ref_axis

    def compute_partials(self, inputs, partials):
        mesh = self.options["mesh"]
        symmetry = self.options["symmetry"]

        # Get mesh parameters and the quarter-chord
        le = mesh[0]
        te = mesh[-1]
        ref_axis = self.ref_axis_pos * te + (1 - self.ref_axis_pos) * le
        x = ref_axis[:, 1]

        # Spanwise(j) index of wing centerline
        n_sym = (len(x) + 1) // 2 - 1

        # Derivative implementation that allows for taper_ratio = 1
        if symmetry:
            # Compute the span
            span = x[-1] - x[0]

            # Distance of each station from left tip(incl. left tip)
            dy = x - x[0]

            # Compute the derivative vector wrt to the taper_ratio
            # Note that this isn't sensitive to the taper_ratio itself,
            # only the span station spacing allowing for taper_ratio = 1
            # This is simply the derivative of the linear interpolation
            # wrt to the end point which is the taper_ratio
            dtaper = np.ones(len(x)) + (-dy / span)
        else:
            # Compute the semi-span considering each semi-span might be
            # perturbed
            span1 = x[n_sym] - x[0]

            # Distance of each left span station from left tip(incl. left tip)
            dy1 = x[: n_sym + 1] - x[0]

            # Compute the left half of the derivative vector wrt to the taper_ratio
            dtaper1 = np.ones(n_sym + 1) + (-dy1 / span1)

            # Compute the semi-span
            span2 = x[-1] - x[n_sym]

            # Distance of each right span station from centerline
            dy2 = x[n_sym + 1 :] - x[n_sym]

            # Compute the right half of the derivative vector wrt to the taper_ratio
            dtaper2 = dy2 / span2

            # Concatinate the two parts of the deritivative vector
            dtaper = np.concatenate([dtaper1, dtaper2])

        # Broadcast d (taper)/ d(taper_ratio) onto each mesh spanwise station in a similar fasion to the compute method
        # This works as only taper is directly sensitive to taper_ratio
        partials["mesh", "taper"] = np.einsum("ijk, j->ijk", mesh - ref_axis, dtaper)


class ScaleX(om.ExplicitComponent):
    """
    OpenMDAO component that manipulates the mesh by modifying the chords along the span of the
    wing by scaling only the x-coord.

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the initial aerodynamic surface.
    chord[ny] : numpy array
        Spanwise distribution of the chord scaler.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh with the new chord lengths.
    """

    def initialize(self):
        """
        Declare options.
        """
        self.options.declare("val", desc="Initial value for chord lengths")
        self.options.declare("mesh_shape", desc="Tuple containing mesh shape (nx, ny).")
        self.options.declare(
            "ref_axis_pos",
            default=0.25,
            desc="Fraction of the chord to use as the reference axis",
        )

    def setup(self):
        mesh_shape = self.options["mesh_shape"]
        val = self.options["val"]
        self.ref_axis_pos = self.options["ref_axis_pos"]
        self.add_input("chord", units=None, val=val)
        self.add_input("in_mesh", shape=mesh_shape, units="m")

        self.add_output("mesh", shape=mesh_shape, units="m")

        # Compute total number of array entries in mesh array
        nx, ny, _ = mesh_shape
        nn = nx * ny * 3

        # Setup the  d mesh/ d chord jacobian

        # All mesh array entries are sensitive to chord
        rows = np.arange(nn)

        # Repeat each spanwise index 3 times since all three coordiantes of each
        # spanwise point is sentive to chord
        # col = np.tile(np.zeros(3), ny) + np.repeat(np.arange(ny), 3)
        # Removed redundant preallocation step
        col = np.repeat(np.arange(ny), 3)

        # At each spanwise station there are nx chorwise point so repeat
        # the pattern nx times
        cols = np.tile(col, nx)

        self.declare_partials("mesh", "chord", rows=rows, cols=cols)

        # Setup the  d mesh/ d in_mesh jacobian

        # Diagonal part of jacobian. Mesh maps directly to in_mesh at first.
        p_rows = np.arange(nn)

        # Off-diagonal part of the jacobian. Off-diagonal part exists as we translate the mesh to its
        # references axis prior to applying the chord distribution. The ref_axis position itself is sensitive
        # to the mesh LE and TE. The LE and TE parts of the mesh are already part of the main diagonal so the off diagonal
        # terms the conver the sensitivies of the remainder of the mesh to the ref_axis.

        # Entries sensitive to trailing edge contribution of ref_axis location. Note that
        # the last row(TE) is dropped here as it's entries are covered as part
        # of the main diagonal.
        te_rows = np.arange(((nx - 1) * ny * 3))

        # Entries sensitive to leading edge contribution of ref_axis location. This is an offset
        # of the te_rows but really it's the entire mesh except the leading edge row as it's entries are covered as part
        # of the main diagonal.
        le_rows = te_rows + ny * 3

        # Incidies of LE row repeated nx-1 times
        le_cols = np.tile(np.arange(3 * ny), nx - 1)

        # Incidies of TE row. Done by offsetting le_cols
        te_cols = le_cols + ny * 3 * (nx - 1)

        # Concactenate rows and cols together
        rows = np.concatenate([p_rows, te_rows, le_rows])
        cols = np.concatenate([p_rows, te_cols, le_cols])

        self.declare_partials("mesh", "in_mesh", rows=rows, cols=cols)

    def compute(self, inputs, outputs):
        mesh = inputs["in_mesh"]
        chord_dist = inputs["chord"]

        # Get trailing edge coordinates (ny, 3)
        te = mesh[-1]
        # Get leading edge coordinates (ny, 3)
        le = mesh[0]
        # Linear interpolation to compute the ref_axis coordinates (ny, 3)
        ref_axis = self.ref_axis_pos * te + (1 - self.ref_axis_pos) * le

        # Modify the mesh based on the chord scaling distribution
        # j - spanwise station index (ny)
        # Broadcast chord_dist array over the mesh along spanwise(j) index multiply it by the x and z coordinates
        outputs["mesh"] = np.einsum("ijk,j->ijk", mesh - ref_axis, chord_dist) + ref_axis

    def compute_partials(self, inputs, partials):
        mesh = inputs["in_mesh"]
        chord_dist = inputs["chord"]

        # Get trailing edge coordinates (ny, 3)
        te = mesh[-1]
        # Get leading edge coordinates (ny, 3)
        le = mesh[0]
        # Linear interpolation to compute the ref_axis coordinates (ny, 3)
        ref_axis = self.ref_axis_pos * te + (1 - self.ref_axis_pos) * le

        # Since we are multiplying the mesh at each spanwise station by chord_dist at that station
        # the deritive with respect to chord is just the mesh itself(offset to ref_axis)
        partials["mesh", "chord"] = (mesh - ref_axis).flatten()

        # Compute total number of array entries in mesh array
        nx, ny, _ = mesh.shape
        nn = nx * ny * 3

        # The diagonol part of the d mesh/ d in_mesh jacobian is just the chord_dist broadcast over the
        # leading edge row of ones then tiled nx times to account for the rest of the mesh rows
        d_mesh = np.einsum("i,ij->ij", chord_dist, np.ones((ny, 3))).flatten()
        partials["mesh", "in_mesh"][:nn] = np.tile(d_mesh, nx)

        # Broadcast (1 - chord_dist) onto a single row of the mesh. Result is needed in all
        # ref_axis related sensitivities.
        d_qc = (np.einsum("ij,i->ij", np.ones((ny, 3)), 1.0 - chord_dist)).flatten()

        # Off-diagonal parts of the jacobian
        nnq = (nx - 1) * ny * 3

        # Sensitivies of non-TE parts of mesh to TE contribution to ref_axis
        partials["mesh", "in_mesh"][nn : nn + nnq] = np.tile(self.ref_axis_pos * d_qc, nx - 1)

        # Sensitivies of non-LE parts of mesh to LE contribution to ref_axis
        partials["mesh", "in_mesh"][nn + nnq :] = np.tile((1 - self.ref_axis_pos) * d_qc, nx - 1)

        # ref_axis related sensitivities have contributions on the main diagonol for the LE and TE
        # themselves.
        nnq = ny * 3

        # Sentivities of TE part of mesh to TE contribution to ref_axis
        partials["mesh", "in_mesh"][nn - nnq : nn] += self.ref_axis_pos * d_qc

        # Sentivities of LE part of mesh to LE contribution to ref_axis
        partials["mesh", "in_mesh"][:nnq] += (1 - self.ref_axis_pos) * d_qc


class Sweep(om.ExplicitComponent):
    """
    OpenMDAO component that manipulates the mesh applying shearing sweep. Positive sweeps back.

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the initial aerodynamic surface.
    sweep : float
        Shearing sweep angle in degrees.
    symmetry : boolean
        Flag set to true if surface is reflected about y=0 plane.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the swept aerodynamic surface.
    """

    def initialize(self):
        """
        Declare options.
        """
        self.options.declare("val", desc="Initial value for x shear.")
        self.options.declare("mesh_shape", desc="Tuple containing mesh shape (nx, ny).")
        self.options.declare(
            "symmetry", default=False, desc="Flag set to true if surface is reflected about y=0 plane."
        )

    def setup(self):
        mesh_shape = self.options["mesh_shape"]
        val = self.options["val"]

        self.add_input("sweep", val=val, units="deg")
        self.add_input("in_mesh", shape=mesh_shape, units="m")

        self.add_output("mesh", shape=mesh_shape, units="m")

        # Declare d mesh/ d sweep jacobian

        # compute total number of points in mesh
        nx, ny, _ = mesh_shape
        nn = nx * ny

        # x-coodinates of entire swept mesh are only sensitive to the scalar sweep(col 0)
        rows = 3 * np.arange(nn)
        cols = np.zeros(nn)

        self.declare_partials("mesh", "sweep", rows=rows, cols=cols)

        # Declare d mesh/ d in_mesh jacobian

        # Diagonal part just passes in_mesh to mesh
        # compute total number of entries in mesh array
        nn = nx * ny * 3
        # Entire swept mesh is sensitive to in_mesh
        n_rows = np.arange(nn)

        # Off-diagonal part to account for distance from symmetry plane part of the sweep calculation
        # This part has sensitive to the y-coordinate of the symmetry plane leading edge and the y-coordinates
        # of the leading edge
        if self.options["symmetry"]:
            # Sensitivity to symmetry plane position
            # y-coodinate index of the symmetry plane leading edge
            y_cp = ny * 3 - 2

            # Fill array with y_cp to cover entire mesh except right tip
            sym_cols = np.tile(y_cp, nx * (ny - 1))

            # x-coordinates indicies of entire mesh except right tip(LE + offset for remainder of mesh)
            sym_rows = np.tile(3 * np.arange(ny - 1), nx) + np.repeat(3 * ny * np.arange(nx), ny - 1)

            # Sensitivity to spanwise station position

            # y-coordinates indices of leading edge except right tip repeated for entire mesh
            span_cols = np.tile(3 * np.arange(ny - 1) + 1, nx)
        else:

            # Sensitivity to symmetry plane position

            # y-coodinate of the center line leading edge
            y_cp = 3 * (ny + 1) // 2 - 2

            # index of center line
            n_sym = (ny - 1) // 2

            # This line generates the x-coordiantes for the leading edge for both spans of the wing
            # Start by tiling the x incicides for the first span twice then adding an offset for the second span.
            # Note the first terms of the repeating addition is 0 so the left span stays as initially generated.
            sym_row = np.tile(3 * np.arange(n_sym), 2) + np.repeat([0, 3 * (n_sym + 1)], n_sym)

            # Repeat this nx times to cover all rows and add the offset so that jacobian covers rest of mesh
            sym_rows = np.tile(sym_row, nx) + np.repeat(3 * ny * np.arange(nx), ny - 1)

            # Repeat y_cp n_sym times
            sym_col = np.tile(y_cp, n_sym)

            # Sensitivity to spanwise station position

            # y-coordinate indicies of left span of mesh
            span_col1 = 3 * np.arange(n_sym) + 1

            # y-coordiante indicies of right span of mesh
            span_col2 = 3 * np.arange(n_sym) + 4 + 3 * n_sym

            # neat trick: swap columns on reflected side so we can assign in just two operations
            # This is performance improving feature that takes advantage of the fact that this part of the
            # jacobian will have either + or - tan(theta) in it. We are simply grouping the cols that will have +
            # entries and - entries together. Ignore the variable naming here as we just want to able to use the
            # same declare partials for both symmetry and no symmetry cases.

            # Group + entries and repeat to cover all chordwise points
            sym_cols = np.tile(np.concatenate([sym_col, span_col2]), nx)

            # Group - entires and repeat to cover all chordwise points
            span_cols = np.tile(np.concatenate([span_col1, sym_col]), nx)

        rows = np.concatenate(([n_rows, sym_rows, sym_rows]))
        cols = np.concatenate(([n_rows, sym_cols, span_cols]))

        self.declare_partials("mesh", "in_mesh", rows=rows, cols=cols)

    def compute(self, inputs, outputs):
        symmetry = self.options["symmetry"]
        sweep_angle = inputs["sweep"][0]
        mesh = inputs["in_mesh"]

        # Get the mesh parameters and desired sweep angle
        nx, ny, _ = mesh.shape
        le = mesh[0]
        p180 = np.pi / 180
        tan_theta = np.tan(p180 * sweep_angle)

        # If symmetric, simply vary the x-coord based on the distance from the
        # center of the wing
        if symmetry:
            y0 = le[-1, 1]
            dx = -(le[:, 1] - y0) * tan_theta

        # Else, vary the x-coord on either side of the wing
        else:
            ny2 = (ny - 1) // 2
            y0 = le[ny2, 1]

            dx_right = (le[ny2:, 1] - y0) * tan_theta
            dx_left = -(le[:ny2, 1] - y0) * tan_theta
            dx = np.hstack((dx_left, dx_right))

        # dx added to mesh x coordinates spanwise.
        outputs["mesh"][:] = mesh
        outputs["mesh"][:, :, 0] += dx

    def compute_partials(self, inputs, partials):
        symmetry = self.options["symmetry"]
        sweep_angle = inputs["sweep"][0]
        mesh = inputs["in_mesh"]

        # Get the mesh parameters and desired sweep angle
        nx, ny, _ = mesh.shape
        le = mesh[0]
        p180 = np.pi / 180
        tan_theta = np.tan(p180 * sweep_angle)

        # Derivative of tan(theta) wrt to theta
        dtan_dtheta = p180 / np.cos(p180 * sweep_angle) ** 2

        # Multiply derivative by distance from center of wing
        if symmetry:
            y0 = le[-1, 1]

            dx_dtheta = -(le[:, 1] - y0) * dtan_dtheta
        else:
            # j index of centerline
            ny2 = (ny - 1) // 2
            # y coordinate of centerline
            y0 = le[ny2, 1]

            dx_dtheta_right = (le[ny2:, 1] - y0) * dtan_dtheta
            dx_dtheta_left = -(le[:ny2, 1] - y0) * dtan_dtheta
            dx_dtheta = np.hstack((dx_dtheta_left, dx_dtheta_right))

        partials["mesh", "sweep"] = np.tile(dx_dtheta, nx)

        # Diagonal part of d mesh/ d in_mesh is just 1 to pass the in_mesh through
        nn = nx * ny * 3
        partials["mesh", "in_mesh"][:nn] = 1.0

        # Assign tan and then -tan to off diagonal parts to account for spanwise station sensitivity
        nn2 = nx * (ny - 1)
        partials["mesh", "in_mesh"][nn : nn + nn2] = tan_theta
        partials["mesh", "in_mesh"][nn + nn2 :] = -tan_theta


class ShearX(om.ExplicitComponent):
    """
    OpenMDAO component that manipulates the mesh by shearing the wing in the x direction
    (distributed sweep).

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

    def initialize(self):
        """
        Declare options.
        """
        self.options.declare("val", desc="Initial value for x shear.")
        self.options.declare("mesh_shape", desc="Tuple containing mesh shape (nx, ny).")

    def setup(self):
        mesh_shape = self.options["mesh_shape"]
        val = self.options["val"]

        self.add_input("xshear", val=val, units="m")
        self.add_input("in_mesh", shape=mesh_shape, units="m")

        self.add_output("mesh", shape=mesh_shape, units="m")

        nx, ny, _ = mesh_shape

        nn = nx * ny

        # Derivative of mesh wrt to xshear vector

        # Vector of all mesh array x entries
        rows = 3 * np.arange(nn)

        # Tile vector of all spanwise stations by number of chordwise panels
        cols = np.tile(np.arange(ny), nx)

        # Jacobian entries pass the columns to the rows one to one
        val = np.ones(nn)

        self.declare_partials("mesh", "xshear", rows=rows, cols=cols, val=val)

        # Derivative of mesh wrt in_mesh is just identity
        nn = nx * ny * 3
        rows = np.arange(nn)
        cols = np.arange(nn)
        val = np.ones(nn)

        self.declare_partials("mesh", "in_mesh", rows=rows, cols=cols, val=val)

    def compute(self, inputs, outputs):
        outputs["mesh"][:] = inputs["in_mesh"]

        # Add the xshear distribution to all x coordinates
        outputs["mesh"][:, :, 0] += inputs["xshear"]


class Stretch(om.ExplicitComponent):
    """
    OpenMDAO component that manipulates the mesh by stretching the mesh in spanwise direction to
    reach specified span

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

    def initialize(self):
        """
        Declare options.
        """
        self.options.declare("val", desc="Initial value for span.")
        self.options.declare("mesh_shape", desc="Tuple containing mesh shape (nx, ny).")
        self.options.declare(
            "symmetry", default=False, desc="Flag set to true if surface is reflected about y=0 plane."
        )
        self.options.declare(
            "ref_axis_pos",
            default=0.25,
            desc="Fraction of the chord to use as the reference axis",
        )

    def setup(self):
        mesh_shape = self.options["mesh_shape"]
        val = self.options["val"]
        self.ref_axis_pos = self.options["ref_axis_pos"]

        self.add_input("span", val=val, units="m")
        self.add_input("in_mesh", shape=mesh_shape, units="m")

        self.add_output("mesh", shape=mesh_shape, units="m")

        # Declare derivative of mesh wrt to span vector
        nx, ny, _ = mesh_shape
        nn = nx * ny

        # All y components of every mesh point
        rows = 3 * np.arange(nn) + 1

        # All mesh points sensitive to the scalar(col 0)
        cols = np.zeros(nn)

        self.declare_partials("mesh", "span", rows=rows, cols=cols)

        # Declare derivative of the mesh wrt to the in_mesh

        # First: x and z on diag is identity.
        # Note this is just the x diag. We will get z by offseting by 2 later.
        nn = nx * ny
        xz_diag = 3 * np.arange(nn)

        # Second: y at the corners of the mesh
        # Four columns at le (tip, root) and te (tip, root)
        i_le0 = 1
        i_le1 = ny * 3 - 2
        i_te0 = (nx - 1) * ny * 3 + 1
        i_te1 = nn * 3 - 2

        # Tile all the y indices of the mesh 4 times for each corner
        rows_4c = np.tile(3 * np.arange(nn) + 1, 4)

        # Tile each corner index nn times and concatenate them together
        cols_4c = np.concatenate([np.tile(i_le0, nn), np.tile(i_le1, nn), np.tile(i_te0, nn), np.tile(i_te1, nn)])

        # Third: y indicies for the rest of the mesh
        # Diagonal stripes

        # y incides of the LE other than corners
        base = 3 * np.arange(1, ny - 1) + 1

        # Tile the base vector nx times to cover rest of mesh and add repeating offset so it covers the y indices only
        row_dg = np.tile(base, nx) + np.repeat(ny * 3 * np.arange(nx), ny - 2)

        # Tile rows_dg twice to account for two contributions to the derivative(LE and TE)
        rows_dg = np.tile(row_dg, 2)

        # Tile the base nx times so its size covers the rest of the mesh other than corners
        col_dg = np.tile(base, nx)

        # Concatenate the result with a version offset by entire mesh minus the trailing edge so that TE is covered
        cols_dg = np.concatenate([col_dg, col_dg + 3 * ny * (nx - 1)])

        # Concatenate all contributions together
        # x diag
        # z diag (x diag offset by 2)
        # 4 corners of mesh
        # diagonals for y indices of remaining mesh
        rows = np.concatenate([xz_diag, xz_diag + 2, rows_4c, rows_dg])
        cols = np.concatenate([xz_diag, xz_diag + 2, cols_4c, cols_dg])

        self.declare_partials("mesh", "in_mesh", rows=rows, cols=cols)

    def compute(self, inputs, outputs):
        symmetry = self.options["symmetry"]
        span = inputs["span"][0]
        mesh = inputs["in_mesh"]

        # Set the span along the quarter-chord line
        le = mesh[0]
        te = mesh[-1]
        ref_axis = self.ref_axis_pos * te + (1 - self.ref_axis_pos) * le

        # The user always deals with the full span, so if they input a specific
        # span value and have symmetry enabled, we divide this value by 2.
        if symmetry:
            span /= 2.0

        # Compute the previous span and determine the scalar needed to reach the
        # desired span
        prev_span = ref_axis[-1, 1] - ref_axis[0, 1]
        s = ref_axis[:, 1] / prev_span

        outputs["mesh"][:] = mesh
        outputs["mesh"][:, :, 1] = s * span

    def compute_partials(self, inputs, partials):
        symmetry = self.options["symmetry"]
        span = inputs["span"][0]
        mesh = inputs["in_mesh"]
        nx, ny, _ = mesh.shape

        # Set the span along the reference axis line
        le = mesh[0]
        te = mesh[-1]
        ref_axis = self.ref_axis_pos * te + (1 - self.ref_axis_pos) * le

        # The user always deals with the full span, so if they input a specific
        # span value and have symmetry enabled, we divide this value by 2.
        if symmetry:
            span /= 2.0

        # Compute the previous span and determine the scalar needed to reach the
        # desired span
        prev_span = ref_axis[-1, 1] - ref_axis[0, 1]
        s = ref_axis[:, 1] / prev_span

        # Compute derivative of mesh wrt to span vector

        if symmetry:
            # Tile half the scalar vector s by nx since we only consider half the span
            partials["mesh", "span"] = np.tile(0.5 * s, nx)
        else:
            # Tile the scalar vector s by nx s
            partials["mesh", "span"] = np.tile(s, nx)

        # Compute the derivative of mesh wrt to in_mesh

        # derivative of s wrt to the prev_span
        d_prev_span = -ref_axis[:, 1] / prev_span**2

        # derivative of s wrt to the ref axis (first and last points)
        d_prev_span_qc0 = np.zeros((ny,))
        d_prev_span_qc1 = np.zeros((ny,))

        # First point and last point only sensitive to ref_axis
        d_prev_span_qc0[0] = d_prev_span_qc1[-1] = 1.0 / prev_span

        # Cover the x and z diagonals with 1s
        nn = nx * ny * 2
        partials["mesh", "in_mesh"][:nn] = 1.0

        # LE tip partials. d mesh / d(le tip position)
        nn2 = nx * ny
        partials["mesh", "in_mesh"][nn : nn + nn2] = np.tile(
            -(1 - self.ref_axis_pos) * span * (d_prev_span - d_prev_span_qc0), nx
        )

        # LE root partials d mesh / d(le root position)
        nn3 = nn + nn2 * 2
        partials["mesh", "in_mesh"][nn + nn2 : nn3] = np.tile(
            (1 - self.ref_axis_pos) * span * (d_prev_span + d_prev_span_qc1), nx
        )

        # TE tip partials d mesh / d(te tip position)
        nn4 = nn3 + nn2
        partials["mesh", "in_mesh"][nn3:nn4] = np.tile(-self.ref_axis_pos * span * (d_prev_span - d_prev_span_qc0), nx)

        # TE root partials d mesh / d(te root position)
        nn5 = nn4 + nn2
        partials["mesh", "in_mesh"][nn4:nn5] = np.tile(self.ref_axis_pos * span * (d_prev_span + d_prev_span_qc1), nx)

        # Non corner LE partials d mesh/ d(le except corners)
        nn6 = nn5 + nx * (ny - 2)
        partials["mesh", "in_mesh"][nn5:nn6] = (1 - self.ref_axis_pos) * span / prev_span

        # Non corner TE partials d mesh/ d(te except corners)
        partials["mesh", "in_mesh"][nn6:] = self.ref_axis_pos * span / prev_span


class ShearY(om.ExplicitComponent):
    """
    OpenMDAO component that manipulates the mesh by shearing the wing in the y direction
    (distributed sweep).

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the initial aerodynamic surface.
    yshear[ny] : numpy array
        Distance to translate wing in y direction.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh with the new chord lengths.
    """

    def initialize(self):
        """
        Declare options.
        """
        self.options.declare("val", desc="Initial value for y shear.")
        self.options.declare("mesh_shape", desc="Tuple containing mesh shape (nx, ny).")

    def setup(self):
        mesh_shape = self.options["mesh_shape"]
        val = self.options["val"]

        self.add_input("yshear", val=val, units="m")
        self.add_input("in_mesh", shape=mesh_shape, units="m")

        self.add_output("mesh", shape=mesh_shape, units="m")

        nx, ny, _ = mesh_shape

        nn = nx * ny
        rows = 3.0 * np.arange(nn) + 1
        cols = np.tile(np.arange(ny), nx)
        val = np.ones(nn)

        self.declare_partials("mesh", "yshear", rows=rows, cols=cols, val=val)

        nn = nx * ny * 3
        rows = np.arange(nn)
        cols = np.arange(nn)
        val = np.ones(nn)

        self.declare_partials("mesh", "in_mesh", rows=rows, cols=cols, val=val)

    def compute(self, inputs, outputs):
        outputs["mesh"][:] = inputs["in_mesh"]
        outputs["mesh"][:, :, 1] += inputs["yshear"]


class Dihedral(om.ExplicitComponent):
    """
    OpenMDAO component that manipulates the mesh by applying dihedral angle. Positive angles up.

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the initial aerodynamic surface.
    dihedral : float
        Dihedral angle in degrees.
    symmetry : boolean
        Flag set to true if surface is reflected about y=0 plane.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the aerodynamic surface with dihedral angle.
    """

    def initialize(self):
        """
        Declare options.
        """
        self.options.declare("val", desc="Initial value for dihedral.")
        self.options.declare("mesh_shape", desc="Tuple containing mesh shape (nx, ny).")
        self.options.declare(
            "symmetry", default=False, desc="Flag set to true if surface is reflected about y=0 plane."
        )

    def setup(self):
        mesh_shape = self.options["mesh_shape"]
        val = self.options["val"]

        self.add_input("dihedral", val=val, units="deg")
        self.add_input("in_mesh", shape=mesh_shape, units="m")

        self.add_output("mesh", shape=mesh_shape, units="m")

        nx, ny, _ = mesh_shape
        nn = nx * ny
        rows = 3 * np.arange(nn) + 2
        cols = np.zeros(nn)

        self.declare_partials("mesh", "dihedral", rows=rows, cols=cols)

        nn = nx * ny * 3
        n_rows = np.arange(nn)

        if self.options["symmetry"]:
            y_cp = ny * 3 - 2
            te_cols = np.tile(y_cp, nx * (ny - 1))
            te_rows = np.tile(3 * np.arange(ny - 1) + 2, nx) + np.repeat(3 * ny * np.arange(nx), ny - 1)
            se_cols = np.tile(3 * np.arange(ny - 1) + 1, nx)
        else:
            y_cp = 3 * (ny + 1) // 2 - 2
            n_sym = (ny - 1) // 2

            te_row = np.tile(3 * np.arange(n_sym) + 2, 2) + np.repeat([0, 3 * (n_sym + 1)], n_sym)
            te_rows = np.tile(te_row, nx) + np.repeat(3 * ny * np.arange(nx), ny - 1)

            te_col = np.tile(y_cp, n_sym)
            se_col1 = 3 * np.arange(n_sym) + 1
            se_col2 = 3 * np.arange(n_sym) + 4 + 3 * n_sym

            # neat trick: swap columns on reflected side so we can assign in just two operations
            te_cols = np.tile(np.concatenate([te_col, se_col2]), nx)
            se_cols = np.tile(np.concatenate([se_col1, te_col]), nx)

        rows = np.concatenate(([n_rows, te_rows, te_rows]))
        cols = np.concatenate(([n_rows, te_cols, se_cols]))

        self.declare_partials("mesh", "in_mesh", rows=rows, cols=cols)

    def compute(self, inputs, outputs):
        symmetry = self.options["symmetry"]
        dihedral_angle = inputs["dihedral"][0]
        mesh = inputs["in_mesh"]

        # Get the mesh parameters and desired sweep angle
        _, ny, _ = mesh.shape
        le = mesh[0]
        p180 = np.pi / 180
        tan_theta = np.tan(p180 * dihedral_angle)

        # If symmetric, simply vary the z-coord based on the distance from the
        # center of the wing
        if symmetry:
            y0 = le[-1, 1]
            dz = -(le[:, 1] - y0) * tan_theta

        else:
            ny2 = (ny - 1) // 2
            y0 = le[ny2, 1]
            dz_right = (le[ny2:, 1] - y0) * tan_theta
            dz_left = -(le[:ny2, 1] - y0) * tan_theta
            dz = np.hstack((dz_left, dz_right))

        # dz added spanwise.
        outputs["mesh"][:] = mesh
        outputs["mesh"][:, :, 2] += dz

    def compute_partials(self, inputs, partials):
        symmetry = self.options["symmetry"]
        dihedral_angle = inputs["dihedral"][0]
        mesh = inputs["in_mesh"]

        # Get the mesh parameters and desired sweep angle
        nx, ny, _ = mesh.shape
        le = mesh[0]
        p180 = np.pi / 180
        tan_theta = np.tan(p180 * dihedral_angle)
        dtan_dangle = p180 / np.cos(p180 * dihedral_angle) ** 2

        # If symmetric, simply vary the z-coord based on the distance from the
        # center of the wing
        if symmetry:
            y0 = le[-1, 1]
            dz_dtheta = -(le[:, 1] - y0) * dtan_dangle

        else:
            ny2 = (ny - 1) // 2
            y0 = le[ny2, 1]

            ddz_right = (le[ny2:, 1] - y0) * dtan_dangle
            ddz_left = -(le[:ny2, 1] - y0) * dtan_dangle
            dz_dtheta = np.hstack((ddz_left, ddz_right))

        # dz added spanwise.
        partials["mesh", "dihedral"] = np.tile(dz_dtheta, nx)

        nn = nx * ny * 3
        partials["mesh", "in_mesh"][:nn] = 1.0

        nn2 = nx * (ny - 1)
        partials["mesh", "in_mesh"][nn : nn + nn2] = tan_theta
        partials["mesh", "in_mesh"][nn + nn2 :] = -tan_theta


class ShearZ(om.ExplicitComponent):
    """
    OpenMDAO component that manipulates the mesh by shearing the wing in the z direction
    (distributed sweep).

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

    def initialize(self):
        """
        Declare options.
        """
        self.options.declare("val", desc="Initial value for z shear.")
        self.options.declare("mesh_shape", desc="Tuple containing mesh shape (nx, ny).")

    def setup(self):
        mesh_shape = self.options["mesh_shape"]
        val = self.options["val"]

        self.add_input("zshear", val=val, units="m")
        self.add_input("in_mesh", shape=mesh_shape, units="m")

        self.add_output("mesh", shape=mesh_shape, units="m")

        nx, ny, _ = mesh_shape

        nn = nx * ny
        rows = 3.0 * np.arange(nn) + 2
        cols = np.tile(np.arange(ny), nx)
        val = np.ones(nn)

        self.declare_partials("mesh", "zshear", rows=rows, cols=cols, val=val)

        nn = nx * ny * 3
        rows = np.arange(nn)
        cols = np.arange(nn)
        val = np.ones(nn)

        self.declare_partials("mesh", "in_mesh", rows=rows, cols=cols, val=val)

    def compute(self, inputs, outputs):
        outputs["mesh"][:] = inputs["in_mesh"]
        outputs["mesh"][:, :, 2] += inputs["zshear"]


class Rotate(om.ExplicitComponent):
    """
    OpenMDAO component that manipulates the mesh by compute rotation matrices given mesh and
    rotation angles in degrees.

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
        applied perpendicular to the wing (say, in the case of a winglet).

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the twisted aerodynamic surface.
    """

    def initialize(self):
        """
        Declare options.
        """
        self.options.declare("val", desc="Initial value for dihedral.")
        self.options.declare("mesh_shape", desc="Tuple containing mesh shape (nx, ny).")
        self.options.declare(
            "symmetry", default=False, desc="Flag set to true if surface is reflected about y=0 plane."
        )
        self.options.declare(
            "rotate_x",
            default=True,
            desc="Flag set to True if the user desires the twist variable to "
            "always be applied perpendicular to the wing (say, in the case of "
            "a winglet).",
        )
        self.options.declare(
            "ref_axis_pos",
            default=0.25,
            desc="Fraction of the chord to use as the reference axis",
        )

    def setup(self):
        mesh_shape = self.options["mesh_shape"]
        val = self.options["val"]
        self.ref_axis_pos = self.options["ref_axis_pos"]

        self.add_input("twist", val=val, units="deg")
        self.add_input("in_mesh", shape=mesh_shape, units="m")

        self.add_output("mesh", shape=mesh_shape, units="m")

        nx, ny, _ = mesh_shape
        nn = nx * ny * 3
        rows = np.arange(nn)
        col = np.tile(np.zeros(3), ny) + np.repeat(np.arange(ny), 3)
        cols = np.tile(col, nx)

        self.declare_partials("mesh", "twist", rows=rows, cols=cols)

        row_base = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
        col_base = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])

        # Diagonal
        nn = nx * ny
        dg_row = np.tile(row_base, nn) + np.repeat(3 * np.arange(nn), 9)
        dg_col = np.tile(col_base, nn) + np.repeat(3 * np.arange(nn), 9)

        # Leading and Trailing edge on diagonal terms.
        row_base_y = np.tile(row_base, ny) + np.repeat(3 * np.arange(ny), 9)
        col_base_y = np.tile(col_base, ny) + np.repeat(3 * np.arange(ny), 9)
        nn2 = 3 * ny
        te_dg_row = np.tile(row_base_y, nx - 1) + np.repeat(nn2 * np.arange(nx - 1), 9 * ny)
        le_dg_col = np.tile(col_base_y, nx - 1)
        le_dg_row = te_dg_row + nn2
        te_dg_col = le_dg_col + 3 * ny * (nx - 1)

        # Leading and Trailing edge off diagonal terms.
        if self.options["symmetry"]:
            row_base_y = np.tile(row_base, ny - 1) + np.repeat(3 * np.arange(ny - 1), 9)
            col_base_y = np.tile(col_base + 3, ny - 1) + np.repeat(3 * np.arange(ny - 1), 9)

            nn2 = 3 * ny
            te_od_row = np.tile(row_base_y, nx) + np.repeat(nn2 * np.arange(nx), 9 * (ny - 1))
            le_od_col = np.tile(col_base_y, nx)
            te_od_col = le_od_col + 3 * ny * (nx - 1)

            rows = np.concatenate([dg_row, le_dg_row, te_dg_row, te_od_row, te_od_row])
            cols = np.concatenate([dg_col, le_dg_col, te_dg_col, le_od_col, te_od_col])

        else:
            n_sym = (ny - 1) // 2

            row_base_y1 = np.tile(row_base, n_sym) + np.repeat(3 * np.arange(n_sym), 9)
            col_base_y1 = np.tile(col_base + 3, n_sym) + np.repeat(3 * np.arange(n_sym), 9)

            row_base_y2 = row_base_y1 + 3 * n_sym + 3
            col_base_y2 = col_base_y1 + 3 * n_sym - 3

            nn2 = 3 * ny

            te_od_row1 = np.tile(row_base_y1, nx) + np.repeat(nn2 * np.arange(nx), 9 * n_sym)
            le_od_col1 = np.tile(col_base_y1, nx)
            te_od_col1 = le_od_col1 + 3 * ny * (nx - 1)
            te_od_row2 = np.tile(row_base_y2, nx) + np.repeat(nn2 * np.arange(nx), 9 * n_sym)
            le_od_col2 = np.tile(col_base_y2, nx)
            te_od_col2 = le_od_col2 + 3 * ny * (nx - 1)

            rows = np.concatenate([dg_row, le_dg_row, te_dg_row, te_od_row1, te_od_row2, te_od_row1, te_od_row2])
            cols = np.concatenate([dg_col, le_dg_col, te_dg_col, le_od_col1, le_od_col2, te_od_col1, te_od_col2])

        self.declare_partials("mesh", "in_mesh", rows=rows, cols=cols)

    def compute(self, inputs, outputs):
        symmetry = self.options["symmetry"]
        rotate_x = self.options["rotate_x"]
        theta_y = inputs["twist"]
        mesh = inputs["in_mesh"]

        te = mesh[-1]
        le = mesh[0]
        ref_axis = self.ref_axis_pos * te + (1 - self.ref_axis_pos) * le

        _, ny, _ = mesh.shape

        if rotate_x:
            # Compute spanwise z displacements along quarter chord
            if symmetry:
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

                # Concatenate thetas
                rad_theta_x = np.concatenate((theta_x_left, np.zeros(1), theta_x_right))

        else:
            rad_theta_x = 0.0

        rad_theta_y = theta_y * np.pi / 180.0

        mats = np.zeros((ny, 3, 3), dtype=type(rad_theta_y[0]))

        cos_rtx = np.cos(rad_theta_x)
        cos_rty = np.cos(rad_theta_y)
        sin_rtx = np.sin(rad_theta_x)
        sin_rty = np.sin(rad_theta_y)

        mats[:, 0, 0] = cos_rty
        mats[:, 0, 2] = sin_rty
        mats[:, 1, 0] = sin_rtx * sin_rty
        mats[:, 1, 1] = cos_rtx
        mats[:, 1, 2] = -sin_rtx * cos_rty
        mats[:, 2, 0] = -cos_rtx * sin_rty
        mats[:, 2, 1] = sin_rtx
        mats[:, 2, 2] = cos_rtx * cos_rty

        outputs["mesh"] = np.einsum("ikj, mij -> mik", mats, mesh - ref_axis) + ref_axis

    def compute_partials(self, inputs, partials):
        symmetry = self.options["symmetry"]
        rotate_x = self.options["rotate_x"]
        theta_y = inputs["twist"]
        mesh = inputs["in_mesh"]

        te = mesh[-1]
        le = mesh[0]
        ref_axis = self.ref_axis_pos * te + (1 - self.ref_axis_pos) * le

        nx, ny, _ = mesh.shape

        if rotate_x:
            # Compute spanwise z displacements along quarter chord
            if symmetry:
                dz_qc = ref_axis[:-1, 2] - ref_axis[1:, 2]
                dy_qc = ref_axis[:-1, 1] - ref_axis[1:, 1]
                theta_x = np.arctan(dz_qc / dy_qc)

                # Prepend with 0 so that root is not rotated
                rad_theta_x = np.append(theta_x, 0.0)

                fact = 1.0 / (1.0 + (dz_qc / dy_qc) ** 2)

                dthx_dq = np.zeros((ny, 3))
                dthx_dq[:-1, 1] = -dz_qc * fact / dy_qc**2
                dthx_dq[:-1, 2] = fact / dy_qc

            else:
                root_index = int((ny - 1) / 2)
                dz_qc_left = ref_axis[:root_index, 2] - ref_axis[1 : root_index + 1, 2]
                dy_qc_left = ref_axis[:root_index, 1] - ref_axis[1 : root_index + 1, 1]
                theta_x_left = np.arctan(dz_qc_left / dy_qc_left)
                dz_qc_right = ref_axis[root_index + 1 :, 2] - ref_axis[root_index:-1, 2]
                dy_qc_right = ref_axis[root_index + 1 :, 1] - ref_axis[root_index:-1, 1]
                theta_x_right = np.arctan(dz_qc_right / dy_qc_right)

                # Concatenate thetas
                rad_theta_x = np.concatenate((theta_x_left, np.zeros(1), theta_x_right))

                fact_left = 1.0 / (1.0 + (dz_qc_left / dy_qc_left) ** 2)
                fact_right = 1.0 / (1.0 + (dz_qc_right / dy_qc_right) ** 2)

                dthx_dq = np.zeros((ny, 3))
                dthx_dq[:root_index, 1] = -dz_qc_left * fact_left / dy_qc_left**2
                dthx_dq[root_index + 1 :, 1] = -dz_qc_right * fact_right / dy_qc_right**2
                dthx_dq[:root_index, 2] = fact_left / dy_qc_left
                dthx_dq[root_index + 1 :, 2] = fact_right / dy_qc_right

        else:
            rad_theta_x = 0.0

        deg2rad = np.pi / 180.0
        rad_theta_y = theta_y * deg2rad

        mats = np.zeros((ny, 3, 3), dtype=type(rad_theta_y[0]))

        cos_rtx = np.cos(rad_theta_x)
        cos_rty = np.cos(rad_theta_y)
        sin_rtx = np.sin(rad_theta_x)
        sin_rty = np.sin(rad_theta_y)

        mats[:, 0, 0] = cos_rty
        mats[:, 0, 2] = sin_rty
        mats[:, 1, 0] = sin_rtx * sin_rty
        mats[:, 1, 1] = cos_rtx
        mats[:, 1, 2] = -sin_rtx * cos_rty
        mats[:, 2, 0] = -cos_rtx * sin_rty
        mats[:, 2, 1] = sin_rtx
        mats[:, 2, 2] = cos_rtx * cos_rty

        dmats_dthy = np.zeros((ny, 3, 3))
        dmats_dthy[:, 0, 0] = -sin_rty * deg2rad
        dmats_dthy[:, 0, 2] = cos_rty * deg2rad
        dmats_dthy[:, 1, 0] = sin_rtx * cos_rty * deg2rad
        dmats_dthy[:, 1, 2] = sin_rtx * sin_rty * deg2rad
        dmats_dthy[:, 2, 0] = -cos_rtx * cos_rty * deg2rad
        dmats_dthy[:, 2, 2] = -cos_rtx * sin_rty * deg2rad

        d_dthetay = np.einsum("ikj, mij -> mik", dmats_dthy, mesh - ref_axis)
        partials["mesh", "twist"] = d_dthetay.flatten()

        nn = nx * ny * 9
        partials["mesh", "in_mesh"][:nn] = np.tile(mats.flatten(), nx)

        # Quarter chord direct contribution.
        eye = np.tile(np.eye(3).flatten(), ny).reshape(ny, 3, 3)
        d_qch = (eye - mats).flatten()

        nqc = ny * 9
        partials["mesh", "in_mesh"][:nqc] += (1 - self.ref_axis_pos) * d_qch
        partials["mesh", "in_mesh"][nn - nqc : nn] += self.ref_axis_pos * d_qch

        if rotate_x:
            dmats_dthx = np.zeros((ny, 3, 3))
            dmats_dthx[:, 1, 0] = cos_rtx * sin_rty
            dmats_dthx[:, 1, 1] = -sin_rtx
            dmats_dthx[:, 1, 2] = -cos_rtx * cos_rty
            dmats_dthx[:, 2, 0] = sin_rtx * sin_rty
            dmats_dthx[:, 2, 1] = cos_rtx
            dmats_dthx[:, 2, 2] = -sin_rtx * cos_rty

            d_dthetax = np.einsum("ikj, mij -> mik", dmats_dthx, mesh - ref_axis)
            d_dq = np.einsum("ijk, jm -> ijkm", d_dthetax, dthx_dq)

            d_dq_flat = d_dq.flatten()

            del_n = nn - 9 * ny
            nn2 = nn + del_n
            nn3 = nn2 + del_n
            partials["mesh", "in_mesh"][nn:nn2] = (1 - self.ref_axis_pos) * d_dq_flat[-del_n:]
            partials["mesh", "in_mesh"][nn2:nn3] = self.ref_axis_pos * d_dq_flat[:del_n]

            # Contribution back to main diagonal.
            del_n = 9 * ny
            partials["mesh", "in_mesh"][:nqc] += (1 - self.ref_axis_pos) * d_dq_flat[:del_n]
            partials["mesh", "in_mesh"][nn - nqc : nn] += self.ref_axis_pos * d_dq_flat[-del_n:]

            # Quarter chord direct contribution.
            d_qch_od = np.tile(d_qch.flatten(), nx - 1)
            partials["mesh", "in_mesh"][nn:nn2] += (1 - self.ref_axis_pos) * d_qch_od
            partials["mesh", "in_mesh"][nn2:nn3] += self.ref_axis_pos * d_qch_od

            # off-off diagonal pieces
            if symmetry:
                d_dq_flat = d_dq[:, :-1, :, :].flatten()

                del_n = nn - 9 * nx
                nn4 = nn3 + del_n
                partials["mesh", "in_mesh"][nn3:nn4] = -(1 - self.ref_axis_pos) * d_dq_flat
                nn5 = nn4 + del_n
                partials["mesh", "in_mesh"][nn4:nn5] = -self.ref_axis_pos * d_dq_flat

            else:
                d_dq_flat1 = d_dq[:, :root_index, :, :].flatten()
                d_dq_flat2 = d_dq[:, root_index + 1 :, :, :].flatten()

                del_n = nx * root_index * 9
                nn4 = nn3 + del_n
                partials["mesh", "in_mesh"][nn3:nn4] = -(1 - self.ref_axis_pos) * d_dq_flat1
                nn5 = nn4 + del_n
                partials["mesh", "in_mesh"][nn4:nn5] = -(1 - self.ref_axis_pos) * d_dq_flat2
                nn6 = nn5 + del_n
                partials["mesh", "in_mesh"][nn5:nn6] = -self.ref_axis_pos * d_dq_flat1
                nn7 = nn6 + del_n
                partials["mesh", "in_mesh"][nn6:nn7] = -self.ref_axis_pos * d_dq_flat2
