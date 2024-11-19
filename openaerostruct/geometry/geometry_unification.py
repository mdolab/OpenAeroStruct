import numpy as np
import openmdao.api as om


def compute_uni_mesh_dims(sections):
    """
    Function computes the dimensions of the unified mesh as well as an index array

    Parameters
    ----------
    sections : list
        List of section OpenAeroStruct surface dictionaries

    Returns
    -------
    uni_mesh_indices : numpy array
        Array of indicies matching the size of the mesh array
    uni_nx : int
        Number of chordwise points in the unified mesh
    uni_ny : int
        Number of spanwise points in the unified mesh
    """
    uni_nx = sections[0]["mesh"].shape[0]
    uni_ny = 0
    for i_sec in range(len(sections)):
        if i_sec == len(sections) - 1:
            uni_ny += sections[i_sec]["mesh"].shape[1]
        else:
            uni_ny += sections[i_sec]["mesh"].shape[1] - 1
    uni_mesh_indices = np.arange(uni_nx * uni_ny * 3).reshape((uni_nx, uni_ny, 3))

    return uni_mesh_indices, uni_nx, uni_ny


def compute_uni_mesh_index_blocks(sections, uni_mesh_indices):
    """
    Function that computes the index block that corresponds to each individual wing section

    Parameters
    ----------
    sections : list
        List of section OpenAeroStruct surface dictionaries
    uni_mesh_indices : numpy array
        Array of indicies matching the size of the mesh array

    Returns
    -------
    blocks : list
        List of numpy arrays containing the incidies corresponding to each section. Basically breaks the uni_mesh_indices
        array up into pieces that correspond to each section.
    """
    blocks = []

    # cursor to track the y position of each section along the unified mesh
    y_curr = 0

    for i_sec in range(len(sections)):
        mesh = sections[i_sec]["mesh"]
        ny = mesh.shape[1]

        if i_sec == len(sections) - 1:
            block = uni_mesh_indices[:, y_curr:, :]
            y_curr += ny
        else:
            block = uni_mesh_indices[:, y_curr : y_curr + (ny - 1), :]
            y_curr += ny - 1

        blocks.append(block.flatten())

    return blocks


def unify_mesh(sections):
    """
    Function that produces a unified mesh from all the individual wing section meshes.

    Parameters
    ----------
    sections : list
        List of section OpenAeroStruct surface dictionaries

    Returns
    -------
    uni_mesh : numpy array
        Unfied surface mesh in OAS format
    """
    for i_sec in np.arange(0, len(sections) - 1):
        mesh = sections[i_sec]["mesh"]

        import copy

        if i_sec == 0:
            uni_mesh = copy.deepcopy(mesh[:, :-1, :])
        else:
            uni_mesh = np.concatenate([uni_mesh, mesh[:, :-1, :]], axis=1)
    # Stitch the results into a singular mesh
    mesh = sections[len(sections) - 1]["mesh"]
    if len(sections) == 1:
        import copy

        uni_mesh = copy.deepcopy(mesh)
    else:
        uni_mesh = np.concatenate([uni_mesh, mesh], axis=1)

    return uni_mesh


class GeomMultiUnification(om.ExplicitComponent):
    """
    OpenMDAO component that combines the meshes associated with each individual section
    of a multi-section wing into a unified mesh. Section 0 as defined by the user(typically
    the most inboard section) will be added to the unified mesh unchanged. The remaining sections
    are then sequencially added outboard to section 0 along the y-direction. The mesh unification
    component operates under the assumption that it is unifying section that are C0 continous along
    the span. The most inboard column of points of section i+1 will be coincident with the
    most outbord column of points of section i. As result, the inboard column of points for each section
    after section 0 are removed prior to the section being added to the unified mesh. The component also
    produces the appropriate sparse Jacobian matrix for the column removal so that the partials can be computed
    appropriately.

    Parameters
    ----------
    sections : list
        A list of section surface dictionaries in OAS format.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the unified surface
    """

    def initialize(self):
        """
        Declare options.
        """
        self.options.declare("sections", types=list, desc="A list of section surface dictionaries to be unified.")
        self.options.declare("surface_name", types=str, desc="The name of the multi-section surface")

    def setup(self):
        sections = self.options["sections"]
        name = self.options["surface_name"]

        # Get the unified mesh size, index array, and block of indicies for each section
        [uni_mesh_indices, uni_nx, uni_ny] = compute_uni_mesh_dims(sections)
        uni_mesh_blocks = compute_uni_mesh_index_blocks(sections, uni_mesh_indices)
        uni_mesh_name = "{}_uni_mesh".format(name)

        # Loop through each section to build the sparsity pattern for the Jacobian. This Jacobian is fixed based on the mesh size so can be declared here
        for i_sec, section in enumerate(sections):
            mesh = section["mesh"]
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = section["name"]

            mesh_name = "{}_def_mesh".format(name)

            self.add_input(mesh_name, shape=(nx, ny, 3), units="m", tags=["mphys_coupling"])

            # Generate index array
            mesh_indices = np.arange(nx * ny * 3).reshape((nx, ny, 3))

            if i_sec == len(sections) - 1:
                cols = mesh_indices.flatten()
            else:
                # If not section N then the most inboard column is disregarded
                cols = mesh_indices[:, :-1, :].flatten()

            # Get data from section block in unified mesh
            rows = uni_mesh_blocks[i_sec]

            # Fill non zero Jacobian entries with ones
            data = np.ones_like(rows)

            self.declare_partials(uni_mesh_name, mesh_name, val=data, rows=rows, cols=cols)

        self.add_output(uni_mesh_name, shape=(uni_nx, uni_ny, 3), units="m")

        # Unify the t/c output of each section if that has been specified
        if "t_over_c_cp" in sections[0].keys():
            uni_tc_name = "{}_uni_t_over_c".format(self.options["surface_name"])

            n_acc = 0
            for section in sections:
                name = section["name"]
                t_over_c_name = "{}_t_over_c".format(name)
                n = int(ny - 1)
                self.add_input(t_over_c_name, shape=(n), tags=["mphys_coupling"])

                self.declare_partials(
                    uni_tc_name, t_over_c_name, rows=np.arange(0, n) + n_acc, cols=np.arange(0, n), val=np.ones(n)
                )

                n_acc += n

            self.add_output(uni_tc_name, shape=(uni_ny - 1))

    def compute(self, inputs, outputs):
        sections = self.options["sections"]
        surface_name = self.options["surface_name"]
        uni_mesh_name = "{}_uni_mesh".format(surface_name)

        # Loop through all sections to unify the mesh
        for i_sec in np.arange(0, len(sections) - 1):
            name = sections[i_sec]["name"]
            mesh_name = "{}_def_mesh".format(name)

            if i_sec == 0:
                uni_mesh = inputs[mesh_name][:, :-1, :]
            else:
                uni_mesh = np.concatenate([uni_mesh, inputs[mesh_name][:, :-1, :]], axis=1)

        mesh_name = "{}_def_mesh".format(sections[len(sections) - 1]["name"])
        if len(sections) == 1:
            uni_mesh = inputs[mesh_name]
        else:
            uni_mesh = np.concatenate([uni_mesh, inputs[mesh_name]], axis=1)

        if "t_over_c_cp" in sections[0].keys():
            uni_tc_name = "{}_uni_t_over_c".format(self.options["surface_name"])
            for i_sec, section in enumerate(sections):
                name = section["name"]
                t_over_c_name = "{}_t_over_c".format(name)
                t_over_c = inputs[t_over_c_name]

                if i_sec == 0:
                    uni_t_over_c = inputs[t_over_c_name]
                else:
                    uni_t_over_c = np.concatenate([uni_t_over_c, t_over_c])
            outputs[uni_tc_name] = uni_t_over_c

        outputs[uni_mesh_name] = uni_mesh
