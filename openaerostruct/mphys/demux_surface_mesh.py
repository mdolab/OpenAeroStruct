import openmdao.api as om
from mphys.core import MPhysVariables

from openaerostruct.mphys.utils import get_number_of_nodes, get_src_indices


class DemuxSurfaceMesh(om.ExplicitComponent):
    """
    Demux surface coordinates from single flattened array.
    Mphys always passes coordinate as a single flattened array,
    but OAS expects them in a series of 3D arrays (one for each surface).
    This component is responsible handling the conversion between the two.

    Parameters
    ----------
    MPhysVariables.Aerodynamics.Surface.COORDINATES[system_size*3] : numpy array
        Flattened aero mesh coordinates for all lifting surfaces.

    Returns
    -------
    def_mesh[nx, ny, 3] : numpy array
        Array defining the nodal coordinates of each lifting surface.
    """

    def initialize(self):
        self.options.declare("surfaces", default=None, desc="oas surface dicts", recordable=False)

    def setup(self):
        self.surfaces = self.options["surfaces"]

        self.nnodes = get_number_of_nodes(self.surfaces)
        self.src_indices = get_src_indices(self.surfaces)

        # OpenMDAO part of setup
        self.add_input(
            MPhysVariables.Aerodynamics.Surface.COORDINATES,
            distributed=False,
            shape=self.nnodes * 3,
            units="m",
            desc="flattened aero mesh coordinates for all oas surfaces",
            tags=["mphys_coupling"],
        )
        for surface in self.surfaces:
            surf_name = surface["name"]
            mesh = surface["mesh"]
            self.add_output(
                f"{surf_name}_def_mesh",
                distributed=False,
                shape=mesh.shape,
                units="m",
                desc="Array defining the nodal coordinates of the lifting surface.",
                tags=["mphys_coupling"],
            )

    def compute(self, inputs, outputs):
        for surface in self.surfaces:
            surf_name = surface["name"]
            outputs[surf_name + "_def_mesh"] = inputs[MPhysVariables.Aerodynamics.Surface.COORDINATES][
                self.src_indices[surf_name]
            ]

    def compute_jacvec_product(self, inputs, d_inputs, d_outputs, mode):
        if mode == "fwd":
            for surface in self.surfaces:
                surf_name = surface["name"]
                if MPhysVariables.Aerodynamics.Surface.COORDINATES in d_inputs and surf_name + "_def_mesh" in d_outputs:
                    d_outputs[surf_name + "_def_mesh"] += d_inputs[MPhysVariables.Aerodynamics.Surface.COORDINATES][
                        self.src_indices[surf_name]
                    ]
        if mode == "rev":
            for surface in self.surfaces:
                surf_name = surface["name"]
                if MPhysVariables.Aerodynamics.Surface.COORDINATES in d_inputs and surf_name + "_def_mesh" in d_outputs:
                    d_inputs[MPhysVariables.Aerodynamics.Surface.COORDINATES][self.src_indices[surf_name]] += d_outputs[
                        surf_name + "_def_mesh"
                    ]
