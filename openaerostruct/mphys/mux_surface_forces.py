import openmdao.api as om
from mphys.core import MPhysVariables

from openaerostruct.mphys.utils import get_number_of_nodes, get_src_indices


class MuxSurfaceForces(om.ExplicitComponent):
    """
    Demux surface coordinates from flattened array.
    Mphys expects forces to be passed as a single flattened array,
    but OAS outputs them as a series of 3D arrays (one for each surface).
    This component is responsible handling the conversion between the two.

    Parameters
    ----------
    mesh_point_forces[nx, ny, 3] : numpy array
        The aeordynamic forces evaluated at the mesh nodes for each lifting surface.
        There is one of these per surface.

    Returns
    -------
    f_aero[system_size*3] : numpy array
        Flattened array of aero nodal forces for all lifting surfaces.
    """

    def initialize(self):
        self.options.declare("surfaces", default=None, desc="oas surface dicts", recordable=False)

    def setup(self):
        self.surfaces = self.options["surfaces"]

        self.nnodes = get_number_of_nodes(self.surfaces)
        self.src_indices = get_src_indices(self.surfaces)

        # OpenMDAO part of setup
        for surface in self.surfaces:
            surf_name = surface["name"]
            mesh = surface["mesh"]
            self.add_input(
                f"{surf_name}_mesh_point_forces",
                distributed=False,
                shape=mesh.shape,
                units="N",
                desc="Array defining the aero forces " "on mesh nodes of the lifting surface.",
                tags=["mphys_coupling"],
            )

        self.add_output(
            MPhysVariables.Aerodynamics.Surface.LOADS,
            distributed=False,
            shape=self.nnodes * 3,
            val=0.0,
            units="N",
            desc="flattened aero forces for all oas surfaces",
            tags=["mphys_coupling"],
        )

    def compute(self, inputs, outputs):
        for surface in self.surfaces:
            surf_name = surface["name"]
            outputs[MPhysVariables.Aerodynamics.Surface.LOADS][self.src_indices[surf_name]] = inputs[
                surf_name + "_mesh_point_forces"
            ]

    def compute_jacvec_product(self, inputs, d_inputs, d_outputs, mode):
        if mode == "fwd":
            for surface in self.surfaces:
                surf_name = surface["name"]
                if (
                    MPhysVariables.Aerodynamics.Surface.LOADS in d_outputs
                    and surf_name + "_mesh_point_forces" in d_inputs
                ):
                    d_outputs[MPhysVariables.Aerodynamics.Surface.LOADS][self.src_indices[surf_name]] += d_inputs[
                        surf_name + "_mesh_point_forces"
                    ]

        if mode == "rev":
            for surface in self.surfaces:
                surf_name = surface["name"]
                if (
                    MPhysVariables.Aerodynamics.Surface.LOADS in d_outputs
                    and surf_name + "_mesh_point_forces" in d_inputs
                ):
                    d_inputs[surf_name + "_mesh_point_forces"] += d_outputs[MPhysVariables.Aerodynamics.Surface.LOADS][
                        self.src_indices[surf_name]
                    ]
