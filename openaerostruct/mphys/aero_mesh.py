import numpy as np
import openmdao.api as om
from mphys.core import MPhysVariables

from openaerostruct.mphys.utils import get_number_of_nodes, get_src_indices


class AeroMesh(om.IndepVarComp):
    """
    Component to read the initial mesh coordinates with OAS.
    Only the root will be responsible for this information.
    The mesh will be broadcasted to all other processors in a following step.
    """

    def initialize(self):
        self.options.declare("surfaces", default=None, desc="oas surface dicts", recordable=False)

    def setup(self):
        if self.comm.rank == 0:
            self.surfaces = self.options["surfaces"]
            nnodes = get_number_of_nodes(self.surfaces)
            src_indices = get_src_indices(self.surfaces)
            xpts = np.zeros(nnodes * 3)
            for surface in self.surfaces:
                surf_name = surface["name"]
                xpts[src_indices[surf_name]] = surface["mesh"]
        else:
            xpts = np.zeros(0)
        self.add_output(
            MPhysVariables.Aerodynamics.Surface.Mesh.COORDINATES,
            distributed=True,
            val=xpts,
            shape=xpts.size,
            units="m",
            desc="aero node coordinates",
            tags=["mphys_coordinates"],
        )
