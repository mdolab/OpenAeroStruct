import openmdao.api as om
from openaerostruct.transfer.displacement_transfer import DisplacementTransfer
from openaerostruct.transfer.compute_transformation_matrix import ComputeTransformationMatrix


class DisplacementTransferGroup(om.Group):
    """
    These components take the displacements and rotations obtained by
    solving the FEM problem and applies them to the aerodynamic mesh
    to produce a deformed aerodynamic mesh.
    """

    def initialize(self):
        self.options.declare("surface", types=dict)

    def setup(self):
        surface = self.options["surface"]

        self.add_subsystem(
            "compute_transformation_matrix", ComputeTransformationMatrix(surface=surface), promotes=["*"]
        )

        self.add_subsystem("displacement_transfer", DisplacementTransfer(surface=surface), promotes=["*"])
