import unittest
import importlib
import numpy as np

from openaerostruct.mphys import AeroBuilder
from openaerostruct.utils.testing import get_default_surfaces

# check if mphys/mpi4py is available
try:
    from mpi4py import MPI

    mpi_flag = True
except ImportError:
    mpi_flag = False

mphys_flag = importlib.util.find_spec("mphys") is not None


@unittest.skipUnless(mphys_flag and mpi_flag, "mphys/mpi4py is required.")
class Test(unittest.TestCase):
    def setUp(self):
        self.surfaces = get_default_surfaces()
        comm = MPI.COMM_WORLD
        # Create mphys builder for aero solver
        self.aero_builder = AeroBuilder(self.surfaces)
        self.aero_builder.initialize(comm)

    def test_tagged_indices(self):
        wing_nnodes = self.surfaces[0]["mesh"].size // 3
        tail_nnodes = self.surfaces[1]["mesh"].size // 3
        with self.subTest(case="wing"):
            wing_inds = self.aero_builder.get_tagged_indices(["wing"])
            np.testing.assert_equal(wing_inds, np.arange(0, wing_nnodes))

        with self.subTest(case="tail"):
            tail_inds = self.aero_builder.get_tagged_indices(["tail"])
            np.testing.assert_equal(tail_inds, np.arange(wing_nnodes, wing_nnodes + tail_nnodes))

        with self.subTest(case="wing+tail"):
            wt_inds = self.aero_builder.get_tagged_indices(["wing", "tail"])
            np.testing.assert_equal(wt_inds, np.arange(0, wing_nnodes + tail_nnodes))

        with self.subTest(case="all"):
            wt_inds = self.aero_builder.get_tagged_indices(-1)
            np.testing.assert_equal(wt_inds, np.arange(0, wing_nnodes + tail_nnodes))


if __name__ == "__main__":
    unittest.main()
