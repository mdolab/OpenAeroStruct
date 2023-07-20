import unittest
import numpy as np

from openaerostruct.mphys import AeroBuilder
from openaerostruct.utils.testing import run_test, get_default_surfaces

# check if mphys/mpi4py is available
try:
    from mpi4py import MPI
    import mphys

    mphys_mpi_flag = True
except ImportError:
    mphys_mpi_flag = False


@unittest.skipUnless(mphys_mpi_flag, "mphys/mpi4py is required.")
class Test(unittest.TestCase):
    def setUp(self):
        surfaces = get_default_surfaces()
        comm = MPI.COMM_WORLD
        # Create mphys builder for aero solver
        self.aero_builder = AeroBuilder(surfaces)
        self.aero_builder.initialize(comm)

    def test_tagged_indices(self):
        with self.subTest(case="wing"):
            wing_inds = self.aero_builder.get_tagged_indices(["wing"])
            np.testing.assert_equal(wing_inds, np.arange(0, 8))

        with self.subTest(case="tail"):
            tail_inds = self.aero_builder.get_tagged_indices(["tail"])
            np.testing.assert_equal(tail_inds, np.arange(8, 23))

        with self.subTest(case="wing+tail"):
            wt_inds = self.aero_builder.get_tagged_indices(["wing", "tail"])
            np.testing.assert_equal(wt_inds, np.arange(0, 23))

        with self.subTest(case="all"):
            wt_inds = self.aero_builder.get_tagged_indices(-1)
            np.testing.assert_equal(wt_inds, np.arange(0, 23))


if __name__ == "__main__":
    unittest.main()
