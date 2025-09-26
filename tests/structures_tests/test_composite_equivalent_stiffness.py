"""Test the utility function to compute the composite's equivalent stiffness."""

import unittest

from openmdao.utils.assert_utils import assert_near_equal
from openaerostruct.structures.utils import compute_composite_stiffness


class TestCompositeEquivalentStiffness(unittest.TestCase):
    """Test the effective stiffness computation for an isotropic material."""

    def test(self):
        # isotropic material stiffness values
        E = 70.0e9
        nu = 0.33
        G = E / (2 * (1 + nu))

        # setup OAS surface dictionary for the composite laminate
        surf_dict = {
            "ply_fractions": [0.25, 0.25, 0.25, 0.25],
            "ply_angles": [0, 45, -45, 90],
            "E1": E,
            "E2": E,
            "nu12": nu,
            "G12": G,
        }

        # compute the effective stiffness
        compute_composite_stiffness(surf_dict)
        E_eff = surf_dict["E"]
        G_eff = surf_dict["G"]

        # make sure the results are close to the expected values
        assert_near_equal(E_eff, E, 1e-6)
        assert_near_equal(G_eff, G, 1e-6)

        # another test with a different setup
        surf_dict = {
            "ply_fractions": [0.1, 0.2, 0.3, 0.25, 0.15],
            "ply_angles": [12, -20, 27, 99, -5],
            "E1": E,
            "E2": E,
            "nu12": nu,
            "G12": G,
        }
        compute_composite_stiffness(surf_dict)
        E_eff = surf_dict["E"]
        G_eff = surf_dict["G"]

        # make sure the results are close to the expected values
        assert_near_equal(E_eff, E, 1e-6)
        assert_near_equal(G_eff, G, 1e-6)


if __name__ == "__main__":
    unittest.main()
