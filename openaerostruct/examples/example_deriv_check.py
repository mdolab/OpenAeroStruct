import unittest

import openmdao.api as om

from numpy.testing import assert_almost_equal

from openaerostruct.structures.section_properties_tube import SectionPropertiesTube
from openaerostruct.utils.testing import get_default_surfaces, run_test, view_mat


class Test(unittest.TestCase):
    """
    This class contains two tests on a single component.
    Each test verifies the derivatives found within the component by using
    OpenMDAO's derivative approximation functions to compare against
    finite-difference or complex-step.
    """

    def test_quick(self):
        """ Short pre-setup test to compare component derivs. """
        surfaces = get_default_surfaces()

        comp = SectionPropertiesTube(surface=surfaces[0])

        run_test(self, comp, complex_flag=True)

    def test_detailed(self):
        """ This is a longer version of the previous method, with plotting. """

        # Load in default lifting surfaces to setup the comparison
        surfaces = get_default_surfaces()

        # Instantiate an OpenMDAO problem and add the component we want to test
        # as asubsystem, giving that component a default lifting surface
        prob = om.Problem()
        prob.model.add_subsystem('tube', SectionPropertiesTube(surface=surfaces[0]))

        # Set up the problem and ensure it uses complex arrays so we can check
        # the derivatives using complex step
        prob.setup(force_alloc_complex=True)

        # Actually run the model, which is just a component in this case, then
        # check the derivatives and store the results in the `check` dict
        prob.run_model()
        check = prob.check_partials(compact_print=False)

        # Loop through this `check` dictionary and visualize the approximated
        # and computed derivatives
        for key, subjac in check[list(check.keys())[0]].items():
            print()
            print(key)
            view_mat(subjac['J_fd'],subjac['J_fwd'],key)

        # Loop through the `check` dictionary and perform assert that the
        # approximated deriv must be very close to the computed deriv
        for key, subjac in check[list(check.keys())[0]].items():
            if subjac['magnitude'].fd > 1e-6:
                assert_almost_equal(
                    subjac['rel error'].forward, 0., err_msg='deriv of %s wrt %s' % key)
                assert_almost_equal(
                    subjac['rel error'].reverse, 0., err_msg='deriv of %s wrt %s' % key)

# Run the tests included in this script
if __name__ == '__main__':
    unittest.main()
