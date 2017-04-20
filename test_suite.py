from __future__ import division, print_function
import sys
from time import time
import unittest
import numpy as np

# Append the parent directory to the system path so we can call those Python
# files. If you have OpenAeroStruct in your PYTHONPATH, this is not necessary.
from os import sys, path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))

from OpenAeroStruct import OASProblem

try:
    import OAS_API
    fortran_flag = True
except:
    fortran_flag = False

class TestAero(unittest.TestCase):

    currentResult = []

    def run(self, result=None):
        self.currentResult.append(result) # remember result for use in tearDown
        unittest.TestCase.run(self, result) # call superclass run method

    def test_aero_analysis_flat(self):
        OAS_prob = OASProblem({'type' : 'aero',
                               'optimize' : False,
                               'record_db' : False})
        OAS_prob.add_surface({'span_cos_spacing' : 0})
        OAS_prob.setup()
        OAS_prob.run()
        prob = OAS_prob.prob

        self.assertAlmostEqual(prob['wing_perf.CL'], .46173591841167, places=5)
        self.assertAlmostEqual(prob['wing_perf.CD'], .005524603647, places=5)

    def test_aero_analysis_flat_multiple(self):
        OAS_prob = OASProblem({'type' : 'aero',
                               'optimize' : False,
                               'record_db' : False})
        OAS_prob.add_surface({'span_cos_spacing' : 0.})
        OAS_prob.add_surface({'name' : 'tail',
                              'span_cos_spacing' : 0.,
                              'offset' : np.array([0., 0., 1000000.])})
        OAS_prob.setup()
        OAS_prob.run()
        prob = OAS_prob.prob
        self.assertAlmostEqual(prob['wing_perf.CL'], .46173591841167, places=5)
        self.assertAlmostEqual(prob['tail_perf.CL'], .46173591841167, places=5)

    def test_aero_analysis_flat_side_by_side(self):
        OAS_prob = OASProblem({'type' : 'aero',
                               'optimize' : False,
                               'record_db' : False})
        OAS_prob.add_surface({'name' : 'wing',
                              'span' : 5.,
                              'num_y' : 3,
                              'span_cos_spacing' : 0.,
                              'symmetry' : False,
                              'offset' : np.array([0., -2.5, 0.])})
        OAS_prob.add_surface({'name' : 'tail',
                              'span' : 5.,
                              'num_y' : 3,
                              'span_cos_spacing' : 0.,
                              'symmetry' : False,
                              'offset' : np.array([0., 2.5, 0.])})
        OAS_prob.setup()
        OAS_prob.run()
        prob = OAS_prob.prob
        self.assertAlmostEqual(prob['wing_perf.CL'], 0.46173591841167183, places=5)
        self.assertAlmostEqual(prob['tail_perf.CL'], 0.46173591841167183, places=5)
        self.assertAlmostEqual(prob['wing_perf.CD'], .005524603647, places=5)
        self.assertAlmostEqual(prob['tail_perf.CD'], .005524603647, places=5)

    def test_aero_analysis_flat_full(self):
        OAS_prob = OASProblem({'type' : 'aero',
                               'optimize' : False,
                               'record_db' : False})
        surf_dict = {'symmetry' : False}
        OAS_prob.add_surface(surf_dict)
        OAS_prob.setup()
        OAS_prob.run()
        prob = OAS_prob.prob
        self.assertAlmostEqual(prob['wing_perf.CL'], .45655138, places=5)
        self.assertAlmostEqual(prob['wing_perf.CD'], 0.0055402121081108589, places=5)

    def test_aero_analysis_flat_viscous_full(self):
        OAS_prob = OASProblem({'type' : 'aero',
                               'optimize' : False,
                               'record_db' : False,
                               'with_viscous' : True})
        surf_dict = {'symmetry' : False}
        OAS_prob.add_surface(surf_dict)
        OAS_prob.setup()
        OAS_prob.run()
        prob = OAS_prob.prob
        self.assertAlmostEqual(prob['wing_perf.CL'], .45655138, places=5)
        self.assertAlmostEqual(prob['wing_perf.CD'], 0.018942466133780547, places=5)

    if fortran_flag:
        def test_aero_optimization(self):
            # Need to use SLSQP here because SNOPT finds a different optimum
            OAS_prob = OASProblem({'type' : 'aero',
                                   'optimize' : True,
                                   'record_db' : False,
                                   'optimizer' : 'SLSQP'})
            OAS_prob.add_surface()

            OAS_prob.add_desvar('wing.twist_cp', lower=-10., upper=15.)
            OAS_prob.add_desvar('wing.sweep', lower=10., upper=30.)
            OAS_prob.add_desvar('wing.dihedral', lower=-10., upper=20.)
            OAS_prob.add_constraint('wing_perf.CL', equals=0.5)
            OAS_prob.add_objective('wing_perf.CD', scaler=1e4)

            OAS_prob.setup()

            OAS_prob.run()
            prob = OAS_prob.prob
            self.assertAlmostEqual(prob['wing_perf.CD'], 0.0049392534859265614, places=5)

    if fortran_flag:
        def test_aero_optimization_fd(self):
            # Need to use SLSQP here because SNOPT finds a different optimum
            OAS_prob = OASProblem({'type' : 'aero',
                                   'optimize' : True,
                                   'record_db' : False,
                                   'optimizer' : 'SLSQP',
                                   'force_fd' : True})
            OAS_prob.add_surface()

            OAS_prob.add_desvar('wing.twist_cp', lower=-10., upper=15.)
            OAS_prob.add_desvar('wing.sweep', lower=10., upper=30.)
            OAS_prob.add_desvar('wing.dihedral', lower=-10., upper=20.)
            OAS_prob.add_desvar('wing.taper', lower=.5, upper=2.)
            OAS_prob.add_constraint('wing_perf.CL', equals=0.5)
            OAS_prob.add_objective('wing_perf.CD', scaler=1e4)

            OAS_prob.setup()

            OAS_prob.run()
            prob = OAS_prob.prob
            self.assertAlmostEqual(prob['wing_perf.CD'], 0.0038513637269919619, places=5)

    if fortran_flag:
        def test_aero_optimization_chord_monotonic(self):
            OAS_prob = OASProblem({'type' : 'aero',
                                   'optimize' : True,
                                   'record_db' : False,
                                   'with_viscous' : False})
            OAS_prob.add_surface({
            'chord_cp' : np.random.random(5),
            'num_y' : 11,
            'monotonic_con' : ['chord'],
            'span_cos_spacing' : 0.,
            })

            OAS_prob.add_desvar('wing.chord_cp', lower=0.1, upper=5.)
            OAS_prob.add_desvar('alpha', lower=-10., upper=10.)
            OAS_prob.add_constraint('wing_perf.CL', equals=0.1)
            OAS_prob.add_constraint('wing.S_ref', equals=20)
            OAS_prob.add_constraint('wing.monotonic_chord', upper=0.)
            OAS_prob.add_objective('wing_perf.CD', scaler=1e4)

            OAS_prob.setup()

            OAS_prob.run()
            prob = OAS_prob.prob
            self.assertAlmostEqual(prob['wing_perf.CD'], 0.00057432581266351113, places=5)
            self.assertAlmostEqual(prob['wing.monotonic_chord'][0], -1.710374671671999, places=4)

    if fortran_flag:
        def test_aero_optimization_chord_monotonic_no_sym(self):
            OAS_prob = OASProblem({'type' : 'aero',
                                   'optimize' : True,
                                   'record_db' : False,
                                   'symmetry' : False,
                                   'with_viscous' : False})
            OAS_prob.add_surface({
            'chord_cp' : np.random.random(5),
            'num_y' : 11,
            'monotonic_con' : ['chord'],
            'span_cos_spacing' : 0.,
            })

            OAS_prob.add_desvar('wing.chord_cp', lower=0.1, upper=5.)
            OAS_prob.add_desvar('alpha', lower=-10., upper=10.)
            OAS_prob.add_constraint('wing_perf.CL', equals=0.1)
            OAS_prob.add_constraint('wing.S_ref', equals=20)
            OAS_prob.add_constraint('wing.monotonic_chord', upper=0.)
            OAS_prob.add_objective('wing_perf.CD', scaler=1e4)

            OAS_prob.setup()

            OAS_prob.run()
            prob = OAS_prob.prob
            self.assertAlmostEqual(prob['wing_perf.CD'], 0.00057432581266351113, places=5)
            self.assertAlmostEqual(prob['wing.monotonic_chord'][0], -1.710374671671999, places=2)

    if fortran_flag:
        def test_aero_viscous_optimization(self):
            OAS_prob = OASProblem({'type' : 'aero',
                                   'optimize' : True,
                                   'record_db' : False,
                                   'with_viscous' : True})
            OAS_prob.add_surface()

            OAS_prob.add_desvar('wing.twist_cp', lower=-10., upper=15.)
            OAS_prob.add_desvar('wing.sweep', lower=10., upper=30.)
            OAS_prob.add_desvar('wing.dihedral', lower=-10., upper=20.)
            OAS_prob.add_constraint('wing_perf.CL', equals=0.5)
            OAS_prob.add_objective('wing_perf.CD', scaler=1e4)

            OAS_prob.setup()

            OAS_prob.run()
            prob = OAS_prob.prob
            self.assertAlmostEqual(prob['wing_perf.CD'], 0.019234984422361764, places=5)

    if fortran_flag:
        def test_aero_viscous_chord_optimization(self):
            # Need to use SLSQP here because SNOPT finds a different optimum
            OAS_prob = OASProblem({'type' : 'aero',
                                   'optimize' : True,
                                   'record_db' : False,
                                   'optimizer' : 'SLSQP',
                                   'with_viscous' : True})
            OAS_prob.add_surface()

            OAS_prob.add_desvar('wing.chord_cp', lower=0.1, upper=3.)
            OAS_prob.add_desvar('alpha', lower=-10., upper=10.)
            OAS_prob.add_constraint('wing_perf.CL', equals=0.5)
            OAS_prob.add_constraint('wing.S_ref', equals=20)
            OAS_prob.add_objective('wing_perf.CD', scaler=1e4)

            OAS_prob.setup()

            OAS_prob.run()
            prob = OAS_prob.prob
            self.assertAlmostEqual(prob['wing_perf.CD'], 0.023570021288019886, places=5)

    if fortran_flag:
        def test_aero_multiple_opt(self):
            OAS_prob = OASProblem({'type' : 'aero',
                                   'optimize' : True,
                                   'record_db' : False})
            surf_dict = { 'name' : 'wing',
                          'span' : 5.,
                          'num_y' : 3,
                          'span_cos_spacing' : 0.}
            OAS_prob.add_surface(surf_dict)
            surf_dict.update({'name' : 'tail',
                           'offset' : np.array([0., 0., 10.])})
            OAS_prob.add_surface(surf_dict)

            OAS_prob.add_desvar('tail.twist_cp', lower=-10., upper=15.)
            OAS_prob.add_desvar('tail.sweep', lower=10., upper=30.)
            OAS_prob.add_desvar('tail.dihedral', lower=-10., upper=20.)
            OAS_prob.add_desvar('tail.taper', lower=.5, upper=2.)
            OAS_prob.add_constraint('tail_perf.CL', equals=0.5)
            OAS_prob.add_objective('tail_perf.CD', scaler=1e4)

            OAS_prob.setup()

            OAS_prob.run()
            prob = OAS_prob.prob
            self.assertAlmostEqual(prob['wing_perf.CL'], 0.41532382375677429, places=4)
            self.assertAlmostEqual(prob['tail_perf.CL'], .5, places=5)
            self.assertAlmostEqual(prob['wing_perf.CD'], .0075400306289957033, places=5)
            self.assertAlmostEqual(prob['tail_perf.CD'], 0.008087914662238814, places=5)


class TestStruct(unittest.TestCase):

    currentResult = []

    def run(self, result=None):
        self.currentResult.append(result) # remember result for use in tearDown
        unittest.TestCase.run(self, result) # call superclass run method

    def test_struct_analysis(self):
        OAS_prob = OASProblem({'type' : 'struct',
                               'optimize' : False,
                               'record_db' : False})
        surf_dict = {'symmetry' : False,
                    't_over_c' : 0.15}
        OAS_prob.add_surface(surf_dict)
        OAS_prob.setup()
        OAS_prob.run()
        prob = OAS_prob.prob
        self.assertAlmostEqual(prob['wing.structural_weight'], 3952.539819242561, places=3)

    def test_struct_analysis_symmetry(self):
        OAS_prob = OASProblem({'type' : 'struct',
                               'optimize' : False,
                               'record_db' : False})
        surf_dict = {'symmetry' : True,
                    't_over_c' : 0.15}
        OAS_prob.add_surface(surf_dict)
        OAS_prob.setup()
        OAS_prob.run()
        prob = OAS_prob.prob
        self.assertAlmostEqual(prob['wing.structural_weight'], 3952.539819242561, places=3)

    if fortran_flag:
        def test_struct_optimization(self):
            OAS_prob = OASProblem({'type' : 'struct',
                                   'optimize' : True,
                                   'record_db' : False})
            OAS_prob.add_surface({'symmetry' : False,
                                't_over_c': 0.15})

            OAS_prob.add_desvar('wing.thickness_cp', lower=0.001, upper=0.25, scaler=1e2)
            OAS_prob.add_constraint('wing.failure', upper=0.)
            OAS_prob.add_constraint('wing.thickness_intersects', upper=0.)
            OAS_prob.add_objective('wing.structural_weight', scaler=1e-3)

            OAS_prob.setup()

            OAS_prob.run()
            prob = OAS_prob.prob

            self.assertAlmostEqual(prob['wing.structural_weight'], 542.94945887080837, places=2)

    if fortran_flag:
        def test_struct_optimization_symmetry(self):
            OAS_prob = OASProblem({'type' : 'struct',
                                   'optimize' : True,
                                   'record_db' : False})
            OAS_prob.add_surface({'t_over_c': 0.15})

            OAS_prob.add_desvar('wing.thickness_cp', lower=0.001, upper=0.25, scaler=1e2)
            OAS_prob.add_constraint('wing.failure', upper=0.)
            OAS_prob.add_constraint('wing.thickness_intersects', upper=0.)
            OAS_prob.add_objective('wing.structural_weight', scaler=1e-3)

            OAS_prob.setup()

            OAS_prob.run()
            prob = OAS_prob.prob
            self.assertAlmostEqual(prob['wing.structural_weight'], 539.67300476678736, places=2)

    if fortran_flag:
        def test_struct_optimization_symmetry_exact(self):
            OAS_prob = OASProblem({'type' : 'struct',
                                   'optimize' : True,
                                   'record_db' : False})
            OAS_prob.add_surface({'exact_failure_constraint' : True,
                                't_over_c': 0.15})

            OAS_prob.add_desvar('wing.thickness_cp', lower=0.001, upper=0.25, scaler=1e2)
            OAS_prob.add_constraint('wing.failure', upper=0.)
            OAS_prob.add_constraint('wing.thickness_intersects', upper=0.)
            OAS_prob.add_objective('wing.structural_weight', scaler=1e-3)

            OAS_prob.setup()

            OAS_prob.run()
            prob = OAS_prob.prob
            self.assertAlmostEqual(prob['wing.structural_weight'], 536.44271005219036, places=2)


class TestAeroStruct(unittest.TestCase):

    currentResult = []

    def run(self, result=None):
        self.currentResult.append(result) # remember result for use in tearDown
        unittest.TestCase.run(self, result) # call superclass run method

    def test_aerostruct_analysis(self):
        OAS_prob = OASProblem({'type' : 'aerostruct',
                               'optimize' : False,
                               'record_db' : False})
        surf_dict = {'num_y' : 13,
                  'num_x' : 2,
                  'wing_type' : 'CRM',
                  'CD0' : 0.015,
                  'symmetry' : False}
        OAS_prob.add_surface(surf_dict)
        OAS_prob.setup()
        OAS_prob.run()
        prob = OAS_prob.prob
        self.assertAlmostEqual(prob['wing_perf.CL'], 0.53470610411546438)
        self.assertAlmostEqual(prob['wing_perf.failure'], -0.88555523062055341, places=5)
        self.assertAlmostEqual(prob['fuelburn'], 99367.43763782426, places=2)

    def test_aerostruct_analysis_symmetry(self):
        OAS_prob = OASProblem({'type' : 'aerostruct',
                               'optimize' : False,
                               'record_db' : False})
        surf_dict = {'symmetry' : True,
                  'num_y' : 13,
                  'num_x' : 2,
                  'wing_type' : 'CRM',
                  'CD0' : 0.015}
        OAS_prob.add_surface(surf_dict)
        OAS_prob.setup()
        OAS_prob.run()
        prob = OAS_prob.prob
        self.assertAlmostEqual(prob['wing_perf.CL'], 0.57096081806728627)
        self.assertAlmostEqual(prob['wing_perf.failure'], -0.8969902215221538, places=5)
        self.assertAlmostEqual(prob['fuelburn'], 105120.8601409321, places=2)

    def test_aerostruct_analysis_symmetry_deriv(self):
        OAS_prob = OASProblem({'type' : 'aerostruct',
                               'optimize' : False,
                               'record_db' : False})
        surf_dict = {'symmetry' : True,
                  'num_y' : 7,
                  'num_x' : 2,
                  'wing_type' : 'CRM',
                  'CD0' : 0.015}
        OAS_prob.add_surface(surf_dict)
        OAS_prob.setup()
        OAS_prob.run()
        prob = OAS_prob.prob

        data = prob.check_partial_derivatives(out_stream=None)

        new_dict = {}
        for key1 in data.keys():
            for key2 in data[key1].keys():
                for key3 in data[key1][key2].keys():
                    if 'rel' in key3:
                        error = np.linalg.norm(data[key1][key2][key3])
                        new_key = key1+'_'+key2[0]+'_'+key2[1]+'_'+key3
                        new_dict.update({new_key : error})

        for key in new_dict.keys():
            error = new_dict[key]
            if not np.isnan(error):

                # The FD check is not valid for these cases
                if 'assembly_forces_Iy' in key or 'assembly_forces_J' in key or \
                'assembly_forces_A' in key or 'assembly_K_loads' in key or \
                'assembly_forces_loads' in key or 'assembly_forces_Iz' in key or \
                'assembly_forces_nodes' in key or 'CM_wing_S_ref' in key or \
                'CM_rho' in key:
                    pass
                elif 'K' in key or 'vonmises' in key:
                    self.assertAlmostEqual(0., error, places=0)
                else:
                    self.assertAlmostEqual(0., error, places=2)

    if fortran_flag:
        def test_aerostruct_optimization(self):
            OAS_prob = OASProblem({'type' : 'aerostruct',
                                   'optimize' : True,
                                   'with_viscous' : True,
                                   'record_db' : False})
            surf_dict = {'num_y' : 7,
                      'num_x' : 2,
                      'wing_type' : 'CRM',
                      'CD0' : 0.015,
                      'symmetry' : False,
                      'num_twist_cp' : 2,
                      'num_thickness_cp' : 2}
            OAS_prob.add_surface(surf_dict)

            OAS_prob.add_desvar('wing.twist_cp', lower=-15., upper=15.)
            OAS_prob.add_desvar('wing.thickness_cp', lower=0.01, upper=0.5, scaler=1e2)
            OAS_prob.add_constraint('wing_perf.failure', upper=0.)
            OAS_prob.add_constraint('wing_perf.thickness_intersects', upper=0.)
            OAS_prob.add_desvar('alpha', lower=-10., upper=10.)
            OAS_prob.add_constraint('eq_con', equals=0.)
            OAS_prob.add_objective('fuelburn', scaler=1e-5)

            OAS_prob.setup()

            OAS_prob.run()
            prob = OAS_prob.prob
            self.assertAlmostEqual(prob['fuelburn'], 90621.476887491794, places=0)
            self.assertAlmostEqual(prob['wing_perf.failure'], 0., places=4)

    if fortran_flag:
        def test_aerostruct_optimization_symmetry(self):
            OAS_prob = OASProblem({'type' : 'aerostruct',
                                   'optimize' : True,
                                   'with_viscous' : True,
                                   'record_db' : False})
            surf_dict = {'symmetry' : True,
                      'num_y' : 7,
                      'num_x' : 3,
                      'wing_type' : 'CRM',
                      'CD0' : 0.015,
                      'num_twist_cp' : 2,
                      'num_thickness_cp' : 2}
            OAS_prob.add_surface(surf_dict)

            OAS_prob.add_desvar('wing.twist_cp', lower=-15., upper=15.)
            OAS_prob.add_desvar('wing.thickness_cp', lower=0.01, upper=0.5, scaler=1e2)
            OAS_prob.add_constraint('wing_perf.failure', upper=0.)
            OAS_prob.add_constraint('wing_perf.thickness_intersects', upper=0.)
            OAS_prob.add_desvar('alpha', lower=-10., upper=10.)
            OAS_prob.add_constraint('eq_con', equals=0.)
            OAS_prob.add_objective('fuelburn', scaler=1e-4)

            OAS_prob.setup()

            OAS_prob.run()
            prob = OAS_prob.prob
            self.assertAlmostEqual(prob['fuelburn'], 91972.867902621932, places=0)
            self.assertAlmostEqual(prob['wing_perf.failure'], 0, places=5)

    if fortran_flag:
        def test_aerostruct_optimization_symmetry_multiple(self):
            OAS_prob = OASProblem({'type' : 'aerostruct',
                                   'optimize' : True,
                                   'with_viscous' : True,
                                   'record_db' : False})
            surf_dict = {'name' : 'wing',
                         'symmetry' : True,
                         'num_y' : 5,
                         'num_x' : 2,
                         'wing_type' : 'CRM',
                         'CD0' : 0.015,
                         'num_twist_cp' : 2,
                         'num_thickness_cp' : 2}
            OAS_prob.add_surface(surf_dict)
            surf_dict.update({'name' : 'tail',
                              'offset':np.array([0., 0., 10.])})
            OAS_prob.add_surface(surf_dict)

            # Add design variables and constraints for both the wing and tail
            OAS_prob.add_desvar('wing.twist_cp', lower=-15., upper=15.)
            OAS_prob.add_desvar('wing.thickness_cp', lower=0.01, upper=0.5, scaler=1e2)
            OAS_prob.add_constraint('wing_perf.failure', upper=0.)
            OAS_prob.add_constraint('wing_perf.thickness_intersects', upper=0.)
            OAS_prob.add_desvar('tail.twist_cp', lower=-15., upper=15.)
            OAS_prob.add_desvar('tail.thickness_cp', lower=0.01, upper=0.5, scaler=1e2)
            OAS_prob.add_constraint('tail_perf.failure', upper=0.)
            OAS_prob.add_constraint('tail_perf.thickness_intersects', upper=0.)

            OAS_prob.add_desvar('alpha', lower=-10., upper=10.)
            OAS_prob.add_constraint('eq_con', equals=0.)
            OAS_prob.add_objective('fuelburn', scaler=1e-5)

            OAS_prob.setup()

            OAS_prob.run()
            prob = OAS_prob.prob
            self.assertAlmostEqual(prob['fuelburn'], 213046.09109747928, places=1)
            self.assertAlmostEqual(prob['wing_perf.failure'], 0, places=5)


if __name__ == "__main__":

    # Get user-supplied argument if provided
    try:
        arg = sys.argv[1]
        arg_provided = True
    except:
        arg_provided = False

    # Based on user input, run one subgroup of tests
    if arg_provided:
        if 'aero' == arg:
            test_classes = [TestAero]
        elif 'struct' == arg:
            test_classes = [TestStruct]
        elif 'aerostruct' == arg:
            test_classes = [TestAeroStruct]
    else:
        arg = 'full'
        test_classes = [TestAero, TestStruct, TestAeroStruct]

    print()
    print('+==================================================+')
    print('             Running ' + arg + ' test suite')
    print('+==================================================+')
    print()

    failures = []
    errors = []
    num_tests = 0

    # Loop through each requested discipline test
    for test_class in test_classes:

        # Set up the test suite and run the tests corresponding to this subgroup
        suite = unittest.TestLoader().loadTestsFromTestCase(test_class)
        unittest.TextTestRunner().run(suite)

        failures.extend(test_class.currentResult[-1].failures)
        errors.extend(test_class.currentResult[-1].errors)
        num_tests += len(test_class.currentResult)

    # Print results and force an exit if an error or failure occurred
    print()
    if len(failures) or len(errors):
        print("There have been errors or failures! Please check the log to " +
              "see which tests failed.\n")
        sys.exit(1)

    else:
        print("Successfully ran {} tests with no errors!\n".format(num_tests))
