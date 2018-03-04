from __future__ import division, print_function
import numpy as np

from openmdao.api import Component, Group, IndepVarComp, Problem, ScipyOptimizer
from scipy.linalg import lu_factor, lu_solve
import OAS_API

nx = 2
ny = 3

class simple_multiplication(Component):

    def __init__(self):
        super(simple_multiplication, self).__init__()

        self.add_param('x', val=np.ones(nx))
        self.add_output('y', val=np.ones(ny))

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['y'] = OAS_API.oas_api.mult(ny, params['x'])

    def linearize(self, params, unknowns, resids):
        jac = self.alloc_jacobian()
        jac['y', 'x'][:] = 0.
        return jac

    def apply_linear(self, params, unknowns, dparams, dunknowns, dresids, mode):
        if mode == 'fwd':
            dresids['y'] += OAS_API.oas_api.mult_d(ny, params['x'], dparams['x'])[1]

        if mode == 'rev':
            dparams['x'] += OAS_API.oas_api.mult_b(params['x'], unknowns['y'], dresids['y'])

root = Group()
root.add('x', IndepVarComp('x', np.linspace(0., 1., nx)))
root.add('simp_mult', simple_multiplication())
root.connect('x.x', 'simp_mult.x')

prob = Problem(root)

prob.driver = ScipyOptimizer()
prob.driver.options['optimizer'] = 'SLSQP'
prob.driver.options['maxiter'] = 11

# prob.print_all_convergence(level=2)


prob.driver.add_desvar('x.x', lower=-50, upper=50)
prob.driver.add_objective('simp_mult.y')

prob.setup()
prob.run()

result = prob['simp_mult.x']
print(result)
