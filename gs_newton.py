from openmdao.solvers.solver_base import NonLinearSolver
from openmdao.solvers.newton import Newton
from openmdao.solvers.nl_gauss_seidel import NLGaussSeidel


class HybridGSNewton(NonLinearSolver):

    def __init__(self):
        super(HybridGSNewton, self).__init__()

        self.nlgs = NLGaussSeidel()
        self.newton = Newton()

        self.nlgs.options['maxiter'] = 5
        self.newton.options['maxiter'] = 1
        self.nlgs.options['atol'] = 1e-10
        self.newton.options['atol'] = 1e-10

    def setup(self, sub):
        """ Initialize sub solvers.

        Args
        ----
        sub: `System`
            System that owns this solver.
        """
        self.nlgs.setup(sub)
        self.newton.setup(sub)

        # Declare to the world we need derivs
        self.ln_solver = self.newton.ln_solver

    def solve(self, params, unknowns, resids, system, metadata=None):

        self.nlgs.solve(params, unknowns, resids, system, metadata)
        self.newton.solve(params, unknowns, resids, system, metadata)
