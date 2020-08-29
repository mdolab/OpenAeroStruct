import numpy as np

import openmdao.api as om


class CreateRHS(om.ExplicitComponent):
    """
    Compute the right-hand-side of the K * u = f linear system to solve for the displacements.
    The RHS is based on the loads. For the aerostructural case, these are
    recomputed at each design point based on the aerodynamic loads.

    Parameters
    ----------
    loads[ny, 6] : numpy array
        Flattened array containing the loads applied on the FEM component,
        computed from the sectional forces.

    Returns
    -------
    forces[6*(ny+1)] : numpy array
        Right-hand-side of the linear system. The loads from the aerodynamic
        analysis or the user-defined loads.
    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        surface = self.options['surface']

        self.ny = surface['mesh'].shape[1]

        self.add_input('total_loads', val=np.zeros((self.ny, 6)), units='N')
        self.add_output('forces', val=np.ones(((self.ny+1)*6)), units='N')

        n = self.ny * 6
        arange = np.arange((n))
        self.declare_partials('forces', 'total_loads', val=1., rows=arange, cols=arange)

    def compute(self, inputs, outputs):
        outputs['forces'][:] = 0.

        # Populate the right-hand side of the linear system using the
        # prescribed or computed loads
        outputs['forces'][:6*self.ny] += inputs['total_loads'].reshape(self.ny*6)

        # Remove extremely small values from the RHS so the linear system
        # can more easily be solved
        outputs['forces'][np.abs(outputs['forces']) < 1e-6] = 0.
