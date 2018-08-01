from __future__ import print_function, division
import numpy as np

from openmdao.api import ExplicitComponent
from openaerostruct.structures.utils import \
    _assemble_system

try:
    from openaerostruct.fortran import OAS_API
    fortran_flag = True
    data_type = float
except:
    fortran_flag = False
    data_type = complex

class AssembleK(ExplicitComponent):
    """
    Compute the displacements and rotations by solving the linear system
    using the structural stiffness matrix.

    Parameters
    ----------
    A[ny-1] : numpy array
        Areas for each FEM element.
    Iy[ny-1] : numpy array
        Area moment of inertia around the y-axis for each FEM element.
    Iz[ny-1] : numpy array
        Area moment of inertia around the z-axis for each FEM element.
    J[ny-1] : numpy array
        Polar moment of inertia for each FEM element.
    nodes[ny, 3] : numpy array
        Flattened array with coordinates for each FEM node.

    Returns
    -------
    K[6*(ny+1), 6*(ny+1)] : numpy array
        Stiffness matrix for the entire FEM system. Used to solve the linear
        system K * u = f to obtain the displacements, u.
    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        surface = self.options['surface']

        self.ny = surface['num_y']

        self.size = size = 6 * self.ny + 6

        self.add_input('A', val=np.ones((self.ny - 1)), units='m**2')#, dtype=data_type))
        self.add_input('Iy', val=np.ones((self.ny - 1)), units='m**4')#, dtype=data_type))
        self.add_input('Iz', val=np.ones((self.ny - 1)), units='m**4')#, dtype=data_type))
        self.add_input('J', val=np.ones((self.ny - 1)), units='m**4')#, dtype=data_type))
        self.add_input('nodes', val=np.ones((self.ny, 3)), units='m')#, dtype=data_type))
        self.add_output('K', val=np.ones((size, size), dtype=data_type), units='N/m')

        # Get material properties from the surface dictionary
        self.E = surface['E']
        self.G = surface['G']

        # Set up arrays to easily set up the K matrix
        self.const2 = np.array([
            [1, -1],
            [-1, 1],
        ], dtype=data_type)
        self.const_y = np.array([
            [12, -6, -12, -6],
            [-6, 4, 6, 2],
            [-12, 6, 12, 6],
            [-6, 2, 6, 4],
        ], dtype=data_type)
        self.const_z = np.array([
            [12, 6, -12, 6],
            [6, 4, -6, 2],
            [-12, -6, 12, -6],
            [6, 2, -6, 4],
        ], dtype=data_type)
        self.x_gl = np.array([1, 0, 0], dtype=data_type)

        self.K_elem = np.zeros((12, 12), dtype=data_type)
        self.T_elem = np.zeros((12, 12), dtype=data_type)
        self.T = np.zeros((3, 3), dtype=data_type)

        self.K = np.zeros((size, size), dtype=data_type)
        self.forces = np.zeros((size), dtype=data_type)

        self.K_a = np.zeros((2, 2), dtype=data_type)
        self.K_t = np.zeros((2, 2), dtype=data_type)
        self.K_y = np.zeros((4, 4), dtype=data_type)
        self.K_z = np.zeros((4, 4), dtype=data_type)

        self.S_a = np.zeros((2, 12), dtype=data_type)
        self.S_a[(0, 1), (0, 6)] = 1.

        self.S_t = np.zeros((2, 12), dtype=data_type)
        self.S_t[(0, 1), (3, 9)] = 1.

        self.S_y = np.zeros((4, 12), dtype=data_type)
        self.S_y[(0, 1, 2, 3), (2, 4, 8, 10)] = 1.

        self.S_z = np.zeros((4, 12), dtype=data_type)
        self.S_z[(0, 1, 2, 3), (1, 5, 7, 11)] = 1.

        self.declare_partials('*', '*')

        if not fortran_flag:
            self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):

        # Find constrained nodes based on closeness to central point
        nodes = inputs['nodes']
        dist = nodes - np.array([5., 0, 0])
        idx = (np.linalg.norm(dist, axis=1)).argmin()
        self.cons = idx

        self.K = \
            _assemble_system(inputs['nodes'],
                             inputs['A'], inputs['J'], inputs['Iy'],
                             inputs['Iz'], self.K_a, self.K_t,
                             self.K_y, self.K_z, self.cons,
                             self.E, self.G, self.x_gl, self.T, self.K_elem,
                             self.S_a, self.S_t, self.S_y, self.S_z,
                             self.T_elem, self.const2, self.const_y,
                             self.const_z, self.ny, self.size,
                             self.K)

        outputs['K'] = self.K

    if fortran_flag:
        def compute_partials(self, inputs, partials):

            for param in inputs:

                d_inputs = {}
                d_inputs[param] = inputs[param].copy()
                d_outputs = {}

                for j, val in enumerate(np.array(d_inputs[param]).flatten()):
                    d_in_b = np.array(d_inputs[param]).flatten()
                    d_in_b[:] = 0.
                    d_in_b[j] = 1.
                    d_inputs[param] = d_in_b.reshape(d_inputs[param].shape)

                    # Find constrained nodes based on closeness to specified cg point
                    nodes = inputs['nodes']
                    dist = nodes - np.array([5., 0, 0])
                    idx = (np.linalg.norm(dist, axis=1)).argmin()
                    self.cons = idx

                    A = inputs['A']
                    J = inputs['J']
                    Iy = inputs['Iy']
                    Iz = inputs['Iz']

                    if 'nodes' not in d_inputs:
                        d_inputs['nodes'] = inputs['nodes'] * 0
                    if 'A' not in d_inputs:
                        d_inputs['A'] = inputs['A'] * 0
                    if 'J' not in d_inputs:
                        d_inputs['J'] = inputs['J'] * 0
                    if 'Iy' not in d_inputs:
                        d_inputs['Iy'] = inputs['Iy'] * 0
                    if 'Iz' not in d_inputs:
                        d_inputs['Iz'] = inputs['Iz'] * 0

                    K, Kd = OAS_API.oas_api.assemblestructmtx_d(nodes, d_inputs['nodes'], A, d_inputs['A'],
                                                 J, d_inputs['J'], Iy, d_inputs['Iy'],
                                                 Iz, d_inputs['Iz'],
                                                 self.K_a, self.K_t, self.K_y, self.K_z,
                                                 self.cons, self.E, self.G, self.x_gl, self.T,
                                                 self.K_elem, self.S_a, self.S_t, self.S_y, self.S_z, self.T_elem,
                                                 self.const2, self.const_y, self.const_z)

                    partials['K', param][:, j] = Kd.flatten()
