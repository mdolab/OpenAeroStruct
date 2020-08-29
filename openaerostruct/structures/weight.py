import numpy as np

import openmdao.api as om
from openaerostruct.structures.utils import norm


class Weight(om.ExplicitComponent):
    """ Compute total weight and center-of-gravity location of the spar elements.

    parameters
    ----------
    A[ny-1] : numpy array
        Areas for each FEM element.
    nodes[ny, 3] : numpy array
        Flattened array with coordinates for each FEM node.

    Returns
    -------
    structural_mass : float
        Weight of the wing structure.
    elmenet_weight[ny-1] : float
        Weight of each element.

    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        self.surface = surface = self.options['surface']

        self.ny = ny = surface['mesh'].shape[1]

        self.add_input('A', val=np.ones((self.ny - 1)), units='m**2')
        self.add_input('nodes', val=np.zeros((self.ny, 3)), units='m')

        self.add_output('structural_mass', val=0., units='kg')
        self.add_output('element_mass', val=np.zeros((self.ny-1)), units='kg')

        self.declare_partials('structural_mass', ['A','nodes'])

        row_col = np.arange(self.ny-1, dtype=int)
        self.declare_partials('element_mass','A', rows=row_col, cols=row_col)

        dimensions = 3
        rows=np.empty((dimensions*2*(ny-1)))
        cols=np.empty((dimensions*2*(ny-1)))
        for i in range (ny-1):
            rows[i*dimensions*2:(i+1)*dimensions*2] = i
            cols[i*dimensions*2:(i+1)*dimensions*2] = np.linspace(i*dimensions,i*dimensions+(dimensions*2-1),dimensions*2)
        self.declare_partials('element_mass','nodes', rows=rows, cols=cols)

        self.set_check_partial_options('*', method='cs', step=1e-40)

    def compute(self, inputs, outputs):
        A = inputs['A']
        nodes = inputs['nodes']
        mrho = self.surface['mrho']
        wwr = self.surface['wing_weight_ratio']

        # Calculate the volume and weight of the structure
        element_volumes = norm(nodes[1:, :] - nodes[:-1, :], axis=1) * A

        # nodes[1:, :] - nodes[:-1, :] this is the delta array of the different between the points
        element_mass = element_volumes * mrho * wwr
        weight = np.sum(element_mass)

        # If the tube is symmetric, double the computed weight
        if self.surface['symmetry']:
            weight *= 2.

        outputs['structural_mass'] = weight
        outputs['element_mass'] = element_mass

    def compute_partials(self, inputs, partials):

        A = inputs['A']
        nodes = inputs['nodes']
        mrho = self.surface['mrho']
        wwr = self.surface['wing_weight_ratio']
        ny = self.ny

        # Calculate the volume and weight of the structure
        const0 = nodes[1:, :] - nodes[:-1, :]
        const1 = np.linalg.norm(const0, axis=1)
        element_volumes = const1 * A
        volume = np.sum(element_volumes)
        const2 = mrho * wwr
        weight = volume * const2

        # First we will solve for dweight_dA
        # Calculate the volume and weight of the total structure
        norms = const1.reshape(1, -1)

        # Multiply by the material density and force of gravity
        dweight_dA = norms * const2

        # Account for symmetry
        if self.surface['symmetry']:
            dweight_dA *= 2.

        # Save the result to the jacobian dictionary
        partials['structural_mass', 'A'] = dweight_dA

        # Next, we will compute the derivative of weight wrt nodes.
        # Here we're using results from AD to compute the derivative
        # Initialize the reverse seeds.
        nodesb = np.zeros(nodes.shape)
        tempb = (const0) * (A / norms).reshape(-1, 1)
        nodesb[1:, :] += tempb
        nodesb[:-1, :] -= tempb

        # Apply the multipliers for material properties and symmetry
        nodesb *= mrho * wwr

        if self.surface['symmetry']:
            nodesb *= 2.

        # Store the flattened array in the jacobian dictionary
        partials['structural_mass', 'nodes'] = nodesb.reshape(1, -1)

        # Element_weight Partials
        partials['element_mass','A'] = const1 * const2

        precalc = np.sum(np.power(const0,2),axis=1)
        d__dprecalc = 0.5 * precalc**(-.5)

        dimensions = 3
        for i in range(ny-1):
            first_part = const0[i,:] * d__dprecalc[i] * 2 * (-1) * A[i] * const2
            second_part = const0[i,:] * d__dprecalc[i] * 2 * A[i] * const2
            partials['element_mass', 'nodes'][i*dimensions*2:(i+1)*dimensions*2] = np.append(first_part,second_part)
