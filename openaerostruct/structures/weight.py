from __future__ import print_function, division
import numpy as np

from openmdao.api import ExplicitComponent

try:
    from openaerostruct.fortran import OAS_API
    fortran_flag = True
    data_type = float
except:
    fortran_flag = False
    data_type = complex

class Weight(ExplicitComponent):
    """ Compute total weight and center-of-gravity location of the spar elements.

    parameters
    ----------
    A[ny-1] : numpy array
        Areas for each FEM element.
    nodes[ny, 3] : numpy array
        Flattened array with coordinates for each FEM node.

    Returns
    -------
    structural_weight : float
        Weight of the structural spar.
    cg_location[3] : numpy array
        Location of the structural spar's cg.

    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        self.surface = surface = self.options['surface']

        self.ny = surface['num_y']

        self.add_input('A', val=np.zeros((self.ny - 1)), units='m**2')#, dtype=data_type))
        self.add_input('nodes', val=np.zeros((self.ny, 3)), units='m')#, dtype=data_type))
        self.add_input('load_factor', val=1.)
        self.add_output('structural_weight', val=0., units='N')
        self.add_output('element_weights', val=np.zeros((self.ny-1)), units='N')

        self.declare_partials('*', '*')

        self.declare_partials('element_weights', '*', method='fd')

    def compute(self, inputs, outputs):
        A = inputs['A']
        nodes = inputs['nodes']

        # Calculate the volume and weight of the structure
        element_volumes = np.linalg.norm(nodes[1:, :] - nodes[:-1, :], axis=1) * A
        element_weights = element_volumes * self.surface['mrho'] * 9.81 * self.surface['wing_weight_ratio'] * inputs['load_factor']
        weight = np.sum(element_weights)
        volume = np.sum(element_volumes)

        # If the tube is symmetric, double the computed weight and set the
        # y-location of the cg to 0, at the symmetry plane
        if self.surface['symmetry']:
            weight *= 2.

        outputs['structural_weight'] = weight
        outputs['element_weights'] = element_weights

    def compute_partials(self, inputs, partials):

        A = inputs['A']
        nodes = inputs['nodes']

        # Calculate the volume and weight of the structure
        element_volumes = np.linalg.norm(nodes[1:, :] - nodes[:-1, :], axis=1) * A
        volume = np.sum(element_volumes)
        weight = volume * self.surface['mrho'] * 9.81 * self.surface['wing_weight_ratio'] * inputs['load_factor']

        # First we will solve for dweight_dA
        # Calculate the volume and weight of the total structure
        norms = np.linalg.norm(nodes[1:, :] - nodes[:-1, :], axis=1).reshape(1, -1)

        # Multiply by the material density and force of gravity
        dweight_dA = norms * self.surface['mrho'] * 9.81 * self.surface['wing_weight_ratio'] * inputs['load_factor']

        # Account for symmetry
        if self.surface['symmetry']:
            dweight_dA *= 2.

        # Save the result to the jacobian dictionary
        partials['structural_weight', 'A'] = dweight_dA

        # Next, we will compute the derivative of weight wrt nodes.
        # Here we're using results from AD to compute the derivative
        # Initialize the reverse seeds.
        nodesb = np.zeros(nodes.shape)
        tempb = (nodes[1:, :] - nodes[:-1, :]) * (A / norms).reshape(-1, 1)
        nodesb[1:, :] += tempb
        nodesb[:-1, :] -= tempb

        # Apply the multipliers for material properties and symmetry
        nodesb *= self.surface['mrho'] * 9.81 * self.surface['wing_weight_ratio'] * inputs['load_factor']

        if self.surface['symmetry']:
            nodesb *= 2.
            partials['structural_weight', 'load_factor'] = weight * 2
        else:
            partials['structural_weight', 'load_factor'] = weight

        # Store the flattened array in the jacobian dictionary
        partials['structural_weight', 'nodes'] = nodesb.reshape(1, -1)
