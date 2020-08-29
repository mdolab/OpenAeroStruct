import numpy as np

import openmdao.api as om
from openaerostruct.structures.utils import norm
from openaerostruct.utils.constants import grav_constant


class ComputePointMassLoads(om.ExplicitComponent):
    """
    Compute the loads on the structure due to point masses.
    The current method adds loads and moments to all of the structural nodes, but
    the nodes closest to the point masses receive proportionally larger loads.

    Parameters
    ----------
    point_mass_locations[n_point_masses, 3] : numpy array
        XYZ location for each point mass in the global frame.
    point_masses[n_point_masses] : numpy array
        Actual magnitude of each point mass, in same order as point_mass_locations.
    nodes[ny, 3] : numpy array
        Flattened array with coordinates for each FEM node.
    load_factor : float
        Load factor for the flight point. Multiplier on the effects of gravity.

    Returns
    -------
    nodal_weightings[n_point_masses, ny] : numpy array
        The normalized weighting factor for each node. The closest nodes have the
        greatest weighting, while farther away nodes have less.
    loads_from_point_masses[ny, 6] : numpy array
        The actual loads array that will be added to the total loads array.
        This is cumulative and includes the forces and moments from all point masses.
    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        self.surface = surface = self.options['surface']
        self.ny = surface['mesh'].shape[1]
        self.n_point_masses = surface['n_point_masses']

        self.add_input('point_mass_locations', val=np.zeros((self.n_point_masses, 3)), units='m')
        self.add_input('point_masses', val=np.zeros((self.n_point_masses)), units='kg')
        self.add_input('nodes', val=np.zeros((self.ny, 3)), units='m')
        self.add_input('load_factor', val=1.0)

        self.add_output('nodal_weightings', val=np.zeros((self.n_point_masses, self.ny)))
        self.add_output('loads_from_point_masses', val=np.zeros((self.ny, 6)), units='N') ## WARNING!!! UNITS ARE A MIXTURE OF N & N*m

        self.set_check_partial_options(wrt='*', method='fd')
        self.declare_partials('*', '*', method='cs')

    def compute(self, inputs, outputs):
        nodes = inputs['nodes']
        nodal_weightings = outputs['nodal_weightings']
        loads_from_point_masses = outputs['loads_from_point_masses']
        loads_from_point_masses[:] = 0.

        # Loop through each point mass location, incrementing idx
        for idx, point_mass_location in enumerate(inputs['point_mass_locations']):

            # Get the vector between the nodes and the point mass location
            xyz_dist = point_mass_location - nodes

            # The y-distance between the nodes and the point mass location
            span_dist = xyz_dist[:, 1]

            # Compute the normalized inverse distance weightings for all of the
            # nodes. These weightings determine the amount of the force and
            # moment that each of the nodes receive.
            # nodal_weightings are scalars for each node.
            dist10 = span_dist**10  # TODO: investigate effect of power here
            inv_dist10 = 1 / (dist10 + 1e-10)
            nodal_weightings[idx, :] = inv_dist10 / np.sum(inv_dist10)

            # Create an array with the inverse distance weighted vectors in the
            # downward direction to account for gravity. Each node has a 3-vector
            # in the downward direction whose magnitude is based on nodal_weightings.
            directional_weightings = np.outer(nodal_weightings[idx, :], np.array([[0., 0., -1.]]))

            # Compute the perceived weight due to the point mass on each node in N
            weight_forces = directional_weightings * grav_constant * inputs['load_factor'] * inputs['point_masses'][idx]

            # Actually add the perceived weights to the load array
            loads_from_point_masses[:, :3] += weight_forces

            # Compute the moments based on the Euclidean distance vectors from
            # the nodes to the point mass, crossed with the perceived weight
            moments = np.cross(xyz_dist, weight_forces)

            # Accumulate the moments to the load array
            loads_from_point_masses[:, 3:] += moments
