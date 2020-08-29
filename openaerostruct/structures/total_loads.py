import numpy as np

import openmdao.api as om


class TotalLoads(om.ExplicitComponent):
    """
    Add the loads from the aerodynamics, structural weight, and fuel weight.

    Parameters
    ----------
    loads[ny, 6] : numpy array
        Array containing the loads applied on the FEM component,
        computed from the sectional forces.
    struct_weight_loads[ny, 6] : numpy array
        Array containing the loads applied on the FEM component,
        computed from the weight of the wing-structure segments.
    fuel_weight_loads[ny, 6] : numpy array
        Array containing the loads applied on the FEM component,
        computed from the weight of the fuel.
    loads_from_point_masses[ny, 6] : numpy array
        The cumulative loads from all point masses.

    Returns
    -------
    total_loads[ny, 6] : numpy array
        Flattened array containing the total loads applied on the FEM.
    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        self.surface = surface = self.options['surface']
        self.ny = surface['mesh'].shape[1]

        self.add_input('loads', val=np.ones((self.ny, 6)), units='N')
        if surface['struct_weight_relief']:
            self.add_input('struct_weight_loads', val=np.zeros((self.ny, 6)), units='N')
        if surface['distributed_fuel_weight']:
            self.add_input('fuel_weight_loads', val=np.zeros((self.ny, 6)), units='N')
        if 'n_point_masses' in surface.keys():
            self.add_input('loads_from_point_masses', val=np.zeros((self.ny, 6)), units='N')
            self.add_input('loads_from_thrusts', val=np.zeros((self.ny, 6)), units='N')

        self.add_output('total_loads', val=np.ones((self.ny, 6)), units='N')

        arange = np.arange(self.ny * 6)

        self.declare_partials('total_loads', 'loads',
            rows=arange, cols=arange, val=1.)

        if self.surface['struct_weight_relief']:
            self.declare_partials('total_loads', 'struct_weight_loads',
                rows=arange, cols=arange, val=1.)

        if self.surface['distributed_fuel_weight']:
            self.declare_partials('total_loads', 'fuel_weight_loads',
                rows=arange, cols=arange, val=1.)

        if 'n_point_masses' in surface.keys():
            self.declare_partials('total_loads', 'loads_from_point_masses',
                rows=arange, cols=arange, val=1.)
            self.declare_partials('total_loads', 'loads_from_thrusts',
                rows=arange, cols=arange, val=1.)

    def compute(self, inputs, outputs):
        outputs['total_loads'] = inputs['loads']

        if self.surface['struct_weight_relief']:
            outputs['total_loads'] += inputs['struct_weight_loads']

        if self.surface['distributed_fuel_weight']:
            outputs['total_loads'] += inputs['fuel_weight_loads']

        if 'n_point_masses' in self.surface.keys():
            outputs['total_loads'] += inputs['loads_from_point_masses']
            outputs['total_loads'] += inputs['loads_from_thrusts']
