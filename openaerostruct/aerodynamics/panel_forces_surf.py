import numpy as np

import openmdao.api as om


class PanelForcesSurf(om.ExplicitComponent):
    """
    Take in the computed panel forces and convert them to sectional forces
    for each surface. Basically just takes the long array that has info
    for all surfaces and creates a new output for each surface with only
    that surface's panel forces.

    Parameters
    ----------
    panel_forces[system_size, 3] : numpy array
        All of the forces acting on all panels in the total system.

    Returns
    -------
    sec_forces[nx-1, ny-1, 3] : numpy array
        Only the panel forces for one individual lifting surface.
        There is one of these per surface.
    """

    def initialize(self):
        self.options.declare('surfaces', types=list)

    def setup(self):
        surfaces = self.options['surfaces']

        system_size = 0

        # Loop through the surfaces to get the total system size
        for surface in surfaces:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]

            system_size += (nx - 1) * (ny - 1)

        arange = np.arange(3 * system_size)

        self.add_input('panel_forces', shape=(system_size, 3), units='N')

        # Loop through the surfaces and add the output of sec_forces based on
        # the size of each surface. Here we keep track of the total indices
        # from panel_forces to make sure the forces go to the correct output
        ind1, ind2 = 0, 0
        for surface in surfaces:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']

            sec_forces_name = '{}_sec_forces'.format(name)

            ind2 += (nx - 1) * (ny - 1) * 3

            self.add_output(sec_forces_name, shape=(nx - 1, ny - 1, 3), units='N')

            rows = np.arange((nx - 1) * (ny - 1) * 3)
            cols = arange[ind1:ind2]
            self.declare_partials(sec_forces_name, 'panel_forces', val=1., rows=rows, cols=cols)

            ind1 += (nx - 1) * (ny - 1) * 3

    def compute(self, inputs, outputs):
        surfaces = self.options['surfaces']

        ind1, ind2 = 0, 0
        for surface in surfaces:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']

            sec_forces_name = '{}_sec_forces'.format(name)

            ind2 += (nx - 1) * (ny - 1)

            # Just pluck out the relevant forces and reshape them
            outputs[sec_forces_name] = inputs['panel_forces'][ind1:ind2].reshape(
                (nx - 1, ny - 1, 3))

            ind1 += (nx - 1) * (ny - 1)
