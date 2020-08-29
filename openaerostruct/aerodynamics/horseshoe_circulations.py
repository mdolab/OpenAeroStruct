import numpy as np
from scipy.sparse import csc_matrix

import openmdao.api as om


class HorseshoeCirculations(om.ExplicitComponent):
    """
    Convert the previously-computed vortex ring circulations into horseshoe
    circulations. Vortex rings and horseshoe vortices produce the same linear
    space, but with a different parameterization. It's easier to compute the
    circulations using a vortex ring approach, but it's easier to compute the
    forces acting on the surface by using the horseshoe circulations.
    That's why we have this component, to convert from one circulation
    space to the other.

    Parameters
    ----------
    circulations[system_size] : numpy array
        The vortex ring circulations obtained by solving the AIC linear system.

    Returns
    -------
    horseshoe_circulations[system_size] : numpy array
        The equivalent horseshoe circulations obtained by intelligently summing
        the vortex ring circulations, accounting for overlaps between rings.
    """

    def initialize(self):
        self.options.declare('surfaces', types=list)

    def setup(self):
        surfaces = self.options['surfaces']

        system_size = 0

        # Loop through all the surfaces to obtain the total system size,
        # which is the number of panels in the total system.
        for surface in surfaces:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]

            system_size += (nx - 1) * (ny - 1)

        self.system_size = system_size

        self.add_input('circulations', shape=system_size, units='m**2/s')
        self.add_output('horseshoe_circulations', shape=system_size, units='m**2/s')

        # To convert between the two circulations, we simply need to set up a
        # matrix that linearly transforms the vortex ring circulations to
        # the horseshoe circulations. Again, because this is a linear
        # transformation, the derivatives are in fact the matrix itself.
        data = [np.ones(system_size)]
        rows = [np.arange(system_size)]
        cols = [np.arange(system_size)]

        ind_1 = 0
        ind_2 = 0
        for surface in surfaces:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            num = (nx - 1) * (ny - 1)

            ind_2 += num

            arange = np.arange(num).reshape((nx - 1), (ny - 1))

            data_ = -np.ones((nx - 2) * (ny - 1))
            rows_ = ind_1 + arange[1:, :].flatten()
            cols_ = ind_1 + arange[:-1, :].flatten()

            data.append(data_)
            rows.append(rows_)
            cols.append(cols_)

            ind_1 += num

        data = np.concatenate(data)
        rows = np.concatenate(rows)
        cols = np.concatenate(cols)

        # Actually create the sparse matrix based on these rows and cols
        self.mtx = csc_matrix((data, (rows, cols)), shape=(system_size, system_size))

        self.declare_partials('horseshoe_circulations', 'circulations', val=data, rows=rows, cols=cols)

    def compute(self, inputs, outputs):
        outputs['horseshoe_circulations'] = self.mtx.dot(inputs['circulations'])
