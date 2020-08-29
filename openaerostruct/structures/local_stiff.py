import numpy as np

import openmdao.api as om


coeffs_2 = np.array([
    [ 1., -1.],
    [-1.,  1.],
])

coeffs_y = np.array([
    [ 12.,  -6., -12.,  -6.],
    [ -6.,   4.,   6.,   2.],
    [-12.,   6.,  12.,   6.],
    [ -6.,   2.,   6.,   4.],
])

coeffs_z = np.array([
    [ 12.,   6., -12.,   6.],
    [  6.,   4.,  -6.,   2.],
    [-12.,  -6.,  12.,  -6.],
    [  6.,   2.,  -6.,   4.],
])


class LocalStiff(om.ExplicitComponent):

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        surface = self.options['surface']

        self.ny = ny = surface['mesh'].shape[1]

        self.add_input('A', shape=ny - 1, units='m**2')
        self.add_input('J', shape=ny - 1, units='m**4')
        self.add_input('Iy', shape=ny - 1, units='m**4')
        self.add_input('Iz', shape=ny - 1, units='m**4')
        self.add_input('element_lengths', shape=ny - 1, units='m')

        self.add_output('local_stiff', shape=(ny - 1, 12, 12))

        rows = np.arange(144 * (ny - 1))
        cols = np.outer(np.arange(ny - 1), np.ones(144, int)).flatten()

        self.declare_partials('local_stiff', 'A', rows=rows, cols=cols)
        self.declare_partials('local_stiff', 'J', rows=rows, cols=cols)
        self.declare_partials('local_stiff', 'Iy', rows=rows, cols=cols)
        self.declare_partials('local_stiff', 'Iz', rows=rows, cols=cols)
        self.declare_partials('local_stiff', 'element_lengths', rows=rows, cols=cols)

    def compute(self, inputs, outputs):
        surface = self.options['surface']

        ny = self.ny
        E = surface['E']
        G = surface['G']

        A  = inputs['A']
        Iy = inputs['Iy']
        Iz = inputs['Iz']
        J  = inputs['J']
        L  = inputs['element_lengths']

        outputs['local_stiff'] = 0.

        for i in range(2):
            for j in range(2):
                outputs['local_stiff'][:, 0 + i, 0 + j] = E * A / L * coeffs_2[i, j]
                outputs['local_stiff'][:, 2 + i, 2 + j] = G * J / L * coeffs_2[i, j]

        for i in range(4):
            for j in range(4):
                outputs['local_stiff'][:, 4 + i, 4 + j] = E * Iy / L ** 3 * coeffs_y[i, j]
                outputs['local_stiff'][:, 8 + i, 8 + j] = E * Iz / L ** 3 * coeffs_z[i, j]

        for i in [1, 3]:
            for j in range(4):
                outputs['local_stiff'][:, 4 + i, 4 + j] *= L
                outputs['local_stiff'][:, 8 + i, 8 + j] *= L
        for i in range(4):
            for j in [1, 3]:
                outputs['local_stiff'][:, 4 + i, 4 + j] *= L
                outputs['local_stiff'][:, 8 + i, 8 + j] *= L

    def compute_partials(self, inputs, partials):
        surface = self.options['surface']
        ny = surface['mesh'].shape[1]
        E = surface['E']
        G = surface['G']

        A  = inputs['A']
        Iy = inputs['Iy']
        Iz = inputs['Iz']
        J  = inputs['J']
        L  = inputs['element_lengths']

        derivs_A = partials['local_stiff', 'A'].reshape((ny - 1, 12, 12))
        derivs_Iy = partials['local_stiff', 'Iy'].reshape((ny - 1, 12, 12))
        derivs_Iz = partials['local_stiff', 'Iz'].reshape((ny - 1, 12, 12))
        derivs_J = partials['local_stiff', 'J'].reshape((ny - 1, 12, 12))
        derivs_L = partials['local_stiff', 'element_lengths'].reshape((ny - 1, 12, 12))

        derivs_A[:] = 0.
        derivs_Iy[:] = 0.
        derivs_Iz[:] = 0.
        derivs_J[:] = 0.
        derivs_L[:] = 0.

        for i in range(2):
            for j in range(2):
                derivs_A[:, 0 + i, 0 + j] = E / L * coeffs_2[i, j]
                derivs_L[:, 0 + i, 0 + j] = -E * A / L ** 2 * coeffs_2[i, j]

                derivs_J[:, 2 + i, 2 + j] = G / L * coeffs_2[i, j]
                derivs_L[:, 2 + i, 2 + j] = -G * J / L ** 2 * coeffs_2[i, j]

        for i in range(4):
            for j in range(4):
                derivs_Iy[:, 4 + i, 4 + j] = E / L ** 3 * coeffs_y[i, j]
                derivs_L [:, 4 + i, 4 + j] = -3 * E * Iy / L ** 4 * coeffs_y[i, j]

                derivs_Iz[:, 8 + i, 8 + j] = E / L ** 3 * coeffs_z[i, j]
                derivs_L [:, 8 + i, 8 + j] = -3 * E * Iz / L ** 4 * coeffs_z[i, j]

        for i in [1, 3]:
            for j in range(4):
                derivs_Iy[:, 4 + i, 4 + j] = E / L ** 2 * coeffs_y[i, j]
                derivs_L[:, 4 + i, 4 + j] = -2 * E * Iy / L ** 3 * coeffs_y[i, j]

                derivs_Iz[:, 8 + i, 8 + j] = E / L ** 2 * coeffs_z[i, j]
                derivs_L[:, 8 + i, 8 + j] = -2 * E * Iz / L ** 3 * coeffs_z[i, j]
        for i in range(4):
            for j in [1, 3]:
                derivs_Iy[:, 4 + i, 4 + j] = E / L ** 2 * coeffs_y[i, j]
                derivs_L[:, 4 + i, 4 + j] = -2 * E * Iy / L ** 3 * coeffs_y[i, j]

                derivs_Iz[:, 8 + i, 8 + j] = E / L ** 2 * coeffs_z[i, j]
                derivs_L[:, 8 + i, 8 + j] = -2 * E * Iz / L ** 3 * coeffs_z[i, j]

        for i in [1, 3]:
            for j in [1, 3]:
                derivs_Iy[:, 4 + i, 4 + j] = E / L * coeffs_y[i, j]
                derivs_L[:, 4 + i, 4 + j] = -E * Iy / L ** 2 * coeffs_y[i, j]

                derivs_Iz[:, 8 + i, 8 + j] = E / L * coeffs_z[i, j]
                derivs_L[:, 8 + i, 8 + j] = -E * Iz / L ** 2 * coeffs_z[i, j]
