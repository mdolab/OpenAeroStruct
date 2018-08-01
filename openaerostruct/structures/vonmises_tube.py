from __future__ import print_function, division
import numpy as np

from openmdao.api import ExplicitComponent

from openaerostruct.structures.utils import norm, unit

try:
    from openaerostruct.fortran import OAS_API
    fortran_flag = True
    data_type = float
except:
    fortran_flag = False
    data_type = complex

class VonMisesTube(ExplicitComponent):
    """ Compute the von Mises stress in each element.

    parameters
    ----------
    nodes[ny, 3] : numpy array
        Flattened array with coordinates for each FEM node.
    radius[ny-1] : numpy array
        Radii for each FEM element.
    disp[ny, 6] : numpy array
        Displacements of each FEM node.

    Returns
    -------
    vonmises[ny-1, 2] : numpy array
        von Mises stress magnitudes for each FEM element.

    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        self.surface = surface = self.options['surface']

        self.ny = surface['num_y']

        self.add_input('nodes', val=np.zeros((self.ny, 3)), units='m')#,  dtype=data_type))
        self.add_input('radius', val=np.zeros((self.ny - 1)), units='m')#,  dtype=data_type))
        self.add_input('disp', val=np.zeros((self.ny, 6)), units='m')#,  dtype=data_type))

        self.add_output('vonmises', val=np.zeros((self.ny-1, 2)), units='N/m**2')#,dtype=data_type))

        self.E = surface['E']
        self.G = surface['G']

        self.T = np.zeros((3, 3), dtype=data_type)
        self.x_gl = np.array([1, 0, 0], dtype=data_type)

        self.declare_partials('*', '*')

        if not fortran_flag:
            self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        radius = inputs['radius']
        disp = inputs['disp']
        nodes = inputs['nodes']
        T = self.T
        E = self.E
        G = self.G
        x_gl = self.x_gl

        if fortran_flag:
            vm = OAS_API.oas_api.calc_vonmises(nodes, radius, disp, E, G, x_gl)
            outputs['vonmises'] = vm

        else:

            num_elems = self.ny - 1
            for ielem in range(self.ny-1):

                P0 = nodes[ielem, :]
                P1 = nodes[ielem+1, :]
                L = norm(P1 - P0)

                x_loc = unit(P1 - P0)
                y_loc = unit(np.cross(x_loc, x_gl))
                z_loc = unit(np.cross(x_loc, y_loc))

                T[0, :] = x_loc
                T[1, :] = y_loc
                T[2, :] = z_loc

                u0x, u0y, u0z = T.dot(disp[ielem, :3])
                r0x, r0y, r0z = T.dot(disp[ielem, 3:])
                u1x, u1y, u1z = T.dot(disp[ielem+1, :3])
                r1x, r1y, r1z = T.dot(disp[ielem+1, 3:])

                tmp = np.sqrt((r1y - r0y)**2 + (r1z - r0z)**2)
                sxx0 = E * (u1x - u0x) / L + E * radius[ielem] / L * tmp
                sxx1 = E * (u0x - u1x) / L + E * radius[ielem] / L * tmp
                sxt = G * radius[ielem] * (r1x - r0x) / L

                outputs['vonmises'][ielem, 0] = np.sqrt(sxx0**2 + 3 * sxt**2)
                outputs['vonmises'][ielem, 1] = np.sqrt(sxx1**2 + 3 * sxt**2)


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

                    radius = inputs['radius']
                    disp = inputs['disp']
                    nodes = inputs['nodes']

                    E = self.E
                    G = self.G
                    x_gl = self.x_gl

                    if 'nodes' not in d_inputs:
                        d_inputs['nodes'] = inputs['nodes'] * 0
                    if 'radius' not in d_inputs:
                        d_inputs['radius'] = inputs['radius'] * 0
                    if 'disp' not in d_inputs:
                        d_inputs['disp'] = inputs['disp'] * 0

                    _, vonmisesd = OAS_API.oas_api.calc_vonmises_d(nodes, d_inputs['nodes'], radius, d_inputs['radius'], disp, d_inputs['disp'], E, G, x_gl)

                    partials['vonmises', param][:, j] = vonmisesd.flatten()
