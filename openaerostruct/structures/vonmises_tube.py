import numpy as np

import openmdao.api as om

from openaerostruct.structures.utils import norm, unit, norm_d, unit_d, cross_d

class VonMisesTube(om.ExplicitComponent):
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

        ny = self.ny = surface['mesh'].shape[1]

        self.add_input('nodes', val=np.zeros((ny, 3)), units='m')
        self.add_input('radius', val=np.zeros((ny - 1)), units='m')
        self.add_input('disp', val=np.zeros((ny, 6)), units='m')

        self.add_output('vonmises', val=np.zeros((ny-1, 2)), units='N/m**2')

        self.E = surface['E']
        self.G = surface['G']

        row = np.concatenate([np.zeros(6), np.ones(6)])
        rows = np.tile(row, ny-1) + np.repeat(2*np.arange(ny-1), 12)
        col = np.tile(np.arange(6), 2)
        cols = np.tile(col, ny-1) + np.repeat(3*np.arange(ny-1), 12)

        self.declare_partials('*', 'nodes', rows=rows, cols=cols)

        rows = np.arange(2 * (ny-1))
        cols = np.repeat(np.arange(ny-1), 2)

        self.declare_partials('*', 'radius', rows=rows, cols=cols)

        row = np.concatenate([np.zeros(12), np.ones(12)])
        rows = np.tile(row, ny-1) + np.repeat(2*np.arange(ny-1), 24)
        col = np.tile(np.arange(12), 2)
        cols = np.tile(col, ny-1) + np.repeat(6*np.arange(ny-1), 24)

        self.declare_partials('*', 'disp', rows=rows, cols=cols)

    def compute(self, inputs, outputs):
        dtype = float
        if self.under_complex_step:
            dtype = complex

        self.T = np.zeros((3, 3),dtype=dtype)
        self.x_gl = np.array([1, 0, 0],dtype=dtype)
        radius = inputs['radius']
        disp = inputs['disp']
        nodes = inputs['nodes']
        T = self.T
        E = self.E
        G = self.G
        x_gl = self.x_gl

        num_elems = self.ny - 1
        for ielem in range(num_elems):

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

    def compute_partials(self, inputs, partials):

        radius = inputs['radius']
        disp = inputs['disp']
        nodes = inputs['nodes']
        T = self.T
        E = self.E
        G = self.G
        x_gl = self.x_gl

        num_elems = self.ny - 1
        for ielem in range(num_elems):

            # Compute the coordinate delta between the two element end points
            P0 = nodes[ielem, :]
            P1 = nodes[ielem+1, :]
            dP = P1 - P0

            # Compute the derivative of element length
            L = norm(dP)
            dLddP = norm_d(dP)

            # unit function converts a vector to a unit vector
            # calculate the transormation to the local element frame.
            # We use x_gl to provide a reference axis to reference
            x_loc = unit(dP)
            dxdP = unit_d(dP)

            y_loc = unit(np.cross(x_loc, x_gl))
            dtmpdx, _ = cross_d(x_loc, x_gl)
            dydtmp = unit_d(np.cross(x_loc, x_gl))
            dydP = dydtmp.dot(dtmpdx).dot(dxdP)

            z_loc = unit(np.cross(x_loc, y_loc))
            dtmpdx, dtmpdy = cross_d(x_loc,y_loc)
            dzdtmp = unit_d(np.cross(x_loc, y_loc))
            dzdP = dzdtmp.dot(dtmpdx).dot(dxdP) + dzdtmp.dot(dtmpdy).dot(dydP)

            T[0, :] = x_loc
            T[1, :] = y_loc
            T[2, :] = z_loc

            u0x = x_loc.dot(disp[ielem, :3])
            r0x, r0y, r0z = T.dot(disp[ielem, 3:])
            u1x = x_loc.dot(disp[ielem+1, :3])
            r1x, r1y, r1z = T.dot(disp[ielem+1, 3:])

            # #$$$$$$$$$$$$$$$$$$$$$$$$$$

            # The derivatives of the above code wrt displacement all boil down to sections of the T matrix
            dxddisp = T[0,:]
            dyddisp = T[1,:]
            dzddisp = T[2,:]

            #The derivatives of the above code wrt T all boil down to sections of the #displacement vector
            du0dloc = disp[ielem, :3]
            dr0dloc = disp[ielem, 3:]
            du1dloc = disp[ielem+1, :3]
            dr1dloc = disp[ielem+1, 3:]

            #$$$$$$$$$$$$$$$$$$$$$$$$$$
            # Original code
            # $$$$$$$$$$$$
            tmp = np.sqrt((r1y - r0y)**2 + (r1z - r0z)**2) + 1e-50 #added eps to avoid 0 disp singularity
            sxx0 = E * (u1x - u0x) / L + E * radius[ielem] / L * tmp
            sxx1 = E * (u0x - u1x) / L + E * radius[ielem] / L * tmp
            sxt = G * radius[ielem] * (r1x - r0x) / L
            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

            dtmpdr0y = 1/tmp * (r1y - r0y)*-1
            dtmpdr1y = 1/tmp * (r1y - r0y)
            dtmpdr0z = 1/tmp * (r1z - r0z)*-1
            dtmpdr1z = 1/tmp * (r1z - r0z)

            # Combine all of the derivtives for tmp
            dtmpdDisp = np.zeros(12)
            dr0xdDisp = np.zeros(12)
            dr1xdDisp = np.zeros(12)

            #r0 term
            dr0xdDisp[3:6] = dxddisp
            dr1xdDisp[9:12] = dxddisp

            dtmpdDisp[3:6] = dtmpdr0y*dyddisp
            dtmpdDisp[3:6] += dtmpdr0z*dzddisp
            dtmpdDisp[9:12] = dtmpdr1y*dyddisp
            dtmpdDisp[9:12] += dtmpdr1z*dzddisp

            # x_loc, y_loc and z_loc terms
            # (dttmpx_loc is zeros, so don't compute with it)
            dtmpdy_loc = dtmpdr0y*dr0dloc + dtmpdr1y*dr1dloc
            dtmpdz_loc = dtmpdr0z*dr0dloc + dtmpdr1z*dr1dloc

            dtmpdP = dtmpdy_loc.dot(dydP) + dtmpdz_loc.dot(dzdP)

            dsxx0dtmp = E * radius[ielem] / L
            dsxx0du0x = -E / L
            dsxx0du1x = E / L
            dsxx0dL =  -E * (u1x - u0x) / (L*L) - E * radius[ielem] / (L*L) * tmp

            dsxx1dtmp = E * radius[ielem] / L
            dsxx1du0x = E / L
            dsxx1du1x = -E / L
            dsxx1dL = -E * (u0x - u1x) / (L*L) - E * radius[ielem] / (L*L) * tmp

            dsxx0dP = dsxx0dtmp * dtmpdP + \
                      dsxx0du0x*du0dloc.dot(dxdP) + dsxx0du1x*du1dloc.dot(dxdP)+\
                      dsxx0dL*dLddP

            dsxx1dP = dsxx1dtmp * dtmpdP + \
                      dsxx1du0x*du0dloc.dot(dxdP)+dsxx1du1x*du1dloc.dot(dxdP)+\
                      dsxx1dL*dLddP

            # Combine sxx0 and sxx1 terms

            # Start with the tmp term
            dsxx0dDisp = dsxx0dtmp * dtmpdDisp
            dsxx1dDisp = dsxx1dtmp * dtmpdDisp

            # Now add the direct u dep
            dsxx0dDisp[0:3] = dsxx0du0x * dxddisp
            dsxx0dDisp[6:9] = dsxx0du1x * dxddisp

            dsxx1dDisp[0:3] = dsxx1du0x * dxddisp
            dsxx1dDisp[6:9] = dsxx1du1x * dxddisp

            # Combine sxt term
            dsxtdr0x = -G * radius[ielem] / L
            dsxtdr1x = G * radius[ielem] / L
            dsxtdL =  - G * radius[ielem] * (r1x - r0x) / (L*L)

            dsxtdP = dsxtdr0x*(dr0dloc.dot(dxdP)) + dsxtdr1x*(dr1dloc.dot(dxdP)) + \
                     dsxtdL*dLddP
            #disp
            dsxtdDisp = dsxtdr0x * dr0xdDisp + dsxtdr1x * dr1xdDisp

            #radius derivatives
            dsxxdrad = E / L * tmp
            dsxtdrad = G * (r1x - r0x)/L

            fact = 1.0 / (np.sqrt(sxx0**2 + 3 * sxt**2))
            dVm0dsxx0 = sxx0 * fact
            dVm0dsxt = 3 * sxt * fact

            fact = 1.0 / (np.sqrt(sxx1**2 + 3 * sxt**2))
            dVm1dsxx1 = sxx1 * fact
            dVm1dsxt = 3 * sxt * fact

            ii = 2 * ielem
            partials['vonmises', 'radius'][ii] = dVm0dsxx0*dsxxdrad + dVm0dsxt*dsxtdrad
            partials['vonmises', 'radius'][ii+1] = dVm1dsxx1*dsxxdrad + dVm1dsxt*dsxtdrad

            ii = 24 * ielem
            partials['vonmises', 'disp'][ii:ii+12] = dVm0dsxx0*dsxx0dDisp + dVm0dsxt*dsxtdDisp
            partials['vonmises', 'disp'][ii+12:ii+24] = dVm1dsxx1*dsxx1dDisp + dVm1dsxt*dsxtdDisp

            # Compute terms for the nodes
            ii = 12 * ielem

            dVm0_dnode = dVm0dsxx0*dsxx0dP + dVm0dsxt*dsxtdP
            partials['vonmises', 'nodes'][ii:ii+3] = -dVm0_dnode
            partials['vonmises', 'nodes'][ii+3:ii+6] = dVm0_dnode

            dVM1_dnode = dVm1dsxx1*dsxx1dP + dVm1dsxt*dsxtdP
            partials['vonmises', 'nodes'][ii+6:ii+9] = -dVM1_dnode
            partials['vonmises', 'nodes'][ii+9:ii+12] = dVM1_dnode
