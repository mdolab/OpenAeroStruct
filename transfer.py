""" Define the transfer components to couple aero and struct analyses. """

from __future__ import division
import numpy

from openmdao.api import Component


class TransferDisplacements(Component):
    """
    Perform displacement transfer.

    Apply the computed displacements on the original mesh to obtain
    the deformed mesh.

    Parameters
    ----------
    mesh : array_like
        Flattened array defining the lifting surfaces.
    disp : array_like
        Flattened array containing displacements on the FEM component.
        Contains displacements for all six degrees of freedom, including
        displacements in the x, y, and z directions, and rotations about the
        x, y, and z axes.

    Returns
    -------
    def_mesh : array_like
        Flattened array defining the lifting surfaces after deformation.

    """

    def __init__(self, aero_ind, fem_ind, fem_origin=0.35):
        super(TransferDisplacements, self).__init__()

        tot_n = numpy.sum(aero_ind[:, 2])
        self.aero_ind = aero_ind
        self.fem_ind = fem_ind
        tot_n_fem = numpy.sum(fem_ind[:, 0])
        self.fem_origin = fem_origin

        self.add_param('mesh', val=numpy.zeros((tot_n, 3), dtype="complex"))
        self.add_param('disp', val=numpy.zeros((tot_n_fem, 6),
                       dtype="complex"))
        self.add_output('def_mesh', val=numpy.zeros((tot_n, 3),
                        dtype="complex"))

        self.deriv_options['type'] = 'cs'
        # self.deriv_options['form'] = 'central'
        #self.deriv_options['extra_check_partials_form'] = "central"

    def solve_nonlinear(self, params, unknowns, resids):
        #
        # from time import time
        # st = time()

        for i_surf, row in enumerate(self.fem_ind):
            nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = \
                self.aero_ind[i_surf, :]
            n_fem, i_fem = row

            mesh = params['mesh'][i:i+n, :].reshape(nx, ny, 3)
            disp = params['disp'][i_fem:i_fem+n_fem]

            w = self.fem_origin
            ref_curve = (1-w) * mesh[0, :, :] + w * mesh[-1, :, :]

            Smesh = numpy.zeros(mesh.shape, dtype="complex")
            for ind in xrange(nx):
                Smesh[ind, :, :] = mesh[ind, :, :] - ref_curve

            def_mesh = numpy.zeros(mesh.shape, dtype="complex")
            cos, sin = numpy.cos, numpy.sin
            for ind in xrange(ny):
                dx, dy, dz, rx, ry, rz = disp[ind, :]

                # 1 eye from the axis rotation matrices
                # -3 eye from subtracting Smesh three times
                T = -2 * numpy.eye(3, dtype="complex")
                T[ 1:,  1:] += [[cos(rx), -sin(rx)], [ sin(rx), cos(rx)]]
                T[::2, ::2] += [[cos(ry),  sin(ry)], [-sin(ry), cos(ry)]]
                T[ :2,  :2] += [[cos(rz), -sin(rz)], [ sin(rz), cos(rz)]]

                def_mesh[:, ind, :] += Smesh[:, ind, :].dot(T)
                def_mesh[:, ind, 0] += dx
                def_mesh[:, ind, 1] += dy
                def_mesh[:, ind, 2] += dz

            unknowns['def_mesh'][i:i+n, :] = \
                (def_mesh + mesh).reshape(n, 3).astype("complex")
            # print time() - st


class TransferLoads(Component):
    """
    Perform aerodynamic load transfer.

    Apply the computed sectional forces on the aerodynamic surfaces to
    obtain the deformed mesh FEM loads.

    Parameters
    ----------
    def_mesh : array_like
        Flattened array defining the lifting surfaces after deformation.
    sec_forces : array_like
        Flattened array containing the sectional forces acting on each panel.
        Stored in Fortran order (only relevant when more than one chordwise
        panel).

    Returns
    -------
    loads : array_like
        Flattened array containing the loads applied on the FEM component,
        computed from the sectional forces.

    """

    def __init__(self, aero_ind, fem_ind, fem_origin=0.35):
        super(TransferLoads, self).__init__()

        tot_n = numpy.sum(aero_ind[:, 2])
        tot_panels = numpy.sum(aero_ind[:, 4])
        self.aero_ind = aero_ind
        self.fem_ind = fem_ind
        tot_n_fem = numpy.sum(fem_ind[:, 0])
        self.fem_origin = fem_origin

        self.add_param('def_mesh', val=numpy.zeros((tot_n, 3)))
        self.add_param('sec_forces', val=numpy.zeros((tot_panels, 3),
                       dtype="complex"))
        self.add_output('loads', val=numpy.zeros((tot_n_fem, 6),
                        dtype="complex"))

        self.deriv_options['type'] = 'cs'
        self.deriv_options['form'] = 'central'
        #self.deriv_options['extra_check_partials_form'] = "central"

    def solve_nonlinear(self, params, unknowns, resids):
        for i_surf, row in enumerate(self.fem_ind):
            nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = \
                self.aero_ind[i_surf, :]
            n_fem, i_fem = row

            mesh = params['def_mesh'][i:i+n, :].reshape(nx, ny, 3)

            sec_forces = params['sec_forces'][i_panels:i_panels+n_panels, :]. \
                reshape(nx-1, ny-1, 3, order='F')
            sec_forces = numpy.sum(sec_forces, axis=0)

            w = 0.25
            a_pts = 0.5 * (1-w) * mesh[:-1, :-1, :] + \
                    0.5 *   w   * mesh[1:, :-1, :] + \
                    0.5 * (1-w) * mesh[:-1,  1:, :] + \
                    0.5 *   w   * mesh[1:,  1:, :]

            w = self.fem_origin
            s_pts = 0.5 * (1-w) * mesh[:-1, :-1, :] + \
                    0.5 *   w   * mesh[1:, :-1, :] + \
                    0.5 * (1-w) * mesh[:-1,  1:, :] + \
                    0.5 *   w   * mesh[1:,  1:, :]

            moment = numpy.zeros((ny - 1, 3), dtype="complex")
            for ind in xrange(ny - 1):
                r = a_pts[0, ind, :] - s_pts[0, ind, :]
                F = sec_forces[ind, :]
                moment[ind, :] = numpy.cross(r, F)

            loads = numpy.zeros((ny, 6), dtype="complex")
            loads[:-1, :3] += 0.5 * sec_forces[:, :]
            loads[ 1:, :3] += 0.5 * sec_forces[:, :]
            loads[:-1, 3:] += 0.5 * moment
            loads[ 1:, 3:] += 0.5 * moment

            unknowns['loads'][i_fem:i_fem+n_fem, :] = loads
