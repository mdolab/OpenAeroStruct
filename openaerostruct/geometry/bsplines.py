from __future__ import print_function, division
import numpy as np
import scipy.sparse

from openmdao.api import ExplicitComponent


def get_bspline_mtx(num_cp, num_pt, order=4):
    """ Create Jacobian to fit a bspline to a set of data.

    Parameters
    ----------
    num_cp : int
        Number of control points.
    num_pt : int
        Number of points.
    order : int, optional
        Order of b-spline fit.

    Returns
    -------
    out : CSR sparse matrix
        Matrix that gives the points vector when multiplied by the control
        points vector.

    """

    knots = np.zeros(num_cp + order)
    knots[order-1:num_cp+1] = np.linspace(0, 1, num_cp - order + 2)
    knots[num_cp+1:] = 1.0
    t_vec = np.linspace(0, 1, num_pt)

    basis = np.zeros(order)
    arange = np.arange(order)
    data = np.zeros((num_pt, order))
    rows = np.zeros((num_pt, order), int)
    cols = np.zeros((num_pt, order), int)

    for ipt in range(num_pt):
        t = t_vec[ipt]

        i0 = -1
        for ind in range(order, num_cp+1):
            if (knots[ind-1] <= t) and (t < knots[ind]):
                i0 = ind - order
        if t == knots[-1]:
            i0 = num_cp - order

        basis[:] = 0.
        basis[-1] = 1.

        for i in range(2, order+1):
            l = i - 1
            j1 = order - l
            j2 = order
            n = i0 + j1
            if knots[n+l] != knots[n]:
                basis[j1-1] = (knots[n+l] - t) / \
                              (knots[n+l] - knots[n]) * basis[j1]
            else:
                basis[j1-1] = 0.
            for j in range(j1+1, j2):
                n = i0 + j
                if knots[n+l-1] != knots[n-1]:
                    basis[j-1] = (t - knots[n-1]) / \
                                (knots[n+l-1] - knots[n-1]) * basis[j-1]
                else:
                    basis[j-1] = 0.
                if knots[n+l] != knots[n]:
                    basis[j-1] += (knots[n+l] - t) / \
                                  (knots[n+l] - knots[n]) * basis[j]
            n = i0 + j2
            if knots[n+l-1] != knots[n-1]:
                basis[j2-1] = (t - knots[n-1]) / \
                              (knots[n+l-1] - knots[n-1]) * basis[j2-1]
            else:
                basis[j2-1] = 0.

        data[ipt, :] = basis
        rows[ipt, :] = ipt
        cols[ipt, :] = i0 + arange

    data, rows, cols = data.flatten(), rows.flatten(), cols.flatten()

    return scipy.sparse.csr_matrix((data, (rows, cols)),
                                   shape=(num_pt, num_cp))


class Bsplines(ExplicitComponent):
    """
    General function to translate from control points to actual points
    using a b-spline representation.

    Parameters
    ----------
    cpname : string
        Name of the OpenMDAO component containing the control point values.
    ptname : string
        Name of the OpenMDAO component that will contain the interpolated
        b-spline values.
    n_input : int
        Number of input control points.
    n_output : int
        Number of outputted interpolated b-spline points.
    """

    def initialize(self):
        self.options.declare('jac', is_valid=lambda jac: len(jac.shape) == 2)
        self.options.declare('num_cp', types=int)
        self.options.declare('num_pt', types=int)
        self.options.declare('in_name', types=str)
        self.options.declare('out_name', types=str)

    def setup(self):
        meta = self.options

        meta['jac'] = get_bspline_mtx(meta['num_cp'], meta['num_pt'], order=min(meta['num_cp'], 4))

        # Add the correct units; all bspline quantities are in m except twist
        if 'twist' in meta['in_name']:
            units = None
        else:
            units = 'm'

        self.add_input(meta['in_name'], val=np.zeros(meta['num_cp']), units=units)
        self.add_output(meta['out_name'], val=np.zeros(meta['num_pt']), units=units)

        meta = self.options

        jac = meta['jac'].tocoo()

        self.declare_partials(
            meta['out_name'], meta['in_name'], val=jac.data, rows=jac.row, cols=jac.col)

    def compute(self, inputs, outputs):
        meta = self.options
        outputs[meta['out_name']] = meta['jac'] * inputs[meta['in_name']]
