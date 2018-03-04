from __future__ import division, print_function
import numpy as np

from openmdao.api import Component

class MaterialsTube(Component):
    """
    Compute geometric properties for a tube element.
    The thicknesses are added to the interior of the element, so the
    'radius' value is the outer radius of the tube.

    Parameters
    ----------
    radius : numpy array
        Outer radii for each FEM element.
    thickness : numpy array
        Tube thickness for each FEM element.

    Returns
    -------
    A : numpy array
        Cross-sectional area for each FEM element.
    Iy : numpy array
        Area moment of inertia around the y-axis for each FEM element.
    Iz : numpy array
        Area moment of inertia around the z-axis for each FEM element.
    J : numpy array
        Polar moment of inertia for each FEM element.
    """

    def __init__(self, surface):
        super(MaterialsTube, self).__init__()

        self.surface = surface

        self.ny = surface['num_y']
        self.mesh = surface['mesh']
        name = surface['name']

        self.add_param('radius', val=np.zeros((self.ny - 1)))
        self.add_param('thickness', val=np.zeros((self.ny - 1)))
        self.add_output('A', val=np.zeros((self.ny - 1)))
        self.add_output('Iy', val=np.zeros((self.ny - 1)))
        self.add_output('Iz', val=np.zeros((self.ny - 1)))
        self.add_output('J', val=np.zeros((self.ny - 1)))

        self.arange = np.arange((self.ny - 1))

    def solve_nonlinear(self, params, unknowns, resids):
        name = self.surface['name']
        pi = np.pi

        # Add thickness to the interior of the radius.
        # The outer radius is the params['radius'] amount.
        r1 = params['radius'] - params['thickness']
        r2 = params['radius']

        # Compute the area, area moments of inertia, and polar moment of inertia
        unknowns['A'] = pi * (r2**2 - r1**2)
        unknowns['Iy'] = pi * (r2**4 - r1**4) / 4.
        unknowns['Iz'] = pi * (r2**4 - r1**4) / 4.
        unknowns['J'] = pi * (r2**4 - r1**4) / 2.

    def linearize(self, params, unknowns, resids):
        name = self.surface['name']
        jac = self.alloc_jacobian()

        pi = np.pi
        radius = params['radius'].real
        t = params['thickness'].real
        r1 = radius - t
        r2 = radius

        dr1_dr = 1.
        dr2_dr = 1.
        dr1_dt = -1.
        dr2_dt =  0.

        r1_3 = r1**3
        r2_3 = r2**3

        a = self.arange
        jac['A', 'radius'][a, a] = 2 * pi * (r2 * dr2_dr - r1 * dr1_dr)
        jac['A', 'thickness'][a, a] = 2 * pi * (r2 * dr2_dt - r1 * dr1_dt)
        jac['Iy', 'radius'][a, a] = pi * (r2_3 * dr2_dr - r1_3 * dr1_dr)
        jac['Iy', 'thickness'][a, a] = pi * (r2_3 * dr2_dt - r1_3 * dr1_dt)
        jac['Iz', 'radius'][a, a] = pi * (r2_3 * dr2_dr - r1_3 * dr1_dr)
        jac['Iz', 'thickness'][a, a] = pi * (r2_3 * dr2_dt - r1_3 * dr1_dt)
        jac['J', 'radius'][a, a] = 2 * pi * (r2_3 * dr2_dr - r1_3 * dr1_dr)
        jac['J', 'thickness'][a, a] = 2 * pi * (r2_3 * dr2_dt - r1_3 * dr1_dt)

        return jac
