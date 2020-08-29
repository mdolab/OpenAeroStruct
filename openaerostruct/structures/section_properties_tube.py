import numpy as np

import openmdao.api as om

class SectionPropertiesTube(om.ExplicitComponent):
    """
    Compute geometric properties for a tube element.
    The thicknesses are added to the interior of the element, so the
    'radius' value is the outer radius of the tube.

    parameters
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

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        self.surface = surface = self.options['surface']

        self.ny = surface['mesh'].shape[1]

        self.add_input('radius', val=np.ones((self.ny - 1)), units='m')
        self.add_input('thickness', val=np.ones((self.ny - 1)) * .1, units='m')
        self.add_output('A', val=np.zeros((self.ny - 1)), units='m**2')
        self.add_output('Iy', val=np.zeros((self.ny - 1)), units='m**4')
        self.add_output('Iz', val=np.zeros((self.ny - 1)), units='m**4')
        self.add_output('J', val=np.zeros((self.ny - 1)), units='m**4')

        a = np.arange((self.ny - 1))
        self.declare_partials('*', '*', rows=a, cols=a)
        self.set_check_partial_options(wrt='*', method='cs')

    def compute(self, inputs, outputs):
        pi = np.pi

        # Add thickness to the interior of the radius.
        # The outer radius is the inputs['radius'] amount.
        r1 = inputs['radius'] - inputs['thickness']
        r2 = inputs['radius']

        # Compute the area, area moments of inertia, and polar moment of inertia
        outputs['A'] = pi * (r2**2 - r1**2)
        outputs['Iy'] = pi * (r2**4 - r1**4) / 4.
        outputs['Iz'] = pi * (r2**4 - r1**4) / 4.
        outputs['J'] = pi * (r2**4 - r1**4) / 2.

    def compute_partials(self, inputs, partials):
        pi = np.pi
        radius = inputs['radius'].real
        t = inputs['thickness'].real
        r1 = radius - t
        r2 = radius

        dr1_dr = 1.
        dr2_dr = 1.
        dr1_dt = -1.
        dr2_dt =  0.

        r1_3 = r1**3
        r2_3 = r2**3

        partials['A', 'radius'] = 2 * pi * (r2 * dr2_dr - r1 * dr1_dr)
        partials['A', 'thickness'] = 2 * pi * (r2 * dr2_dt - r1 * dr1_dt)
        partials['Iy', 'radius'] = pi * (r2_3 * dr2_dr - r1_3 * dr1_dr)
        partials['Iy', 'thickness'] = pi * (r2_3 * dr2_dt - r1_3 * dr1_dt)
        partials['Iz', 'radius'] = pi * (r2_3 * dr2_dr - r1_3 * dr1_dr)
        partials['Iz', 'thickness'] = pi * (r2_3 * dr2_dt - r1_3 * dr1_dt)
        partials['J', 'radius'] = 2 * pi * (r2_3 * dr2_dr - r1_3 * dr1_dr)
        partials['J', 'thickness'] = 2 * pi * (r2_3 * dr2_dt - r1_3 * dr1_dt)
