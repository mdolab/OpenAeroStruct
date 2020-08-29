import numpy as np

import openmdao.api as om

from openaerostruct.functionals.moment_coefficient import MomentCoefficient
from openaerostruct.functionals.total_lift_drag import TotalLiftDrag
from openaerostruct.functionals.sum_areas import SumAreas

class TotalAeroPerformance(om.Group):
    """
    Group to contain the total aerodynamic performance components.
    """

    def initialize(self):
        self.options.declare('surfaces', types=list)
        self.options.declare('user_specified_Sref', types=bool)

    def setup(self):
        surfaces = self.options['surfaces']

        if not self.options['user_specified_Sref']:
            self.add_subsystem('sum_areas',
                SumAreas(surfaces=surfaces),
                promotes_inputs=['*S_ref'],
                promotes_outputs=['S_ref_total'])

        self.add_subsystem('CL_CD',
                 TotalLiftDrag(surfaces=surfaces),
                 promotes_inputs=['*CL', '*CD', '*S_ref', 'S_ref_total'],
                 promotes_outputs=['CL', 'CD'])

        self.add_subsystem('moment',
                 MomentCoefficient(surfaces=surfaces),
                 promotes_inputs=['v', 'cg', 'rho', '*S_ref', '*b_pts', '*widths', '*chords', '*sec_forces', 'S_ref_total'],
                 promotes_outputs=['CM'])
