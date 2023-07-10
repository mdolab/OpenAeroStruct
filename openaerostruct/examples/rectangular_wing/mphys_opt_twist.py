"""
Use Mphys/OpenVSP interface to perform inviscid drag minimization of
a rectangular wing with respect to spanwise twist, subject to a lift constraint.
The expected result from lifting line theory should produce an elliptical lift distrbution.
Check output directory for Tecplot solution files.
"""

import os
import openmdao.api as om
import numpy as np

from mphys import Multipoint
from mphys.scenario_aerodynamic import ScenarioAerodynamic
from pygeo.mphys import OM_DVGEOCOMP

from openaerostruct.geometry.utils import generate_vsp_surfaces
from openaerostruct.mphys import AeroBuilder


class Top(Multipoint):
    def setup(self):
        # VSP model
        vsp_file = os.path.join(os.path.dirname(__file__), "rect_wing.vsp3")

        # Generate half-body mesh of rectangular wing
        surfaces = generate_vsp_surfaces(vsp_file, symmetry=True, include=["WingGeom"])

        # Define input surface dictionary for our wing
        surf_options = {
            "type": "aero",
            "S_ref_type": "wetted",  # how we compute the wing area,
            # can be 'wetted' or 'projected'
            "twist_cp": np.zeros(3),  # Define twist using 3 B-spline cp's
            "ref_axis_pos": 0.25,  # Define the reference axis position. 0 is the leading edge, 1 is the trailing edge.
            # Aerodynamic performance of the lifting surface at
            # an angle of attack of 0 (alpha=0).
            # These CL0 and CD0 values are added to the CL and CD
            # obtained from aerodynamic analysis of the surface to get
            # the total CL and CD.
            # These CL0 and CD0 values do not vary wrt alpha.
            "CL0": 0.0,  # CL of the surface at alpha=0
            "CD0": 0.0,  # CD of the surface at alpha=0
            # Airfoil properties for viscous drag calculation
            "k_lam": 0.05,  # percentage of chord with laminar
            # flow, used for viscous drag
            "t_over_c": 0.12,  # thickness over chord ratio (NACA0015)
            "c_max_t": 0.303,  # chordwise location of maximum (NACA0015)
            # thickness
            "with_viscous": False,  # if true, compute viscous drag,
            "with_wave": False,
        }  # end of surface dictionary

        # Update each surface with default options
        for surface in surfaces:
            surface.update(surf_options)

        mach = 0.0
        aoa = 5.0
        beta = 0.0
        rho = 0.38
        vel = 248.136
        re = 1e6

        dvs = self.add_subsystem("dvs", om.IndepVarComp(), promotes=["*"])
        dvs.add_output("aoa", val=aoa, units="deg")
        dvs.add_output("yaw", val=beta, units="deg")
        dvs.add_output("rho", val=rho, units="kg/m**3")
        dvs.add_output("mach", mach)
        dvs.add_output("v", vel, units="m/s")
        dvs.add_output("reynolds", re, units="1/m")
        dvs.add_output("cg", val=np.zeros((3)), units="m")

        # Create mphys builder for aero solver
        aero_builder = AeroBuilder(surfaces)
        aero_builder.initialize(self.comm)

        # Create mesh component and connect with solver
        self.add_subsystem("mesh", aero_builder.get_mesh_coordinate_subsystem())

        # add the geometry component, we dont need a builder because we do it here.
        self.add_subsystem("geometry", OM_DVGEOCOMP(file=vsp_file, type="vsp"))
        # add pointset
        self.geometry.nom_add_discipline_coords("aero")

        self.mphys_add_scenario("cruise", ScenarioAerodynamic(aero_builder=aero_builder))
        self.connect("mesh.x_aero0", "geometry.x_aero_in")
        self.connect("geometry.x_aero0", "cruise.x_aero")

        # Connect dv ivc's to solver
        for dv in ["aoa", "yaw", "rho", "mach", "v", "reynolds", "cg"]:
            self.connect(dv, f"cruise.{dv}")

    def configure(self):
        # create geometric DV setup
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_0", "Twist", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_1", "Twist", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_2", "Twist", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_3", "Twist", scaledStep=False)
        self.geometry.nom_addVSPVariable("WingGeom", "XSec_4", "Twist", scaledStep=False)


prob = om.Problem()
prob.model = Top()

# Set optimizer as model driver
prob.driver = om.ScipyOptimizeDriver()
prob.driver.options["debug_print"] = ["nl_cons", "objs", "desvars"]

# Setup problem and add design variables, constraint, and objective
prob.model.add_design_var("geometry.WingGeom:XSec_0:Twist", lower=-10.0, upper=10.0)
prob.model.add_design_var("geometry.WingGeom:XSec_1:Twist", lower=-10.0, upper=10.0)
prob.model.add_design_var("geometry.WingGeom:XSec_2:Twist", lower=-10.0, upper=10.0)
prob.model.add_design_var("geometry.WingGeom:XSec_3:Twist", lower=-10.0, upper=10.0)
prob.model.add_design_var("geometry.WingGeom:XSec_4:Twist", lower=-10.0, upper=10.0)
prob.model.add_constraint("cruise.WingGeom.CL", equals=0.5)
prob.model.add_objective("cruise.WingGeom.CD", scaler=1e4)

prob.setup()

# Create a n^2 diagram for user to view model connections
# om.n2(prob)
prob.run_model()
prob.run_driver()

print("CL", prob["cruise.WingGeom.CL"][0])
print("CD", prob["cruise.WingGeom.CD"][0])
