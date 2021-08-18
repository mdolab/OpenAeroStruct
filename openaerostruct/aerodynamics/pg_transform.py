import openmdao.api as om

from openaerostruct.aerodynamics.pg_wind_rotation import RotateToWindFrame, RotateFromWindFrame
from openaerostruct.aerodynamics.pg_scale import ScaleToPrandtlGlauert, ScaleFromPrandtlGlauert


class PGTransform(om.Group):
    """
    This group is responsible for transforming the VLM geometries from
    physical coordinates to Prandtl-Glauert coordinates. This allows the
    compressible aerodynamic problem to be solved as an equivelent
    incompressible problem.

    The transform can be broken down into two steps:

        1. Rotate the geometry from the body frame to the wind frame so
        x-axis is parallel to freestream velocity.

        2. Scale wind frame coordinates by Prandtl-Glauert factor to retrieve
        equivalent Prandtl-Glauert geometry.
    """

    def initialize(self):
        self.options.declare("surfaces", types=list)
        self.options.declare(
            "rotational", False, types=bool, desc="Set to True to turn on support for computing angular velocities"
        )

    def setup(self):
        surfaces = self.options["surfaces"]
        rotational = self.options["rotational"]

        # Create component to rotate mesh geometry to wind frame,
        # s.t. the freestream velocity is along the x-axis
        self.add_subsystem(
            "rotate",
            RotateToWindFrame(surfaces=surfaces, rotational=rotational),
            promotes_inputs=["*"],
            promotes_outputs=["*"],
        )

        # Scale y and z direction by Prandtl Glauert factor
        self.add_subsystem(
            "scale",
            ScaleToPrandtlGlauert(surfaces=surfaces, rotational=rotational),
            promotes_inputs=["*"],
            promotes_outputs=["*"],
        )


class InversePGTransform(om.Group):
    """
    This group is responsible for transforming the solved incompressible
    forces in the Prandtl-Glauert domain to the compressible forces in the
    physical aerodynamic frame. This is the reverse procedure used in the
    PGTransform group.

    The inverse transform can be broken down into two steps:

        1. Scale Prandtl-Glauert force vectors to physical force vectors in wind
        frame.

        2. Rotate physical forces from wind frame to the aerodynamic frame.
    """

    def initialize(self):
        self.options.declare("surfaces", types=list)

    def setup(self):
        surfaces = self.options["surfaces"]

        # Scale force vectors back to compressible values
        scale_sys = ScaleFromPrandtlGlauert(surfaces=surfaces)
        self.add_subsystem("scale", scale_sys, promotes_inputs=["*"], promotes_outputs=["*"])

        # Rotate forces back to aerodynamic frame
        rot_sys = RotateFromWindFrame(surfaces=surfaces)
        self.add_subsystem("rotate", rot_sys, promotes_inputs=["*"], promotes_outputs=["*"])
