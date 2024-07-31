import openmdao.api as om

# from openaerostruct.structures.energy import Energy
# from openaerostruct.structures.weight import Weight
# from openaerostruct.structures.spar_within_wing import SparWithinWing
from openaerostruct.structures.vonmises_tube import VonMisesTube
from openaerostruct.structures.vonmises_wingbox import VonMisesWingbox
from openaerostruct.structures.tsaiwu_wingbox import TsaiWuWingbox
from openaerostruct.structures.non_intersecting_thickness import NonIntersectingThickness
from openaerostruct.structures.failure_exact import FailureExact
from openaerostruct.structures.failure_ks import FailureKS


class SpatialBeamFunctionals(om.Group):
    """Group that contains the spatial beam functionals used to evaluate
    performance."""

    def initialize(self):
        self.options.declare("surface", types=dict)

    def setup(self):
        surface = self.options["surface"]

        # Commented out energy for now since we haven't ever used its output
        # self.add_subsystem('energy',
        #          Energy(surface=surface),
        #          promotes=['*'])

        if surface["fem_model_type"] == "tube":
            self.add_subsystem(
                "thicknessconstraint",
                NonIntersectingThickness(surface=surface),
                promotes_inputs=["thickness", "radius"],
                promotes_outputs=["thickness_intersects"],
            )

            self.add_subsystem(
                "vonmises",
                VonMisesTube(surface=surface),
                promotes_inputs=["radius", "nodes", "disp"],
                promotes_outputs=["vonmises"],
            )

        elif surface["fem_model_type"] == "wingbox":

            if "useComposite" in surface.keys() and surface["useComposite"]:  # using the Composite wingbox
                subsystemname = "tsaiwu_sr"
                promotedoutput = "tsaiwu_sr"
            else:  # using the Isotropic wingbox
                subsystemname = "vonmises"
                promotedoutput = "vonmises"

            self.add_subsystem(
                subsystemname,
                TsaiWuWingbox(surface=surface),
                promotes_inputs=[
                    "Qz",
                    "J",
                    "A_enc",
                    "spar_thickness",
                    "htop",
                    "hbottom",
                    "hfront",
                    "hrear",
                    "nodes",
                    "disp",
                ],
                promotes_outputs=[promotedoutput],
            )
        else:
            raise NameError("Please select a valid `fem_model_type` from either `tube` or `wingbox`.")

        # The following component has not been fully tested so we leave it
        # commented out for now. Use at your own risk.
        # self.add_subsystem('sparconstraint',
        #          SparWithinWing(surface=surface),
        #          promotes=['*'])

        if surface["exact_failure_constraint"]:
            if "useComposite" in surface.keys() and surface["useComposite"]:  # using the Composite wingbox
                self.add_subsystem(
                    "failure",
                    FailureExact(surface=surface),
                    promotes_inputs=["tsaiwu_sr"],
                    promotes_outputs=["failure"],
                )
            else:  # using the Isotropic structures
                self.add_subsystem(
                    "failure", FailureExact(surface=surface), promotes_inputs=["vonmises"], promotes_outputs=["failure"]
                )
        else:
            if "useComposite" in surface.keys() and surface["useComposite"]:  # using the Composite wingbox
                self.add_subsystem(
                    "failure", FailureKS(surface=surface), promotes_inputs=["tsaiwu_sr"], promotes_outputs=["failure"]
                )
            else:
                self.add_subsystem(
                    "failure", FailureKS(surface=surface), promotes_inputs=["vonmises"], promotes_outputs=["failure"]
                )
