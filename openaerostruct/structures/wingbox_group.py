import numpy as np
import openmdao.api as om
from openaerostruct.structures.section_properties_wingbox import SectionPropertiesWingbox
from openaerostruct.structures.wingbox_geometry import WingboxGeometry


class WingboxGroup(om.Group):
    """Group that contains everything needed for a structural-only problem."""

    def initialize(self):
        self.options.declare("surface", types=dict)

    def setup(self):
        surface = self.options["surface"]
        ny = surface["mesh"].shape[1]

        if "spar_thickness_cp" in surface.keys():
            n_cp = len(surface["spar_thickness_cp"])
            # Add bspline components for active bspline geometric variables.
            x_interp = np.linspace(0.0, 1.0, int(ny - 1))
            comp = self.add_subsystem(
                "spar_thickness_bsp",
                om.SplineComp(
                    method="bsplines", x_interp_val=x_interp, num_cp=n_cp, interp_options={"order": min(n_cp, 4)}
                ),
                promotes_inputs=["spar_thickness_cp"],
                promotes_outputs=["spar_thickness"],
            )
            comp.add_spline(y_cp_name="spar_thickness_cp", y_interp_name="spar_thickness", y_units="m")
            self.set_input_defaults("spar_thickness_cp", val=surface["spar_thickness_cp"], units="m")

        if "skin_thickness_cp" in surface.keys():
            n_cp = len(surface["skin_thickness_cp"])
            # Add bspline components for active bspline geometric variables.
            x_interp = np.linspace(0.0, 1.0, int(ny - 1))
            comp = self.add_subsystem(
                "skin_thickness_bsp",
                om.SplineComp(
                    method="bsplines", x_interp_val=x_interp, num_cp=n_cp, interp_options={"order": min(n_cp, 4)}
                ),
                promotes_inputs=["skin_thickness_cp"],
                promotes_outputs=["skin_thickness"],
            )
            comp.add_spline(y_cp_name="skin_thickness_cp", y_interp_name="skin_thickness", y_units="m")
            self.set_input_defaults("skin_thickness_cp", val=surface["skin_thickness_cp"], units="m")

        self.add_subsystem(
            "wingbox_geometry",
            WingboxGeometry(surface=surface),
            promotes_inputs=["mesh"],
            promotes_outputs=["fem_chords", "fem_twists", "streamwise_chords"],
        )

        self.add_subsystem(
            "wingbox",
            SectionPropertiesWingbox(surface=surface),
            promotes_inputs=[
                "spar_thickness",
                "skin_thickness",
                "t_over_c",
                "fem_chords",
                "fem_twists",
                "streamwise_chords",
            ],
            promotes_outputs=["A", "Iy", "Qz", "Iz", "J", "A_enc", "A_int", "htop", "hbottom", "hfront", "hrear"],
        )
