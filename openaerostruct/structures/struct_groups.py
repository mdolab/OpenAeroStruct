import openmdao.api as om
from openaerostruct.geometry.geometry_group import Geometry
from openaerostruct.structures.spatial_beam_states import SpatialBeamStates
from openaerostruct.structures.spatial_beam_functionals import SpatialBeamFunctionals
from openaerostruct.structures.spatial_beam_setup import SpatialBeamSetup
from openaerostruct.structures.tube_group import TubeGroup
from openaerostruct.structures.wingbox_group import WingboxGroup


class SpatialBeamAlone(om.Group):
    """Group that contains everything needed for a structural-only problem."""

    def initialize(self):
        self.options.declare("surface", types=dict)

    def setup(self):
        surface = self.options["surface"]

        self.add_subsystem(
            "geometry", Geometry(surface=surface), promotes_inputs=[], promotes_outputs=["mesh", "t_over_c"]
        )

        if surface["fem_model_type"] == "tube":
            tube_promotes_input = []
            tube_promotes_output = ["A", "Iy", "Iz", "J", "radius", "thickness"]
            if "thickness_cp" in surface.keys():
                tube_promotes_input.append("thickness_cp")
            if "radius_cp" not in surface.keys():
                tube_promotes_input = tube_promotes_input + ["mesh", "t_over_c"]

            self.add_subsystem(
                "tube_group",
                TubeGroup(surface=surface),
                promotes_inputs=tube_promotes_input,
                promotes_outputs=tube_promotes_output,
            )
        elif surface["fem_model_type"] == "wingbox":
            wingbox_promotes_in = ["mesh", "t_over_c"]
            wingbox_promotes_out = ["A", "Iy", "Iz", "J", "Qz", "A_enc", "A_int", "htop", "hbottom", "hfront", "hrear"]
            if "skin_thickness_cp" in surface.keys() and "spar_thickness_cp" in surface.keys():
                wingbox_promotes_in.append("skin_thickness_cp")
                wingbox_promotes_in.append("spar_thickness_cp")
                wingbox_promotes_out.append("skin_thickness")
                wingbox_promotes_out.append("spar_thickness")
            elif "skin_thickness_cp" in surface.keys() or "spar_thickness_cp" in surface.keys():
                raise NameError("Please have both skin and spar thickness as design variables, not one or the other.")

            self.add_subsystem(
                "wingbox_group",
                WingboxGroup(surface=surface),
                promotes_inputs=wingbox_promotes_in,
                promotes_outputs=wingbox_promotes_out,
            )
        else:
            raise NameError("Please select a valid `fem_model_type` from either `tube` or `wingbox`.")

        if surface["fem_model_type"] == "tube":
            self.add_subsystem(
                "struct_setup",
                SpatialBeamSetup(surface=surface),
                promotes_inputs=["mesh", "A", "Iy", "Iz", "J"],
                promotes_outputs=["nodes", "local_stiff_transformed", "structural_mass", "cg_location", "element_mass"],
            )
        else:
            self.add_subsystem(
                "struct_setup",
                SpatialBeamSetup(surface=surface),
                promotes_inputs=["mesh", "A", "Iy", "Iz", "J", "A_int"],
                promotes_outputs=[
                    "nodes",
                    "local_stiff_transformed",
                    "structural_mass",
                    "cg_location",
                    "element_mass",
                ],
            )

        promotes = []
        if surface["struct_weight_relief"]:
            promotes = promotes + list(set(["nodes", "element_mass", "load_factor"]))
        if surface["distributed_fuel_weight"]:
            promotes = promotes + list(set(["nodes", "load_factor"]))
        if "n_point_masses" in surface.keys():
            promotes = promotes + list(
                set(["point_mass_locations", "point_masses", "nodes", "load_factor", "engine_thrusts"])
            )

        self.add_subsystem(
            "struct_states",
            SpatialBeamStates(surface=surface),
            promotes_inputs=["local_stiff_transformed", "forces", "loads"] + promotes,
            promotes_outputs=["disp"],
        )

        if surface["fem_model_type"] == "tube":
            self.add_subsystem(
                "struct_funcs",
                SpatialBeamFunctionals(surface=surface),
                promotes_inputs=["thickness", "radius", "nodes", "disp"],
                promotes_outputs=["thickness_intersects", "vonmises", "failure"],
            )
        else:
            self.add_subsystem(
                "struct_funcs",
                SpatialBeamFunctionals(surface=surface),
                promotes_inputs=[
                    "spar_thickness",
                    "disp",
                    "Qz",
                    "J",
                    "A_enc",
                    "htop",
                    "hbottom",
                    "hfront",
                    "hrear",
                    "nodes",
                ],
                promotes_outputs=["vonmises", "failure"],
            )
