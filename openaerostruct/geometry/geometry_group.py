import numpy as np

import openmdao.api as om


class Geometry(om.Group):
    """
    Group that contains all components needed for any type of OAS problem.

    Because we use this general group, there's some logic to figure out which
    components to add and which connections to make.
    This is especially true for all of the geometric manipulation types, such
    as twist, sweep, etc., in that we handle the creation of these parameters
    differently if the user wants to have them vary in the optimization problem.
    """

    def initialize(self):
        self.options.declare("surface", types=dict)
        self.options.declare("DVGeo", default=None)
        self.options.declare("connect_geom_DVs", default=True)
        # The option "connect_geom_DVs" is no longer necessary, but we still keep it to be backward compatible.

    def setup(self):
        surface = self.options["surface"]

        # Get the surface name and create a group to contain components
        # only for this surface
        ny = surface["mesh"].shape[1]

        if self.options["DVGeo"]:
            from openaerostruct.geometry.ffd_component import GeometryMesh

            self.set_input_defaults("shape", val=np.zeros((surface["mx"], surface["my"])), units="m")

            if "t_over_c_cp" in surface.keys():
                n_cp = len(surface["t_over_c_cp"])
                # Add bspline components for active bspline geometric variables.
                x_interp = np.linspace(0.0, 1.0, int(ny - 1))
                comp = self.add_subsystem(
                    "t_over_c_bsp",
                    om.SplineComp(
                        method="bsplines", x_interp_val=x_interp, num_cp=n_cp, interp_options={"order": min(n_cp, 4)}
                    ),
                    promotes_inputs=["t_over_c_cp"],
                    promotes_outputs=["t_over_c"],
                )
                comp.add_spline(y_cp_name="t_over_c_cp", y_interp_name="t_over_c")
                if surface.get("t_over_c_cp_dv", True):
                    self.set_input_defaults("t_over_c_cp", val=surface["t_over_c_cp"])

            self.add_subsystem(
                "mesh",
                GeometryMesh(surface=surface, DVGeo=self.options["DVGeo"]),
                promotes_inputs=["shape"],
                promotes_outputs=["mesh"],
            )

        else:
            from openaerostruct.geometry.geometry_mesh import GeometryMesh

            bsp_inputs = []

            if "twist_cp" in surface.keys():
                n_cp = len(surface["twist_cp"])
                # Add bspline components for active bspline geometric variables.
                x_interp = np.linspace(0.0, 1.0, int(ny))
                comp = self.add_subsystem(
                    "twist_bsp",
                    om.SplineComp(
                        method="bsplines", x_interp_val=x_interp, num_cp=n_cp, interp_options={"order": min(n_cp, 4)}
                    ),
                    promotes_inputs=["twist_cp"],
                    promotes_outputs=["twist"],
                )
                comp.add_spline(y_cp_name="twist_cp", y_interp_name="twist", y_units="deg")
                bsp_inputs.append("twist")

                # Since default assumption is that we want tail rotation as a design variable, add this to allow for trimmed drag polar where the tail rotation should not be a design variable
                if surface.get("twist_cp_dv", True):
                    self.set_input_defaults("twist_cp", val=surface["twist_cp"], units="deg")

            if "chord_cp" in surface.keys():
                n_cp = len(surface["chord_cp"])
                # Add bspline components for active bspline geometric variables.
                x_interp = np.linspace(0.0, 1.0, int(ny))
                comp = self.add_subsystem(
                    "chord_bsp",
                    om.SplineComp(
                        method="bsplines", x_interp_val=x_interp, num_cp=n_cp, interp_options={"order": min(n_cp, 4)}
                    ),
                    promotes_inputs=["chord_cp"],
                    promotes_outputs=["chord"],
                )
                comp.add_spline(y_cp_name="chord_cp", y_interp_name="chord", y_units="m")
                bsp_inputs.append("chord")
                if surface.get("chord_cp_dv", True):
                    self.set_input_defaults("chord_cp", val=surface["chord_cp"], units="m")

            if "t_over_c_cp" in surface.keys():
                n_cp = len(surface["t_over_c_cp"])
                # Add bspline components for active bspline geometric variables.
                x_interp = np.linspace(0.0, 1.0, int(ny - 1))
                comp = self.add_subsystem(
                    "t_over_c_bsp",
                    om.SplineComp(
                        method="bsplines", x_interp_val=x_interp, num_cp=n_cp, interp_options={"order": min(n_cp, 4)}
                    ),
                    promotes_inputs=["t_over_c_cp"],
                    promotes_outputs=["t_over_c"],
                )
                comp.add_spline(y_cp_name="t_over_c_cp", y_interp_name="t_over_c")
                if surface.get("t_over_c_cp_dv", True):
                    self.set_input_defaults("t_over_c_cp", val=surface["t_over_c_cp"])

            if "xshear_cp" in surface.keys():
                n_cp = len(surface["xshear_cp"])
                # Add bspline components for active bspline geometric variables.
                x_interp = np.linspace(0.0, 1.0, int(ny))
                comp = self.add_subsystem(
                    "xshear_bsp",
                    om.SplineComp(
                        method="bsplines", x_interp_val=x_interp, num_cp=n_cp, interp_options={"order": min(n_cp, 4)}
                    ),
                    promotes_inputs=["xshear_cp"],
                    promotes_outputs=["xshear"],
                )
                comp.add_spline(y_cp_name="xshear_cp", y_interp_name="xshear", y_units="m")
                bsp_inputs.append("xshear")
                if surface.get("xshear_cp_dv", True):
                    self.set_input_defaults("xshear_cp", val=surface["xshear_cp"], units="m")

            if "yshear_cp" in surface.keys():
                n_cp = len(surface["yshear_cp"])
                # Add bspline components for active bspline geometric variables.
                x_interp = np.linspace(0.0, 1.0, int(ny))
                comp = self.add_subsystem(
                    "yshear_bsp",
                    om.SplineComp(
                        method="bsplines", x_interp_val=x_interp, num_cp=n_cp, interp_options={"order": min(n_cp, 4)}
                    ),
                    promotes_inputs=["yshear_cp"],
                    promotes_outputs=["yshear"],
                )
                comp.add_spline(y_cp_name="yshear_cp", y_interp_name="yshear", y_units="m")
                bsp_inputs.append("yshear")
                if surface.get("yshear_cp_dv", True):
                    self.set_input_defaults("yshear_cp", val=surface["yshear_cp"], units="m")

            if "zshear_cp" in surface.keys():
                n_cp = len(surface["zshear_cp"])
                # Add bspline components for active bspline geometric variables.
                x_interp = np.linspace(0.0, 1.0, int(ny))
                comp = self.add_subsystem(
                    "zshear_bsp",
                    om.SplineComp(
                        method="bsplines", x_interp_val=x_interp, num_cp=n_cp, interp_options={"order": min(n_cp, 4)}
                    ),
                    promotes_inputs=["zshear_cp"],
                    promotes_outputs=["zshear"],
                )
                comp.add_spline(y_cp_name="zshear_cp", y_interp_name="zshear", y_units="m")
                bsp_inputs.append("zshear")
                if surface.get("zshear_cp_dv", True):
                    self.set_input_defaults("zshear_cp", val=surface["zshear_cp"], units="m")

            if "sweep" in surface.keys():
                bsp_inputs.append("sweep")
                if surface.get("sweep_dv", True):
                    self.set_input_defaults("sweep", val=surface["sweep"], units="deg")

            if "span" in surface.keys():
                bsp_inputs.append("span")
                if surface.get("span_dv", True):
                    self.set_input_defaults("span", val=surface["span"], units="m")

            if "dihedral" in surface.keys():
                bsp_inputs.append("dihedral")
                if surface.get("dihedral_dv", True):
                    self.set_input_defaults("dihedral", val=surface["dihedral"], units="deg")

            if "taper" in surface.keys():
                bsp_inputs.append("taper")
                if surface.get("taper_dv", True):
                    self.set_input_defaults("taper", val=surface["taper"])

            self.add_subsystem(
                "mesh", GeometryMesh(surface=surface), promotes_inputs=bsp_inputs, promotes_outputs=["mesh"]
            )
