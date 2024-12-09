import numpy as np

import openmdao.api as om
from openaerostruct.utils.check_surface_dict import check_surface_dict_keys
import openaerostruct.geometry.geometry_mesh_gen as meshGen
from openaerostruct.geometry.geometry_unification import GeomMultiUnification
from openaerostruct.geometry.geometry_multi_join import GeomMultiJoin
from openaerostruct.utils.interpolation import get_normalized_span_coords


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

        # key validation of the surface dict
        check_surface_dict_keys(surface)

        # Get the surface name and create a group to contain components
        # only for this surface

        if self.options["DVGeo"]:
            from openaerostruct.geometry.ffd_component import GeometryMesh

            self.set_input_defaults("shape", val=np.zeros((surface["mx"], surface["my"])), units="m")

            if "t_over_c_cp" in surface.keys():
                n_cp = len(surface["t_over_c_cp"])
                # Add bspline components for active bspline geometric variables.
                x_interp = get_normalized_span_coords(surface, mid_panel=True)
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
                x_interp = get_normalized_span_coords(surface)
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
                x_interp = get_normalized_span_coords(surface)
                comp = self.add_subsystem(
                    "chord_bsp",
                    om.SplineComp(
                        method="bsplines", x_interp_val=x_interp, num_cp=n_cp, interp_options={"order": min(n_cp, 4)}
                    ),
                    promotes_inputs=["chord_cp"],
                    promotes_outputs=["chord"],
                )
                comp.add_spline(y_cp_name="chord_cp", y_interp_name="chord", y_units=None)
                bsp_inputs.append("chord")
                if surface.get("chord_cp_dv", True):
                    self.set_input_defaults("chord_cp", val=surface["chord_cp"], units=None)

            if "t_over_c_cp" in surface.keys():
                n_cp = len(surface["t_over_c_cp"])
                # Add bspline components for active bspline geometric variables.
                x_interp = get_normalized_span_coords(surface, mid_panel=True)
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
                x_interp = get_normalized_span_coords(surface)
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
                x_interp = get_normalized_span_coords(surface)
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
                x_interp = get_normalized_span_coords(surface)
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


# Function that constructs the individual section surface data dictionaries
def build_sections(surface):
    """This function returns an OpenMDAO Independent Variable Component with an output vector appropriately
    named and sized to function as an unified B-spline that joins multiple sections by construction.

    Parameters
    ----------
    surface: dict
        OpenAeroStruct multi-section surface dictionary

    Returns
    -------
    section_surfaces : list
        List of OpenAeroStruct surface dictionaries for each individual surface

    """
    # Get number of sections
    num_sections = surface["num_sections"]

    if surface["meshes"] == "gen-meshes":
        # Verify that all required inputs for automatic mesh generation are provided for each section
        if len(surface["ny"]) != num_sections:
            raise ValueError("Number of spanwise points needs to be provided for each section")
        if len(surface["taper"]) != num_sections:
            raise ValueError("Taper needs to be provided for each section")
        if len(surface["span"]) != num_sections:
            raise ValueError("Span needs to be provided for each section")
        if len(surface["sweep"]) != num_sections:
            raise ValueError("Sweep needs to be provided for each section")

        # Generate unified and individual section meshes
        mesh, sec_meshes = meshGen.generate_mesh(surface)
    else:
        # Allow user to provide mesh for each section
        if len(surface["meshes"]) != num_sections:
            raise ValueError("A mesh needs to be provided for each section.")
        sec_meshes = surface["meshes"]

    if len(surface["sec_name"]) != num_sections:
        raise ValueError("A name needs to be provided for each section.")

    # List of support keys for multi-section wings
    # NOTE: make sure this is consistent to the documentation's surface dict page
    target_keys = [
        # Essential Info
        "num_section",
        "symmetry",
        "S_ref_type",
        "ref_axis_pos",
        # wing definition
        "span",
        "taper",
        "sweep",
        "dihedral",
        "twist_cp",
        "chord_cp",
        "xshear_cp",
        "yshear_cp",
        "zshear_cp",
        # aerodynamics
        "CL0",
        "CD0",
        "with_viscous",
        "with_wave",
        "groundplane",
        "k_lam",
        "t_over_c_cp",
        "c_max_t",
    ]

    # Constructs a list of section dictionaries and adds the specified supported keys and values from the mult-section surface dictionary.
    surface_sections = []
    num_sections = surface["num_sections"]

    for i in range(num_sections):
        section = {}
        for k in set(surface).intersection(target_keys):
            if type(surface[k]) is list:
                section[k] = surface[k][i]
            else:
                section[k] = surface[k]
        section["mesh"] = sec_meshes[i]
        section["name"] = surface["sec_name"][i]
        surface_sections.append(section)
    return surface_sections


class MultiSecGeometry(om.Group):
    """
    Group that contains the section geometery groups for the multi-section surface


    This group handles the creation of each section geometry group based on parameters
    supplied in the multi-section surface dictionary. Meshes for each section can be
    provided by the user or automatically generated based on parameters supplied in the
    surface dictionary. The group also adds a mesh unification component that combines the
    individual section for each mesh into a singular unified mesh for use in aero components.
    Optionally, the joining component can be added that computes the edge distances between sections.
    This information can be used to set a distance constraint along the specified axes if needed.
    """

    def initialize(self):
        self.options.declare("surface", types=dict)  # Multi-section surface dictionary
        self.options.declare(
            "joining_comp", types=bool, default=False
        )  # Specify if a distance computation component should be added
        self.options.declare(
            "dim_constr", types=list, default=[]
        )  # List of arrays corresponding to each shared edge between section along the surface. Each array inidicates along which axes the distance constarint is applied([x y z])

    def setup(self):
        surface = self.options["surface"]
        joining_comp = self.options["joining_comp"]
        dc = self.options["dim_constr"]

        # key validation of the surface dict
        check_surface_dict_keys(surface)

        sec_dicts = build_sections(surface)

        section_names = []
        for sec in sec_dicts:
            geom_group = Geometry(surface=sec)
            self.add_subsystem(sec["name"], geom_group)
            section_names.append(sec["name"])

        # Add the mesh unification component
        unification_name = "{}_unification".format(surface["name"])

        uni_mesh = GeomMultiUnification(sections=sec_dicts, surface_name=surface["name"])
        self.add_subsystem(unification_name, uni_mesh)

        # Connect each section mesh to mesh unification component inputs
        for sec_name in section_names:
            self.connect("{}.mesh".format(sec_name), "{}.{}_def_mesh".format(unification_name, sec_name))

        # Connect each section t over c B-spline to t over c unification component if needed
        if "t_over_c_cp" in surface.keys():
            for sec_name in section_names:
                self.connect("{}.t_over_c".format(sec_name), "{}.{}_t_over_c".format(unification_name, sec_name))

        if joining_comp:
            # Add section joining component to output edge distances
            joining_name = "{}_joining".format(surface["name"])

            join = GeomMultiJoin(sections=sec_dicts, dim_constr=dc)
            self.add_subsystem(joining_name, join)

            for sec_name in section_names:
                self.connect("{}.mesh".format(sec_name), "{}.{}_join_mesh".format(joining_name, sec_name))
