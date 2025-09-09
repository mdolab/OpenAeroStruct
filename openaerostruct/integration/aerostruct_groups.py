from openaerostruct.aerodynamics.geometry import VLMGeometry
from openaerostruct.geometry.geometry_group import Geometry
from openaerostruct.geometry.utils import build_section_dicts
from openaerostruct.geometry.geometry_unification import GeomMultiUnification
from openaerostruct.geometry.geometry_multi_join import GeomMultiJoin
from openaerostruct.transfer.displacement_transfer_group import DisplacementTransferGroup
from openaerostruct.structures.spatial_beam_setup import SpatialBeamSetup
from openaerostruct.structures.spatial_beam_states import SpatialBeamStates
from openaerostruct.aerodynamics.functionals import VLMFunctionals
from openaerostruct.structures.spatial_beam_functionals import SpatialBeamFunctionals
from openaerostruct.functionals.total_performance import TotalPerformance
from openaerostruct.transfer.load_transfer import LoadTransfer
from openaerostruct.aerodynamics.states import VLMStates
from openaerostruct.aerodynamics.compressible_states import CompressibleVLMStates
from openaerostruct.structures.tube_group import TubeGroup
from openaerostruct.structures.wingbox_group import WingboxGroup
from openaerostruct.utils.check_surface_dict import check_surface_dict_keys
from openaerostruct.structures.utils import compute_composite_stiffness
import openmdao.api as om


class AerostructGeometry(om.Group):
    def initialize(self):
        self.options.declare("surface", types=dict)
        self.options.declare("DVGeo", default=None)
        self.options.declare("connect_geom_DVs", default=True)

    def setup(self):
        surface = self.options["surface"]
        DVGeo = self.options["DVGeo"]
        connect_geom_DVs = self.options["connect_geom_DVs"]

        # key validation of the surface dict
        check_surface_dict_keys(surface)

        geom_promotes_in = []
        geom_promotes_out = ["mesh"]

        if connect_geom_DVs:
            # If connect_geom_DVs is true, then we promote all of the geometric design variables.
            # If it's false, then we do not promote them, which means that the geometry at each AeroStruct point is independent,
            # and the user can provide different values at each point.
            # This is useful when you want to have morphing DVs, such as twist or span, that are different at each point in a multipoint scheme.
            if "twist_cp" in surface.keys():
                geom_promotes_in.append("twist_cp")
            if "t_over_c_cp" in surface.keys():
                geom_promotes_out.append("t_over_c")
            if "sweep" in surface.keys():
                geom_promotes_in.append("sweep")
            if "taper" in surface.keys():
                geom_promotes_in.append("taper")
            if "mx" in surface.keys():
                geom_promotes_in.append("shape")

        self.add_subsystem(
            "geometry",
            Geometry(surface=surface, DVGeo=DVGeo),
            promotes_inputs=geom_promotes_in,
            promotes_outputs=geom_promotes_out,
        )

        if surface["fem_model_type"].lower() == "tube":
            tube_promotes_input = []
            tube_promotes_output = ["A", "Iy", "Iz", "J", "radius", "thickness"]
            if "thickness_cp" in surface.keys() and connect_geom_DVs:
                tube_promotes_input.append("thickness_cp")
            if "radius_cp" not in surface.keys():
                tube_promotes_input = tube_promotes_input + ["mesh", "t_over_c"]

            self.add_subsystem(
                "tube_group",
                TubeGroup(surface=surface),
                promotes_inputs=tube_promotes_input,
                promotes_outputs=tube_promotes_output,
            )
        elif (
            surface["fem_model_type"].lower() == "wingbox"
        ):  # connections and nomenclature remains the same for both isotropic and composite wingbox
            wingbox_promotes_in = ["mesh", "t_over_c"]
            wingbox_promotes_out = [
                "A",
                "Iy",
                "Iz",
                "J",
                "Qz",
                "A_enc",
                "A_int",
                "htop",
                "hbottom",
                "hfront",
                "hrear",
            ]
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

        if surface["fem_model_type"].lower() == "wingbox":  # same for both isotropic and composite wingbox
            promotes = ["A_int"]
        else:
            promotes = []

        self.add_subsystem(
            "struct_setup",
            SpatialBeamSetup(surface=surface),
            promotes_inputs=["mesh", "A", "Iy", "Iz", "J"] + promotes,
            promotes_outputs=["nodes", "local_stiff_transformed", "structural_mass", "cg_location", "element_mass"],
        )


class MultiSecAerostructGeometry(om.Group):
    """
    Group that contains the section geometery groups and structural geometry group
    for the multi-section aerostruct surface


    This group handles the creation of each section geometry group based on parameters
    supplied in the multi-section surface dictionary. The group also adds the appropriate
    structural group for aerostructural analysis. Meshes for each section can be
    provided by the user or automatically generated based on parameters supplied in the
    surface dictionary. The group also adds a mesh unification component that combines the
    individual section for each mesh into a singular unified mesh for use in aero components.
    Note that the structural nodes will based on the unified mesh. Optionally, the joining component can be added that computes the edge distances between sections.
    This information can be used to set a distance constraint along the specified axes if needed.
    """

    def initialize(self):
        self.options.declare("surface", types=dict)  # Multi-section surface dictionary
        self.options.declare("connect_geom_DVs", default=True)
        self.options.declare(
            "joining_comp", types=bool, default=False
        )  # Specify if a distance computation component should be added
        self.options.declare(
            "dim_constr", types=list, default=[]
        )  # List of arrays corresponding to each shared edge between section along the surface. Each array inidicates along which axes the distance constarint is applied([x y z])
        self.options.declare("shift_uni_mesh", types=bool, default=True)  # Flag to apply mesh shifting or not

    def setup(self):
        surface = self.options["surface"]
        connect_geom_DVs = self.options["connect_geom_DVs"]
        joining_comp = self.options["joining_comp"]
        dc = self.options["dim_constr"]
        shift_uni_mesh = self.options["shift_uni_mesh"]

        # key validation of the surface dict
        check_surface_dict_keys(surface)

        """
        ### Multi-section surface geometry setup ###
        """

        sec_dicts = build_section_dicts(surface)

        section_names = []
        for sec in sec_dicts:
            geom_group = Geometry(surface=sec)
            self.add_subsystem(sec["name"], geom_group)
            section_names.append(sec["name"])

        # Add the mesh unification component
        unification_name = "{}_unification".format(surface["name"])

        promotes_outputs = [("{}_uni_mesh".format(surface["name"]), "mesh")]
        if "t_over_c_cp" in surface.keys():
            promotes_outputs += [("{}_uni_t_over_c".format(surface["name"]), "t_over_c")]

        uni_mesh = GeomMultiUnification(sections=sec_dicts, surface_name=surface["name"], shift_uni_mesh=shift_uni_mesh)
        self.add_subsystem(unification_name, uni_mesh, promotes_outputs=promotes_outputs)

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

        """
        ### Structural geometry setup ###
        """
        if surface["fem_model_type"] == "tube":
            tube_promotes_input = []
            tube_promotes_output = ["A", "Iy", "Iz", "J", "radius", "thickness"]
            if "thickness_cp" in surface.keys() and connect_geom_DVs:
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

        if surface["fem_model_type"] == "wingbox":
            promotes = ["A_int"]
        else:
            promotes = []

        self.add_subsystem(
            "struct_setup",
            SpatialBeamSetup(surface=surface),
            promotes_inputs=["mesh", "A", "Iy", "Iz", "J"] + promotes,
            promotes_outputs=["nodes", "local_stiff_transformed", "structural_mass", "cg_location", "element_mass"],
        )


class CoupledAS(om.Group):
    def initialize(self):
        self.options.declare("surface", types=dict)

    def setup(self):
        surface = self.options["surface"]

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

        self.add_subsystem(
            "def_mesh",
            DisplacementTransferGroup(surface=surface),
            promotes_inputs=["nodes", "mesh", "disp"],
            promotes_outputs=["def_mesh"],
        )

        self.add_subsystem(
            "aero_geom",
            VLMGeometry(surface=surface),
            promotes_inputs=["def_mesh"],
            promotes_outputs=["b_pts", "widths", "lengths_spanwise", "lengths", "chords", "normals", "S_ref"],
        )

        self.linear_solver = om.LinearRunOnce()


class CoupledPerformance(om.Group):
    def initialize(self):
        self.options.declare("surface", types=dict)

    def setup(self):
        surface = self.options["surface"]

        self.add_subsystem(
            "aero_funcs",
            VLMFunctionals(surface=surface),
            promotes_inputs=[
                "v",
                "alpha",
                "beta",
                "Mach_number",
                "re",
                "rho",
                "widths",
                "lengths_spanwise",
                "lengths",
                "S_ref",
                "sec_forces",
                "t_over_c",
            ],
            promotes_outputs=["CDv", "CDw", "L", "D", "CL1", "CDi", "CD", "CL", "Cl"],
        )

        if surface["fem_model_type"].lower() == "tube":
            self.add_subsystem(
                "struct_funcs",
                SpatialBeamFunctionals(surface=surface),
                promotes_inputs=["thickness", "radius", "nodes", "disp"],
                promotes_outputs=["thickness_intersects", "vonmises", "failure"],
            )

        elif surface["fem_model_type"].lower() == "wingbox":
            if "useComposite" in surface.keys() and surface["useComposite"]:  # using the Composite Wing Box
                promotedoutput = "tsaiwu_sr"
            else:  # using the isotropic Wing Box
                promotedoutput = "vonmises"

            self.add_subsystem(
                "struct_funcs",
                SpatialBeamFunctionals(surface=surface),
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
                promotes_outputs=[promotedoutput, "failure"],
            )
        else:
            raise NameError("Please select a valid `fem_model_type` from either `tube` or `wingbox`.")


class AerostructPoint(om.Group):
    def initialize(self):
        self.options.declare("surfaces", types=list)
        self.options.declare("user_specified_Sref", types=bool, default=False)
        self.options.declare("internally_connect_fuelburn", types=bool, default=True)
        self.options.declare(
            "compressible",
            types=bool,
            default=False,
            desc="Turns on compressibility correction for moderate Mach number " "flows. Defaults to False.",
        )
        self.options.declare(
            "rotational", False, types=bool, desc="Set to True to turn on support for computing angular velocities"
        )

    def setup(self):
        surfaces = self.options["surfaces"]
        rotational = self.options["rotational"]

        # Check for multi-section surfaces and create suitable surface dictionaries for them
        for i, surface in enumerate(surfaces):
            # If multisection mesh then build a single surface with the unified mesh data
            if "is_multi_section" in surface.keys():
                import copy

                target_keys = [
                    # Essential Info
                    "name",
                    "symmetry",
                    "S_ref_type",
                    "ref_axis_pos",
                    "mesh",
                    # aerodynamics
                    "CL0",
                    "CD0",
                    "with_viscous",
                    "with_wave",
                    "groundplane",
                    "k_lam",
                    "t_over_c_cp",
                    "c_max_t",
                    # structures
                    "fem_model_type",
                    "E",
                    "G",
                    "yield",
                    "mrho",
                    "fem_origin",
                    "wing_weight_ratio",
                    "exact_failure_constraint",
                    "struct_weight_relief",
                    "distributed_fuel_weight",
                    "fuel_density",
                    "Wf_reserve",
                    "n_point_masses",
                    # structural parameterization tube
                    "thickness_cp",
                    "radius_cp",
                    # structural parameterization wingbox
                    "spar_thickness_cp",
                    "skin_thickness_cp",
                    "original_wingbox_airfoil_t_over_c",
                    "strength_factor_for_upper_skin",
                    "data_x_upper",
                    "data_y_upper",
                    "data_x_lower",
                    "data_y_lower",
                ]

                # Constructs a surface dictionary and adds the specified supported keys and values from the mult-section surface dictionary.
                aeroStructSurface = {}
                for k in set(surface).intersection(target_keys):
                    aeroStructSurface[k] = surface[k]
                surfaces[i] = copy.deepcopy(aeroStructSurface)

        coupled = om.Group()

        for surface in surfaces:
            name = surface["name"]

            # if useComposite is enabled, compute the effective E and G values for the composite material
            useComposite = "useComposite" in surface.keys() and surface["useComposite"]
            if useComposite:
                compute_composite_stiffness(surface)

            # Connect the output of the loads component with the FEM
            # displacement parameter. This links the coupling within the coupled
            # group that necessitates the subgroup solver.
            coupled.connect(name + "_loads.loads", name + ".loads")

            # Perform the connections with the modified names within the
            # 'aero_states' group.
            coupled.connect(name + ".normals", "aero_states." + name + "_normals")
            coupled.connect(name + ".def_mesh", "aero_states." + name + "_def_mesh")

            # Connect the results from 'coupled' to the performance groups
            coupled.connect(name + ".def_mesh", name + "_loads.def_mesh")
            coupled.connect("aero_states." + name + "_sec_forces", name + "_loads.sec_forces")

            # Connect the results from 'aero_states' to the performance groups
            self.connect("coupled.aero_states." + name + "_sec_forces", name + "_perf" + ".sec_forces")

            # Connection performance functional variables
            self.connect(name + "_perf.CL", "total_perf." + name + "_CL")
            self.connect(name + "_perf.CD", "total_perf." + name + "_CD")
            self.connect("coupled.aero_states." + name + "_sec_forces", "total_perf." + name + "_sec_forces")
            self.connect("coupled." + name + ".chords", name + "_perf.aero_funcs.chords")

            # Connect parameters from the 'coupled' group to the performance
            # groups for the individual surfaces.
            self.connect("coupled." + name + ".disp", name + "_perf.disp")
            self.connect("coupled." + name + ".S_ref", name + "_perf.S_ref")
            self.connect("coupled." + name + ".widths", name + "_perf.widths")
            # self.connect('coupled.' + name + '.chords', name + '_perf.chords')
            self.connect("coupled." + name + ".lengths", name + "_perf.lengths")
            self.connect("coupled." + name + ".lengths_spanwise", name + "_perf.lengths_spanwise")

            # Connect parameters from the 'coupled' group to the total performance group.
            self.connect("coupled." + name + ".S_ref", "total_perf." + name + "_S_ref")
            self.connect("coupled." + name + ".widths", "total_perf." + name + "_widths")
            self.connect("coupled." + name + ".chords", "total_perf." + name + "_chords")
            self.connect("coupled." + name + ".b_pts", "total_perf." + name + "_b_pts")

            # Add components to the 'coupled' group for each surface.
            # The 'coupled' group must contain all components and parameters
            # needed to converge the aerostructural system.
            coupled_AS_group = CoupledAS(surface=surface)

            if (
                surface["distributed_fuel_weight"]
                or "n_point_masses" in surface.keys()
                or surface["struct_weight_relief"]
            ):
                prom_in = ["load_factor"]
            else:
                prom_in = []

            coupled.add_subsystem(name, coupled_AS_group, promotes_inputs=prom_in)

        # check for ground effect and if so, promote
        ground_effect = False
        for surface in surfaces:
            if surface.get("groundplane", False):
                ground_effect = True

        if self.options["compressible"] is True:
            aero_states = CompressibleVLMStates(surfaces=surfaces, rotational=rotational)
            prom_in = ["v", "alpha", "beta", "rho", "Mach_number"]
        else:
            aero_states = VLMStates(surfaces=surfaces, rotational=rotational)
            prom_in = ["v", "alpha", "beta", "rho"]
        if ground_effect:
            prom_in.append("height_agl")

        if rotational:
            prom_in.extend(["omega", ("cg", "empty_cg")])

        # Add a single 'aero_states' component for the whole system within the
        # coupled group.
        coupled.add_subsystem("aero_states", aero_states, promotes_inputs=prom_in)

        # Explicitly connect parameters from each surface's group and the common
        # 'aero_states' group.
        for surface in surfaces:
            name = surface["name"]

            # Add a loads component to the coupled group
            coupled.add_subsystem(name + "_loads", LoadTransfer(surface=surface))

        """
        ### Change the solver settings here ###
        """

        # Set solver properties for the coupled group
        # coupled.linear_solver = ScipyKrylov()
        # coupled.linear_solver.precon = om.LinearRunOnce()

        coupled.nonlinear_solver = om.NonlinearBlockGS(use_aitken=True)
        coupled.nonlinear_solver.options["maxiter"] = 100
        coupled.nonlinear_solver.options["atol"] = 1e-7
        coupled.nonlinear_solver.options["rtol"] = 1e-30
        coupled.nonlinear_solver.options["iprint"] = 2
        coupled.nonlinear_solver.options["err_on_non_converge"] = True

        # coupled.linear_solver = om.DirectSolver()

        coupled.linear_solver = om.DirectSolver(assemble_jac=True)
        coupled.options["assembled_jac_type"] = "csc"

        # coupled.nonlinear_solver = om.NewtonSolver(solve_subsystems=True)
        # coupled.nonlinear_solver.options['maxiter'] = 50

        """
        ### End change of solver settings ###
        """
        prom_in = ["v", "alpha", "beta", "rho"]
        if self.options["compressible"] is True:
            prom_in.append("Mach_number")
        if ground_effect:
            prom_in.append("height_agl")
        if rotational:
            prom_in.extend(["omega", "empty_cg"])

        # Add the coupled group to the model problem
        self.add_subsystem("coupled", coupled, promotes_inputs=prom_in)

        for surface in surfaces:
            name = surface["name"]

            # Add a performance group which evaluates the data after solving
            # the coupled system
            perf_group = CoupledPerformance(surface=surface)

            self.add_subsystem(
                name + "_perf", perf_group, promotes_inputs=["rho", "v", "alpha", "beta", "re", "Mach_number"]
            )

        # Add functionals to evaluate performance of the system.
        # Note that only the interesting results are promoted here; not all
        # of the parameters.
        self.add_subsystem(
            "total_perf",
            TotalPerformance(
                surfaces=surfaces,
                user_specified_Sref=self.options["user_specified_Sref"],
                internally_connect_fuelburn=self.options["internally_connect_fuelburn"],
            ),
            promotes_inputs=[
                "v",
                "rho",
                "empty_cg",
                "total_weight",
                "CT",
                "speed_of_sound",
                "R",
                "Mach_number",
                "W0",
                "load_factor",
                "S_ref_total",
            ],
            promotes_outputs=["L_equals_W", "fuelburn", "CL", "CD", "CM", "cg"],
        )
