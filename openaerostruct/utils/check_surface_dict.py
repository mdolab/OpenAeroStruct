import warnings


def check_surface_dict_keys(surface):
    """
    Key valication function for the OAS surface dict.
    Shows a warning if a user provided a key that is (likely) not implemented in OAS.

    Parameters
    ----------
    surface : dict
        User-defined surface dict
    """

    # NOTE: make sure this is consistent to the documentation's surface dict page
    keys_implemented = [
        # wing definition
        "name",
        "symmetry",
        "S_ref_type",
        "mesh",
        "span",
        "taper",
        "sweep",
        "dihedral",
        "twist_cp",
        "chord_cp",
        "xshear_cp",
        "yshear_cp",
        "zshear_cp",
        "ref_axis_pos",
        # aerodynamics
        "CL0",
        "CD0",
        "with_viscous",
        "with_wave",
        "groundplane",
        "k_lam",
        "t_over_c_cp",
        "c_max_t",
        # structure
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
        # tube structure
        "thickness_cp",
        "radius_cp",
        # wingbox structure
        "spar_thickness_cp",
        "skin_thickness_cp",
        "original_wingbox_airfoil_t_over_c",
        "strength_factor_for_upper_skin",
        "data_x_upper",
        "data_y_upper",
        "data_x_lower",
        "data_y_lower",
        # tsaiwu_wingbox structure
        "useComposite",
        # FFD
        "mx",
        "my",
    ]
    # keys that are required when useComposite is True
    compositeInputs = [
        "composite_safetyfactor",
        "plyangles",
        "plyfractions",
        "E1",
        "E2",
        "nu12",
        "G12",
        "sigma_t1",
        "sigma_c1",
        "sigma_t2",
        "sigma_c2",
        "sigma_12max",
    ]

    keys_implemented += compositeInputs

    for key in surface.keys():
        if key not in keys_implemented:
            warnings.warn(
                "Key `{}` in surface dict is (likely) not supported in OAS and will be ignored".format(key),
                category=RuntimeWarning,
                stacklevel=2,
            )

    # adding checks for using the composite failure model
    # check1: if useComposite is True, then the following keys must be present
    if surface.get("useComposite", False):
        for key in compositeInputs:
            if key not in surface.keys():
                raise ValueError(
                    "Key `{}` is required when `useComposite` is True".format(key),
                )

    # check2: if useComposite is True, then 'fem_model_type' must be 'wingbox'
    if surface.get("useComposite", False) and surface.get("fem_model_type", "") != "wingbox":
        raise ValueError(
            "`fem_model_type` must be 'wingbox' when `useComposite` is True",
        )

    # check3: if useComposite is True, then length of plyangles and plyfractions must be equal
    if surface.get("useComposite", False):
        if len(surface["plyangles"]) != len(surface["plyfractions"]):
            raise ValueError(
                "Length of `plyangles` and `plyfractions` must be equal",
            )
