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
        # FFD
        "mx",
        "my",
    ]

    for key in surface.keys():
        if key not in keys_implemented:
            warnings.warn(
                "Key `{}` in surface dict is (likely) not supported in OAS and will be ignored".format(key),
                category=RuntimeWarning,
                stacklevel=2,
            )
