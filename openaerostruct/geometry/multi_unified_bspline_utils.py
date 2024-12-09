import numpy as np
import openmdao.api as om


def build_multi_spline(out_name, num_sections, control_points):
    """This function returns an OpenMDAO Independent Variable Component with an output vector appropriately
    named and sized to function as an unified B-spline that joins multiple sections by construction.

    Parameters
    ----------
    out_name: string
        Name of the output to assign to the B-spline
    num_sections : int
        Number of sections
    control_points: list
        List of B-spline control point arrays corresponding to each section

    Returns
    -------
    spline_control : OpenMDAO component object
        The unified B-spline indpendent variable component

    """
    if len(control_points) != num_sections:
        raise Exception("Target sections need to match with control points!")

    single_sections = len([cp for cp in control_points if len(cp) == 1])

    control_poin_vec = np.ones(len(np.concatenate(control_points)) - (num_sections - 1 - single_sections))

    spline_control = om.IndepVarComp()
    spline_control.add_output("{}_spline".format(out_name), val=control_poin_vec)

    return spline_control


def connect_multi_spline(prob, section_surfaces, sec_cp, out_name, comp_name, return_bind_inds=False):
    """This function connects the the unified B-spline component with the individual B-splines
    of each section. There is a point of overlap at each section so that each edge control point control the edge
    controls points of each section's B-spline. This is how section joining by consturction is acheived.
    An issue occurs however when a B-spline in a particular section only has one control point. In this case the one
    section control point is bound to the left edge B-spline component control point. As result, there is nothing to
    maintain C0 continuity with the next section. As result a constraint will need to be manually set. To facilitate this,
    the array bind_inds will contain a list of the B-spline control point indicies that will need to be manually constrained to
    their previous sections.


    Parameters
    ----------
    prob : OpenMDAO problem object
        The OpenAeroStruct problem object with the unified B-spline component added.
    section_surfaces : list
        List of the surface dictionaries for each section.
    sec_cp : list
        List of B-spline control point arrays for each section.
    out_name: string
        Name of the unified B-spline component output to connect from
    comp_name: string
        Name of the unified B-spline component added to the problem object
    return_bind_inds: bool
        Return list of unjoined unified B-spline inidices. Default is False.

    Returns
    -------
    bind_inds : list
        List of unified B-spline control point indicies not connected due to the presence of a single control point section.(Only if return bind_inds specified)

    """
    acc = 0
    bind_inds = []
    for i, section in enumerate(section_surfaces):
        point_count = len(sec_cp[i])
        src_inds = np.arange(acc, acc + point_count)
        acc += point_count - 1
        if point_count == 1:
            acc += 1
            bind_inds.append(acc)
        prob.model.connect(
            "{}.{}".format(comp_name, out_name) + "_spline",
            "surface." + section["name"] + ".{}".format(out_name),
            src_indices=src_inds,
        )

    if return_bind_inds:
        return bind_inds
