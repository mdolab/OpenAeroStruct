import numpy as np
import openmdao.api as om



def build_multi_spline(name, num_sections, control_points):

    if len(control_points) != num_sections:
        raise Exception("Target sections need to match with control points!")
    
    single_sections = len([cp for cp in control_points if len(cp)==1])
    
    control_poin_vec = np.ones(len(np.concatenate(control_points))-(num_sections-1-single_sections))
    
    spline_control = om.IndepVarComp()
    spline_control.add_output("{}_spline".format(name), val=control_poin_vec)

    return spline_control




