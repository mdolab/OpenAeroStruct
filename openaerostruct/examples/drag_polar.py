import numpy as np
import matplotlib.pylab as plt

import openmdao.api as om
from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.geometry.geometry_group import Geometry
from openaerostruct.aerodynamics.aero_groups import AeroPoint


def compute_drag_polar(Mach, alphas, surfaces, trimmed=False):

    if isinstance(surfaces, dict):
        surfaces = [
            surfaces,
        ]

    # Create the OpenMDAO problem
    prob = om.Problem()
    # Create an independent variable component that will supply the flow
    # conditions to the problem.
    indep_var_comp = om.IndepVarComp()
    indep_var_comp.add_output('v', val=248.136, units='m/s')
    indep_var_comp.add_output('alpha', val=0.0, units='deg')
    indep_var_comp.add_output('Mach_number', val=Mach)
    indep_var_comp.add_output('re', val=1.0e6, units='1/m')
    indep_var_comp.add_output('rho', val=0.38, units='kg/m**3')
    indep_var_comp.add_output('cg', val=np.zeros((3)), units='m')
    # Add this IndepVarComp to the problem model
    prob.model.add_subsystem('prob_vars', indep_var_comp, promotes=['*'])

    for surface in surfaces:
        name = surface['name']
        # Create and add a group that handles the geometry for the
        # aerodynamic lifting surface
        geom_group = Geometry(surface=surface)
        prob.model.add_subsystem(name, geom_group)

        # Connect the mesh from the geometry component to the analysis point
        prob.model.connect(name + '.mesh', 'aero.' + name + '.def_mesh')
        # Perform the connections with the modified names within the
        # 'aero_states' group.
        prob.model.connect(name + '.mesh', 'aero.aero_states.' + name + '_def_mesh')

    # Create the aero point group, which contains the actual aerodynamic
    # analyses
    point_name = 'aero'
    aero_group = AeroPoint(surfaces=surfaces)
    prob.model.add_subsystem(point_name, aero_group, promotes_inputs=['v', 'alpha', 'Mach_number', 're', 'rho', 'cg'])

    # For trimmed polar, setup balance component
    if trimmed == True:
        bal = om.BalanceComp()
        bal.add_balance(name='tail_rotation', rhs_val=0.0, units='deg')
        prob.model.add_subsystem('balance', bal, promotes_outputs=['tail_rotation'])
        prob.model.connect('aero.CM', 'balance.lhs:tail_rotation', src_indices=[1])
        prob.model.connect('tail_rotation', 'tail.twist_cp')

        prob.model.nonlinear_solver = om.NonlinearBlockGS(use_aitken=True)

        prob.model.nonlinear_solver.options['iprint'] = 2
        prob.model.nonlinear_solver.options['maxiter'] = 100
        prob.model.linear_solver = om.DirectSolver()

    prob.setup()

    # prob['tail_rotation'] = -0.75

    prob.run_model()
    # prob.check_partials(compact_print = True)
    # prob.model.list_outputs(prom_name = True)

    prob.model.list_outputs(residuals=True)

    CLs = []
    CDs = []
    CMs = []

    for a in alphas:
        prob['alpha'] = a
        prob.run_model()
        CLs.append(prob['aero.CL'][0])
        CDs.append(prob['aero.CD'][0])
        CMs.append(prob['aero.CM'][1])  # Take only the longitudinal CM
        # print(a, prob['aero.CL'], prob['aero.CD'], prob['aero.CM'][1])

    # Plot CL vs alpha and drag polar
    fig, axes = plt.subplots(nrows=3)
    axes[0].plot(alphas, CLs)
    axes[1].plot(alphas, CMs)
    axes[2].plot(CLs, CDs)
    fig.savefig('drag_polar.pdf')
    # plt.show()

    return CLs, CDs, CMs


if __name__ == '__main__':
    # Create a dictionary to store options about the mesh
    mesh_dict = {'num_y': 7, 'num_x': 2, 'wing_type': 'CRM', 'symmetry': True, 'num_twist_cp': 5}

    # Generate the aerodynamic mesh based on the previous dictionary
    mesh, twist_cp = generate_mesh(mesh_dict)

    # Create a dictionary with info and options about the wing
    wing_surface = {
        # Wing definition
        'name': 'wing',  # name of the surface
        'symmetry': True,  # if true, model one half of wing
        # reflected across the plane y = 0
        'S_ref_type': 'wetted',  # how we compute the wing area,
        # can be 'wetted' or 'projected'
        'fem_model_type': 'tube',
        'twist_cp': twist_cp,
        'mesh': mesh,
        # Aerodynamic performance of the lifting surface at
        # an angle of attack of 0 (alpha=0).
        # These CL0 and CD0 values are added to the CL and CD
        # obtained from aerodynamic analysis of the surface to get
        # the total CL and CD.
        # These CL0 and CD0 values do not vary wrt alpha.
        'CL0': 0.0,  # CL of the surface at alpha=0
        'CD0': 0.015,  # CD of the surface at alpha=0
        # Airfoil properties for viscous drag calculation
        'k_lam': 0.05,  # percentage of chord with laminar
        # flow, used for viscous drag
        't_over_c': 0.15,  # thickness over chord ratio (NACA0015)
        'c_max_t': 0.303,  # chordwise location of maximum (NACA0015)
        # thickness
        'with_viscous': True,  # if true, compute viscous drag
        'with_wave': False,
    }

    # Create a dictionary to store options about the tail surface
    mesh_dict = {'num_y': 7, 'num_x': 2, 'wing_type': 'rect', 'symmetry': True, 'offset': np.array([50, 0.0, 0.0])}

    mesh = generate_mesh(mesh_dict)

    tail_surface = {
        # Wing definition
        'name': 'tail',  # name of the surface
        'symmetry': True,  # if true, model one half of wing
        # reflected across the plane y = 0
        'S_ref_type': 'wetted',  # how we compute the wing area,
        # can be 'wetted' or 'projected'
        'twist_cp': np.zeros((1)),
        'twist_cp_dv': False,
        'mesh': mesh,
        # Aerodynamic performance of the lifting surface at
        # an angle of attack of 0 (alpha=0).
        # These CL0 and CD0 values are added to the CL and CD
        # obtained from aerodynamic analysis of the surface to get
        # the total CL and CD.
        # These CL0 and CD0 values do not vary wrt alpha.
        'CL0': 0.0,  # CL of the surface at alpha=0
        'CD0': 0.0,  # CD of the surface at alpha=0
        'fem_origin': 0.35,
        # Airfoil properties for viscous drag calculation
        'k_lam': 0.05,  # percentage of chord with laminar
        # flow, used for viscous drag
        't_over_c': 0.15,  # thickness over chord ratio (NACA0015)
        'c_max_t': 0.303,  # chordwise location of maximum (NACA0015)
        # thickness
        'with_viscous': True,  # if true, compute viscous drag
        'with_wave': False,
    }

    surfaces = [wing_surface, tail_surface]

    Mach = 0.82
    alphas = np.linspace(-10, 15, 25)
    # alphas = [0.]

    CL, CD, CM = compute_drag_polar(Mach, alphas, surfaces, trimmed=True)
