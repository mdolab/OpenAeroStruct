import numpy as np
import matplotlib.pylab as plt

import openmdao.api as om
from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.geometry.geometry_group import Geometry
from openaerostruct.aerodynamics.aero_groups import AeroPoint


def compute_drag_polar_ground_effect(Mach, alphas, heights, surfaces, trimmed=False, visualize=False):
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
    indep_var_comp.add_output('height_agl', val=8000, units='m')
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
    prob.model.add_subsystem(point_name, aero_group, promotes_inputs=['v', 'alpha', 'Mach_number', 're', 'rho', 'cg', 'height_agl'])

    # For trimmed polar, setup balance component
    if trimmed == True:
        bal = om.BalanceComp()
        bal.add_balance(name='tail_rotation', rhs_val=0.0, units='deg')
        prob.model.add_subsystem('balance', bal, promotes_outputs=['tail_rotation'])
        prob.model.connect('aero.CM', 'balance.lhs:tail_rotation', src_indices=[1])
        prob.model.connect('tail_rotation', 'tail.twist_cp')

        prob.model.nonlinear_solver = om.NewtonSolver(solve_subsystems=True)

        prob.model.nonlinear_solver.options['iprint'] = 2
        prob.model.nonlinear_solver.options['maxiter'] = 10
        prob.model.linear_solver = om.DirectSolver()
    if visualize:
        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options['tol'] = 1e-9
        recorder = om.SqliteRecorder("polar_ground_effect.db")
        prob.driver.add_recorder(recorder)
        prob.driver.recording_options['record_derivatives'] = True
        prob.driver.recording_options['includes'] = ['*']
        prob.driver.options['maxiter'] = 1
        # Setup problem and add design variables, constraint, and objective
        prob.model.add_design_var('cg', lower=-0.0, upper=0.0)
        prob.model.add_objective('aero.CL', scaler=1e4)

    prob.setup()

    # prob['tail_rotation'] = -0.75
    span = prob['wing.mesh.stretch.span']

    prob.run_model()
    # prob.check_partials(compact_print = True)
    # prob.model.list_outputs(prom_name = True)

    prob.model.list_outputs(residuals=True)

    CLs = []
    CDs = []
    CMs = []
    for height in heights:
        CLs_h = []
        CDs_h = []
        CMs_h = []
        for a in alphas:
            prob['alpha'] = a
            prob['height_agl'] = height
            prob.run_model()
            CLs_h.append(prob['aero.CL'][0])
            CDs_h.append(prob['aero.CD'][0])
            CMs_h.append(prob['aero.CM'][1])  # Take only the longitudinal CM
            # print(a, prob['aero.CL'], prob['aero.CD'], prob['aero.CM'][1])
        CLs.append(CLs_h)
        CDs.append(CDs_h)
        CMs.append(CMs_h)
    # Plot the drag polar
    if visualize:
        prob.run_driver()
    fig, ax = plt.subplots(nrows=1)
    for ih, height in enumerate(heights):
        ax.plot(CLs[ih], CDs[ih], label='h/b='+'%.1f' % (height/span)[0])
        ax.set_ylabel('CDi')
        ax.set_xlabel('CL')
    ax.legend()
    fig, ax = plt.subplots(nrows=1)
    # compute ground effect correction factor and plot
    CLs_arr = np.array(CLs)
    CDs_arr = np.array(CDs)
    k = CDs_arr[:,-1] / CLs_arr[:,-1] ** 2
    k = k / k[0]
    hob = np.array(heights) / span
    ax.plot(hob, k)
    ax.set_xscale('log')
    ax.set_xlabel('h/b')
    ax.set_ylabel('$C_{Di}/C_{Di,\infty}$')
    plt.show()
    return CLs, CDs, CMs


if __name__ == '__main__':
    # Create a dictionary to store options about the mesh

    nx = 3
    ny = 7
    mesh_dict = {
        'num_y': ny,
        'num_x': nx,
        'wing_type': 'rect',
        'symmetry': True,
        'span': 12.0,
        'chord': 1,
        'span_cos_spacing': 1.0,
        'chord_cos_spacing': 0.0,
    }

    # Generate half-wing mesh of rectangular wing
    mesh = generate_mesh(mesh_dict)

 # Create a dictionary with info and options about the wing
    wing_surface = {
        # Wing definition
        'name': 'wing',  # name of the surface
        'symmetry': True,  # if true, model one half of wing
        'groundplane': True,
        # reflected across the plane y = 0
        'S_ref_type': 'wetted',  # how we compute the wing area,
        # can be 'wetted' or 'projected'
        'fem_model_type': 'tube',
        # 'twist_cp': twist_cp,
        'mesh': mesh,
        # Aerodynamic performance of the lifting surface at
        # an angle of attack of 0 (alpha=0).
        # These CL0 and CD0 values are added to the CL and CD
        # obtained from aerodynamic analysis of the surface to get
        # the total CL and CD.
        # These CL0 and CD0 values do not vary wrt alpha.
        'CL0': 0.0,  # CL of the surface at alpha=0
        'CD0': 0.0,  # CD of the surface at alpha=0
        # Airfoil properties for viscous drag calculation
        'k_lam': 0.05,  # percentage of chord with laminar
        # flow, used for viscous drag
        't_over_c': 0.15,  # thickness over chord ratio (NACA0015)
        'c_max_t': 0.303,  # chordwise location of maximum (NACA0015)
        # thickness
        'with_viscous': False,  # if true, compute viscous drag
        'with_wave': False,
    }
    surfaces = [wing_surface]

    Mach = 0.3
    alphas = np.linspace(-5, 5, 10)
    heights = np.array([1200.,24.,12.,6.,3.,1.2])
    # alphas = [0.]

    CL, CD, CM = compute_drag_polar_ground_effect(Mach, alphas, heights, surfaces, trimmed=False, visualize=True)