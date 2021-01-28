import openmdao.api as om

from openmdao.utils.assert_utils import assert_check_partials
import numpy as np
from openaerostruct.geometry.utils import generate_mesh


def view_mat(mat1, mat2=None, key='Title', tol=1e-10):  # pragma: no cover
    """
    Helper function used to visually examine matrices. It plots mat1 and mat2 side by side,
    and shows the difference between the two.

    Parameters
    ----------
    mat1 : numpy array
        The Jacobian approximated by openMDAO
    mat2 : numpy array
        The Jacobian computed by compute_partials
    key : str
        The name of the tuple (of, wrt) for which the Jacobian is computed
    tol : float (Optional)
        The tolerance, below which the two numbers are considered the same for
        plotting purposes.

    """
    import matplotlib.pyplot as plt
    if len(mat1.shape) > 2:
        mat1 = np.sum(mat1, axis=2)
    if mat2 is not None:
        if len(mat2.shape) > 2:
            mat2 = np.sum(mat2, axis=2)
        vmin = np.amin(np.hstack((mat1.flatten(),mat2.flatten())))
        vmax = np.amax(np.hstack((mat1.flatten(),mat2.flatten())))
    else:
        vmin = np.amin(np.hstack((mat1.flatten())))
        vmax = np.amax(np.hstack((mat1.flatten())))
    if vmax-vmin < tol: # add small difference for plotting if both values are the same
        vmin = vmin - tol
        vmax = vmax + tol

    if mat2 is not None:
        fig, ax = plt.subplots(ncols=3,figsize=(12,6))
        ax[0].imshow(mat1.real, interpolation='none',vmin=vmin,vmax=vmax)
        ax[0].set_title('Approximated Jacobian')

        im = ax[1].imshow(mat2.real, interpolation='none',vmin=vmin,vmax=vmax)
        fig.colorbar(im, orientation='horizontal',ax=ax[0:2].ravel().tolist())
        ax[1].set_title('User-Defined Jacobian')

        diff = mat2.real - mat1.real
        if np.max(np.abs(diff).flatten()) < tol: # add small difference for plotting if diff is small
            vmin = -1*tol
            vmax = tol
        im2 = ax[2].imshow(diff, interpolation='none', vmin=vmin,vmax=vmax)
        fig.colorbar(im2, orientation='horizontal',ax=ax[2],aspect=10)
        ax[2].set_title('Difference')

    else:
        mtx = np.hstack((mat1.flatten()))
        vmin = np.nanmin(mtx[mtx != -np.inf])
        vmax = np.nanmax(mtx[mtx != np.inf])
        print(vmin, vmax)
        fig = plt.figure(figsize=(12,6))
        ax = plt.gca()
        im = plt.imshow(mat1.real, interpolation='none',vmin=vmin,vmax=vmax)
        fig.colorbar(im, orientation='horizontal', ax=ax, aspect=10)
        plt.title('Jacobian')

    plt.suptitle(key)
    plt.show()

def run_test(test_obj, comp, complex_flag=False, compact_print=True, method='fd', step=1e-6, atol=1e-5, rtol=1e-5, view=False):
    prob = om.Problem()
    prob.model.add_subsystem('comp', comp)
    prob.setup(force_alloc_complex=complex_flag)

    prob.run_model()

    if method=='cs':
        step = 1e-40

    check = prob.check_partials(compact_print=compact_print, method=method, step=step)

    if view:
        # Loop through this `check` dictionary and visualize the approximated
        # and computed derivatives
        for key, subjac in check[list(check.keys())[0]].items():
            view_mat(subjac['J_fd'],subjac['J_fwd'],key)

    assert_check_partials(check, atol=atol, rtol=rtol)

    return prob

def get_default_surfaces():
    # Create a dictionary to store options about the mesh
    mesh_dict = {'num_y' : 7,
                 'num_x' : 2,
                 'wing_type' : 'CRM',
                 'symmetry' : True,
                 'num_twist_cp' : 5}

    # Generate the aerodynamic mesh based on the previous dictionary
    mesh, twist_cp = generate_mesh(mesh_dict)

    wing_dict = {'name' : 'wing',
                 'num_y' : 4,
                 'num_x' : 2,
                 'symmetry' : True,
                 'groundplane' : False,
                 'S_ref_type' : 'wetted',
                 'CL0' : 0.1,
                 'CD0' : 0.1,
                 'mesh' : mesh,

                 # Airfoil properties for viscous drag calculation
                 'k_lam' : 0.05,         # percentage of chord with laminar
                                         # flow, used for viscous drag
                 't_over_c_cp' : np.array([0.15]),      # thickness over chord ratio (NACA0015)
                 'c_max_t' : .303,       # chordwise location of maximum (NACA0015)
                                         # thickness
                 'with_viscous' : True,  # if true, compute viscous drag
                 'with_wave' : False, # if true, computes wave drag
                 'fem_model_type' : 'tube',

                 # Structural values are based on aluminum 7075
                 'E' : 70.e9,            # [Pa] Young's modulus of the spar
                 'G' : 30.e9,            # [Pa] shear modulus of the spar
                 'yield' : 500.e6 / 2.5, # [Pa] yield stress divided by 2.5 for limiting case
                 'mrho' : 3.e3,          # [kg/m^3] material density
                 'fem_origin' : 0.35,    # normalized chordwise location of the spar
                 'wing_weight_ratio' : 2.,
                 'struct_weight_relief' : False,    # True to add the weight of the structure to the loads on the structure
                 'distributed_fuel_weight' : False,    # True to add the weight of the structure to the loads on the structure
                 'Wf_reserve' : 10000.,

                 }

    # Create a dictionary to store options about the mesh
    mesh_dict = {'num_y' : 5,
                 'num_x' : 3,
                 'wing_type' : 'rect',
                 'symmetry' : False}

    # Generate the aerodynamic mesh based on the previous dictionary
    mesh = generate_mesh(mesh_dict)

    tail_dict = {'name' : 'tail',
                 'num_y' : 5,
                 'num_x' : 3,
                 'symmetry' : False,
                 'mesh' : mesh}

    surfaces = [wing_dict, tail_dict]

    return surfaces


def get_ground_effect_surfaces():
    # Create a dictionary to store options about the mesh
    mesh_dict = {'num_y' : 7,
                 'num_x' : 2,
                 'wing_type' : 'CRM',
                 'symmetry' : True,
                 'num_twist_cp' : 5}

    # Generate the aerodynamic mesh based on the previous dictionary
    mesh, twist_cp = generate_mesh(mesh_dict)

    wing_dict = {'name' : 'wing',
                 'num_y' : 4,
                 'num_x' : 2,
                 'symmetry' : True,
                 'groundplane' : True,
                 'S_ref_type' : 'wetted',
                 'CL0' : 0.1,
                 'CD0' : 0.1,
                 'mesh' : mesh,

                 # Airfoil properties for viscous drag calculation
                 'k_lam' : 0.05,         # percentage of chord with laminar
                                         # flow, used for viscous drag
                 't_over_c_cp' : np.array([0.15]),      # thickness over chord ratio (NACA0015)
                 'c_max_t' : .303,       # chordwise location of maximum (NACA0015)
                                         # thickness
                 'with_viscous' : True,  # if true, compute viscous drag
                 'with_wave' : False, # if true, computes wave drag
                 'fem_model_type' : 'tube',

                 # Structural values are based on aluminum 7075
                 'E' : 70.e9,            # [Pa] Young's modulus of the spar
                 'G' : 30.e9,            # [Pa] shear modulus of the spar
                 'yield' : 500.e6 / 2.5, # [Pa] yield stress divided by 2.5 for limiting case
                 'mrho' : 3.e3,          # [kg/m^3] material density
                 'fem_origin' : 0.35,    # normalized chordwise location of the spar
                 'wing_weight_ratio' : 2.,
                 'struct_weight_relief' : False,    # True to add the weight of the structure to the loads on the structure
                 'distributed_fuel_weight' : False,    # True to add the weight of the structure to the loads on the structure
                 'Wf_reserve' : 10000.,

                 }

    # Create a dictionary to store options about the mesh
    mesh_dict = {'num_y' : 5,
                 'num_x' : 3,
                 'wing_type' : 'rect',
                 'symmetry' : True}

    # Generate the aerodynamic mesh based on the previous dictionary
    mesh = generate_mesh(mesh_dict)

    tail_dict = {'name' : 'tail',
                 'num_y' : 3,
                 'num_x' : 3,
                 'symmetry' : True,
                 'groundplane': True,
                 'mesh' : mesh}

    surfaces = [wing_dict, tail_dict]

    return surfaces

def assert_check_totals(totals, atol=1e-6, rtol=1e-6):
    for of_wrt_pair in totals.keys():
        of, wrt = of_wrt_pair
        total_dict = totals[of_wrt_pair]
        abs_error = total_dict['abs error']
        rel_error = total_dict['rel error']
        for i_mode, mode in enumerate(['forward','forward_reverse','reverse']):
            if rel_error[i_mode] is not None:
                if np.abs(rel_error[i_mode])-rtol > 0.0:
                    err_str = ('Total deriv (mode ' + mode + ') of ' + of +
                          ' with respect to ' + wrt + ' had rel error ' +
                          str(rel_error[i_mode]) + ' which is greater than tol ' + str(rtol))
                    raise ValueError(err_str)
            if abs_error[i_mode] is not None:
                if np.abs(abs_error[i_mode])-atol > 0.0:
                    err_str = ('Total deriv (mode ' + mode + ') of ' + of +
                          ' with respect to ' + wrt + ' had abs error ' +
                          str(abs_error[i_mode]) + ' which is greater than tol ' + str(atol))
                    raise ValueError(err_str)