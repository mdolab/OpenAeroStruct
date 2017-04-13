"""
Define the structural analysis component using spatial beam theory.

Each FEM element has 6-DOF; translation in the x, y, and z direction and
rotation about the x, y, and z-axes.

"""

from __future__ import division, print_function
import numpy as np

from openmdao.api import Component, Group
from scipy.linalg import lu_factor, lu_solve

try:
    import OAS_API
    fortran_flag = True
    data_type = float
except:
    fortran_flag = False
    data_type = complex


def norm(vec):
    return np.sqrt(np.sum(vec**2))


def unit(vec):
    return vec / norm(vec)


def radii(mesh, t_c=0.15):
    """ Obtain the radii of the FEM component based on chord. """
    vectors = mesh[-1, :, :] - mesh[0, :, :]
    chords = np.sqrt(np.sum(vectors**2, axis=1))
    chords = 0.5 * chords[:-1] + 0.5 * chords[1:]
    return t_c * chords


def _assemble_system(nodes, A, J, Iy, Iz, loads,
                     K_a, K_t, K_y, K_z,
                     cons, E, G, x_gl, T,
                     K_elem, S_a, S_t, S_y, S_z, T_elem,
                     const2, const_y, const_z, n, size, K, forces):

    """
    Assemble the structural stiffness matrix based on 6 degrees of freedom
    per element.

    Can be run in Fortran or Python code depending on the flags used.
    """

    # Populate the right-hand side of the linear system using the
    # prescribed or computed loads
    forces[:] = 0.0
    forces[:6*n] = loads.reshape(n*6)
    forces[np.abs(forces) < 1e-6] = 0.

    # Fortran
    if fortran_flag:
        K = OAS_API.oas_api.assemblestructmtx(nodes, A, J, Iy, Iz,
                                     K_a, K_t, K_y, K_z,
                                     cons, E, G, x_gl, T,
                                     K_elem, S_a, S_t, S_y, S_z, T_elem,
                                     const2, const_y, const_z)

    # Python
    else:

        K[:] = 0.

        # Loop over each element
        for ielem in range(n-1):

            # Obtain the element nodes
            P0 = nodes[ielem, :]
            P1 = nodes[ielem+1, :]

            x_loc = unit(P1 - P0)
            y_loc = unit(np.cross(x_loc, x_gl))
            z_loc = unit(np.cross(x_loc, y_loc))

            T[0, :] = x_loc
            T[1, :] = y_loc
            T[2, :] = z_loc

            for ind in range(4):
                T_elem[3*ind:3*ind+3, 3*ind:3*ind+3] = T

            L = norm(P1 - P0)
            EA_L = E * A[ielem] / L
            GJ_L = G * J[ielem] / L
            EIy_L3 = E * Iy[ielem] / L**3
            EIz_L3 = E * Iz[ielem] / L**3

            K_a[:, :] = EA_L * const2
            K_t[:, :] = GJ_L * const2

            K_y[:, :] = EIy_L3 * const_y
            K_y[1, :] *= L
            K_y[3, :] *= L
            K_y[:, 1] *= L
            K_y[:, 3] *= L

            K_z[:, :] = EIz_L3 * const_z
            K_z[1, :] *= L
            K_z[3, :] *= L
            K_z[:, 1] *= L
            K_z[:, 3] *= L

            K_elem[:] = 0
            K_elem += S_a.T.dot(K_a).dot(S_a)
            K_elem += S_t.T.dot(K_t).dot(S_t)
            K_elem += S_y.T.dot(K_y).dot(S_y)
            K_elem += S_z.T.dot(K_z).dot(S_z)

            res = T_elem.T.dot(K_elem).dot(T_elem)

            in0, in1 = ielem, ielem+1

            # Populate the full matrix with stiffness
            # contributions from each node
            K[6*in0:6*in0+6, 6*in0:6*in0+6] += res[:6, :6]
            K[6*in1:6*in1+6, 6*in0:6*in0+6] += res[6:, :6]
            K[6*in0:6*in0+6, 6*in1:6*in1+6] += res[:6, 6:]
            K[6*in1:6*in1+6, 6*in1:6*in1+6] += res[6:, 6:]

        # Include a scaled identity matrix in the rows and columns
        # corresponding to the structural constraints.
        # Hardcoded 1 constraint for now.
        for ind in range(1):
            for k in range(6):
                K[-6+k, 6*cons+k] = 1.e9
                K[6*cons+k, -6+k] = 1.e9

    return K, forces


class AssembleK(Component):
    """
    Compute the displacements and rotations by solving the linear system
    using the structural stiffness matrix.

    Parameters
    ----------
    A[ny-1] : numpy array
        Areas for each FEM element.
    Iy[ny-1] : numpy array
        Mass moment of inertia around the y-axis for each FEM element.
    Iz[ny-1] : numpy array
        Mass moment of inertia around the z-axis for each FEM element.
    J[ny-1] : numpy array
        Polar moment of inertia for each FEM element.
    nodes[ny, 3] : numpy array
        Flattened array with coordinates for each FEM node.
    loads[ny, 6] : numpy array
        Flattened array containing the loads applied on the FEM component,
        computed from the sectional forces.

    Returns
    -------
    K[(nx-1)*(ny-1), (nx-1)*(ny-1)] : numpy array
        Stiffness matrix for the entire FEM system. Used to solve the linear
        system K * u = f to obtain the displacements, u.
    forces[(nx-1)*(ny-1)] : numpy array
        Right-hand-side of the linear system. The loads from the aerodynamic
        analysis or the user-defined loads.
    """

    def __init__(self, surface, cg_x=5):
        super(AssembleK, self).__init__()

        self.ny = surface['num_y']

        self.size = size = 6 * self.ny + 6

        self.add_param('A', val=np.zeros((self.ny - 1), dtype=data_type))
        self.add_param('Iy', val=np.zeros((self.ny - 1), dtype=data_type))
        self.add_param('Iz', val=np.zeros((self.ny - 1), dtype=data_type))
        self.add_param('J', val=np.zeros((self.ny - 1), dtype=data_type))
        self.add_param('nodes', val=np.zeros((self.ny, 3), dtype=data_type))
        self.add_param('loads', val=np.zeros((self.ny, 6), dtype=data_type))

        self.add_output('K', val=np.zeros((size, size), dtype=data_type))
        self.add_output('forces', val=np.zeros((size), dtype=data_type))

        self.E = surface['E']
        self.G = surface['G']
        self.cg_x = cg_x

        self.const2 = np.array([
            [1, -1],
            [-1, 1],
        ], dtype=data_type)
        self.const_y = np.array([
            [12, -6, -12, -6],
            [-6, 4, 6, 2],
            [-12, 6, 12, 6],
            [-6, 2, 6, 4],
        ], dtype=data_type)
        self.const_z = np.array([
            [12, 6, -12, 6],
            [6, 4, -6, 2],
            [-12, -6, 12, -6],
            [6, 2, -6, 4],
        ], dtype=data_type)
        self.x_gl = np.array([1, 0, 0], dtype=data_type)

        self.K_elem = np.zeros((12, 12), dtype=data_type)
        self.T_elem = np.zeros((12, 12), dtype=data_type)
        self.T = np.zeros((3, 3), dtype=data_type)

        self.K = np.zeros((size, size), dtype=data_type)
        self.forces = np.zeros((size), dtype=data_type)

        self.K_a = np.zeros((2, 2), dtype=data_type)
        self.K_t = np.zeros((2, 2), dtype=data_type)
        self.K_y = np.zeros((4, 4), dtype=data_type)
        self.K_z = np.zeros((4, 4), dtype=data_type)

        self.S_a = np.zeros((2, 12), dtype=data_type)
        self.S_a[(0, 1), (0, 6)] = 1.

        self.S_t = np.zeros((2, 12), dtype=data_type)
        self.S_t[(0, 1), (3, 9)] = 1.

        self.S_y = np.zeros((4, 12), dtype=data_type)
        self.S_y[(0, 1, 2, 3), (2, 4, 8, 10)] = 1.

        self.S_z = np.zeros((4, 12), dtype=data_type)
        self.S_z[(0, 1, 2, 3), (1, 5, 7, 11)] = 1.

        if not fortran_flag:
            self.deriv_options['type'] = 'cs'
            self.deriv_options['form'] = 'central'

    def solve_nonlinear(self, params, unknowns, resids):

        # Find constrained nodes based on closeness to specified cg point
        nodes = params['nodes']
        dist = nodes - np.array([self.cg_x, 0, 0])
        idx = (np.linalg.norm(dist, axis=1)).argmin()
        self.cons = idx

        loads = params['loads']

        self.K, self.forces = \
            _assemble_system(params['nodes'],
                             params['A'], params['J'], params['Iy'],
                             params['Iz'], loads, self.K_a, self.K_t,
                             self.K_y, self.K_z, self.cons,
                             self.E, self.G, self.x_gl, self.T, self.K_elem,
                             self.S_a, self.S_t, self.S_y, self.S_z,
                             self.T_elem, self.const2, self.const_y,
                             self.const_z, self.ny, self.size,
                             self.K, self.forces)

        unknowns['K'] = self.K
        unknowns['forces'] = self.forces

    def apply_linear(self, params, unknowns, dparams, dunknowns, dresids, mode):

        # Find constrained nodes based on closeness to specified cg point
        nodes = params['nodes']
        dist = nodes - np.array([self.cg_x, 0, 0])
        idx = (np.linalg.norm(dist, axis=1)).argmin()
        self.cons = idx

        loads = params['loads']

        A = params['A']
        J = params['J']
        Iy = params['Iy']
        Iz = params['Iz']

        if mode == 'fwd':
            K, Kd = OAS_API.oas_api.assemblestructmtx_d(nodes, dparams['nodes'], A, dparams['A'],
                                         J, dparams['J'], Iy, dparams['Iy'],
                                         Iz, dparams['Iz'],
                                         self.K_a, self.K_t, self.K_y, self.K_z,
                                         self.cons, self.E, self.G, self.x_gl, self.T,
                                         self.K_elem, self.S_a, self.S_t, self.S_y, self.S_z, self.T_elem,
                                         self.const2, self.const_y, self.const_z)

            dresids['K'] += Kd
            dresids['forces'][:-6] += dparams['loads'].reshape(-1)

        if mode == 'rev':
            nodesb, Ab, Jb, Iyb, Izb = OAS_API.oas_api.assemblestructmtx_b(nodes, A, J, Iy, Iz,
                                self.K_a, self.K_t, self.K_y, self.K_z,
                                self.cons, self.E, self.G, self.x_gl, self.T,
                                self.K_elem, self.S_a, self.S_t, self.S_y, self.S_z, self.T_elem,
                                self.const2, self.const_y, self.const_z, self.K, dresids['K'])

            dparams['nodes'] += nodesb
            dparams['A'] += Ab
            dparams['J'] += Jb
            dparams['Iy'] += Iyb
            dparams['Iz'] += Izb

            dparams['loads'] += dresids['forces'][:-6].reshape(-1, 6)


class SpatialBeamFEM(Component):
    """
    Compute the displacements and rotations by solving the linear system
    using the structural stiffness matrix.
    This component is copied from OpenMDAO's LinearSystem component with the
    names of the parameters and outputs changed to match our problem formulation.

    Parameters
    ----------
    K[6*(ny+1), 6*(ny+1)] : numpy array
        Stiffness matrix for the entire FEM system. Used to solve the linear
        system K * u = f to obtain the displacements, u.
    forces[6*(ny+1)] : numpy array
        Right-hand-side of the linear system. The loads from the aerodynamic
        analysis or the user-defined loads.

    Returns
    -------
    disp_aug[6*(ny+1)] : numpy array
        Augmented displacement array. Obtained by solving the system
        K * u = f, where f is a flattened version of loads.

    """

    def __init__(self, size):
        super(SpatialBeamFEM, self).__init__()

        self.add_param('K', val=np.zeros((size, size), dtype=data_type))
        self.add_param('forces', val=np.zeros((size), dtype=data_type))
        self.add_state('disp_aug', val=np.zeros((size), dtype=data_type))

        self.size = size

        # cache
        self.lup = None
        self.forces_cache = None

    def solve_nonlinear(self, params, unknowns, resids):
        """ Use np to solve Ax=b for x.
        """

        # lu factorization for use with solve_linear
        self.lup = lu_factor(params['K'])

        unknowns['disp_aug'] = lu_solve(self.lup, params['forces'])
        resids['disp_aug'] = params['K'].dot(unknowns['disp_aug']) - params['forces']

    def apply_nonlinear(self, params, unknowns, resids):
        """Evaluating residual for given state."""

        resids['disp_aug'] = params['K'].dot(unknowns['disp_aug']) - params['forces']

    def apply_linear(self, params, unknowns, dparams, dunknowns, dresids, mode):
        """ Apply the derivative of state variable with respect to
        everything."""

        if mode == 'fwd':

            if 'disp_aug' in dunknowns:
                dresids['disp_aug'] += params['K'].dot(dunknowns['disp_aug'])
            if 'K' in dparams:
                dresids['disp_aug'] += dparams['K'].dot(unknowns['disp_aug'])
            if 'forces' in dparams:
                dresids['disp_aug'] -= dparams['forces']

        elif mode == 'rev':

            if 'disp_aug' in dunknowns:
                dunknowns['disp_aug'] += params['K'].T.dot(dresids['disp_aug'])
            if 'K' in dparams:
                dparams['K'] += np.outer(unknowns['disp_aug'], dresids['disp_aug']).T
            if 'forces' in dparams:
                dparams['forces'] -= dresids['disp_aug']

    def solve_linear(self, dumat, drmat, vois, mode=None):
        """ LU backsubstitution to solve the derivatives of the linear system."""

        if mode == 'fwd':
            sol_vec, forces_vec = self.dumat, self.drmat
            t=0
        else:
            sol_vec, forces_vec = self.drmat, self.dumat
            t=1

        if self.forces_cache is None:
            self.forces_cache = np.zeros((self.size, ))
        forces = self.forces_cache

        for voi in vois:
            forces[:] = forces_vec[voi]['disp_aug']

            sol = lu_solve(self.lup, forces, trans=t)

            sol_vec[voi]['disp_aug'] = sol[:]


class SpatialBeamDisp(Component):
    """
    Reshape the flattened displacements from the linear system solution into
    a 2D array so we can more easily use the results.

    The solution to the linear system has additional results due to the
    constraints on the FEM model. The displacements from this portion of
    the linear system are not needed, so we select only the relevant
    portion of the displacements for further calculations.

    Parameters
    ----------
    disp_aug[6*(ny+1)] : numpy array
        Augmented displacement array. Obtained by solving the system
        K * disp_aug = forces, where forces is a flattened version of loads.

    Returns
    -------
    disp[6*ny] : numpy array
        Actual displacement array formed by truncating disp_aug.

    """

    def __init__(self, surface):
        super(SpatialBeamDisp, self).__init__()

        self.ny = surface['num_y']

        self.add_param('disp_aug', val=np.zeros(((self.ny+1)*6), dtype=data_type))
        self.add_output('disp', val=np.zeros((self.ny, 6), dtype=data_type))

    def solve_nonlinear(self, params, unknowns, resids):
        # Obtain the relevant portions of disp_aug and store the reshaped
        # displacements in disp
        unknowns['disp'] = params['disp_aug'][:-6].reshape((-1, 6))

    def linearize(self, params, unknowns, resids):
        jac = self.alloc_jacobian()
        n = self.ny * 6
        jac['disp', 'disp_aug'][:n, :n] = np.eye((n))
        return jac


class ComputeNodes(Component):
    """
    Compute FEM nodes based on aerodynamic mesh.

    The FEM nodes are placed at fem_origin * chord,
    with the default fem_origin = 0.35.

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Array defining the nodal points of the lifting surface.

    Returns
    -------
    nodes[ny, 3] : numpy array
        Flattened array with coordinates for each FEM node.

    """

    def __init__(self, surface):
        super(ComputeNodes, self).__init__()

        self.ny = surface['num_y']
        self.nx = surface['num_x']
        self.fem_origin = surface['fem_origin']

        self.add_param('mesh', val=np.zeros((self.nx, self.ny, 3), dtype=data_type))
        self.add_output('nodes', val=np.zeros((self.ny, 3), dtype=data_type))

    def solve_nonlinear(self, params, unknowns, resids):
        w = self.fem_origin
        mesh = params['mesh']

        unknowns['nodes'] = (1-w) * mesh[0, :, :] + w * mesh[-1, :, :]

    def linearize(self, params, unknowns, resids):
        jac = self.alloc_jacobian()
        w = self.fem_origin

        n = self.ny * 3
        jac['nodes', 'mesh'][:n, :n] = np.eye(n) * (1-w)
        jac['nodes', 'mesh'][:n, -n:] = np.eye(n) * w

        return jac


class SpatialBeamEnergy(Component):
    """ Compute strain energy.

    Parameters
    ----------
    disp[ny, 6] : numpy array
        Actual displacement array formed by truncating disp_aug.
    loads[ny, 6] : numpy array
        Array containing the loads applied on the FEM component,
        computed from the sectional forces.

    Returns
    -------
    energy : float
        Total strain energy of the structural component.

    """

    def __init__(self, surface):
        super(SpatialBeamEnergy, self).__init__()

        ny = surface['num_y']

        self.add_param('disp', val=np.zeros((ny, 6), dtype=data_type))
        self.add_param('loads', val=np.zeros((ny, 6), dtype=data_type))
        self.add_output('energy', val=0.)

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['energy'] = np.sum(params['disp'] * params['loads'])

    def linearize(self, params, unknowns, resids):
        jac = self.alloc_jacobian()
        jac['energy', 'disp'][0, :] = params['loads'].real.flatten()
        jac['energy', 'loads'][0, :] = params['disp'].real.flatten()
        return jac

class SpatialBeamWeight(Component):
    """ Compute total weight.

    Parameters
    ----------
    A[ny-1] : numpy array
        Areas for each FEM element.
    nodes[ny, 3] : numpy array
        Flattened array with coordinates for each FEM node.

    Returns
    -------
    weight : float
        Total weight of the structural spar.
    """

    def __init__(self, surface):
        super(SpatialBeamWeight, self).__init__()

        self.surface = surface

        self.ny = surface['num_y']

        self.add_param('A', val=np.zeros((self.ny - 1), dtype=data_type))
        self.add_param('nodes', val=np.zeros((self.ny, 3), dtype=data_type))
        self.add_output('structural_weight', val=0.)

    def solve_nonlinear(self, params, unknowns, resids):
        A = params['A']
        nodes = params['nodes']

        # Calculate the volume and weight of the total structure
        volume = np.sum(np.linalg.norm(nodes[1:, :] - nodes[:-1, :], axis=1) * A)

        weight = volume * self.surface['mrho'] * 9.81

        if self.surface['symmetry']:
            weight *= 2.

        unknowns['structural_weight'] = weight

    def linearize(self, params, unknowns, resids):
        jac = self.alloc_jacobian()
        A = params['A']
        nodes = params['nodes']

        # First we will solve for dweight_dA
        # Calculate the volume and weight of the total structure
        norms = np.linalg.norm(nodes[1:, :] - nodes[:-1, :], axis=1).reshape(1, -1)

        # Multiply by the material density and force of gravity
        dweight_dA = norms * self.surface['mrho'] * 9.81

        # Account for symmetry
        if self.surface['symmetry']:
            dweight_dA *= 2.

        # Save the result to the jacobian dictionary
        jac['structural_weight', 'A'] = dweight_dA

        # Next, we will compute the derivative of weight wrt nodes.
        # Here we're using results from AD to compute the derivative
        # Initialize the reverse seeds.
        nodesb = np.zeros(nodes.shape)
        tempb = (nodes[1:, :] - nodes[:-1, :]) * (A / norms).reshape(-1, 1)
        nodesb[1:, :] += tempb
        nodesb[:-1, :] -= tempb

        # Apply the multipliers for material properties and symmetry
        nodesb *= self.surface['mrho'] * 9.81

        if self.surface['symmetry']:
            nodesb *= 2.

        # Store the flattened array in the jacobian dictionary
        jac['structural_weight', 'nodes'] = nodesb.reshape(1, -1)

        return jac

class SpatialBeamVonMisesTube(Component):
    """ Compute the von Mises stress in each element.

    Parameters
    ----------
    radius[ny-1] : numpy array
        Radii for each FEM element.
    nodes[ny, 3] : numpy array
        Flattened array with coordinates for each FEM node.
    disp[ny, 6] : numpy array
        Displacements of each FEM node.

    Returns
    -------
    vonmises[ny-1, 2] : numpy array
        von Mises stress magnitudes for each FEM element.

    """

    def __init__(self, surface):
        super(SpatialBeamVonMisesTube, self).__init__()

        self.ny = surface['num_y']

        self.add_param('nodes', val=np.zeros((self.ny, 3),
                       dtype=data_type))
        self.add_param('radius', val=np.zeros((self.ny - 1),
                       dtype=data_type))
        self.add_param('disp', val=np.zeros((self.ny, 6),
                       dtype=data_type))

        self.add_output('vonmises', val=np.zeros((self.ny-1, 2),
                        dtype=data_type))

        if not fortran_flag:
            self.deriv_options['type'] = 'cs'
            self.deriv_options['form'] = 'central'

        self.E = surface['E']
        self.G = surface['G']

        self.T = np.zeros((3, 3), dtype=data_type)
        self.x_gl = np.array([1, 0, 0], dtype=data_type)
        self.t = 0

    def solve_nonlinear(self, params, unknowns, resids):
        radius = params['radius']
        disp = params['disp']
        nodes = params['nodes']
        vonmises = unknowns['vonmises']
        T = self.T
        E = self.E
        G = self.G
        x_gl = self.x_gl

        if fortran_flag:
            vm = OAS_API.oas_api.calc_vonmises(nodes, radius, disp, E, G, x_gl)
            unknowns['vonmises'] = vm

        else:

            num_elems = self.ny - 1
            for ielem in range(self.ny-1):

                P0 = nodes[ielem, :]
                P1 = nodes[ielem+1, :]
                L = norm(P1 - P0)

                x_loc = unit(P1 - P0)
                y_loc = unit(np.cross(x_loc, x_gl))
                z_loc = unit(np.cross(x_loc, y_loc))

                T[0, :] = x_loc
                T[1, :] = y_loc
                T[2, :] = z_loc

                u0x, u0y, u0z = T.dot(disp[ielem, :3])
                r0x, r0y, r0z = T.dot(disp[ielem, 3:])
                u1x, u1y, u1z = T.dot(disp[ielem+1, :3])
                r1x, r1y, r1z = T.dot(disp[ielem+1, 3:])

                tmp = np.sqrt((r1y - r0y)**2 + (r1z - r0z)**2)
                sxx0 = E * (u1x - u0x) / L + E * radius[ielem] / L * tmp
                sxx1 = E * (u0x - u1x) / L + E * radius[ielem] / L * tmp
                sxt = G * radius[ielem] * (r1x - r0x) / L

                vonmises[ielem, 0] = np.sqrt(sxx0**2 + sxt**2)
                vonmises[ielem, 1] = np.sqrt(sxx1**2 + sxt**2)

    def apply_linear(self, params, unknowns, dparams, dunknowns, dresids, mode):

        radius = params['radius'].real
        disp = params['disp'].real
        nodes = params['nodes'].real
        vonmises = unknowns['vonmises'].real
        E = self.E
        G = self.G
        x_gl = self.x_gl

        if mode == 'fwd':
            _, vonmisesd = OAS_API.oas_api.calc_vonmises_d(nodes, dparams['nodes'], radius, dparams['radius'], disp, dparams['disp'], E, G, x_gl)
            dresids['vonmises'] += vonmisesd

        if mode == 'rev':
            nodesb, radiusb, dispb = OAS_API.oas_api.calc_vonmises_b(nodes, radius, disp, E, G, x_gl, vonmises, dresids['vonmises'])
            dparams['nodes'] += nodesb
            dparams['radius'] += radiusb
            dparams['disp'] += dispb

class SpatialBeamFailureKS(Component):
    """
    Aggregate failure constraints from the structure.

    To simplify the optimization problem, we aggregate the individual
    elemental failure constraints using a Kreisselmeier-Steinhauser (KS)
    function.

    The KS function produces a smoother constraint than using a max() function
    to find the maximum point of failure, which produces a better-posed
    optimization problem.

    The rho parameter controls how conservatively the KS function aggregates
    the failure constraints. A lower value is more conservative while a greater
    value is more aggressive (closer approximation to the max() function).

    Parameters
    ----------
    vonmises[ny-1, 2] : numpy array
        von Mises stress magnitudes for each FEM element.

    Returns
    -------
    failure : float
        KS aggregation quantity obtained by combining the failure criteria
        for each FEM node. Used to simplify the optimization problem by
        reducing the number of constraints.

    """

    def __init__(self, surface, rho=100):
        super(SpatialBeamFailureKS, self).__init__()

        self.ny = surface['num_y']

        self.add_param('vonmises', val=np.zeros((self.ny-1, 2), dtype=data_type))
        self.add_output('failure', val=0.)

        self.sigma = surface['stress']
        self.rho = rho

    def solve_nonlinear(self, params, unknowns, resids):
        sigma = self.sigma
        rho = self.rho
        vonmises = params['vonmises']

        fmax = np.max(vonmises/sigma - 1)

        nlog, nsum, nexp = np.log, np.sum, np.exp
        ks = 1 / rho * nlog(nsum(nexp(rho * (vonmises/sigma - 1 - fmax))))
        unknowns['failure'] = fmax + ks

    def linearize(self, params, unknowns, resids):
        jac = self.alloc_jacobian()
        vonmises = params['vonmises']
        sigma = self.sigma
        rho = self.rho

        # Find the location of the max stress constraint
        fmax = np.max(vonmises / sigma - 1)
        i, j = np.where((vonmises/sigma - 1)==fmax)
        i, j = i[0], j[0]

        # Set incoming seed as 1 so we simply get the jacobian entries
        ksb = 1.

        # Use results from the AD code to compute the jacobian entries
        tempb0 = ksb / (rho * np.sum(np.exp(rho * (vonmises/sigma - fmax - 1))))
        tempb = np.exp(rho*(vonmises/sigma-fmax-1))*rho*tempb0
        fmaxb = ksb - np.sum(tempb)

        # Populate the entries
        derivs = tempb / sigma
        derivs[i, j] += fmaxb / sigma

        # Reshape and save them to the jac dict
        jac['failure', 'vonmises'] = derivs.reshape(1, -1)

        return jac

class SpatialBeamFailureExact(Component):
    """
    Outputs individual failure constraints on each FEM element.

    Parameters
    ----------
    vonmises[ny-1, 2] : numpy array
        von Mises stress magnitudes for each FEM element.

    Returns
    -------
    failure[ny-1, 2] : numpy array
        Array of failure conditions. Positive if element has failed.

    """

    def __init__(self, surface):
        super(SpatialBeamFailureExact, self).__init__()

        self.ny = surface['num_y']

        self.add_param('vonmises', val=np.zeros((self.ny-1, 2), dtype=data_type))
        self.add_output('failure', val=np.zeros((self.ny-1, 2), dtype=data_type))

        self.sigma = surface['stress']

    def solve_nonlinear(self, params, unknowns, resids):
        sigma = self.sigma
        vonmises = params['vonmises']

        unknowns['failure'] = vonmises/sigma - 1

    def linearize(self, params, unknowns, resids):
        return {('failure', 'vonmises') : np.eye(((self.ny-1)*2)) / self.sigma}

class SparWithinWing(Component):
    """

    .. warning::
        This component has not been extensively tested.
        It may require additional coding to work as intended.

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Array defining the nodal points of the lifting surface.
    radius[ny-1] : numpy array
        Radius of each element of the FEM spar.

    Returns
    -------
    spar_within_wing[ny-1] : numpy array
        If all the values are negative, each element is within the wing,
        based on the surface's t_over_c value and the current chord.
        Set a constraint with
        `OASProblem.add_constraint('spar_within_wing', upper=0.)`
    """

    def __init__(self, surface):
        super(SparWithinWing, self).__init__()

        self.surface = surface
        self.ny = surface['num_y']
        self.nx = surface['num_x']

        self.add_param('mesh', val=np.zeros((self.nx, self.ny, 3), dtype=data_type))
        self.add_param('radius', val=np.zeros((self.ny-1), dtype=data_type))
        self.add_output('spar_within_wing', val=np.zeros((self.ny-1), dtype=data_type))

        self.deriv_options['type'] = 'fd'
        self.deriv_options['form'] = 'central'

    def solve_nonlinear(self, params, unknowns, resids):
        mesh = params['mesh']
        max_radius = radii(mesh, self.surface['t_over_c'])
        unknowns['spar_within_wing'] = params['radius'] - max_radius

    def linearize(self, params, unknowns, resids):
        jac = {}
        jac['spar_within_wing', 'radius'] = -np.eye(self.ny-1)
        fd_jac = self.fd_jacobian(params, unknowns, resids,
                                        fd_params=['mesh'],
                                        fd_unknowns=['spar_within_wing'],
                                        fd_states=[])
        jac.update(fd_jac)
        return jac

class NonIntersectingThickness(Component):
    """

    Parameters
    ----------
    thickness[ny-1] : numpy array
        Thickness of each element of the FEM spar.
    radius[ny-1] : numpy array
        Radius of each element of the FEM spar.

    Returns
    -------
    thickness_intersects[ny-1] : numpy array
        If all the values are negative, each element does not intersect itself.
        If a value is positive, then the thickness within the hollow spar
        intersects itself and presents an impossible design.
        Add a constraint as
        `OASProblem.add_constraint('thickness_intersects', upper=0.)`
    """

    def __init__(self, surface):
        super(NonIntersectingThickness, self).__init__()

        self.ny = surface['num_y']
        self.nx = surface['num_x']

        self.add_param('thickness', val=np.zeros((self.ny-1)))
        self.add_param('radius', val=np.zeros((self.ny-1)))
        self.add_output('thickness_intersects', val=np.zeros((self.ny-1)))

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['thickness_intersects'] = 2 * params['thickness'] - params['radius']

    def linearize(self, params, unknowns, resids):
        jac = {}
        jac['thickness_intersects', 'thickness'] = 2 * np.eye(self.ny-1)
        jac['thickness_intersects', 'radius'] = -np.eye(self.ny-1)
        return jac


class SpatialBeamStates(Group):
    """ Group that contains the spatial beam states. """

    def __init__(self, surface):
        super(SpatialBeamStates, self).__init__()

        size = 6 * surface['num_y'] + 6

        self.add('nodes',
                 ComputeNodes(surface),
                 promotes=['*'])
        self.add('assembly',
                 AssembleK(surface),
                 promotes=['*'])
        self.add('fem',
                 SpatialBeamFEM(size),
                 promotes=['*'])
        self.add('disp',
                 SpatialBeamDisp(surface),
                 promotes=['*'])


class SpatialBeamFunctionals(Group):
    """ Group that contains the spatial beam functionals used to evaluate
    performance. """

    def __init__(self, surface):
        super(SpatialBeamFunctionals, self).__init__()

        self.add('energy',
                 SpatialBeamEnergy(surface),
                 promotes=['*'])
        self.add('structural_weight',
                 SpatialBeamWeight(surface),
                 promotes=['*'])
        self.add('vonmises',
                 SpatialBeamVonMisesTube(surface),
                 promotes=['*'])
        self.add('thicknessconstraint',
                 NonIntersectingThickness(surface),
                 promotes=['*'])
        # The following component has not been fully tested so we leave it
        # commented out for now. Use at own risk.
        # self.add('sparconstraint',
        #          SparWithinWing(surface),
        #          promotes=['*'])

        if surface['exact_failure_constraint']:
            self.add('failure',
                     SpatialBeamFailureExact(surface),
                     promotes=['*'])
        else:
            self.add('failure',
                    SpatialBeamFailureKS(surface),
                    promotes=['*'])
