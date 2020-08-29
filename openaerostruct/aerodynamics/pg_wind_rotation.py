import numpy as np

import openmdao.api as om


class RotateToWindFrame(om.ExplicitComponent):
    """
    Rotate the VLM Geometry from the standard aerodynamic to the wind frame.
    In the wind frame the freestream will be along the x-axis.

    This transformation is given by the following rotation matrix:

    | x_wind |   | cosb*cosa, -sinb, cosb*sina |   | x_aero |
    | y_wind | = | sinb*cosa,  cosb, sinb*sina | . | y_aero |
    | z_wind |   | -sina,         0,      cosa |   | z_aero |


    Where "a" is the angle of attack and "b" is the sideslip angle.

    Parameters
    ----------
    def_mesh[nx, ny, 3] : numpy array
        Array defining the nodal coordinates of the lifting surface in aero
        frame.
    bound_vecs[num_eval_points, 3] : numpy array
        The vectors representing the bound vortices for each panel in the
        problem.
        This array contains points for all lifting surfaces in the problem.
    coll_pts[num_eval_points, 3] : numpy array
        The xyz coordinates of the collocation points used in the VLM analysis.
        This array contains points for all lifting surfaces in the problem.
    force_pts[num_eval_points, 3] : numpy array
        The xyz coordinates of the force points used in the VLM analysis.
        We evaluate the velocity of the air at these points to get the sectional
        forces acting on the panel. This includes both the freestream and the
        induced velocity acting at these points.
        This array contains points for all lifting surfaces in the problem.
    normals[nx-1, ny-1, 3] : numpy array
        The normal vector for each panel in aero frame, computed as the cross of
        the two diagonals from the mesh points.
    rotational_velocities[num_eval_points, 3] : numpy array
        The rotated freestream velocities at each evaluation point for all
        lifting surfaces.
        This array contains points for all lifting surfaces in the problem.
    alpha : float
        Angle of attack in degrees.
    beta : float
        Sideslip angle in degrees.

    Returns
    -------
    def_mesh_w_frame[nx, ny, 3] : numpy array
        Array defining the nodal coordinates of the lifting surface in wind
        frame.
    bound_vecs_w_frame[num_eval_points, 3] : numpy array
        Bound points for the horseshoe vortices in wind frame.
    coll_pts_w_frame[num_eval_points, 3] : numpy array
        Collocation points on the 3/4 chord line where the flow tangency
        condition is satisfed in wind frame.
    force_pts_w_frame[num_eval_points, 3] : numpy array
        Force pts in wind frame.
    normals_w_frame[nx-1, ny-1, 3] : numpy array
        The normal vector for each panel in wind frame.
    rotational_velocities_w_frame[num_eval_points, 3] : numpy array
        Velocity component at collecation points due to rotational velocity in
        wind frame.
    """

    def initialize(self):
        self.options.declare('surfaces', types=list)
        self.options.declare('rotational', False, types=bool,
                             desc="Set to True to turn on support for computing angular velocities")

    def setup(self):
        surfaces = self.options['surfaces']
        rotational = self.options['rotational']

        # Loop through all the surfaces to determine the total number
        # of evaluation points.
        num_eval_points = 0
        for surface in surfaces:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]

            num_eval_points += (nx - 1) * (ny - 1)

        self.add_input('alpha', val=0., units='rad')
        self.add_input('beta', val=0., units='rad')
        self.add_input('coll_pts', shape=(num_eval_points, 3), units='m')
        self.add_input('force_pts', shape=(num_eval_points, 3), units='m')
        self.add_input('bound_vecs', shape=(num_eval_points, 3), units='m')

        self.add_output('coll_pts_w_frame', shape=(num_eval_points, 3), units='m')
        self.add_output('force_pts_w_frame', shape=(num_eval_points, 3), units='m')
        self.add_output('bound_vecs_w_frame', shape=(num_eval_points, 3), units='m')

        if rotational:
            self.add_input('rotational_velocities', shape=(num_eval_points, 3), units='m/s')
            self.add_output('rotational_velocities_w_frame', shape=(num_eval_points, 3), units='m/s')

        for surface in surfaces:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']

            mesh_name = '{}_def_mesh'.format(name)
            self.add_input(mesh_name, shape=(nx, ny, 3), units='m')

            normals_name = '{}_normals'.format(name)
            self.add_input(normals_name, shape=(nx - 1, ny - 1, 3))

            mesh_name = '{}_def_mesh_w_frame'.format(name)
            self.add_output(mesh_name, shape=(nx, ny, 3), units='m')

            normals_name = '{}_normals_w_frame'.format(name)
            self.add_output(normals_name, shape=(nx - 1, ny - 1, 3))

        # We'll compute all of sensitivities associated with angle of attack and
        # sideslip number through complex-step. Since it's a scalar this is
        # pretty cheap.
        self.declare_partials('*', 'alpha', method='cs')
        self.declare_partials('*', 'beta', method='cs')

        row = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
        col = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])

        rows = np.tile(row, num_eval_points) + np.repeat(3*np.arange(num_eval_points), 9)
        cols = np.tile(col, num_eval_points) + np.repeat(3*np.arange(num_eval_points), 9)

        self.declare_partials('coll_pts_w_frame', 'coll_pts', rows=rows, cols=cols)
        self.declare_partials('force_pts_w_frame', 'force_pts', rows=rows, cols=cols)
        self.declare_partials('bound_vecs_w_frame', 'bound_vecs', rows=rows, cols=cols)

        if rotational:
            self.declare_partials('rotational_velocities_w_frame', 'rotational_velocities', rows=rows, cols=cols)

        for surface in surfaces:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']

            nn = (nx-1) * (ny-1)
            rows = np.tile(row, nn) + np.repeat(3*np.arange(nn), 9)
            cols = np.tile(col, nn) + np.repeat(3*np.arange(nn), 9)

            wrt_name = '{}_normals'.format(name)
            of_name = '{}_normals_w_frame'.format(name)
            self.declare_partials(of_name, wrt_name, rows=rows, cols=cols)

            nn = nx * ny
            rows = np.tile(row, nn) + np.repeat(3*np.arange(nn), 9)
            cols = np.tile(col, nn) + np.repeat(3*np.arange(nn), 9)

            wrt_name = '{}_def_mesh'.format(name)
            of_name = '{}_def_mesh_w_frame'.format(name)
            self.declare_partials(of_name, wrt_name, rows=rows, cols=cols)

    def compute(self, inputs, outputs):
        rotational = self.options['rotational']

        alpha = inputs['alpha']
        beta = inputs['beta']

        cosa = np.cos(alpha)
        sina = np.sin(alpha)
        cosb = np.cos(beta)
        sinb = np.sin(beta)

        # Define aero->wind rotation matrix
        Tw = np.array([[cosb*cosa, -sinb, cosb*sina],
                       [sinb*cosa,  cosb, sinb*sina],
                       [-sina,         0, cosa     ]], alpha.dtype)

        # Use einsum for fast vectorized matrix multiplication
        outputs['bound_vecs_w_frame'] = np.einsum('lk,jk->jl', Tw, inputs['bound_vecs'])
        outputs['coll_pts_w_frame'] = np.einsum('lk,jk->jl', Tw, inputs['coll_pts'])
        outputs['force_pts_w_frame'] = np.einsum('lk,jk->jl', Tw, inputs['force_pts'])

        if rotational:
            outputs['rotational_velocities_w_frame'] = \
                np.einsum('lk,jk->jl', Tw, inputs['rotational_velocities'])

        for surface in self.options['surfaces']:
            name = surface['name']

            wrt_name = '{}_def_mesh'.format(name)
            of_name = '{}_def_mesh_w_frame'.format(name)
            outputs[of_name] = np.einsum('lk,ijk->ijl', Tw, inputs[wrt_name])

            wrt_name = '{}_normals'.format(name)
            of_name = '{}_normals_w_frame'.format(name)
            outputs[of_name] = np.einsum('lk,ijk->ijl', Tw, inputs[wrt_name])

    def compute_partials(self, inputs, partials):
        rotational = self.options['rotational']
        alpha = inputs['alpha']
        beta = inputs['beta']

        cosa = np.cos(alpha)
        sina = np.sin(alpha)
        cosb = np.cos(beta)
        sinb = np.sin(beta)

        num_eval_pts = inputs['bound_vecs'].shape[0]

        # Define aero->wind rotation matrix
        Tw = np.array([[cosb*cosa, -sinb, cosb*sina],
                       [sinb*cosa,  cosb, sinb*sina],
                       [-sina,         0, cosa     ]], alpha.dtype).flatten()


        partials['coll_pts_w_frame', 'coll_pts'] = np.tile(Tw, num_eval_pts)
        partials['force_pts_w_frame', 'force_pts'] = np.tile(Tw, num_eval_pts)
        partials['bound_vecs_w_frame', 'bound_vecs'] = np.tile(Tw, num_eval_pts)

        if rotational:
            partials['rotational_velocities_w_frame', 'rotational_velocities'] = np.tile(Tw, num_eval_pts)

        for surface in self.options['surfaces']:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']

            wrt_name = '{}_normals'.format(name)
            of_name = '{}_normals_w_frame'.format(name)
            nn = (nx-1) * (ny-1)
            partials[of_name, wrt_name] = np.tile(Tw, nn)

            wrt_name = '{}_def_mesh'.format(name)
            of_name = '{}_def_mesh_w_frame'.format(name)
            nn = nx * ny
            partials[of_name, wrt_name] = np.tile(Tw, nn)


class RotateFromWindFrame(om.ExplicitComponent):
    """
    Rotate the aerodynamic sectional and nodal force vectors from the wind to
    the standard aerodynamic frame. This is the reverse operation of the
    RotateToWindFrame component.

    This transformation is given by the following rotation matrix:

    | F_x_aero |   | cosb*cosa, sinb*cosa, -sina |   | F_x_wind |
    | F_y_aero | = | -sinb,          cosb,     0 | . | F_y_wind |
    | F_z_aero |   | cosb*sina, sinb*sina,  cosa |   | F_z_wind |


    Where "a" is the angle of attack and "b" is the sideslip angle.

    Parameters
    ----------
    sec_forces_w_frame[nx-1, ny-1, 3] : numpy array
        Force vectors on each panel (lattice) in wind frame.

    alpha : float
        Angle of attack in degrees.
    beta : float
        Sideslip angle in degrees.

    Returns
    -------
    sec_forces[nx-1, ny-1, 3] : numpy array
        Force vectors on each panel (lattice) in aero frame.
    """

    def initialize(self):
        self.options.declare('surfaces', types=list)

    def setup(self):
        self.add_input('alpha', val=0., units='rad')
        self.add_input('beta', val=0., units='rad')

        # We'll compute all of sensitivities associated with angle of attack and
        # sideslip number through complex-step. Since it's a scalar this is
        # pretty cheap.
        self.declare_partials('*', 'alpha', method='cs')
        self.declare_partials('*', 'beta', method='cs')

        for surface in self.options['surfaces']:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']

            wrt_name = '{}_sec_forces_w_frame'.format(name)
            of_name = '{}_sec_forces'.format(name)

            self.add_input(wrt_name, val=np.zeros((nx-1, ny-1, 3)), units='N')
            self.add_output(of_name, val=np.zeros((nx-1, ny-1, 3)), units='N')

            row = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
            col = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])

            nn = (nx-1) * (ny-1)
            rows = np.tile(row, nn) + np.repeat(3*np.arange(nn), 9)
            cols = np.tile(col, nn) + np.repeat(3*np.arange(nn), 9)

            self.declare_partials(of_name, wrt_name, rows=rows, cols=cols)

    def compute(self, inputs, outputs):
        alpha = inputs['alpha']
        beta = inputs['beta']

        cosa = np.cos(alpha)
        sina = np.sin(alpha)
        cosb = np.cos(beta)
        sinb = np.sin(beta)

        # Define aero->wind rotation matrix
        # wind->aero rotation matrix is given by transpose
        Tw = np.array([[cosb*cosa, -sinb, cosb*sina],
                       [sinb*cosa,  cosb, sinb*sina],
                       [-sina,         0, cosa     ]], alpha.dtype).T

        # Use einsum for fast vectorized matrix multiplication
        for surface in self.options['surfaces']:
            name = surface['name']
            wrt_name = '{}_sec_forces_w_frame'.format(name)
            of_name = '{}_sec_forces'.format(name)

            outputs[of_name] = np.einsum('lk,ijk->ijl', Tw, inputs[wrt_name])

    def compute_partials(self, inputs, partials):
        alpha = inputs['alpha']
        beta = inputs['beta']

        cosa = np.cos(alpha)
        sina = np.sin(alpha)
        cosb = np.cos(beta)
        sinb = np.sin(beta)

        # Define aero->wind rotation matrix
        Tw = np.array([[cosb*cosa, -sinb, cosb*sina],
                       [sinb*cosa,  cosb, sinb*sina],
                       [-sina,         0, cosa     ]], alpha.dtype).T.flatten()

        for surface in self.options['surfaces']:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']

            wrt_name = '{}_sec_forces_w_frame'.format(name)
            of_name = '{}_sec_forces'.format(name)

            nn = (nx-1) * (ny-1)
            partials[of_name, wrt_name] = np.tile(Tw, nn)
