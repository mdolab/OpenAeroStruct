import numpy as np

import openmdao.api as om


class ScaleToPrandtlGlauert(om.ExplicitComponent):
    """
    Scale the wind frame coordinates to get the Prandtl-Glauert transformed
    geometry.

    The Prandtl glauert transformation is defined as below:
        Coordinates
        x_pg = x_wind
        y_pg = B*y_wind
        z_pg = B*z_wind

        Normals
        n_x_pg = B*n_x_wind
        n_y_pg = n_y_wind
        n_z_pg = n_z_wind

        Perturbation velocities
        v_x_pg = B^2*v_x_wind
        v_y_pg = B*v_y_wind
        v_z_pg = B*v_z_wind

    where B = sqrt(1 - M^2).
    Note: The freestream velocity remains untransformed and therefore is not
    included in this component.

    Parameters
    ----------
    def_mesh_w_frame[nx, ny, 3] : numpy array
        Array defining the nodal coordinates of the lifting surface in aero
        frame.
    bound_vecs_w_frame[num_eval_points, 3] : numpy array
        The vectors representing the bound vortices for each panel in the
        problem.
        This array contains points for all lifting surfaces in the problem.
    coll_pts_w_frame[num_eval_points, 3] : numpy array
        The xyz coordinates of the collocation points used in the VLM analysis.
        This array contains points for all lifting surfaces in the problem.
    force_pts_w_frame[num_eval_points, 3] : numpy array
        The xyz coordinates of the force points used in the VLM analysis.
        We evaluate the velocity of the air at these points to get the sectional
        forces acting on the panel. This includes both the freestream and the
        induced velocity acting at these points.
        This array contains points for all lifting surfaces in the problem.
    normals_w_frame[nx-1, ny-1, 3] : numpy array
        The normal vector for each panel in aero frame, computed as the cross of
        the two diagonals from the mesh points.
    rotational_velocities_w_frame[num_eval_points, 3] : numpy array
        The rotated freestream velocities at each evaluation point for all
        lifting surfaces.
        This array contains points for all lifting surfaces in the problem.

    M : float
        Freestream Mach number.

    Returns
    -------
    def_mesh_pg[nx, ny, 3] : numpy array
        Array defining the nodal coordinates of the lifting surface in PG frame.
    bound_vecs_pg[num_eval_points, 3] : numpy array
        Bound points in PG frame.
    coll_pts_pg[num_eval_points, 3] : numpy array
        Collocation points in PG frame.
    force_pts_pg[num_eval_points, 3] : numpy array
        Force points in PG frame.
    normals_pg[nx-1, ny-1, 3] : numpy array
        The normal vector for each panel in PG frame.
    rotational_velocities_pg[num_eval_points, 3] : numpy array
        Velocity component at collocation points due to rotational velocity in PG frame.
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

        self.add_input('Mach_number', val=0.)

        self.add_input('coll_pts_w_frame', shape=(num_eval_points, 3), units='m')
        self.add_input('force_pts_w_frame', shape=(num_eval_points, 3), units='m')
        self.add_input('bound_vecs_w_frame', shape=(num_eval_points, 3), units='m')

        self.add_output('coll_pts_pg', shape=(num_eval_points, 3), units='m')
        self.add_output('force_pts_pg', shape=(num_eval_points, 3), units='m')
        self.add_output('bound_vecs_pg', shape=(num_eval_points, 3), units='m')

        if rotational:
            self.add_input('rotational_velocities_w_frame', shape=(num_eval_points, 3), units='m/s')
            self.add_output('rotational_velocities_pg', shape=(num_eval_points, 3), units='m/s')

        for surface in surfaces:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']

            mesh_name = '{}_def_mesh_w_frame'.format(name)
            self.add_input(mesh_name, shape=(nx, ny, 3), units='m')

            normals_name = '{}_normals_w_frame'.format(name)
            self.add_input(normals_name, shape=(nx - 1, ny - 1, 3))

            mesh_name = '{}_def_mesh_pg'.format(name)
            self.add_output(mesh_name, shape=(nx, ny, 3), units='m')

            normals_name = '{}_normals_pg'.format(name)
            self.add_output(normals_name, shape=(nx - 1, ny - 1, 3))

        # We'll compute all of sensitivities associated with Mach number through
        # complex-step. Since it's a scalar this is pretty cheap.
        self.declare_partials('*', 'Mach_number', method='cs')

        row_col = np.arange(num_eval_points*3)

        self.declare_partials('coll_pts_pg', 'coll_pts_w_frame', rows=row_col, cols=row_col)
        self.declare_partials('force_pts_pg', 'force_pts_w_frame', rows=row_col, cols=row_col)
        self.declare_partials('bound_vecs_pg', 'bound_vecs_w_frame', rows=row_col, cols=row_col)

        if rotational:
            self.declare_partials('rotational_velocities_pg', 'rotational_velocities_w_frame', rows=row_col, cols=row_col)

        for surface in surfaces:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']

            nn = (nx-1) * (ny-1)
            row_col = np.arange(3*nn)

            wrt_name = '{}_normals_w_frame'.format(name)
            of_name = '{}_normals_pg'.format(name)
            self.declare_partials(of_name, wrt_name, rows=row_col, cols=row_col)

            nn = nx * ny
            row_col = np.arange(3*nn)

            wrt_name = '{}_def_mesh_w_frame'.format(name)
            of_name = '{}_def_mesh_pg'.format(name)
            self.declare_partials(of_name, wrt_name, rows=row_col, cols=row_col)

    def compute(self, inputs, outputs):
        rotational = self.options['rotational']

        M = inputs['Mach_number']
        betaPG = np.sqrt(1 - M**2)

        outputs['bound_vecs_pg'] = inputs['bound_vecs_w_frame']
        outputs['bound_vecs_pg'][:, 1] *= betaPG
        outputs['bound_vecs_pg'][:, 2] *= betaPG

        outputs['coll_pts_pg'] = inputs['coll_pts_w_frame']
        outputs['coll_pts_pg'][:, 1] *= betaPG
        outputs['coll_pts_pg'][:, 2] *= betaPG

        outputs['force_pts_pg'] = inputs['force_pts_w_frame']
        outputs['force_pts_pg'][:, 1] *= betaPG
        outputs['force_pts_pg'][:, 2] *= betaPG

        if rotational:
            outputs['rotational_velocities_pg'] = inputs['rotational_velocities_w_frame']
            outputs['rotational_velocities_pg'][:, 0] *= betaPG**2
            outputs['rotational_velocities_pg'][:, 1] *= betaPG
            outputs['rotational_velocities_pg'][:, 2] *= betaPG

        for surface in self.options['surfaces']:
            name = surface['name']

            wrt_name = '{}_def_mesh_w_frame'.format(name)
            of_name = '{}_def_mesh_pg'.format(name)

            outputs[of_name] = inputs[wrt_name]
            outputs[of_name][:, :, 1] *= betaPG
            outputs[of_name][:, :, 2] *= betaPG

            wrt_name = '{}_normals_w_frame'.format(name)
            of_name = '{}_normals_pg'.format(name)

            outputs[of_name] = inputs[wrt_name]
            outputs[of_name][:, :, 0] *= betaPG

    def compute_partials(self, inputs, partials):
        rotational = self.options['rotational']

        M = inputs['Mach_number']
        betaPG = np.sqrt(1 - M**2)
        fact = np.array([1.0, betaPG, betaPG], M.dtype).flatten()
        fact_norm = np.array([betaPG, 1.0, 1.0], M.dtype).flatten()
        num_eval_pts = inputs['bound_vecs_w_frame'].shape[0]

        partials['bound_vecs_pg', 'bound_vecs_w_frame'] = np.tile(fact, num_eval_pts)
        partials['coll_pts_pg', 'coll_pts_w_frame'] = np.tile(fact, num_eval_pts)
        partials['force_pts_pg', 'force_pts_w_frame'] = np.tile(fact, num_eval_pts)

        if rotational:
            fact_rot = np.array([betaPG**2, betaPG, betaPG], M.dtype).flatten()
            partials['rotational_velocities_pg', 'rotational_velocities_w_frame'] = np.tile(fact_rot, num_eval_pts)

        for surface in self.options['surfaces']:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']

            wrt_name = '{}_def_mesh_w_frame'.format(name)
            of_name = '{}_def_mesh_pg'.format(name)
            nn = nx * ny
            partials[of_name, wrt_name] = np.tile(fact, nn)

            wrt_name = '{}_normals_w_frame'.format(name)
            of_name = '{}_normals_pg'.format(name)
            nn = (nx-1) * (ny-1)
            partials[of_name, wrt_name] = np.tile(fact_norm, nn)


class ScaleFromPrandtlGlauert(om.ExplicitComponent):
    """
    Scale the Prandtl-Glauert transformed forces to get the physical forces
    Prandtl-Glauert transformed geometry.

    The inverse Prandtl-Glauert transformation for forces is defined as below:
        F_x_wind = F_x_pg/B^4
        F_y_wind = F_y_pg/B^3
        F_z_wind = F_z_pg/B^3

    where B = sqrt(1 - M^2).

    Parameters
    ----------
    sec_forces_pg[nx-1, ny-1, 3] : numpy array
        Force vectors on each panel (lattice) in PG domain.

    M : float
        Freestream Mach number.

    Returns
    -------
    sec_forces_w_frame[nx-1, ny-1, 3] : numpy array
        Force vectors on each panel (lattice) in wind frame.
    """

    def initialize(self):
        self.options.declare('surfaces', types=list)

    def setup(self):

        self.add_input('Mach_number', val=0.)

        # We'll compute all of sensitivities associated with Mach number through
        # complex-step. Since it's a scalar this is pretty cheap.
        self.declare_partials('*', 'Mach_number', method='cs')

        for surface in self.options['surfaces']:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']

            wrt_name = '{}_sec_forces_pg'.format(name)
            of_name = '{}_sec_forces_w_frame'.format(name)

            self.add_input(wrt_name, val=np.zeros((nx-1, ny-1, 3)), units='N')
            self.add_output(of_name, val=np.zeros((nx-1, ny-1, 3)), units='N')

            nn = (nx-1) * (ny-1)
            row_col = np.arange(3*nn)

            self.declare_partials(of_name, wrt_name, rows=row_col, cols=row_col)

    def compute(self, inputs, outputs):
        M = inputs['Mach_number']
        betaPG = np.sqrt(1 - M**2)

        for surface in self.options['surfaces']:
            name = surface['name']
            wrt_name = '{}_sec_forces_pg'.format(name)
            of_name = '{}_sec_forces_w_frame'.format(name)

            outputs[of_name] = inputs[wrt_name]
            outputs[of_name][:, :, 0] *= (1.0/betaPG**4)
            outputs[of_name][:, :, 1] *= (1.0/betaPG**3)
            outputs[of_name][:, :, 2] *= (1.0/betaPG**3)

    def compute_partials(self, inputs, partials):
        M = inputs['Mach_number']
        betaPG = np.sqrt(1 - M**2)

        term1 = 1.0/betaPG**4
        term2 = 1.0/betaPG**3
        fact = np.array([term1, term2, term2], M.dtype).flatten()

        for surface in self.options['surfaces']:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']

            wrt_name = '{}_sec_forces_pg'.format(name)
            of_name = '{}_sec_forces_w_frame'.format(name)

            nn = (nx-1) * (ny-1)
            partials[of_name, wrt_name] = np.tile(fact, nn)
