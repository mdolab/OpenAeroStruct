import numpy as np

import openmdao.api as om


class RotationalVelocity(om.ExplicitComponent):
    """
    Compute the velocity due to rigid body rotation.

    Parameters
    ----------
    omega[3] : ndarray
        Angular velocity vector for each surface about center of gravity.
        Only available if the rotational options is set to True.
    cg[3] : ndarray
        The x, y, z coordinates of the center of gravity for the entire aircraft.
        Only available if the rotational options is set to True.
    coll_pts[num_eval_points, 3] : ndarray
        The xyz coordinates of the collocation points used in the VLM analysis.
        This array contains points for all lifting surfaces in the problem.
        Only available if the rotational options is set to True.

    Returns
    -------
    rotational_velocities[num_eval_points, 3] : numpy array
        The rotated freestream velocities at each evaluation point for all
        lifting surfaces.
        This array contains points for all lifting surfaces in the problem.
    """

    def initialize(self):
        self.options.declare('surfaces', types=list)

    def setup(self):
        surfaces = self.options['surfaces']

        system_size = 0
        sizes = []

        # Loop through each surface and cumulatively add the number of panels
        # to obtain system_size.
        for surface in surfaces:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            size = (nx - 1) * (ny - 1)
            system_size += size
            sizes.append(size)

        self.system_size = system_size

        self.add_input('coll_pts', shape=(system_size, 3), units='m')
        self.add_input('omega', val=np.zeros((3, )), units='rad/s')
        self.add_input('cg', val=np.ones((3, )), units='m')

        self.add_output('rotational_velocities', shape=(system_size, 3), units='m/s')

        # First Half of cross product
        row = np.array([1, 2, 0])
        col = np.array([2, 0, 1])

        rows1 = np.tile(row, system_size) + np.repeat(3*np.arange(system_size), 3)
        cols1 = np.tile(col, system_size)

        # Second Half of cross product
        rows2 = np.tile(col, system_size) + np.repeat(3*np.arange(system_size), 3)
        cols2 = np.tile(row, system_size)

        rows = np.concatenate([rows1, rows2])
        cols = np.concatenate([cols1, cols2])

        self.declare_partials('rotational_velocities', 'cg', rows=rows, cols=cols)
        self.declare_partials('rotational_velocities', 'omega', rows=rows, cols=cols)

        cols1 = np.tile(col, system_size) + np.repeat(3*np.arange(system_size), 3)
        cols2 = np.tile(row, system_size) + np.repeat(3*np.arange(system_size), 3)
        cols = np.concatenate([cols1, cols2])

        self.declare_partials('rotational_velocities', 'coll_pts', rows=rows, cols=cols)

    def compute(self, inputs, outputs):
        # Angular velocity term
        cg = inputs['cg']
        omega = inputs['omega']
        c_pts = inputs['coll_pts']

        for j in np.arange(c_pts.shape[0]):
            r = c_pts[j, :] - cg

            outputs['rotational_velocities'][j, :] = np.cross(omega, r)

    def compute_partials(self, inputs, J):
        cg = inputs['cg']
        omega = inputs['omega']
        c_pts = inputs['coll_pts']

        surfaces = self.options['surfaces']
        idx = jdx = 0
        ii = self.system_size * 3
        for j, surface in enumerate(surfaces):
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            size = (nx - 1) * (ny - 1)

            r = c_pts[jdx:jdx+size, :] - cg

            # Cross product derivatives organized so we can tile a variable directly into slices

            J['rotational_velocities', 'cg'][idx:idx+size*3] = np.tile(omega, size)
            J['rotational_velocities', 'cg'][idx+ii:idx+ii+size*3] = -np.tile(omega, size)

            J['rotational_velocities', 'coll_pts'][idx:idx+size*3] = -np.tile(omega, size)
            J['rotational_velocities', 'coll_pts'][idx+ii:idx+ii+size*3] = np.tile(omega, size)

            J['rotational_velocities', 'omega'][idx:idx+size*3] = r.flatten()
            J['rotational_velocities', 'omega'][idx+ii:idx+ii+size*3] = -r.flatten()

            idx += 3*size
            jdx += size
