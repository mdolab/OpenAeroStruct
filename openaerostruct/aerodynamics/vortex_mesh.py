import numpy as np

import openmdao.api as om

class VortexMesh(om.ExplicitComponent):
    """
    Compute the vortex mesh based on the deformed aerodynamic mesh.

    Parameters
    ----------
    def_mesh[nx, ny, 3] : numpy array
        We have a mesh for each lifting surface in the problem.
        That is, if we have both a wing and a tail surface, we will have both
        `wing_def_mesh` and `tail_def_mesh` as inputs.
    height_agl : scalar
        If ground effect is turned on, this input defines the height above
        the groud plane (defined from the origin 0,0,0)
    alpha : scalar
        If ground effect is turned on, this input defines the angular
        rotation of the ground plane

    Returns
    -------
    vortex_mesh[nx, ny, 3] : numpy array
        The actual aerodynamic mesh used in VLM calculations, where we look
        at the rings of the panels instead of the panels themselves. That is,
        this mesh coincides with the quarter-chord panel line, except for the
        final row, where it lines up with the trailing edge.
    """

    def initialize(self):
        self.options.declare('surfaces', types=list)

    def setup(self):
        surfaces = self.options['surfaces']

        # Because the vortex_mesh always comes from the deformed mesh in the
        # same way, the Jacobian is fully linear and can be set here instead
        # of doing compute_partials.
        # We do have to account for symmetry here to create a ghost mesh
        # by mirroring the symmetric mesh.

        any_ground_effect = False

        for surface in surfaces:
            mesh=surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']

            mesh_name = '{}_def_mesh'.format(name)
            vortex_mesh_name = '{}_vortex_mesh'.format(name)

            self.add_input(mesh_name, shape=(nx, ny, 3), units='m')

            ground_effect = surface.get('groundplane', False)


            if ground_effect:
                if not any_ground_effect:
                    # only need to add the extra inputs once
                    any_ground_effect = True
                    self._cached_constant_partial_vals = dict()
                    self.add_input('height_agl', val=8000.0, units='m')
                    self.add_input('alpha', val=0.*np.pi/180, units='rad')

            if surface['symmetry']:
                if ground_effect:
                    self.add_output(vortex_mesh_name, shape=(2*nx, ny*2-1, 3), units='m')
                    # these are cheaper to just do with CS
                    self.declare_partials(vortex_mesh_name, ['alpha', 'height_agl'], method='cs')
                    mesh_indices = np.arange(nx * ny * 3).reshape((nx, ny, 3))
                    vor_indices = np.arange(2*nx * (2*ny-1) * 3).reshape((2*nx, (2*ny-1), 3))
                    quadrant_1_indices = vor_indices[:nx,:ny,:]
                    quadrant_2_indices = vor_indices[:nx,ny:,:]
                    quadrant_3_indices = vor_indices[nx:,:ny,:]
                    quadrant_4_indices = vor_indices[nx:,ny:,:]
                else:
                    # no groundplane
                    self.add_output(vortex_mesh_name, shape=(nx, ny*2-1, 3), units='m')
                    mesh_indices = np.arange(nx * ny * 3).reshape((nx, ny, 3))
                    vor_indices = np.arange(nx * (2*ny-1) * 3).reshape((nx, (2*ny-1), 3))
                    quadrant_1_indices = vor_indices[:nx,:ny,:]
                    quadrant_2_indices = vor_indices[:nx,ny:,:]


                # quadrant 1 is just the baseline mesh
                rows = np.tile(quadrant_1_indices[:-1, :, :].flatten(), 2)
                rows = np.hstack((rows, quadrant_1_indices[-1, :, :].flatten()))
                cols = np.concatenate([
                    mesh_indices[:-1, :, :].flatten(),
                    mesh_indices[1:  , :, :].flatten(),
                    mesh_indices[-1  , :, :].flatten(),
                ])

                data = np.concatenate([
                    0.75 * np.ones((nx-1) * ny * 3),
                    0.25 * np.ones((nx-1) * ny * 3),
                    np.ones(ny * 3), ]) # back row,

                # quadrant 2 is the reflection of the baseline across the midline
                # need to build these piecewise xyz because of the midline reflection
                for dim3 in range(3):
                    rows = np.hstack((rows, np.tile(quadrant_2_indices[:-1, :, dim3].flatten(), 2)))
                    rows = np.hstack((rows, quadrant_2_indices[-1  , :, dim3].flatten()))
                    cols = np.concatenate([cols,
                        mesh_indices[:-1, :-1, dim3][:,::-1].flatten(),
                        mesh_indices[1:  , :-1, dim3][:,::-1].flatten(),
                        mesh_indices[-1  , :-1, dim3][::-1].flatten(),
                    ])

                data = np.concatenate([
                    data,
                    0.75 * np.ones((nx-1) * (ny - 1)),
                    0.25 * np.ones((nx-1) * (ny - 1)),
                    np.ones(ny-1),

                    -0.75 * np.ones((nx-1) * (ny - 1)),
                    -0.25 * np.ones((nx-1) * (ny - 1)),
                    -np.ones(ny-1),

                    0.75 * np.ones((nx-1) * (ny - 1)),
                    0.25 * np.ones((nx-1) * (ny - 1)),
                    np.ones(ny-1),
                ])

                if ground_effect:
                    # these reflections (across the groundplane) are more complex because of the alpha rotation
                    # which means that the x and z points of the reflected mesh depend on BOTH the x and z points of the initial mesh
                    # y only depends on y as usual

                    # third quadrant dependencies (x on x, y on y, z on z, x on z, z on x)
                    list_of_deps = [(0,0),(1,1),(2,2),(0,2),(2,0)]
                    for dep_of, dep_on in list_of_deps:
                        rows = np.hstack((rows, np.tile(quadrant_3_indices[:-1, :, dep_of].flatten(), 2)))
                        rows = np.hstack((rows, quadrant_3_indices[-1  , :, dep_of].flatten()))
                        cols = np.concatenate([cols,
                            mesh_indices[:-1, :, dep_on].flatten(),
                            mesh_indices[1:  , :, dep_on].flatten(),
                            mesh_indices[-1  , :, dep_on].flatten(),
                        ])

                    # fourth quadrant dependencies (x on x, y on y, z on z, x on z, z on x)
                    for dep_of, dep_on in list_of_deps:
                        rows = np.hstack((rows, np.tile(quadrant_4_indices[:-1, :, dep_of].flatten(), 2)))
                        rows = np.hstack((rows, quadrant_4_indices[-1  , :, dep_of].flatten()))
                        cols = np.concatenate([cols,
                            mesh_indices[:-1, :-1, dep_on][:,::-1].flatten(),
                            mesh_indices[1:  , :-1, dep_on][:,::-1].flatten(),
                            mesh_indices[-1  , :-1, dep_on][::-1].flatten(),
                        ])

                    # can't declare constant partials because these depend on alpha (and h?)
                    self.declare_partials(vortex_mesh_name, mesh_name, rows=rows, cols=cols)
                    self._cached_constant_partial_vals[name] = data.copy()

                else:
                    # no groundplane, constant partial values
                    self.declare_partials(vortex_mesh_name, mesh_name, val=data, rows=rows, cols=cols)

            else:
                if ground_effect:
                    raise ValueError('Ground effect is not supported without symmetry turned on')

                self.add_output(vortex_mesh_name, shape=(nx, ny, 3), units='m')

                mesh_indices = np.arange(nx * ny * 3).reshape(
                    (nx, ny, 3))

                rows = np.tile(mesh_indices[:(nx-1), :, :].flatten(), 2)
                rows = np.hstack((rows, mesh_indices[-1  , :, :].flatten()))
                cols = np.concatenate([
                    mesh_indices[:-1, :, :].flatten(),
                    mesh_indices[1:  , :, :].flatten(),
                    mesh_indices[-1  , :, :].flatten(),
                ])

                data = np.concatenate([
                    0.75 * np.ones((nx-1) * ny * 3),
                    0.25 * np.ones((nx-1) * ny * 3),
                    np.ones(ny * 3),  # back row
                ])

                self.declare_partials(vortex_mesh_name, mesh_name, val=data, rows=rows, cols=cols)

    def compute(self, inputs, outputs):
        surfaces = self.options['surfaces']

        for surface in surfaces:
            nx = surface['mesh'].shape[0]
            ny = surface['mesh'].shape[1]
            name = surface['name']
            ground_effect = surface.get('groundplane', False)

            mesh_name = '{}_def_mesh'.format(name)
            vortex_mesh_name = '{}_vortex_mesh'.format(name)
            if not ground_effect:
                if surface['symmetry']:
                    mesh = np.zeros((nx, ny*2-1, 3), dtype=type(inputs[mesh_name][0, 0, 0]))
                    mesh[:, :ny, :] = inputs[mesh_name]
                    # indices are numbered from tip to centerline
                    #  reflection is             all but midpoint  in rev order
                    mesh[:, ny:, :] = inputs[mesh_name][:, :-1, :][:, ::-1, :]
                    mesh[:, ny:, 1] *= -1.
                else:
                    mesh = inputs[mesh_name]

                # all but the last station are moved to the quarterchord point
                outputs[vortex_mesh_name][:-1, :, :] = 0.75 * mesh[:-1, :, :] + 0.25 * mesh[1:, :, :]
                # the last one is coincident
                outputs[vortex_mesh_name][-1, :, :] = mesh[-1, :, :]
            else:
                # symmetric in y plus ground plane using the first dimension
                mesh = np.zeros((2*nx, ny*2-1, 3), dtype=type(inputs[mesh_name][0, 0, 0]))

                # regular image
                mesh[:nx, :ny, :] = inputs[mesh_name]
                # indices are numbered from tip to centerline
                #  reflection is             all but midpoint  in rev order
                mesh[:nx, ny:, :] = inputs[mesh_name][:, :-1, :][:, ::-1, :]
                mesh[:nx, ny:, 1] *= -1.

                # alpha = 5*np.pi/180
                alpha = inputs['alpha'][0]
                plane_normal = np.array([np.sin(alpha), 0.0, -np.cos(alpha)]).reshape((1,1,3))
                plane_point = np.zeros((1,1,3)) + plane_normal * inputs['height_agl']

                # reflect about the ground plane
                # plane is defined parallel to the free stream and height_agl from the origin 0 0 0
                v = mesh[:nx,:,:] - plane_point
                temp = np.inner(v,plane_normal).squeeze()[:,:,np.newaxis]
                v_par = temp*plane_normal
                mesh[nx:,:,:] = mesh[:nx,:,:] - 2*v_par

                outputs[vortex_mesh_name][:nx-1, :, :] = 0.75 * mesh[:nx-1, :, :] + 0.25 * mesh[1:nx, :, :]
                outputs[vortex_mesh_name][nx-1, :, :] = mesh[nx-1, :, :]
                outputs[vortex_mesh_name][nx:-1, :, :] = 0.75 * mesh[nx:-1, :, :] + 0.25 * mesh[nx+1:, :, :]
                outputs[vortex_mesh_name][-1, :, :] = mesh[-1, :, :]


    def compute_partials(self, inputs, J):
        surfaces = self.options['surfaces']
        for surface in surfaces:
            mesh=surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']
            ground_effect = surface.get('groundplane', False)

            mesh_name = '{}_def_mesh'.format(name)
            vortex_mesh_name = '{}_vortex_mesh'.format(name)
            if not ground_effect:
                # if ground effect is not enabled the derivatives are constant
                # and this method need nto be called
                pass
            else:
                data = self._cached_constant_partial_vals[name]
                # we've already figured out the partials for quadrants 1 and 2
                # quandrants 3 and 4 are the ground plane reflections which
                # depend on angle of attack so they need to be computed each time

                # first comes quadrant 3
                # x on x, y on y, z on z, x on z, z on x is the order
                alpha = inputs['alpha']
                x_on_x_const = 1 - 2*np.sin(alpha)**2
                z_on_z_const = 1 - 2*np.cos(alpha)**2
                x_on_z_const = 2*np.sin(alpha)*np.cos(alpha)
                z_on_x_const = 2*np.sin(alpha)*np.cos(alpha)

                data = np.concatenate([
                    data,
                    # x on x
                    x_on_x_const * 0.75 * np.ones((nx-1) * ny),
                    x_on_x_const * 0.25 * np.ones((nx-1) * ny),
                    x_on_x_const * np.ones(ny),

                    # y on y
                    0.75 * np.ones((nx-1) * ny),
                    0.25 * np.ones((nx-1) * ny),
                    np.ones(ny),

                    # z on z
                    z_on_z_const * 0.75 * np.ones((nx-1) * ny),
                    z_on_z_const * 0.25 * np.ones((nx-1) * ny),
                    z_on_z_const * np.ones(ny),

                    # x on z
                    x_on_z_const * 0.75 * np.ones((nx-1) * ny),
                    x_on_z_const * 0.25 * np.ones((nx-1) * ny),
                    x_on_z_const * np.ones(ny),

                    # z on x
                    z_on_x_const * 0.75 * np.ones((nx-1) * ny),
                    z_on_x_const * 0.25 * np.ones((nx-1) * ny),
                    z_on_x_const * np.ones(ny),
                ])

                # now quadrant 4 with different dims and reflected y coords

                data = np.concatenate([
                    data,
                    # x on x
                    x_on_x_const * 0.75 * np.ones((nx-1) * (ny - 1)),
                    x_on_x_const * 0.25 * np.ones((nx-1) * (ny - 1)),
                    x_on_x_const * np.ones((ny - 1)),

                    # y on y
                    -0.75 * np.ones((nx-1) * (ny - 1)),
                    -0.25 * np.ones((nx-1) * (ny - 1)),
                    -np.ones((ny - 1)),

                    # z on z
                    z_on_z_const * 0.75 * np.ones((nx-1) * (ny - 1)),
                    z_on_z_const * 0.25 * np.ones((nx-1) * (ny - 1)),
                    z_on_z_const * np.ones((ny - 1)),

                    # x on z
                    x_on_z_const * 0.75 * np.ones((nx-1) * (ny - 1)),
                    x_on_z_const * 0.25 * np.ones((nx-1) * (ny - 1)),
                    x_on_z_const * np.ones((ny - 1)),

                    # z on x
                    z_on_x_const * 0.75 * np.ones((nx-1) * (ny - 1)),
                    z_on_x_const * 0.25 * np.ones((nx-1) * (ny - 1)),
                    z_on_x_const * np.ones((ny - 1)),
                ])

                J[vortex_mesh_name, mesh_name] = data