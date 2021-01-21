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

        # TODO BB mirror the groundplane here also
        # negative Z point
        # do I allow dual axis symm?
        groundplane = False

        for surface in surfaces:
            mesh=surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']

            mesh_name = '{}_def_mesh'.format(name)
            vortex_mesh_name = '{}_vortex_mesh'.format(name)

            self.add_input(mesh_name, shape=(nx, ny, 3), units='m')

            # TODO come back and redo
            if surface['groundplane']:
                if not groundplane:
                    # only need to do this once
                    groundplane = True
                    self.add_input('height_agl', val=10., units='m')

            if surface['symmetry']:
                # TODO BB come back and fix the derivatives which are now very wrong
                if surface['groundplane']:
                    self.add_output(vortex_mesh_name, shape=(2*nx, ny*2-1, 3), units='m')
                else:
                    self.add_output(vortex_mesh_name, shape=(nx, ny*2-1, 3), units='m')

                mesh_indices = np.arange(nx * ny * 3).reshape((nx, ny, 3))
                vor_indices = np.arange(nx * (2*ny-1) * 3).reshape((nx, (2*ny-1), 3))

                # original y side, two passes for the quarter-chording
                rows = np.tile(vor_indices[:(nx-1), :ny, :].flatten(), 2)
                # original y side, last row, one pass (on the TE)
                rows = np.hstack((rows, vor_indices[-1  , :ny, :].flatten()))

                # symmetry y side, x and z coordinates only, two passes for quarter-chording, reversed order in y
                rows = np.hstack((rows, np.tile(vor_indices[:(nx-1), ny:, [0, 2]][:, ::-1, :].flatten(), 2)))
                # last row, x and z coordinates only, one pass, reverse order in y
                rows = np.hstack((rows, vor_indices[-1, ny:, [0, 2]].flatten()[::-1]))

                # symmetry y side, y coordinate only, two passes for quarter-chording, reversed order in y
                rows = np.hstack((rows, np.tile(vor_indices[:(nx-1), ny:, 1][:, ::-1].flatten(), 2)))
                # last row, y coord only, one pass,
                # TODO does this need to be in reversed order?
                rows = np.hstack((rows, vor_indices[-1, ny:, 1].flatten()))

                cols = np.concatenate([
                    mesh_indices[:-1, :, :].flatten(),
                    mesh_indices[1:  , :, :].flatten(),
                    mesh_indices[-1  , :, :].flatten(),

                    mesh_indices[:-1, :-1, [0, 2]].flatten(),
                    mesh_indices[1:  , :-1, [0, 2]].flatten(),
                    mesh_indices[-1  , :-1, [0, 2]][::-1, :].flatten(),

                    mesh_indices[:-1, :-1, 1].flatten(),
                    mesh_indices[1:  , :-1, 1].flatten(),
                    mesh_indices[-1  , :-1, 1][::-1].flatten(),
                ])

                data = np.concatenate([
                    0.75 * np.ones((nx-1) * ny * 3),
                    0.25 * np.ones((nx-1) * ny * 3),
                    np.ones(ny * 3),  # back row

                    0.75 * np.ones((nx-1) * (ny-1) * 2),
                    0.25 * np.ones((nx-1) * (ny-1) * 2),
                    np.ones((ny-1) * 2),  # back row

                    -0.75 * np.ones((nx-1) * (ny-1)),
                    -.25  * np.ones((nx-1) * (ny-1)),
                    -np.ones((ny-1)),  # back row
                ])

                self.declare_partials(vortex_mesh_name, mesh_name, val=data, rows=rows, cols=cols)

            else:
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

            mesh_name = '{}_def_mesh'.format(name)
            vortex_mesh_name = '{}_vortex_mesh'.format(name)
            if not surface['groundplane']:
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
                if not surface['symmetry']:
                    raise ValueError('Symmetry and groundplane must be on at the same time')
                # symmetric in y plus ground plane using the first dimension
                mesh = np.zeros((2*nx, ny*2-1, 3), dtype=type(inputs[mesh_name][0, 0, 0]))

                # regular image
                mesh[:nx, :ny, :] = inputs[mesh_name]
                # indices are numbered from tip to centerline
                #  reflection is             all but midpoint  in rev order
                mesh[:nx, ny:, :] = inputs[mesh_name][:, :-1, :][:, ::-1, :]
                mesh[:nx, ny:, 1] *= -1.

                # ground plane image
                mesh[nx:,:,:] = mesh[:nx,:,:]
                mesh[nx:,:,2] *= -1.
                mesh[nx:,:,2] -= 2*inputs['height_agl']

                # from matplotlib import pyplot as plt

                # for i in [ny*2 - 1]:
                #     xs = mesh[:,:i+1,0].flatten()
                #     ys = mesh[:,:i+1,1].flatten()
                #     zs = mesh[:,:i+1,2].flatten()
                #     fig = plt.figure()
                #     ax = fig.add_subplot(111, projection='3d')
                #     ax.scatter(xs, ys, zs)
                #     plt.show()


                outputs[vortex_mesh_name][:nx-1, :, :] = 0.75 * mesh[:nx-1, :, :] + 0.25 * mesh[1:nx, :, :]
                outputs[vortex_mesh_name][nx-1, :, :] = mesh[nx-1, :, :]
                outputs[vortex_mesh_name][nx:-1, :, :] = 0.75 * mesh[nx:-1, :, :] + 0.25 * mesh[nx+1:, :, :]
                outputs[vortex_mesh_name][-1, :, :] = mesh[-1, :, :]
