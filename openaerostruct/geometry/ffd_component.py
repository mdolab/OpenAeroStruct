""" Manipulate geometry mesh based on high-level design parameters. """

import numpy as np

import openmdao.api as om
try:
    from pygeo import DVGeometry
except:
    pass


class GeometryMesh(om.ExplicitComponent):
    """
    OpenMDAO component that performs mesh manipulation functions. It reads in
    the initial mesh from the surface dictionary and outputs the altered
    mesh based on the geometric design variables.

    Depending on the design variables selected or the supplied geometry information,
    only some of the follow parameters will actually be given to this component.
    If parameters are not active (they do not deform the mesh), then
    they will not be given to this component.

    Parameters
    ----------

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Modified mesh based on the initial mesh in the surface dictionary and
        the geometric design variables.
    """

    def initialize(self):
        self.options.declare('surface', types=dict)
        self.options.declare('DVGeo', types=DVGeometry)

    def setup(self):
        self.surface = surface = self.options['surface']
        self.mx, self.my = self.surface['mx'], self.surface['my']

        self.DVGeo = self.options['DVGeo']

        # Associate a 'reference axis' for large-scale manipulation
        self.DVGeo.addRefAxis('wing_axis', xFraction=0.25, alignIndex='i')

        # Now add local (shape) variables
        self.DVGeo.addGeoDVLocal('shape', lower=-0.5, upper=0.5, axis='z')

        pts = surface['mesh'].reshape(-1, 3)

        self.DVGeo.addPointSet(pts, 'surface')

        coords = self.DVGeo.getLocalIndex(0)
        self.inds = coords[:, 0, :]
        self.inds2 = coords[:, 1, :]

        self.add_input('shape', val=np.zeros((self.mx, self.my)), units='m')

        self.add_output('mesh', val=surface['mesh'], units='m')

        self.declare_partials('*', '*')

    def compute(self, inputs, outputs):
        surface = self.surface

        dvs = self.DVGeo.getValues()

        for i, row in enumerate(self.inds):
            for j, ind in enumerate(row):
                ind2 = self.inds2[i, j]
                dvs['shape'][ind] = inputs['shape'][i, j]
                dvs['shape'][ind2] = inputs['shape'][i, j]

        self.DVGeo.setDesignVars(dvs)
        coords = self.DVGeo.update('surface')

        mesh = coords.copy()
        mesh = mesh.reshape(surface['mesh'].shape)
        outputs['mesh'] = mesh

    def compute_partials(self, inputs, partials):
        self.DVGeo.computeTotalJacobian('surface')
        jac = self.DVGeo.JT['surface'].toarray().T
        my_jac = partials['mesh', 'shape']
        my_jac[:, :] = 0.

        for i, ind in enumerate(self.inds.flatten()):
            ind2 = self.inds2.flatten()[i]
            my_jac[:, i] += jac[:, ind]
            my_jac[:, i] += jac[:, ind2]

        partials['mesh', 'shape'] = my_jac
