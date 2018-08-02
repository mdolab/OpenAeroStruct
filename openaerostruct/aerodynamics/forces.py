from __future__ import print_function, division
import numpy as np

from openmdao.api import ExplicitComponent

from openaerostruct.aerodynamics.utils import _assemble_AIC_mtx, _assemble_AIC_mtx_b, _assemble_AIC_mtx_d

try:
    from openaerostruct.fortran import OAS_API
    fortran_flag = True
    data_type = float
except:
    fortran_flag = False
    data_type = complex

compressible = False

def view_mat(mat):
    """ Helper function used to visually examine matrices. """
    import matplotlib.pyplot as plt
    if len(mat.shape) > 2:
        mat = np.sum(mat, axis=2)
    im = plt.imshow(mat.real, interpolation='none')
    plt.colorbar(im, orientation='horizontal')
    plt.show()

class Forces(ExplicitComponent):
    """ Compute aerodynamic forces acting on each section.

    Note that the first two inputs and the output have the surface name
    prepended on it. E.g., 'def_mesh' on a surface called 'wing' would be
    'wing.def_mesh', etc.

    Parameters
    ----------
    def_mesh[nx, ny, 3] : numpy array
        Array defining the nodal coordinates of the lifting surface.
    b_pts[nx-1, ny, 3] : numpy array
        Bound points for the horseshoe vortices, found along the 1/4 chord.

    circulations[(nx-1)*(ny-1)] : numpy array
        Flattened vector of horseshoe vortex strengths calculated by solving
        the linear system of AIC_mtx * circulations = rhs, where rhs is
        based on the air velocity at each collocation point.
    alpha : float
        Angle of attack in degrees.
    v : float
        Freestream air velocity in m/s.
    rho : float
        Air density in kg/m^3.

    Returns
    -------
    sec_forces[nx-1, ny-1, 3] : numpy array
        Contains the sectional forces acting on each panel.
        Stored in Fortran order (only relevant with more than one chordwise
        panel).

    """

    def initialize(self):
        self.options.declare('surfaces', types=list)

    def setup(self):
        self.surfaces = surfaces = self.options['surfaces']

        tot_panels = 0
        for surface in surfaces:
            name = surface['name']
            ny = surface['num_y']
            nx = surface['num_x']
            tot_panels += (nx - 1) * (ny - 1)

            self.add_input(name + '_def_mesh', val=np.zeros((nx, ny, 3)), units='m')#, dtype=data_type))
            self.add_input(name + '_b_pts', val=np.zeros((nx-1, ny, 3)), units='m')#, dtype=data_type))
            self.add_input(name + '_cos_sweep', val=np.zeros((ny-1)), units='m')#, dtype=data_type))
            self.add_input(name + '_widths', val=np.ones((ny-1)), units='m')#, dtype=data_type))

            self.add_output(name + '_sec_forces', val=np.zeros((nx-1, ny-1, 3)), units='N')#, dtype=data_type))

        self.tot_panels = tot_panels

        self.add_input('circulations', val=np.zeros((tot_panels)), units='m**2/s')
        self.add_input('alpha', val=3.)
        self.add_input('M', val=0.1)
        self.add_input('v', val=10., units='m/s')
        self.add_input('rho', val=3., units='kg/m**3')

        self.mtx = np.zeros((tot_panels, tot_panels, 3), dtype=data_type)
        self.v = np.zeros((tot_panels, 3), dtype=data_type)

        self.declare_partials('*', '*')

        if not fortran_flag:
            self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        circ = inputs['circulations']
        alpha = inputs['alpha'] * np.pi / 180.
        rho = inputs['rho']
        cosa = np.cos(alpha)
        sina = np.sin(alpha)

        # that the collocation points used here are the midpoints of each
        # bound vortex filament, not the collocation points from above
        _assemble_AIC_mtx(self.mtx, inputs, self.surfaces, skip=True)

        # Compute the induced velocities at the midpoints of the
        # bound vortex filaments
        for ind in range(3):
            self.v[:, ind] = self.mtx[:, :, ind].dot(circ)

        # Add the freestream velocity to the induced velocity so that
        # self.v is the total velocity seen at the point
        self.v[:, 0] += cosa * inputs['v']
        self.v[:, 2] += sina * inputs['v']

        i = 0
        for surface in self.surfaces:
            name = surface['name']
            nx = surface['num_x']
            ny = surface['num_y']
            num_panels = (nx - 1) * (ny - 1)

            b_pts = inputs[name + '_b_pts']

            if fortran_flag:
                sec_forces = OAS_API.oas_api.forcecalc(self.v[i:i+num_panels, :], circ[i:i+num_panels], rho, b_pts)
            else:

                # Get the vectors for each bound vortex of the horseshoe vortices
                bound = b_pts[:, 1:, :] - b_pts[:, :-1, :]

                # Cross the obtained velocities with the bound vortex filament
                # vectors
                cross = np.cross(self.v[i:i+num_panels],
                                    bound.reshape(-1, bound.shape[-1], order='F'))

                sec_forces = np.zeros((num_panels, 3), dtype=data_type)
                # Compute the sectional forces acting on each panel
                for ind in range(3):
                    sec_forces[:, ind] = \
                        (inputs['rho'] * circ[i:i+num_panels] * cross[:, ind])

            # Reshape the forces into the expected form
            forces = sec_forces.reshape((nx-1, ny-1, 3), order='F')

            sweep_angle = inputs[name + '_cos_sweep'] / inputs[name + '_widths']
            beta = np.sqrt(1 - inputs['M']**2 * sweep_angle**2)

            if not compressible:
                beta[:] = 1.

            for j, B in enumerate(beta):
                outputs[name + '_sec_forces'][:, j, :] = forces[:, j, :] / B

            i += num_panels


    if fortran_flag:
        if 0:
            def compute_jacvec_product(self, inputs, outputs, d_inputs, d_outputs, mode):
                circ = inputs['circulations']
                alpha = inputs['alpha'] * np.pi / 180.
                rho = inputs['rho']
                cosa = np.cos(alpha)
                sina = np.sin(alpha)

                # Assemble a different matrix here than the AIC_mtx from above; Note
                # that the collocation points used here are the midpoints of each
                # bound vortex filament, not the collocation points from above
                _assemble_AIC_mtx(self.mtx, inputs, self.surfaces, skip=True)

                # Compute the induced velocities at the midpoints of the
                # bound vortex filaments
                for ind in range(3):
                    self.v[:, ind] = self.mtx[:, :, ind].dot(circ)

                # Add the freestream velocity to the induced velocity so that
                # self.v is the total velocity seen at the point
                self.v[:, 0] += cosa * inputs['v']
                self.v[:, 2] += sina * inputs['v']

                if mode == 'fwd':

                    circ = inputs['circulations']
                    alpha = inputs['alpha'] * np.pi / 180.
                    if 'alpha' in d_inputs:
                        alphad = d_inputs['alpha'] * np.pi / 180.
                    else:
                        alphad = 0.

                    if 'circulations' in d_inputs:
                        circ_d = d_inputs['circulations']
                    else:
                        circ_d = np.zeros(circ.shape)
                    cosa = np.cos(alpha)
                    sina = np.sin(alpha)
                    cosad = -sina * alphad
                    sinad = cosa * alphad
                    rho = inputs['rho']
                    v = inputs['v']

                    mtxd = np.zeros(self.mtx.shape)

                    # Actually assemble the AIC matrix
                    _assemble_AIC_mtx_d(mtxd, inputs, d_inputs, self.surfaces, skip=True)

                    vd = np.zeros(self.v.shape)

                    # Compute the induced velocities at the midpoints of the
                    # bound vortex filaments
                    for ind in range(3):
                        vd[:, ind] += mtxd[:, :, ind].dot(circ)
                        vd[:, ind] += self.mtx[:, :, ind].real.dot(circ_d)

                    # Add the freestream velocity to the induced velocity so that
                    # self.v is the total velocity seen at the point
                    if 'v' in d_inputs:
                        v_d = d_inputs['v']
                    else:
                        v_d = 0.
                    vd[:, 0] += cosa * v_d
                    vd[:, 2] += sina * v_d
                    vd[:, 0] += cosad * v
                    vd[:, 2] += sinad * v

                    if 'rho' in d_inputs:
                        rho_d = d_inputs['rho']
                    else:
                        rho_d = 0.

                    i = 0
                    rho = inputs['rho'].real
                    for surface in self.surfaces:
                        name = surface['name']
                        nx = surface['num_x']
                        ny = surface['num_y']

                        num_panels = (nx - 1) * (ny - 1)

                        b_pts = inputs[name + '_b_pts']
                        if name+'_b_pts' in d_inputs:
                            b_pts_d = d_inputs[name + '_b_pts']
                        else:
                            b_pts_d = np.zeros(b_pts.shape)

                        self.compute(inputs, outputs)

                        sec_forces = outputs[name + '_sec_forces'].real

                        sec_forces, sec_forcesd = OAS_API.oas_api.forcecalc_d(self.v[i:i+num_panels, :], vd[i:i+num_panels],
                                                    circ[i:i+num_panels], circ_d[i:i+num_panels],
                                                    rho, rho_d,
                                                    b_pts, b_pts_d)

                        d_outputs[name + '_sec_forces'] += sec_forcesd.reshape((nx-1, ny-1, 3), order='F')
                        i += num_panels


                if mode == 'rev':

                    circ = inputs['circulations']
                    alpha = inputs['alpha'] * np.pi / 180.
                    cosa = np.cos(alpha)
                    sina = np.sin(alpha)

                    i = 0
                    rho = inputs['rho'].real
                    v = inputs['v']
                    vb = np.zeros(self.v.shape)

                    for surface in self.surfaces:
                        name = surface['name']
                        nx = surface['num_x']
                        ny = surface['num_y']
                        num_panels = (nx - 1) * (ny - 1)

                        b_pts = inputs[name + '_b_pts']

                        sec_forcesb = d_outputs[name + '_sec_forces'].reshape((num_panels, 3), order='F')

                        v_b, circb, rhob, bptsb, _ = OAS_API.oas_api.forcecalc_b(self.v[i:i+num_panels, :], circ[i:i+num_panels], rho, b_pts, sec_forcesb)

                        if 'circulations' in d_inputs:
                            d_inputs['circulations'][i:i+num_panels] += circb
                        vb[i:i+num_panels] = v_b
                        if 'rho' in d_inputs:
                            d_inputs['rho'] += rhob
                        if name + '_b_pts' in d_inputs:
                            d_inputs[name + '_b_pts'] += bptsb

                        i += num_panels

                    sinab = inputs['v'] * np.sum(vb[:, 2])
                    if 'v' in d_inputs:
                        d_inputs['v'] += cosa * np.sum(vb[:, 0]) + sina * np.sum(vb[:, 2])
                    cosab = inputs['v'] * np.sum(vb[:, 0])
                    ab = np.cos(alpha) * sinab - np.sin(alpha) * cosab
                    if 'alpha' in d_inputs:
                        d_inputs['alpha'] += np.pi * ab / 180.

                    mtxb = np.zeros(self.mtx.shape)
                    circb = np.zeros(circ.shape)
                    for i in range(3):
                        for j in range(self.tot_panels):
                            mtxb[j, :, i] += circ * vb[j, i]
                            circb += self.mtx[j, :, i].real * vb[j, i]

                    if 'circulations' in d_inputs:
                        d_inputs['circulations'] += circb

                    _assemble_AIC_mtx_b(mtxb, inputs, d_inputs, self.surfaces, skip=True)

        else:
            def compute_partials(self, inputs, partials):
                circ = inputs['circulations']
                alpha = inputs['alpha'] * np.pi / 180.
                rho = inputs['rho']
                cosa = np.cos(alpha)
                sina = np.sin(alpha)

                # Assemble a different matrix here than the AIC_mtx from above; Note
                # that the collocation points used here are the midpoints of each
                # bound vortex filament, not the collocation points from above
                _assemble_AIC_mtx(self.mtx, inputs, self.surfaces, skip=True)

                # Compute the induced velocities at the midpoints of the
                # bound vortex filaments
                for ind in range(3):
                    self.v[:, ind] = self.mtx[:, :, ind].dot(circ)

                # Add the freestream velocity to the induced velocity so that
                # self.v is the total velocity seen at the point
                self.v[:, 0] += cosa * inputs['v']
                self.v[:, 2] += sina * inputs['v']

                not_real_outputs = {}

                i = 0
                for surface in self.surfaces:
                    name = surface['name']
                    nx = surface['num_x']
                    ny = surface['num_y']
                    num_panels = (nx - 1) * (ny - 1)

                    b_pts = inputs[name + '_b_pts']

                    sec_forces = OAS_API.oas_api.forcecalc(self.v[i:i+num_panels, :], circ[i:i+num_panels], rho, b_pts)

                    # Reshape the forces into the expected form
                    not_real_outputs[name + '_sec_forces'] = sec_forces.reshape((nx-1, ny-1, 3), order='F')

                    i += num_panels

                for surface in self.surfaces:


                    name = surface['name']
                    dS = partials[name + '_sec_forces', name + '_cos_sweep'].copy()
                    d_inputs = {}
                    sec_forcesb = np.zeros((surface['num_x'] - 1, surface['num_y'] - 1, 3))

                    sweep_angle = inputs[name + '_cos_sweep'] / inputs[name + '_widths']
                    beta = np.sqrt(1 - inputs['M']**2 * sweep_angle**2)

                    if not compressible:
                        beta[:] = 1.

                    for k, val in enumerate(sec_forcesb.flatten()):
                        for key in inputs:
                            d_inputs[key] = inputs[key].copy()
                            d_inputs[key][:] = 0.

                        sec_forcesb[:] = 0.
                        sec_forcesb = sec_forcesb.flatten()
                        sec_forcesb[k] = 1.
                        sec_forcesb = sec_forcesb.reshape(surface['num_x'] - 1, surface['num_y'] - 1, 3)

                        for i, B in enumerate(beta):
                            sec_forcesb[:, i, :] /= B

                        sec_forcesb = sec_forcesb.reshape((-1, 3), order='F')

                        circ = inputs['circulations']
                        alpha = inputs['alpha'] * np.pi / 180.
                        cosa = np.cos(alpha)
                        sina = np.sin(alpha)

                        ind = 0
                        rho = inputs['rho'].real
                        v = inputs['v']
                        vb = np.zeros(self.v.shape)

                        for surface_ in self.surfaces:
                            name_ = surface_['name']
                            nx_ = surface_['num_x']
                            ny_ = surface_['num_y']
                            num_panels_ = (nx_ - 1) * (ny_ - 1)

                            if name == name_:
                                b_pts = inputs[name_ + '_b_pts']

                                v_b, circb, rhob, bptsb, _ = OAS_API.oas_api.forcecalc_b(self.v[ind:ind+num_panels_, :], circ[ind:ind+num_panels_], rho, b_pts, sec_forcesb)

                                if 'circulations' in d_inputs:
                                    d_inputs['circulations'][ind:ind+num_panels_] += circb
                                vb[ind:ind+num_panels_] = v_b
                                if 'rho' in d_inputs:
                                    d_inputs['rho'] += rhob
                                if name + '_b_pts' in d_inputs:
                                    d_inputs[name_ + '_b_pts'] += bptsb

                            ind += num_panels_

                        sinab = inputs['v'] * np.sum(vb[:, 2])
                        if 'v' in d_inputs:
                            d_inputs['v'] += cosa * np.sum(vb[:, 0]) + sina * np.sum(vb[:, 2])
                        cosab = inputs['v'] * np.sum(vb[:, 0])
                        ab = np.cos(alpha) * sinab - np.sin(alpha) * cosab
                        if 'alpha' in d_inputs:
                            d_inputs['alpha'] += np.pi * ab / 180.

                        mtxb = np.zeros(self.mtx.shape)
                        circb = np.zeros(circ.shape)
                        for i in range(3):
                            for j in range(self.tot_panels):
                                mtxb[j, :, i] += circ * vb[j, i]
                                circb += self.mtx[j, :, i].real * vb[j, i]

                        if 'circulations' in d_inputs:
                            d_inputs['circulations'] += circb

                        _assemble_AIC_mtx_b(mtxb, inputs, d_inputs, self.surfaces, skip=True)

                        for key in d_inputs:
                            partials[name + '_sec_forces', key][k, :] = d_inputs[key].flatten()

                    if compressible:

                        M = inputs['M']
                        S = inputs[name + '_cos_sweep']
                        X = not_real_outputs[name + '_sec_forces']
                        W = inputs[name + '_widths']
                        nx = surface['num_x']
                        ny = surface['num_y']
                        d_M = np.zeros((nx-1, ny-1, 3))
                        d_S = np.zeros((nx-1, ny-1, 3))
                        d_W = np.zeros((nx-1, ny-1, 3))
                        for ix in range(nx-1):
                            for ind in range(3):
                                d_M[ix, :, ind] = (M * S**2 * X[ix, :, ind]) / (W**2 * (1 - (M**2 * S**2) / W**2)**1.5)
                                d_S[ix, :, ind] = (M**2 * S * X[ix, :, ind]) / (W**2 * (1 - (M**2 * S**2) / W**2)**1.5)
                                d_W[ix, :, ind] = (M**2 * S**2 * X[ix, :, ind]) / (W**3 * (1 - (M**2 * S**2) / W**2)**1.5)

                        partials[name + '_sec_forces', 'M'] = d_M.flatten()

                        ny3 = (ny-1) * 3
                        for i in range(ny-1):
                            for ix in range(nx-1):
                                for ind in range(3):
                                    partials[name + '_sec_forces', name + '_cos_sweep'][3*i+ix*ny3:3*(i+1)+ix*ny3, i][ind] = d_S[ix, i, ind]
                                    partials[name + '_sec_forces', name + '_widths'][3*i+ix*ny3:3*(i+1)+ix*ny3, i][ind] = -d_W[ix, i, ind]
