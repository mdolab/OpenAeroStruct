""" Script to plot geometry for aero, struct, or aerostruct cases.

Usage is
`python view_geometry.py aero.db` for aero only,
`python view_geometry.py struct.db` for struct only,
`python view_geometry.py aerostruct.db` for aerostruct, or
`python view_geometry.py __name__` for user-named database.

You can select a certain zoom factor for the 3d view by adding a number as a
last keyword.
The larger the number, the closer the view. Floats or ints are accepted.

Ex: `python view_geometry.py aero.db 1` a wider view than `python view_geometry.py aero.db 5`.

"""


from __future__ import division, print_function
import sys
major_python_version = sys.version_info[0]

from six import iteritems
import numpy as np

try:
    import matplotlib
    matplotlib.use('TkAgg')
    matplotlib.rcParams['lines.linewidth'] = 2
    matplotlib.rcParams['axes.edgecolor'] = 'gray'
    matplotlib.rcParams['axes.linewidth'] = 0.5
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,\
        NavigationToolbar2TkAgg
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    import matplotlib.animation as manimation
    import sqlitedict
except:
    print()
    print("Correct plotting modules not available; please consult import list")
    print()

#####################
# User-set parameters
#####################

db_name = sys.argv[1]

try:
    zoom_scale = sys.argv[2]
except:
    zoom_scale = 2.8

class Display(object):
    def __init__(self, db_name):

        self.f = plt.figure(dpi=100, figsize=(14, 14), facecolor='white')

        self.ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4,
                                   colspan=4, projection='3d')

        self.num_iters = 0
        self.db_name = db_name
        self.show_wing = True
        self.show_tube = True
        self.curr_pos = 0
        self.old_n = 0
        self.aerostruct = False

        self.load_db()

        self.plot_wing()

    def load_db(self):
        self.db = sqlitedict.SqliteDict(self.db_name, 'iterations')

        self.twist = []
        self.mesh = []
        self.def_mesh = []
        self.radius = []
        self.thickness = []
        sec_forces = []
        normals = []
        widths = []
        self.lift = []
        self.lift_ell = []
        self.vonmises = []
        alpha = []
        rho = []
        v = []
        self.CL = []
        self.AR = []
        self.S_ref = []
        self.obj = []

        meta_db = sqlitedict.SqliteDict(self.db_name, 'metadata')
        self.opt = False
        for item in meta_db['Unknowns']:
            if 'is_objective' in meta_db['Unknowns'][item].keys():
                self.obj_key = item
                if major_python_version == 3:
                    keys_length = sum(1 for _ in self.db.keys())
                else:
                    keys_length = len(self.db.keys())
                if keys_length > 2:
                    self.opt = True

        deriv_keys = sqlitedict.SqliteDict(self.db_name, 'derivs').keys()
        deriv_keys = [int(key.split('|')[-1]) for key in deriv_keys]

        for i, (case_name, case_data) in enumerate(iteritems(self.db)):

            if i == 0:
                pass
            elif i not in deriv_keys:
                if deriv_keys:
                    continue # don't plot these cases

            if self.opt:
                self.obj.append(case_data['Unknowns'][self.obj_key])

            names = []
            for key in case_data['Unknowns'].keys():

                # Aerostructural
                if 'coupled' in key and 'loads' in key:
                    self.aerostruct = True
                    names.append(key.split('_')[:-1][0])

                # Aero only
                elif 'def_mesh' in key and 'coupled' not in key:
                    names.append(key.split('.')[0])

                # Structural only
                elif 'disp_aug' in key and 'coupled' not in key:
                    names.append(key.split('.')[0])

            self.names = names
            n_names = len(names)

            self.twist_included = False

            # Loop through each of the surfaces
            for name in names:

                # Check if this is an aerostructual case; treat differently
                # due to the way the problem is organized
                if not self.aerostruct:

                    # A mesh exists for all types of cases
                    self.mesh.append(case_data['Unknowns'][name+'.mesh'])

                    try:
                        self.radius.append(case_data['Unknowns'][name+'.radius'])
                        self.thickness.append(case_data['Unknowns'][name+'.thickness'])
                        self.vonmises.append(
                            np.max(case_data['Unknowns'][name+'.vonmises'], axis=1))
                        self.show_tube = True
                    except:
                        self.show_tube = False
                    try:
                        self.def_mesh.append(case_data['Unknowns'][name+'.def_mesh'])
                        normals.append(case_data['Unknowns'][name+'.normals'])
                        widths.append(case_data['Unknowns'][name+'.widths'])
                        sec_forces.append(case_data['Unknowns']['aero_states.' + name + '_sec_forces'])
                        self.CL.append(case_data['Unknowns'][name+'_perf.CL1'])
                        self.S_ref.append(case_data['Unknowns'][name+'.S_ref'])
                        self.show_wing = True

                        # Not the best solution for now, but this will ensure
                        # that this plots corectly even if twist isn't a desvar
                        try:
                            self.twist.append(case_data['Unknowns'][name+'.twist'])
                            self.twist_included = True
                        except:
                            pass
                    except:
                        self.show_wing = False
                else:
                    self.show_wing, self.show_tube = True, True
                    short_name = name.split('.')[1:][0]

                    self.mesh.append(case_data['Unknowns'][short_name+'.mesh'])
                    self.radius.append(case_data['Unknowns'][short_name+'.radius'])
                    self.thickness.append(case_data['Unknowns'][short_name+'.thickness'])
                    self.vonmises.append(
                        np.max(case_data['Unknowns'][short_name+'_perf.vonmises'], axis=1))
                    self.def_mesh.append(case_data['Unknowns'][name+'.def_mesh'])
                    normals.append(case_data['Unknowns'][name+'.normals'])
                    widths.append(case_data['Unknowns'][name+'.widths'])
                    sec_forces.append(case_data['Unknowns']['coupled.aero_states.' + short_name + '_sec_forces'])
                    self.CL.append(case_data['Unknowns'][short_name+'_perf.CL1'])
                    self.S_ref.append(case_data['Unknowns'][name+'.S_ref'])

                    # Not the best solution for now, but this will ensure
                    # that this plots corectly even if twist isn't a desvar
                    try:
                        self.twist.append(case_data['Unknowns'][short_name+'.twist'])
                        self.twist_included = True
                    except:
                        pass

                if not self.twist_included:
                    ny = self.mesh[0].shape[1]
                    self.twist.append(np.zeros(ny))

            if self.show_wing:
                alpha.append(case_data['Unknowns']['alpha'] * np.pi / 180.)
                rho.append(case_data['Unknowns']['rho'])
                v.append(case_data['Unknowns']['v'])

        if self.opt:
            self.num_iters = np.max([int(len(self.mesh) / n_names) - 1, 1])
        else:
            self.num_iters = 0

        symm_count = 0
        for mesh in self.mesh:
            if np.all(mesh[:, :, 1] >= -1e-8) or np.all(mesh[:, :, 1] <= 1e-8):
                symm_count += 1
        if symm_count == len(self.mesh):
            self.symmetry = True
        else:
            self.symmetry = False

        if self.symmetry:

            new_mesh = []
            if self.show_tube:
                new_r = []
                new_thickness = []
                new_vonmises = []
            if self.show_wing:
                new_twist = []
                new_sec_forces = []
                new_def_mesh = []
                new_widths = []
                new_normals = []

            for i in range(self.num_iters + 1):
                for j, name in enumerate(names):
                    mirror_mesh = self.mesh[i*n_names+j].copy()
                    mirror_mesh[:, :, 1] *= -1.
                    mirror_mesh = mirror_mesh[:, ::-1, :][:, 1:, :]
                    new_mesh.append(np.hstack((self.mesh[i*n_names+j], mirror_mesh)))

                    if self.show_tube:
                        thickness = self.thickness[i*n_names+j]
                        new_thickness.append(np.hstack((thickness, thickness[::-1])))
                        r = self.radius[i*n_names+j]
                        new_r.append(np.hstack((r, r[::-1])))
                        vonmises = self.vonmises[i*n_names+j]
                        new_vonmises.append(np.hstack((vonmises, vonmises[::-1])))

                    if self.show_wing:
                        mirror_mesh = self.def_mesh[i*n_names+j].copy()
                        mirror_mesh[:, :, 1] *= -1.
                        mirror_mesh = mirror_mesh[:, ::-1, :][:, 1:, :]
                        new_def_mesh.append(np.hstack((self.def_mesh[i*n_names+j], mirror_mesh)))

                        mirror_normals = normals[i*n_names+j].copy()
                        mirror_normals = mirror_normals[:, ::-1, :][:, 1:, :]
                        new_normals.append(np.hstack((normals[i*n_names+j], mirror_normals)))

                        mirror_forces = sec_forces[i*n_names+j].copy()
                        mirror_forces = mirror_forces[:, ::-1, :]
                        new_sec_forces.append(np.hstack((sec_forces[i*n_names+j], mirror_forces)))

                        new_widths.append(np.hstack((widths[i*n_names+j], widths[i*n_names+j][::-1])))
                        twist = self.twist[i*n_names+j]
                        new_twist.append(np.hstack((twist, twist[::-1][1:])))

            self.mesh = new_mesh
            if self.show_tube:
                self.thickness = new_thickness
                self.radius = new_r
                self.vonmises = new_vonmises
            if self.show_wing:
                self.def_mesh = new_def_mesh
                self.twist = new_twist
                widths = new_widths
                normals = new_normals
                sec_forces = new_sec_forces

        if self.show_wing:
            for i in range(self.num_iters + 1):
                for j, name in enumerate(names):
                    m_vals = self.mesh[i*n_names+j].copy()
                    cvec = m_vals[0, :, :] - m_vals[-1, :, :]
                    chords = np.sqrt(np.sum(cvec**2, axis=1))
                    chords = 0.5 * (chords[1:] + chords[:-1])
                    a = alpha[i]
                    cosa = np.cos(a)
                    sina = np.sin(a)

                    forces = np.sum(sec_forces[i*n_names+j], axis=0)

                    lift = (-forces[:, 0] * sina + forces[:, 2] * cosa) / \
                        widths[i*n_names+j]/0.5/rho[i]/v[i]**2

                    span = (m_vals[0, :, 1] / (m_vals[0, -1, 1] - m_vals[0, 0, 1]))
                    span = span - (span[0] + .5)

                    lift_area = np.sum(lift * (span[1:] - span[:-1]))

                    lift_ell = 4 * lift_area / np.pi * np.sqrt(1 - (2*span)**2)

                    self.lift.append(lift)
                    self.lift_ell.append(lift_ell)

                    wingspan = np.abs(m_vals[0, -1, 1] - m_vals[0, 0, 1])
                    self.AR.append(wingspan**2 / self.S_ref[i*n_names+j])

            # recenter def_mesh points for better viewing
            for i in range(self.num_iters + 1):
                center = np.zeros((3))
                for j in range(n_names):
                    center += np.mean(self.def_mesh[i*n_names+j], axis=(0,1))
                for j in range(n_names):
                    self.def_mesh[i*n_names+j] -= center / n_names

        # recenter mesh points for better viewing
        for i in range(self.num_iters + 1):
            center = np.zeros((3))
            for j in range(n_names):
                center += np.mean(self.mesh[i*n_names+j], axis=(0,1))
            for j in range(n_names):
                self.mesh[i*n_names+j] -= center / n_names

        if self.show_wing:
            self.min_twist, self.max_twist = self.get_list_limits(self.twist)
            diff = (self.max_twist - self.min_twist) * 0.05
            self.min_twist -= diff
            self.max_twist += diff
            self.min_l, self.max_l = self.get_list_limits(self.lift)
            self.min_le, self.max_le = self.get_list_limits(self.lift_ell)
            self.min_l, self.max_l = min(self.min_l, self.min_le), max(self.max_l, self.max_le)
            diff = (self.max_l - self.min_l) * 0.05
            self.min_l -= diff
            self.max_l += diff
        if self.show_tube:
            self.min_t, self.max_t = self.get_list_limits(self.thickness)
            diff = (self.max_t - self.min_t) * 0.05
            self.min_t -= diff
            self.max_t += diff
            self.min_vm, self.max_vm = self.get_list_limits(self.vonmises)
            diff = (self.max_vm - self.min_vm) * 0.05
            self.min_vm -= diff
            self.max_vm += diff

    def plot_wing(self):

        n_names = len(self.names)
        self.ax.cla()
        az = self.ax.azim
        el = self.ax.elev
        dist = self.ax.dist

        for j, name in enumerate(self.names):
            mesh0 = self.mesh[self.curr_pos*n_names+j].copy()

            self.ax.set_axis_off()

            if self.show_wing:
                def_mesh0 = self.def_mesh[self.curr_pos*n_names+j]
                x = mesh0[:, :, 0]
                y = mesh0[:, :, 1]
                z = mesh0[:, :, 2]

                try:  # show deformed mesh option may not be available
                    if self.show_def_mesh.get():
                        x_def = def_mesh0[:, :, 0]
                        y_def = def_mesh0[:, :, 1]
                        z_def = def_mesh0[:, :, 2]

                        self.c2.grid(row=0, column=3, padx=5, sticky=Tk.W)
                        if self.ex_def.get():
                            z_def = (z_def - z) * 10 + z_def
                            def_mesh0 = (def_mesh0 - mesh0) * 30 + def_mesh0
                        else:
                            def_mesh0 = (def_mesh0 - mesh0) * 2 + def_mesh0
                        self.ax.plot_wireframe(x_def, y_def, z_def, rstride=1, cstride=1, color='k')
                        self.ax.plot_wireframe(x, y, z, rstride=1, cstride=1, color='k', alpha=.3)
                    else:
                        self.ax.plot_wireframe(x, y, z, rstride=1, cstride=1, color='k')
                        self.c2.grid_forget()
                except:
                    self.ax.plot_wireframe(x, y, z, rstride=1, cstride=1, color='k')

            if self.show_tube:
                r0 = self.radius[self.curr_pos*n_names+j]
                t0 = self.thickness[self.curr_pos*n_names+j]
                colors = t0
                colors = colors / np.max(colors)
                num_circ = 12
                fem_origin = 0.35
                n = mesh0.shape[1]
                p = np.linspace(0, 2*np.pi, num_circ)

                for i, thick in enumerate(t0):
                    r = np.array((r0[i], r0[i]))
                    R, P = np.meshgrid(r, p)
                    X, Z = R*np.cos(P), R*np.sin(P)
                    chords = mesh0[-1, :, 0] - mesh0[0, :, 0]
                    comp = fem_origin * chords + mesh0[0, :, 0]
                    X[:, 0] += comp[i]
                    X[:, 1] += comp[i+1]
                    Z[:, 0] += fem_origin * (mesh0[-1, i, 2] - mesh0[0, i, 2]) + mesh0[0, i, 2]
                    Z[:, 1] += fem_origin * (mesh0[-1, i+1, 2] - mesh0[0, i+1, 2]) + mesh0[0, i+1, 2]
                    Y = np.empty(X.shape)
                    Y[:] = np.linspace(mesh0[0, i, 1], mesh0[0, i+1, 1], 2)
                    col = np.zeros(X.shape)
                    col[:] = colors[i]
                    try:
                        self.ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                            facecolors=cm.viridis(col), linewidth=0)
                    except:
                        self.ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                            facecolors=cm.coolwarm(col), linewidth=0)

        lim = 0.
        for j in range(n_names):
            ma = np.max(self.mesh[self.curr_pos*n_names+j], axis=(0,1,2))
            if ma > lim:
                lim = ma
        lim /= float(zoom_scale)
        self.ax.auto_scale_xyz([-lim, lim], [-lim, lim], [-lim, lim])

        self.ax.view_init(elev=el, azim=az)  # Reproduce view
        self.ax.dist = dist


    def check_length(self):
        # Load the current sqlitedict
        db = sqlitedict.SqliteDict(self.db_name, 'iterations')

        # Get the number of current iterations
        # Minus one because OpenMDAO uses 1-indexing
        self.num_iters = int(db.keys()[-1].split('|')[-1])

    def get_list_limits(self, input_list):
        list_min = 1.e20
        list_max = -1.e20
        for list_ in input_list:
            mi = np.min(list_)
            if mi < list_min:
                list_min = mi
            ma = np.max(list_)
            if ma > list_max:
                list_max = ma

        return list_min, list_max


def disp_plot(db_name):
    disp = Display(db_name)
    plt.show()

if __name__ == '__main__':
    disp_plot(db_name)
