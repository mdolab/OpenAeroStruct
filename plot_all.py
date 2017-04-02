""" Script to plot results from aero, struct, or aerostruct optimization.

Usage is
`python plot_all.py a` for aero only,
`python plot_all.py s` for struct only,
`python plot_all.py as` for aerostruct, or
`python plot_all.py __name__` for user-named database.

The script automatically appends '.db' to the provided name.
Ex: `python plot_all.py example` opens 'example.db'.

"""


from __future__ import division
import tkFont
import Tkinter as Tk
import sys

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

import numpy
import sqlitedict
import aluminum

#####################
# User-set parameters
#####################

if sys.argv[1] == 'as':
    filename = 'aerostruct'
elif sys.argv[1] == 'a':
    filename = 'vlm'
elif sys.argv[1] == 's':
    filename = 'spatialbeam'
else:
    filename = sys.argv[1]

db_name = filename + '.db'


class Display(object):
    def __init__(self, db_name):

        self.root = Tk.Tk()
        self.root.wm_title("Viewer")

        self.f = plt.figure(dpi=100, figsize=(12, 6), facecolor='white')
        self.canvas = FigureCanvasTkAgg(self.f, master=self.root)
        self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        self.options_frame = Tk.Frame(self.root)
        self.options_frame.pack()

        toolbar = NavigationToolbar2TkAgg(self.canvas, self.root)
        toolbar.update()
        self.canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        self.ax = plt.subplot2grid((4, 8), (0, 0), rowspan=4,
                                   colspan=4, projection='3d')

        self.num_iters = 0
        self.db_name = db_name
        self.show_wing = True
        self.show_tube = True
        self.curr_pos = 0
        self.old_n = 0

        self.load_db()

        if self.show_wing and not self.show_tube:
            self.ax2 = plt.subplot2grid((4, 8), (0, 4), rowspan=2, colspan=4)
            self.ax3 = plt.subplot2grid((4, 8), (2, 4), rowspan=2, colspan=4)
        if self.show_tube and not self.show_wing:
            self.ax4 = plt.subplot2grid((4, 8), (0, 4), rowspan=2, colspan=4)
            self.ax5 = plt.subplot2grid((4, 8), (2, 4), rowspan=2, colspan=4)
        if self.show_wing and self.show_tube:
            self.ax2 = plt.subplot2grid((4, 8), (0, 4), colspan=4)
            self.ax3 = plt.subplot2grid((4, 8), (1, 4), colspan=4)
            self.ax4 = plt.subplot2grid((4, 8), (2, 4), colspan=4)
            self.ax5 = plt.subplot2grid((4, 8), (3, 4), colspan=4)

    def load_db(self):
        self.db = sqlitedict.SqliteDict(self.db_name, 'openmdao')
        self.twist = []
        self.mesh = []
        self.def_mesh = []
        self.r = []
        self.t = []
        sec_forces = []
        normals = []
        widths = []
        self.lift = []
        self.lift_ell = []
        self.vonmises = []
        alpha = []
        rho = []
        v = []
        self.aero_ind = []
        self.fem_ind = []
        self.CL = []
        self.AR = []
        self.S_ref = []
        self.obj = []

        for tag in self.db['metadata']:
            for item in self.db['metadata'][tag]:
                for flag in self.db['metadata'][tag][item]:
                    if 'is_objective' in flag:
                        self.obj_key = item
        for case_name, case_data in self.db.iteritems():
            if "metadata" in case_name or "derivs" in case_name:
                continue  # don't plot these cases
            try:
                self.db[case_name + '/derivs']
            except:
                continue

            self.mesh.append(case_data['Unknowns']['mesh'])
            self.aero_ind.append(case_data['Unknowns']['aero_ind'])
            self.obj.append(case_data['Unknowns'][self.obj_key])

            try:
                self.r.append(case_data['Unknowns']['r'])
                self.t.append(case_data['Unknowns']['thickness'])
                self.fem_ind.append(case_data['Unknowns']['fem_ind'])
                self.vonmises.append(
                    numpy.max(case_data['Unknowns']['vonmises'], axis=1))
                self.show_tube = True
            except:
                self.show_tube = False
                pass
            try:
                def_mesh = case_data['Unknowns']['def_mesh']
                self.def_mesh.append(def_mesh)
                nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = self.aero_ind[0][0, :]
                def_mesh = def_mesh[:n, :].reshape(nx, ny, 3)

                self.twist.append(case_data['Unknowns']['twist'])

                normals.append(case_data['Unknowns']['normals'])
                widths.append(case_data['Unknowns']['widths'])
                sec_forces.append(case_data['Unknowns']['sec_forces'])
                alpha.append(case_data['Unknowns']['alpha'] * numpy.pi / 180.)
                rho.append(case_data['Unknowns']['rho'])
                v.append(case_data['Unknowns']['v'])
                self.CL.append(case_data['Unknowns']['CL1'])
                self.S_ref.append(case_data['Unknowns']['S_ref'])
                self.show_wing = True
            except:
                self.show_wing = False
                pass

        self.num_iters = numpy.max([len(self.mesh) - 1, 1])

        if self.show_wing:

            nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = self.aero_ind[0][0, :]

            for i in range(self.num_iters + 1):
                m_vals = self.mesh[i][:n, :].reshape(nx, ny, 3)
                cvec = m_vals[0, :, :] - m_vals[-1, :, :]
                chords = numpy.sqrt(numpy.sum(cvec**2, axis=1))
                chords = 0.5 * (chords[1:] + chords[:-1])
                a = alpha[i]
                cosa = numpy.cos(a)
                sina = numpy.sin(a)
                forces = numpy.sum(sec_forces[i][:n_panels, :].reshape(nx-1, ny-1, 3, order='F'), axis=0)
                widths_ = widths[i][:ny-1]

                lift = (-forces[:, 0] * sina + forces[:, 2] * cosa) / \
                    widths_/0.5/rho[i]/v[i]**2
                # lift = (-forces[:, 0] * sina + forces[:, 2] * cosa)/chords/0.5/rho[i]/v[i]**2
                # lift = (-forces[:, 0] * sina + forces[:, 2] * cosa)*chords/0.5/rho[i]/v[i]**2

                span = (m_vals[0, :, 1] / (m_vals[0, -1, 1] - m_vals[0, 0, 1]))
                span = span - (span[0] + .5)

                lift_area = numpy.sum(lift * (span[1:] - span[:-1]))

                lift_ell = 4 * lift_area / numpy.pi * \
                    numpy.sqrt(1 - (2*span)**2)

                self.lift.append(lift)
                self.lift_ell.append(lift_ell)

                wingspan = numpy.abs(m_vals[0, -1, 1] - m_vals[0, 0, 1])
                self.AR.append(wingspan**2 / self.S_ref[i])

            # recenter def_mesh points for better viewing
            for i in range(self.num_iters + 1):
                center = numpy.mean(self.mesh[i], axis=0)
                self.def_mesh[i] = self.def_mesh[i] - center

        # recenter mesh points for better viewing
        for i in range(self.num_iters + 1):
            # center defined as the average of all nodal points
            center = numpy.mean(self.mesh[i], axis=0)
            # center defined as the mean of the min and max in each direction
            # center = (numpy.max(self.mesh[i], axis=0) + numpy.min(self.mesh[i], axis=0)) / 2
            self.mesh[i] = self.mesh[i] - center


        if self.show_wing:
            self.min_twist, self.max_twist = numpy.min(self.twist), numpy.max(self.twist)
            diff = (self.max_twist - self.min_twist) * 0.05
            self.min_twist -= diff
            self.max_twist += diff
            self.min_l, self.max_l = numpy.min(self.lift), numpy.max(self.lift)
            self.min_le, self.max_le = numpy.min(self.lift_ell), numpy.max(self.lift_ell)
            self.min_l, self.max_l = min(self.min_l, self.min_le), max(self.max_l, self.max_le)
            diff = (self.max_l - self.min_l) * 0.05
            self.min_l -= diff
            self.max_l += diff
        if self.show_tube:
            self.min_t, self.max_t = numpy.min(self.t), numpy.max(self.t)
            diff = (self.max_t - self.min_t) * 0.05
            self.min_t -= diff
            self.max_t += diff
            self.min_vm, self.max_vm = numpy.min(self.vonmises), numpy.max(self.vonmises)
            diff = (self.max_vm - self.min_vm) * 0.05
            self.min_vm -= diff
            self.max_vm += diff

    def plot_sides(self):
        nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = self.aero_ind[0][0, :]
        m_vals = self.mesh[self.curr_pos][:n, :].reshape(nx, ny, 3).copy()
        span = m_vals[0, -1, 1] - m_vals[0, 0, 1]
        rel_span = (m_vals[0, :, 1] - m_vals[0, 0, 1]) * 2 / span - 1
        span_diff = ((m_vals[0, :-1, 1] + m_vals[0, 1:, 1]) / 2 - m_vals[0, 0, 1]) * 2 / span - 1

        if self.show_wing:
            self.ax2.cla()
            self.ax3.cla()
            t_vals = self.twist[self.curr_pos]
            l_vals = self.lift[self.curr_pos]
            le_vals = self.lift_ell[self.curr_pos]

            self.ax2.plot(rel_span, t_vals, lw=2, c='b')
            self.ax2.locator_params(axis='y',nbins=5)
            self.ax2.locator_params(axis='x',nbins=3)
            self.ax2.set_ylim([self.min_twist, self.max_twist])
            self.ax2.set_xlim([-1, 1])
            self.ax2.set_ylabel('twist', rotation="horizontal", ha="right")

            self.ax3.plot(span_diff, l_vals, lw=2, c='b')
            self.ax3.plot(rel_span, le_vals, '--', lw=2, c='g')
            self.ax3.text(0.05, 0.8, 'elliptical',
                transform=self.ax3.transAxes, color='g')
            self.ax3.locator_params(axis='y',nbins=4)
            self.ax3.locator_params(axis='x',nbins=3)
            self.ax3.set_ylim([self.min_l, self.max_l])
            self.ax3.set_xlim([-1, 1])
            self.ax3.set_ylabel('lift', rotation="horizontal", ha="right")

        if self.show_tube:
            n_fem, i_fem = self.fem_ind[0][0, :]
            self.ax4.cla()
            self.ax5.cla()
            thick_vals = self.t[self.curr_pos][i_fem:i_fem+n_fem-1]
            vm_vals = self.vonmises[self.curr_pos][i_fem:i_fem+n_fem-1]

            self.ax4.plot(span_diff, thick_vals, lw=2, c='b')
            self.ax4.locator_params(axis='y',nbins=4)
            self.ax4.locator_params(axis='x',nbins=3)
            self.ax4.set_ylim([self.min_t, self.max_t])
            self.ax4.set_xlim([-1, 1])
            self.ax4.set_ylabel('thickness', rotation="horizontal", ha="right")

            self.ax5.plot(span_diff, vm_vals, lw=2, c='b')
            self.ax5.locator_params(axis='y',nbins=4)
            self.ax5.locator_params(axis='x',nbins=3)
            self.ax5.set_ylim([self.min_vm, self.max_vm])
            self.ax5.set_ylim([0, 25e6])
            self.ax5.set_xlim([-1, 1])
            self.ax5.set_ylabel('von mises', rotation="horizontal", ha="right")
            self.ax5.axhline(aluminum.stress, c='r', lw=2, ls='--')
            self.ax5.text(0.05, 0.85, 'failure limit',
                transform=self.ax5.transAxes, color='r')

    def plot_wing(self):
        self.ax.cla()
        az = self.ax.azim
        el = self.ax.elev
        dist = self.ax.dist
        nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = self.aero_ind[0][0, :]
        mesh0 = self.mesh[self.curr_pos][:n, :].reshape(nx, ny, 3).copy()

        self.ax.set_axis_off()

        if self.show_wing:
            for i_surf, row in enumerate(self.aero_ind[0]):
                nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = row

                mesh0 = self.mesh[self.curr_pos][i:i+n, :].reshape(nx, ny, 3).copy()
                def_mesh0 = self.def_mesh[self.curr_pos][i:i+n, :].reshape(nx, ny, 3)
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
            num_surf = self.fem_ind[0].shape[0]
            for i_surf, row in enumerate(self.fem_ind[0]):
                n_fem, i_fem = row
                nx, ny, n, n_bpts, n_panels, i, i_bpts, i_panels = self.aero_ind[0][i_surf, :]
                mesh0 = self.mesh[self.curr_pos][i:i+n, :].reshape(nx, ny, 3).copy()

                r0 = self.r[self.curr_pos][i_fem-i_surf:i_fem-i_surf+n_fem-1]
                t0 = self.t[self.curr_pos][i_fem-i_surf:i_fem-i_surf+n_fem-1]
                colors = t0
                colors = colors / numpy.max(colors)
                num_circ = 12
                fem_origin = 0.35
                n = mesh0.shape[1]
                p = numpy.linspace(0, 2*numpy.pi, num_circ)
                if self.show_wing:
                    if self.show_def_mesh.get():
                        mesh0[:, :, 2] = def_mesh0[:, :, 2]
                for i, thick in enumerate(t0):
                    r = numpy.array((r0[i], r0[i]))
                    R, P = numpy.meshgrid(r, p)
                    X, Z = R*numpy.cos(P), R*numpy.sin(P)
                    chords = mesh0[-1, :, 0] - mesh0[0, :, 0]
                    comp = fem_origin * chords + mesh0[0, :, 0]
                    X[:, 0] += comp[i]
                    X[:, 1] += comp[i+1]
                    Z[:, 0] += fem_origin * mesh0[-1, i, 2]
                    Z[:, 1] += fem_origin * mesh0[-1, i+1, 2]
                    Y = numpy.empty(X.shape)
                    Y[:] = numpy.linspace(mesh0[0, i, 1], mesh0[0, i+1, 1], 2)
                    col = numpy.zeros(X.shape)
                    col[:] = colors[i]
                    try:
                        self.ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                            facecolors=cm.viridis(col), linewidth=0)
                    except:
                        self.ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                            facecolors=cm.coolwarm(col), linewidth=0)

        lim = numpy.max(numpy.max(mesh0)) / 2.8
        self.ax.auto_scale_xyz([-lim, lim], [-lim, lim], [-lim, lim])
        self.ax.set_title("Major Iteration: {}".format(self.curr_pos))
        round_to_n = lambda x, n: round(x, -int(numpy.floor(numpy.log10(abs(x)))) + (n - 1))
        obj_val = round_to_n(self.obj[self.curr_pos], 7)
        self.ax.text2D(.55, .05, self.obj_key + ': {}'.format(obj_val),
            transform=self.ax.transAxes, color='k')
        if self.show_wing and not self.show_tube:
            span_eff = self.CL[self.curr_pos]**2 / numpy.pi / self.AR[self.curr_pos] / obj_val
            self.ax.text2D(.55, .0, 'e: {}'.format(round_to_n(span_eff[0], 7)),
                transform=self.ax.transAxes, color='k')

        self.ax.view_init(elev=el, azim=az)  # Reproduce view
        self.ax.dist = dist

    def save_video(self):
        FFMpegWriter = manimation.writers['ffmpeg']
        metadata = dict(title='Movie', artist='Matplotlib')
        writer = FFMpegWriter(fps=5, metadata=metadata, bitrate=3000)

        with writer.saving(self.f, "movie.mp4", 100):
            self.curr_pos = 0
            self.update_graphs()
            self.f.canvas.draw()
            plt.draw()
            for i in range(10):
                writer.grab_frame()

            for i in range(self.num_iters):
                self.curr_pos = i
                self.update_graphs()
                self.f.canvas.draw()
                plt.draw()
                writer.grab_frame()

            self.curr_pos = self.num_iters
            self.update_graphs()
            self.f.canvas.draw()
            plt.draw()
            for i in range(20):
                writer.grab_frame()

    def update_graphs(self, e=None):
        if e is not None:
            self.curr_pos = int(e)
            self.curr_pos = self.curr_pos % (self.num_iters + 1)

        self.plot_wing()
        self.plot_sides()
        self.canvas.show()

    def check_length(self):
        db = sqlitedict.SqliteDict(self.db_name, 'openmdao')
        n = 0
        for case_name, case_data in db.iteritems():
            if "metadata" in case_name or "derivs" in case_name:
                continue  # don't plot these cases
            try:
                db[case_name + '/derivs']
            except:
                continue
            n += 1
        self.num_iters = n

    def auto_ref(self):
        """
        Automatically refreshes the history file, which is
        useful if examining a running optimization.
        """
        if self.var_ref.get():
            self.root.after(800, self.auto_ref)
            self.check_length()
            self.update_graphs()

            if self.num_iters > self.old_n:
                self.load_db()
                self.old_n = self.num_iters
                self.draw_slider()

    def save_image(self):
        fname = 'fig' + '.png'
        plt.savefig(fname)

    def quit(self):
        """
        Destroy GUI window cleanly if quit button pressed.
        """
        self.root.quit()
        self.root.destroy()

    def draw_slider(self):
        # scale to choose iteration to view
        self.w = Tk.Scale(
            self.options_frame,
            from_=0, to=self.num_iters,
            orient=Tk.HORIZONTAL,
            resolution=1,
            font=tkFont.Font(family="Helvetica", size=10),
            command=self.update_graphs,
            length=200)

        if self.curr_pos == self.num_iters - 1 or self.curr_pos == 0:
            self.curr_pos = self.num_iters
        self.w.set(self.curr_pos)
        self.w.grid(row=0, column=1, padx=5, sticky=Tk.W)

    def draw_GUI(self):
        """
        Create the frames and widgets in the bottom section of the canvas.
        """
        font = tkFont.Font(family="Helvetica", size=10)

        lab_font = Tk.Label(
            self.options_frame,
            text="Iteration number:",
            font=font)
        lab_font.grid(row=0, column=0, sticky=Tk.S)

        self.draw_slider()

        if self.show_wing and self.show_tube:
            # checkbox to show deformed mesh
            self.show_def_mesh = Tk.IntVar()
            c1 = Tk.Checkbutton(
                self.options_frame,
                text="Show deformed mesh",
                variable=self.show_def_mesh,
                command=self.update_graphs,
                font=font)
            c1.grid(row=0, column=2, padx=5, sticky=Tk.W)

            # checkbox to exaggerate deformed mesh
            self.ex_def = Tk.IntVar()
            self.c2 = Tk.Checkbutton(
                self.options_frame,
                text="Exaggerate deformations",
                variable=self.ex_def,
                command=self.update_graphs,
                font=font)
            self.c2.grid(row=0, column=3, padx=5, sticky=Tk.W)

        # Option to automatically refresh history file
        # especially useful for currently running optimizations
        self.var_ref = Tk.IntVar()
        # self.var_ref.set(1)
        c11 = Tk.Checkbutton(
            self.options_frame,
            text="Automatically refresh",
            variable=self.var_ref,
            command=self.auto_ref,
            font=font)
        c11.grid(row=0, column=4, sticky=Tk.W, pady=6)

        button = Tk.Button(
            self.options_frame,
            text='Save video',
            command=self.save_video,
            font=font)
        button.grid(row=0, column=5, padx=5, sticky=Tk.W)

        button4 = Tk.Button(
            self.options_frame,
            text='Save image',
            command=self.save_image,
            font=font)
        button4.grid(row=0, column=6, padx=5, sticky=Tk.W)

        button5 = Tk.Button(
            self.options_frame,
            text='Quit',
            command=self.quit,
            font=font)
        button5.grid(row=0, column=7, padx=5, sticky=Tk.W)

        self.auto_ref()

def disp_plot(db_name):
    disp = Display(db_name)
    disp.draw_GUI()
    plt.tight_layout()
    disp.root.protocol("WM_DELETE_WINDOW", disp.quit)
    Tk.mainloop()

if __name__ == '__main__':
    disp_plot(db_name)
