import numpy as np

from openmdao.api import ExplicitComponent
import os


class LiftDistribution(ExplicitComponent):
    """This is a post-processing component for writing out the normalized lift
    distribution over all lifting surfaces in a Tecplot format.
    The purpose of this component is to write visualization files, as such it
    has no output variables or sensitivities.

    Parameters
    ----------
    def_mesh[nx, ny, 3] : numpy array
        Array defining the nodal coordinates of the lifting surface.
    sec_forces[nx-1, ny-1, 3] : numpy array
        Contains the sectional forces acting on each panel.
        Stored in Fortran order (only relevant with more than one chordwise
        panel).
    circulations[(nx-1)*(ny-1)] : numpy array
        Flattened vector of horseshoe vortex strengths calculated by solving
        the linear system of AIC_mtx * circulations = rhs, where rhs is
        based on the air velocity at each collocation point.
    alpha : float
        Angle of attack in degrees.
    """

    def initialize(self):
        self.options.declare("surfaces", types=list)
        self.options.declare("Nspan", types=int, default=60)
        self.options.declare("base_name", types=str, default="solution")
        self.options.declare("output_dir", types=str, default="./results")

    def setup(self):
        self.N = self.options["Nspan"]
        # Inputs
        self.add_input("alpha", val=3.0, units="deg", tags=["mphys_input"])
        for surface in self.options["surfaces"]:
            name = surface["name"]
            mesh = surface["mesh"]
            nx, ny, _ = mesh.shape
            self.add_input(name + "_def_mesh", val=mesh, units="m", tags=["mphys_coupling"])
            self.add_input(name + "_sec_forces", val=np.ones((nx - 1, ny - 1, 3)), units="N", tags=["mphys_coupling"])
        # Prevents OpenMDAO from suppressing this component when using group_by_pre_opt_post feature
        self.options["always_opt"] = True

    def compute(self, inputs, outputs):
        """
        Take the lift distribution over each surface and combine it into one
        curve over the whole aircraft.
        """

        # Input parameters
        alpha = inputs["alpha"] * np.pi / 180.0
        cosa = np.cos(alpha)
        sina = np.sin(alpha)

        # Compute the minimum and maximum span coordinate over all surfaces
        y_min = np.inf
        y_max = -np.inf
        for surface in self.options["surfaces"]:
            name = surface["name"]
            def_mesh = inputs[name + "_def_mesh"]
            y_station = def_mesh[0, :, 1]
            if y_min > min(y_station):
                y_min = min(y_station)
            if y_max < max(y_station):
                y_max = max(y_station)

        # Now compute the total lift distribution from y_min to y_max
        y = np.linspace(y_min, y_max, self.N)
        total_lift_dist = np.zeros(self.N)
        # Here we loop through each lifting surface and add its portion of the
        # lift distribution over its portion of the span
        for surface in self.options["surfaces"]:
            name = surface["name"]
            def_mesh = inputs[name + "_def_mesh"]
            sec_forces = inputs[name + "_sec_forces"]

            # Define the y stations we will be interpolating between as the
            # mid-points of the panels
            y_station = np.zeros([def_mesh.shape[1] + 1])
            y_station[1:-1] = (def_mesh[0, 1:, 1] + def_mesh[0, :-1, 1]) / 2
            # Add the first and last node along the span for extrapolating
            # values near the root and stip of the surface
            y_station[0] = def_mesh[0, 0, 1]
            y_station[-1] = def_mesh[0, -1, 1]

            # Now compute the panel widths in the span (y) direction
            # Compute the widths of each panel at the quarter-chord line
            quarter_chord = 0.25 * def_mesh[-1] + 0.75 * def_mesh[0]
            widths = np.abs(quarter_chord[1:, 1] - quarter_chord[:-1, 1])

            # Now compute the lift distribution at the mid-point of each panel
            # Lift distribution: dimensional l(y) = -Fx(y) sin(alpha) + Fz(y) cos(alpha) / widths(y)
            forces = np.sum(sec_forces, axis=0)  # sum section forces in the chordwise x-direction: forces(ny,3)
            lift_dist = (-forces[:, 0] * sina + forces[:, 2] * cosa) / widths[:]

            # Now define the lift distribution of each mid-point from the
            # input variable
            surface_lift_dist = np.zeros_like(y_station)
            surface_lift_dist[1:-1] = lift_dist[:]
            # Extrapolate a constant value to the first and last node
            surface_lift_dist[0] = surface_lift_dist[1]
            surface_lift_dist[-1] = surface_lift_dist[-2]
            # Add the surfaces contribution to the total lift distribution
            total_lift_dist += np.interp(y, y_station, surface_lift_dist, left=0, right=0)

        # Compute the normalized lift distribution by integrating
        total_lift = np.trapz(total_lift_dist, y)
        span = y_max - y_min
        # Normalize the lift distribution so that the area under the curve is
        # unity
        norm_lift_dist = total_lift_dist / total_lift * span
        # Finaly write result to tecplot file
        self.writeToTecplot(y, norm_lift_dist)

    def writeToTecplot(self, y, ld):
        """
        Write lift distribution as Tecplot .dat file.
        """
        if self.comm.rank == 0:
            # Compute the normalized span coordinate
            eta = y / max(abs(y))
            # Compute an elliptical lift distribution w/ unit area for reference
            elliptical = np.sqrt(1 - eta**2) * (4 / np.pi)

            # Now write out tecplot file
            fileName = self.options["base_name"] + "_%.3d_lift.dat" % (self.iter_count)
            filePath = os.path.join(self.options["output_dir"], fileName)
            file_handle = open(filePath, "w")

            file_handle.write('TITLE = "OpenAeroStruct: Lift Distribution Plot"\n')
            file_handle.write('VARIABLES = "eta", "Y", "Normalized Lift", "elliptical"\n')
            file_handle.write('ZONE T="Lift Dist", I=%d, F=POINT\n' % (len(y)))
            for i in range(len(y)):
                file_handle.write("%E %E %E %E\n" % (eta[i], y[i], ld[i], elliptical[i]))
            file_handle.close()
