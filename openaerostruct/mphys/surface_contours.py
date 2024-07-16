import numpy as np

from openmdao.api import ExplicitComponent
import os


class SurfaceContour(ExplicitComponent):
    """This is a post-processing component for writing out the aerodynamic
    solution of all lifting surfaces in a Tecplot format.
    The purpose of this component is to write visualization files, as such it
    has no output variables or sensitivities.

    Parameters
    ----------
    v : float
        Freestream air velocity in m/s.
    rho : float
        Air density in kg/m^3.
    def_mesh[nx, ny, 3] : numpy array
        Array defining the nodal coordinates of the lifting surface.
    circulations[(nx-1)*(ny-1)] : numpy array
        Flattened vector of horseshoe vortex strengths calculated by solving
        the linear system of AIC_mtx * circulations = rhs, where rhs is
        based on the air velocity at each collocation point.
    sec_forces[nx-1, ny-1, 3] : numpy array
        Contains the sectional forces acting on each panel.
    node_forces[nx, ny, 3] : numpy array
        Contains the aerodynamic forces acting on the nodes of the surface.
    """

    def initialize(self):
        self.options.declare("surfaces", types=list)
        self.options.declare("base_name", types=str, default="solution")
        self.options.declare("output_dir", types=str, default="./results")

    def setup(self):
        self.add_input("v", val=1.0, units="m/s", tags=["mphys_input"])
        self.add_input("rho", val=1.0, units="kg/m**3", tags=["mphys_input"])
        tot_panels = 0
        for surface in self.options["surfaces"]:
            name = surface["name"]
            mesh = surface["mesh"]
            nx, ny, _ = mesh.shape
            tot_panels += (nx - 1) * (ny - 1)
            self.add_input(name + "_def_mesh", val=mesh, units="m", tags=["mphys_coupling"])
            self.add_input(name + "_mesh_point_forces", val=np.ones([nx, ny, 3]), units="N", tags=["mphys_coupling"])
            self.add_input(name + "_sec_forces", val=np.ones([nx - 1, ny - 1, 3]), units="N", tags=["mphys_coupling"])
        self.add_input("circulations", val=np.zeros([tot_panels]), units="m**2/s", tags=["mphys_coupling"])
        self.solution_counter = 0
        # Prevents OpenMDAO from suppressing this component when using group_by_pre_opt_post feature
        self.options["always_opt"] = True

    def compute(self, inputs, outputs):
        """
        Compute the circulation and delta-C_P contours over each VLM surface and
        write to Tecplot file.
        """
        # We will store results in dictionary for each surface,
        # this will be used when we write out Tecplot file
        surf_circs = {}
        surf_mesh = {}
        surf_dCp = {}
        surf_forces = {}

        circ = inputs["circulations"]
        i = 0
        # Compute freestream dynamic pressure
        q = 0.5 * inputs["rho"] * inputs["v"] ** 2

        # Loop through each surface and get the circulation and delta_Cp
        for surface in self.options["surfaces"]:
            name = surface["name"]
            nx, ny, _ = surface["mesh"].shape
            num_panels = (nx - 1) * (ny - 1)

            # Unflatten the portion of the circulation vector corresponding to this surface
            surf_circs[name] = np.reshape(circ[i : i + num_panels], (nx - 1, ny - 1), order="C")
            surf_mesh[name] = inputs[name + "_def_mesh"]
            panel_forces = inputs[name + "_sec_forces"].reshape(nx - 1, ny - 1, 3)

            # Compute the normal vector of each panel by crossing the diagonals
            norm_vecs = np.cross(
                surf_mesh[name][:-1, 1:, :] - surf_mesh[name][1:, :-1, :],
                surf_mesh[name][:-1, :-1, :] - surf_mesh[name][1:, 1:, :],
                axis=2,
            )

            # Compute the area of each panel using the cross product magnitude
            panel_areas = np.sqrt(np.sum(norm_vecs**2, axis=2)) * 0.5

            # Normalize the normal vectors
            for j in range(3):
                norm_vecs[:, :, j] /= panel_areas * 2.0

            # Compute the difference in pressure between the upper and lower wing surface
            delta_p = np.zeros([nx - 1, ny - 1])
            for j in range(3):
                delta_p[:, :] += np.real(norm_vecs[:, :, j] * panel_forces[:, :, j] / panel_areas)

            # Normalize pressure change
            surf_dCp[name] = delta_p / q

            # Store forces at nodes
            surf_forces[name] = inputs[name + "_mesh_point_forces"]

            i += num_panels

        # Write out values to tecplot .plt file
        self.write_to_tecplot(surf_mesh, surf_circs, surf_dCp, surf_forces)
        self.solution_counter += 1

    def write_to_tecplot(self, meshes, circs, dCps, panel_forces):
        """
        Write circulation distribution as Tecplot .plt file.
        """
        if self.comm.rank == 0:
            # Now write out tecplot file
            file_name = self.options["base_name"] + "_%.3d_panel.plt" % (self.solution_counter)
            file_path = os.path.join(self.options["output_dir"], file_name)
            file_handle = open(file_path, "w")

            file_handle.write('TITLE = "OpenAeroStruct: Aerodynamic Solution"\n')
            file_handle.write('VARIABLES = "X", "Y", "Z", "Circulation", "delta_Cp", "FX", "FY", "FZ"\n')
            for name in meshes:
                mesh = meshes[name]
                circulations = circs[name]
                delta_Cps = dCps[name]
                surf_forces = panel_forces[name]
                nx = mesh.shape[0]
                ny = mesh.shape[1]

                file_handle.write(f"Zone T={name} I={nx}, J={ny},")
                file_handle.write("DATAPACKING=BLOCK, VARLOCATION=([4,5]=CELLCENTERED)\n")
                # Write out node locations
                for k in range(3):
                    for j in range(ny):
                        for i in range(nx):
                            file_handle.write("%f " % (mesh[i, j, k]))
                    file_handle.write("\n")
                # Write out panel circulations
                for j in range(ny - 1):
                    for i in range(nx - 1):
                        file_handle.write("%f " % (circulations[i, j]))
                file_handle.write("\n")
                # Write out panel delta_Cps
                for j in range(ny - 1):
                    for i in range(nx - 1):
                        file_handle.write("%f " % (delta_Cps[i, j]))
                file_handle.write("\n")
                # Write out panel forces
                for k in range(3):
                    for j in range(ny):
                        for i in range(nx):
                            file_handle.write("%f " % (surf_forces[i, j, k]))
                    file_handle.write("\n")
            file_handle.close()
