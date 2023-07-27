"""
Class definition for the Mphys builder for the aero solver.
"""

import copy
import warnings

import openmdao.api as om

from openaerostruct.mphys.utils import get_number_of_nodes, get_node_indices
from openaerostruct.mphys.aero_mesh import AeroMesh
from openaerostruct.mphys.demux_surface_mesh import DemuxSurfaceMesh
from openaerostruct.mphys.mux_surface_forces import MuxSurfaceForces
from openaerostruct.mphys.aero_solver_group import AeroSolverGroup
from openaerostruct.mphys.aero_funcs_group import AeroFuncsGroup

try:
    from mphys.builder import Builder
    from mphys.distributed_converter import DistributedConverter, DistributedVariableDescription

    mphys_found = True
except ImportError:
    mphys_found = False
    Builder = object


class AeroCouplingGroup(om.Group):
    """
    Group that wraps the aerodynamic states into the Mphys's broader coupling group.

    This is done in four steps:

        1. The deformed aero coordinates are read in as a distributed flattened array
        and split up into multiple 3D serial arrays (one per surface).

        2. The VLM problem is then solved based on the deformed mesh.

        3. The aerodynamic nodal forces for each surface produced by the VLM solver
        are concatonated into a flattened array.

        4. The serial force vector is converted to a distributed array and
        provided as output tothe rest of the Mphys coupling groups.
    """

    def initialize(self):
        self.options.declare("surfaces", default=None, desc="oas surface dicts", recordable=False)
        self.options.declare("compressible", default=True, desc="prandtl glauert compressibiity flag", recordable=True)

    def setup(self):
        self.surfaces = self.options["surfaces"]
        self.compressible = self.options["compressible"]

        self.set_input_defaults("aoa", val=0.0, units="deg")
        self.set_input_defaults("yaw", val=0.0, units="deg")
        if self.compressible:
            self.set_input_defaults("mach", val=0.0)

        nnodes = get_number_of_nodes(self.surfaces)

        # Convert distributed mphys mesh input into a serial vector OAS can use
        in_vars = [DistributedVariableDescription(name="x_aero", shape=(nnodes * 3), tags=["mphys_coordinates"])]

        self.add_subsystem("collector", DistributedConverter(distributed_inputs=in_vars), promotes_inputs=["x_aero"])
        self.connect("collector.x_aero_serial", "demuxer.x_aero")

        # Demux flattened surface mesh vector into seperate vectors for each surface
        self.add_subsystem(
            "demuxer",
            DemuxSurfaceMesh(surfaces=self.surfaces),
            promotes_outputs=["*_def_mesh"],
        )

        # OAS aero states group
        self.add_subsystem(
            "states",
            AeroSolverGroup(surfaces=self.surfaces, compressible=self.compressible),
            promotes_inputs=["*"],
            promotes_outputs=["*"],
        )

        # Mux all surface forces into one flattened array
        self.add_subsystem(
            "muxer",
            MuxSurfaceForces(surfaces=self.surfaces),
            promotes_inputs=["*_mesh_point_forces"],
        )

        # Convert serial force vector to distributed, like mphys expects
        out_vars = [DistributedVariableDescription(name="f_aero", shape=(nnodes * 3), tags=["mphys_coupling"])]

        self.add_subsystem(
            "distributor", DistributedConverter(distributed_outputs=out_vars), promotes_outputs=["f_aero"]
        )
        self.connect("muxer.f_aero", "distributor.f_aero_serial")


class AeroBuilder(Builder):
    """
    Mphys builder class responsible for setting up components of OAS's aerodynamic solver.
    """

    def_options = {"user_specified_Sref": False, "compressible": True, "output_dir": "./", "write_solution": True}

    def __init__(self, surfaces, options=None):
        if not mphys_found:
            raise ImportError(
                "MPhys is required in order to use the OpenAeroStruct mphys module. "
                + "Ensure MPhys is installed properly and can be found on your path."
            )
        self.surfaces = surfaces
        # Copy default options
        self.options = copy.deepcopy(self.def_options)
        # Update with user-defined options
        if options:
            self.options.update(options)

    def initialize(self, comm):
        self.comm = comm
        self.nnodes = get_number_of_nodes(self.surfaces)

    def get_coupling_group_subsystem(self, scenario_name=None):
        return AeroCouplingGroup(surfaces=self.surfaces, compressible=self.options["compressible"])

    def get_mesh_coordinate_subsystem(self, scenario_name=None):
        return AeroMesh(surfaces=self.surfaces)

    def get_post_coupling_subsystem(self, scenario_name=None):
        user_specified_Sref = self.options["user_specified_Sref"]
        return AeroFuncsGroup(
            surfaces=self.surfaces,
            write_solution=self.options["write_solution"],
            output_dir=self.options["output_dir"],
            user_specified_Sref=user_specified_Sref,
            scenario_name=scenario_name,
        )

    def get_ndof(self):
        """
        Tells Mphys this is a 3D problem.
        """
        return 3

    def get_number_of_nodes(self):
        """
        Get the number of nodes on root proc
        """
        if self.comm.rank == 0:
            return self.nnodes
        return 0

    def get_tagged_indices(self, tags):
        """
        Method that returns grid IDs for a list of body/boundary tags.

        Parameters
        ----------
        tags : list[str]
            list of surface names from which to pull the gridIDs

        Returns
        -------
        grid_ids : list[int]
            list of grid IDs that correspond to given body/boundary tags
        """
        if self.comm.rank == 0:
            if tags == -1 or tags == [-1]:
                return list(range(self.nnodes))
            else:
                all_surf_indices = get_node_indices(self.surfaces)
                tagged_indices = []
                for tag in tags:
                    if tag in all_surf_indices:
                        surf_indices = all_surf_indices[tag].flatten()
                        tagged_indices.extend(list(surf_indices))
                    else:
                        warnings.warn(
                            f'Tag name "{tag}" not found in list of added surfaces. Skipping tag.',
                            category=RuntimeWarning,
                            stacklevel=2,
                        )
                return tagged_indices
        return []
