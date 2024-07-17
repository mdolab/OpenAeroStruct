import numpy as np
import openmdao.api as om

#Function computes the dimensions of the unified mesh as well as an index array
def compute_unimesh_dims(sections):
    uni_nx = sections[0]["mesh"].shape[0]
    uni_ny = 0
    for iSec in range(len(sections)):
        if iSec == 0:
            uni_ny += sections[iSec]["mesh"].shape[1]
        elif sections[iSec]["symmetry"] and iSec != 0:
            uni_ny += (sections[iSec]["mesh"].shape[1] - 1)
        else:
            uni_ny += (sections[iSec]["mesh"].shape[1] - 2)
    uni_mesh_indices = np.arange(uni_nx * uni_ny * 3).reshape((uni_nx, uni_ny, 3))

    return uni_mesh_indices,uni_nx,uni_ny


#Function that computes the index block that corresponds to each individual wing section
def compute_unimesh_index_blocks(sections,uni_mesh_indices):
    blocks = []

    #cursor to track the y position of each section along the unified mesh
    ycurr = -1
    for iSec in range(len(sections)):
        mesh = sections[iSec]["mesh"]
        nx = mesh.shape[0]
        ny = mesh.shape[1]

        if sections[iSec]['symmetry']:
            if iSec == 0:
                block = uni_mesh_indices[:,ycurr:ycurr-ny:-1,:][:,::-1,:]
                ycurr -= ny
            else:
                block = uni_mesh_indices[:,ycurr:ycurr-(ny-1):-1,:][:,::-1,:]
                ycurr -= ny-1
        else:
            raise ValueError("Mesh unification not supported without symmetry turned on. For now.")
            if iSec == 0:
                block = uni_mesh_indices[:,-1:-1-ny:-1,:]
            else:
                block = uni_mesh_indices[:,-1-(iSec*ny):-1-(iSec*ny)-(ny-1):-1,:]

        blocks.append(block.flatten())

    return blocks

#Function that produces a unified mesh from all the individual wing sections meshes
def unify_mesh(sections):
    for isec,section in enumerate(sections):
            mesh = section["mesh"]
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = section["name"]
            import copy
            # Stitch the results into a singular mesh
            if isec == 0:
                uniMesh = copy.deepcopy(mesh)
            else:
                if section["symmetry"]:
                    uniMesh = np.concatenate((mesh[:,0 : ny-1,:], uniMesh),axis=1)
                else:
                    uniMesh = np.concatenate((mesh[:,0 : ny/2-1,:], uniMesh, mesh[:,ny/2+1:,:]),axis=1)
    return uniMesh


class GeomMultiUnification(om.ExplicitComponent):
    """
    OpenMDAO component that combines the meshes associated with each individual section
    of a multi-section wing into a unified mesh. Section 0 as defined by the user(typically
    the most inboard section) will be added to the unified mesh unchanged. The remaining sections
    are then sequencially added outboard to section 0 along the y-direction. The mesh unification
    component operates under the assumption that it is unifying section that are C0 continous along
    the span. The most inboard column of points of section i+1 will be coincident with the 
    most outbord column of points of section i. As result, the inboard column of points for each section
    after section 0 are removed prior to the section being added to the unified mesh. The component also 
    produces the appropriate sparse Jacobian matrix for the column removal so that the partials can be computed
    appropriately.

    Parameters
    ----------
    sections : list
        A list of section surface dictionaries in OAS format.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the unified surface
    """
    def initialize(self):
        """
        Declare options.
        """
        self.options.declare("sections", types=list, desc="A list of section surface dictionaries to be unified.")
        self.options.declare("surface_name", types=str, desc="The name of the multi-section surface")

    def setup(self):
        sections = self.options["sections"]
        name = self.options["surface_name"]

        #Get the unified mesh size, index array, and block of indicies for each section
        [uni_mesh_indices, uni_nx, uni_ny] = compute_unimesh_dims(sections)
        uni_mesh_blocks = compute_unimesh_index_blocks(sections,uni_mesh_indices)
        uni_mesh_name = '{}_uni_mesh'.format(name)

        #Loop through each section to build the sparsity pattern for the Jacobian. This Jacobian is fixed based on the mesh size so can be declared here
        for iSec,section in enumerate(sections):
            mesh = section["mesh"]
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = section["name"]

            mesh_name = "{}_def_mesh".format(name)

            self.add_input(mesh_name, shape=(nx, ny, 3), units="m", tags=["mphys_coupling"])
            
            if section["symmetry"]:
                left_wing = abs(section["mesh"][0, 0, 1]) > abs(section["mesh"][0, -1, 1])
                
                #Generate index array
                mesh_indices = np.arange(nx * ny * 3).reshape((nx, ny, 3))

                if iSec == 0:
                    cols = mesh_indices.flatten()
                else:
                    #If not section 0 then the most inboard column is disregarded
                    if left_wing:
                        cols = mesh_indices[:,:-1,:].flatten()
                    else:
                        cols = mesh_indices[:,1:,:].flatten()

                #Get data from section block in unified mesh
                rows = uni_mesh_blocks[iSec]

                #Fill non zero Jacobian entries with ones
                data = np.ones_like(rows)

                self.declare_partials(uni_mesh_name, mesh_name, val=data, rows=rows, cols=cols)
                #self.declare_partials(uni_mesh_name, mesh_name, method='cs')
            else:
                
                mesh_indices = np.arange(nx * ny * 3).reshape((nx, ny, 3))
                if iSec == 0:
                    cols = mesh_indices.flatten()
                else:
                    cols = np.concatenate((mesh_indices[:,:ny/2-1,:].flatten(),mesh_indices[:,ny/2+1:,:].flatten()))

                self.add_output(uni_mesh_name, shape=(uni_nx, uni_ny * 2, 3), units="m")


                raise ValueError("Mesh unification not supported without symmetry turned on.")
        
        self.add_output(uni_mesh_name, shape=(uni_nx, uni_ny, 3), units="m")


    def compute(self, inputs, outputs):
        sections = self.options["sections"]
        surface_name = self.options["surface_name"]
        uni_mesh_name = '{}_uni_mesh'.format(surface_name)
        
        #Loop through all sections to unify the mesh
        for isec,section in enumerate(sections):
            mesh = section["mesh"]
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = section["name"]
            mesh_name = "{}_def_mesh".format(name)

            #Build the unified mesh
            if isec == 0:
                uniMesh = inputs[mesh_name]
            else:
                #If not section 0 then the most inboard column is disregarded
                if section["symmetry"]:
                    left_wing = abs(inputs[mesh_name][0, 0, 1]) > abs(inputs[mesh_name][0, -1, 1])
                    if left_wing:
                        uniMesh = np.concatenate((inputs[mesh_name][:,0 : ny-1,:], uniMesh),axis=1)
                    else:
                        uniMesh = np.concatenate((uniMesh,inputs[mesh_name][:,1 : ny,:]),axis=1)
                else:
                    uniMesh = np.concatenate((inputs[mesh_name][:,0 : ny/2-1,:], uniMesh, inputs[mesh_name][:,ny/2+1:,:]),axis=1)
        outputs[uni_mesh_name] = uniMesh

    def compute_partials(self, inputs, J):
        #Dummy code. Not needed for now.
        sections = self.options["sections"]
        for section in sections:
            mesh = section["mesh"]
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = section["name"]

            mesh_name = "{}_def_mesh".format(name)
            uni_mesh_name = "{}_uni_mesh".format(name)
            continue
