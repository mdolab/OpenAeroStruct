import numpy as np
import openmdao.api as om
from openaerostruct.utils.vector_algebra import compute_norm
from openaerostruct.utils.vector_algebra import compute_norm_deriv



def get_section_edge_left(mesh,v=np.ones(3),edge_cur=0,edges_all_constraints=np.ones(3)):
    nx = mesh.shape[0]
    ny = mesh.shape[1]
    le_index = 0
    te_index  = np.ravel_multi_index((nx-1,0,0),mesh.shape)
    mask = np.array(v, dtype='bool')

    rows = np.arange(0,2*np.sum(v)) + 2*int(np.sum(edges_all_constraints[:edge_cur]))
    cols = np.concatenate([np.arange(le_index,le_index+3)[mask],np.arange(te_index,te_index+3)[mask]])

    return mesh[[0,-1],0][:,np.arange(0,3)[mask]], rows, cols


def get_section_edge_right(mesh,v=np.ones(3),edge_cur=0,edges_all_constraints=np.ones(3)):
    nx = mesh.shape[0]
    ny = mesh.shape[1]
    le_index = np.ravel_multi_index((0,ny-1,0),mesh.shape)
    te_index = np.ravel_multi_index((nx-1,ny-1,0),mesh.shape)
    mask = np.array(v, dtype='bool')

    rows = np.arange(0,2*np.sum(v)) + 2*int(np.sum(edges_all_constraints[:edge_cur]))
    cols = np.concatenate([np.arange(le_index,le_index+3)[mask],np.arange(te_index,te_index+3)[mask]])
   
    return mesh[[0,-1],-1][:,np.arange(0,3)[mask]], rows, cols



class GeomMultiJoin(om.ExplicitComponent):
    """

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
        self.options.declare("sections", types=list, desc="A list of section surface dictionaries to be joined.")
        self.options.declare("dim_constr",types=list, default=[np.ones(3)],desc="A list of vectors of length three corresponding to each edge. Entries corresponding the dimension([x,y,z]) the user wishes to constraint should be set to 1. Remaining entries should be zero.")

    def setup(self):
        sections = self.options["sections"]
        self.num_sections = len(sections)
        self.dim_constr = self.options["dim_constr"]
        edgeTotal = self.num_sections-1
        
        if len(self.dim_constr) != (edgeTotal):
            self.dim_constr = [np.array([1,0,0]) for i in range(edgeTotal)]

        constr_size = 2*np.count_nonzero(np.concatenate(self.dim_constr))
        self.add_output("section_separation", val=np.zeros(constr_size))

        edge_cur = 0
        for iSec, section in enumerate(sections):
            mesh = section["mesh"]
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = section["name"]

            mesh_name = "{}_join_mesh".format(name)
            self.add_input(mesh_name, shape=(nx, ny, 3), units="m")

            if section["symmetry"]:
                left_wing = abs(mesh[0, 0, 1]) > abs(mesh[0, -1, 1])
                if left_wing:
                    if iSec == 0:
                        rows,cols = get_section_edge_left(mesh,self.dim_constr[iSec],edge_cur,self.dim_constr)[1:]
                        vals = -1*np.ones_like(rows)
                    elif iSec < len(sections) - 1:
                        rows1,cols1 = get_section_edge_right(mesh,self.dim_constr[iSec-1],edge_cur,self.dim_constr)[1:]
                        vals1 = np.ones_like(rows1)

                        edge_cur += 1
                        rows2, cols2 = get_section_edge_left(mesh,self.dim_constr[iSec],edge_cur,self.dim_constr)[1:]
                        vals2 = -1*np.ones_like(rows2)

                        rows = np.concatenate([rows1,rows2])
                        cols = np.concatenate([cols1,cols2])
                        vals = np.concatenate([vals1,vals2])
                    else:
                        rows, cols = get_section_edge_right(mesh,self.dim_constr[iSec-1],edge_cur,self.dim_constr)[1:]
                        vals = np.ones_like(rows)
                else:
                    if iSec == 0:
                        rows,cols = get_section_edge_right(mesh,self.dim_constr[iSec],edge_cur,self.dim_constr)[1:]
                        vals = -1*np.ones_like(rows)
                    elif iSec < len(sections) - 1:
                        rows,cols = get_section_edge_left(mesh,self.dim_constr[iSec-1],edge_cur,self.dim_constr)[1:]
                        vals1 = np.ones_like(rows1)

                        edge_cur += 1
                        rows,cols = get_section_edge_right(mesh,self.dim_constr[iSec],edge_cur,self.dim_constr)[1:]
                        vals2 = -1*np.ones_like(rows2)

                        rows = np.concatenate([rows1,rows2])
                        cols = np.concatenate([cols1,cols2])
                        vals = np.concatenate([vals1,vals2])
                    else:
                        rows, cols = get_section_edge_left(mesh,self.dim_constr[iSec-1],edge_cur,self.dim_constr)[1:]
                        vals = np.ones_like(rows)
                        
            
            self.declare_partials("section_separation",mesh_name,rows=rows,cols=cols,val=vals)
            #self.declare_partials("section_separation", mesh_name, method='cs')
        #self.declare_partials(of="*", wrt="*", method="cs")

            
        
    def compute(self, inputs, outputs):
        sections = self.options["sections"]
        edges = []
        edge_constraints = []
        for iSec, section in enumerate(sections):
            name = section["name"]
            mesh_name = "{}_join_mesh".format(name)
            if section["symmetry"]:
                left_wing = abs(inputs[mesh_name][0, 0, 1]) > abs(inputs[mesh_name][0, -1, 1])
                if left_wing:
                    if iSec == 0:
                        edges.append(get_section_edge_left(inputs[mesh_name],self.dim_constr[iSec])[0])
                    elif iSec < len(sections) - 1:
                        edges.append(get_section_edge_right(inputs[mesh_name],self.dim_constr[iSec-1])[0])
                        edges.append(get_section_edge_left(inputs[mesh_name],self.dim_constr[iSec])[0])
                    else:
                        edges.append(get_section_edge_right(inputs[mesh_name],self.dim_constr[iSec-1])[0])
                else:
                    if iSec == 0:
                        edges.append(get_section_edge_right(inputs[mesh_name],self.dim_constr[iSec])[0])
                    elif iSec < len(sections) - 1:
                        edges.append(get_section_edge_left(inputs[mesh_name],self.dim_constr[iSec-1])[0])
                        edges.append(get_section_edge_right(inputs[mesh_name],self.dim_constr[iSec])[0])
                    else:
                        edges.append(get_section_edge_left(inputs[mesh_name],self.dim_constr[iSec-1])[0])
            else:
                raise Exception("Joining requires symmetry on for now.")
            
        for i in range(self.num_sections - 1):
            edge_constraints.append((edges[2*i+1] - edges[2*i]).flatten())

        outputs["section_separation"] = np.array(edge_constraints).flatten()