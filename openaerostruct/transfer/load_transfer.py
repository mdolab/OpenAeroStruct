from __future__ import division, print_function
import numpy as np
from scipy.sparse import coo_matrix, hstack, identity

from openmdao.api import ExplicitComponent

data_type = complex

def _skew(vector):
    out = np.array([[0, -vector[2], vector[1]],\
                    [vector[2], 0, -vector[0]],\
                    [-vector[1], vector[0], 0]])
    return out
        
class LoadTransfer(ExplicitComponent):
    """
    Perform aerodynamic load transfer.

    Apply the computed sectional forces on the aerodynamic surfaces to
    obtain the deformed mesh FEM loads.

    Parameters
    ----------
    def_mesh[nx, ny, 3] : numpy array
        Flattened array defining the lifting surfaces after deformation.
    sec_forces[nx-1, ny-1, 3] : numpy array
        Flattened array containing the sectional forces acting on each panel.
        Stored in Fortran order (only relevant when more than one chordwise
        panel).

    Returns
    -------
    LoadsA[ny,3] : numpy array  = loads[:,:3]
    LoadsB[ny,3]: nympy array   = loads[:,3:]
        Flattened array containing the loads applied on the FEM component,
        computed from the sectional forces.
        
    Loads[ny, 6] : numpy array = [LoadsA,LoadsB]
    
    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        self.surface = surface = self.options['surface']

        self.ny = ny = surface['num_y']
        self.nx = nx = surface['num_x']


        if surface['fem_model_type'] == 'tube':
            self.fem_origin = surface['fem_origin']
        else:
            y_upper = surface['data_y_upper']
            x_upper = surface['data_x_upper']
            y_lower = surface['data_y_lower']

            self.fem_origin = (x_upper[0]  * (y_upper[0]  - y_lower[0]) +
                               x_upper[-1] * (y_upper[-1] - y_lower[-1])) / \
                             ((y_upper[0]  -  y_lower[0]) + (y_upper[-1] - y_lower[-1]))
        self.w1 = w1 = 0.25
        self.w2 = w2 = self.fem_origin

        self.add_input('def_mesh', val=np.random.random((nx, ny, 3)), units='m')
        self.add_input('sec_forces', val=np.random.random((nx-1, ny-1, 3)), units='N')

        self.add_output('loadsA', val=np.zeros((ny, 3)), units='N')
        self.add_output('loadsB', val=np.zeros((ny, 3)), units='N*m')
        self.add_output('loads', val=np.zeros((self.ny,6)), units='N') ## WARNING!!! UNITS ARE A MIXTURE OF N & N*m
        # Well, technically the units of this load array are mixed.
        # The first 3 indices are N and the last 3 are N*m.
        
        
        # Below we create each individual partial derivative
        
        # --------------------------------dloadsA__dsec_forces-------------------------------------
        # this derivative is linear
        rowsA=np.zeros((3*(ny-1)*(nx-1)*2))  
        colsA=np.zeros((3*(ny-1)*(nx-1)*2))
        for x_idx in range(nx-1):
            offset = x_idx*(3*(ny-1))
            for col_idx in range(3*(ny-1)):
                rowsA[2*offset+col_idx*2] = col_idx            
                rowsA[2*offset+col_idx*2+1] = col_idx+3 
                colsA[2*offset+col_idx*2] = offset+col_idx
                colsA[2*offset+col_idx*2+1] = offset+col_idx
        dloadsA__dsec_forces = coo_matrix((0.5*np.ones(len(rowsA)), (rowsA, colsA)), shape=(3*(ny), 3*(ny-1)*(nx-1)))   
        self.declare_partials(of='loadsA', wrt='sec_forces', rows=rowsA, cols=colsA, val=0.5)
        
        # --------------------------------dloadsB__dsec_forces-------------------------------------
        #dloadsB__dmoment
        rowsA=np.zeros((3*(ny-1)*2))  
        colsA=np.zeros((3*(ny-1)*2))
        data=0.5*np.ones((3*(ny-1)*2))
        for col_idx in range(3*(ny-1)):
            rowsA[col_idx*2] = col_idx            
            rowsA[col_idx*2+1] = col_idx+3 
            colsA[col_idx*2] = col_idx
            colsA[col_idx*2+1] = col_idx
        self.dloadsB__dmoment = coo_matrix((data, (rowsA, colsA)), shape=(3*ny,3*(ny-1)))
        
        # dmoment__dsec_forces
        diff = np.ones((nx-1, ny-1, 3))
        rows = np.zeros((9*(ny-1)*(nx-1)))
        cols = np.zeros((9*(ny-1)*(nx-1)))
        data = np.ones((9*(ny-1)*(nx-1)))
        tmp_row = np.array([0,0,0,1,1,1,2,2,2])
        tmp_col = np.array([0,1,2,0,1,2,0,1,2])
        # setup rows and cols
        for nxitr in range(nx-1):
            for nyitr in range(ny-1):
                index = 9*(nyitr+(ny-1)*nxitr)     
                rows[index:index+9] = tmp_row+nyitr*3.    
                cols[index:index+9] = tmp_col+nyitr*3.+3*(ny-1)*nxitr 
                data[index:index+9] = -_skew(diff[nxitr,nyitr,:]).flatten()
        dmoment__dsec_forces = coo_matrix((data, (rows,cols)), shape=(3*(ny-1), 3*(ny-1)*(nx-1)))
        
        # dloadsB__dsec_forces
        dloadsB__dsec_forces = self.dloadsB__dmoment * dmoment__dsec_forces
        dloadsB__dsec_forces = dloadsB__dsec_forces.tocoo() # force the system to be a coo_matrix     
        self.declare_partials(of='loadsB', wrt='sec_forces', rows=dloadsB__dsec_forces.row, cols=dloadsB__dsec_forces.col)

        # --------------------------------dloadsB__ddef_mesh-------------------------------------
        # setup sparce partials for dloadsB__ddef_mesh where:
        # dloadsB__ddef_mesh = dloadsB__dmoment * dmoment__ddiff * ddiff__ddef_mesh

        # dloadsB__dmoment
        rows = np.zeros((3*(ny-1)*2))
        cols = np.zeros((3*(ny-1)*2))
        data = np.ones((3*(ny-1)*2))
        tmp_row = np.array([0,3])
        tmp_col = np.array([0,0])
        for i in range(3*(ny-1)):
            rows[i*2:(1+i)*2] = i+tmp_row
            cols[i*2:(1+i)*2] = i+tmp_col
        self.dloadsB__dmoment = 0.5*coo_matrix((data, (rows,cols)), shape=(3*ny, 3*(ny-1)))
        
        # dmoment__ddiff
        # this will need to be recalculated in def_partials b/c it's a function of sec_forces
        rows = np.zeros((9*(ny-1)*(nx-1)))
        cols = np.zeros((9*(ny-1)*(nx-1)))
        data = np.zeros((9*(ny-1)*(nx-1)))
        tmp_row = np.array([0,0,0,1,1,1,2,2,2])
        tmp_col = np.array([0,1,2,0,1,2,0,1,2])
        tmp_sec_forces = np.ones((nx-1,ny-1,3))
        # setup rows and cols
        for nxitr in range(nx-1):
            for nyitr in range(ny-1):
                index = 9*(nyitr+(ny-1)*nxitr)     
                rows[index:index+9] = tmp_row+nyitr*3.    
                cols[index:index+9] = tmp_col+nyitr*3.+3*(ny-1)*nxitr 
                data[index:index+9] = -_skew(tmp_sec_forces[nxitr,nyitr,:]).flatten()
        dmoment__ddiff = coo_matrix((data, (rows,cols)), shape=(3*(ny-1), 3*(ny-1)*(nx-1)))
        
        # ddiff__ds_pts
        tmp_ddiff__ds_pts = np.zeros((3*(ny-1)*(nx-1),3*(ny-1)))
        for xiter in range(nx-1):
            tmp_ddiff__ds_pts[xiter*3*(ny-1):(1+xiter)*3*(ny-1),0:3*(ny-1)] = np.eye((3*(ny-1)))
        ddiff__ds_pts = coo_matrix(tmp_ddiff__ds_pts)
        
        # ds_pts__ddef_mesh
        tmp_ds_pts__ddef_mesh = np.zeros((3*(ny-1),3*(ny)*(nx)))
        block = np.zeros((3*(ny-1),3*ny))
        block[:,0:-3] = np.eye((3*(ny-1)))
        block[:,3:] = block[:,3:] + np.eye((3*(ny-1)))
        tmp_ds_pts__ddef_mesh[:,0:3*ny] = 0.5*(1-w2)*block
        tmp_ds_pts__ddef_mesh[:,3*ny*(nx-1):] = tmp_ds_pts__ddef_mesh[:,3*ny*(nx-1):] + 0.5 * w2 * block
        ds_pts__ddef_mesh = coo_matrix(tmp_ds_pts__ddef_mesh)

        # ddiff__da_pts
        ddiff__da_pts = identity((3*(ny-1)*(nx-1)))

        # da_pts__ddef_mesh
        tmp_da_pts__ddef_mesh = np.zeros((3*(ny-1)*(nx-1),3*ny*nx))
        block2 = np.zeros((3*(ny-1),6*ny))
        block2[:,:-3*ny] = 0.5 * (1-w1) * block
        block2[:,3*ny:] = 0.5 * w1 * block
        down = 3*(ny-1)
        over = 3*ny
        for xiter in range(nx-1):
            tmp_da_pts__ddef_mesh[xiter*down:xiter*down+3*(ny-1),over*xiter:over*xiter+6*ny] = block2
        da_pts__ddef_mesh = coo_matrix(tmp_da_pts__ddef_mesh)
        
        # to shorten future calculation times we calculate two other partials
        # ddiff__ddef_mesh
        self.ddiff__ddef_mesh = ddiff__da_pts * da_pts__ddef_mesh - ddiff__ds_pts * ds_pts__ddef_mesh
        self.ddiff__ddef_mesh = self.ddiff__ddef_mesh.tocoo()
        
        # dloadsB__ddef_mesh
        dloadsB__ddef_mesh = self.dloadsB__dmoment * dmoment__ddiff * self.ddiff__ddef_mesh
        #dloadsB__ddef_mesh = dloadsB__dmoment * dmoment__ddiff * ddiff__ddef_mesh
        dloadsB__ddef_mesh = dloadsB__ddef_mesh.tocoo()
        rows = dloadsB__ddef_mesh.row
        cols = dloadsB__ddef_mesh.col
        self.declare_partials(of='loadsB', wrt='def_mesh', rows=rows, cols=cols)
        
        #---------------------------------dloads__ddef_mesh
        #dloads__dloadsA
        blk = np.array([0,1,2])
        rows = np.zeros((ny*3))
        data = np.ones((ny*3))
        cols = np.linspace(0,ny*3-1, ny*3)
        for yiter in range(ny):
            rows[yiter*3:(1+yiter)*3] = blk+6*yiter
        dloads__dloadsA = coo_matrix((data, (rows, cols)), shape=(ny*6, ny*3))
        
        #dloads__dloadsB
        rows = rows + 3
        self.dloads__dloadsB = coo_matrix((data, (rows, cols)), shape=(ny*6, ny*3))
        
        dloads__ddef_mesh = self.dloads__dloadsB * dloadsB__ddef_mesh # + self.dloads__dloadsA * dloadsA__ddef_mesh but this last term is zero
        dloads__ddef_mesh = dloads__ddef_mesh.tocoo()
        rows = dloads__ddef_mesh.row
        cols = dloads__ddef_mesh.col
        self.declare_partials(of='loads', wrt='def_mesh', rows=rows, cols=cols)
        
        
        #---------------------------------dloads__dsec_forces------------------------------------
        dloads__dsec_forces = dloads__dloadsA * dloadsA__dsec_forces + self.dloads__dloadsB * dloadsB__dsec_forces
        dloads__dsec_forces = dloads__dsec_forces.tocoo()
        rows = dloads__dsec_forces.row
        cols = dloads__dsec_forces.col
        self.declare_partials(of='loads', wrt='sec_forces', rows=rows, cols=cols)
        
        
        # -------------------------------- Check Partial Options-------------------------------------
        self.set_check_partial_options('*', method='cs', step=1e-40)


        #---------------------------------- dloads__dloadsA__dsec_forces
        self.dloads__dloadsA__dsec_forces = dloads__dloadsA * dloadsA__dsec_forces

    def compute(self, inputs, outputs):
        mesh = inputs['def_mesh'] #[nx, ny, 3]
        sec_forces = inputs['sec_forces']
        
        # Compute the aerodynamic centers at the quarter-chord point of each panel
        # a_pts [nx-1, ny-1, 3]
        a_pts = 0.5 * (1-self.w1) * mesh[:-1, :-1, :] + \
                0.5 *   self.w1   * mesh[1:, :-1, :] + \
                0.5 * (1-self.w1) * mesh[:-1,  1:, :] + \
                0.5 *   self.w1   * mesh[1:,  1:, :]
        # Compute the structural midpoints based on the fem_origin location        
        # s_pts [ny-1, 3]
        s_pts = 0.5 * (1-self.w2) * mesh[0, :-1, :] + \
                0.5 *   self.w2   * mesh[-1, :-1, :] + \
                0.5 * (1-self.w2) * mesh[0,  1:, :] + \
                0.5 *   self.w2   * mesh[-1,  1:, :]
        # Find the moment arm between the aerodynamic centers of each panel
        # and the FEM elements
        # diff [nx-1, ny-1, 3]
        diff = a_pts - s_pts
        moment = np.zeros((self.ny - 1, 3), dtype=np.complex128)
        for ind in range(self.nx-1):
            moment = moment + np.cross(diff[ind, :, :], sec_forces[ind, :, :], axis=1)

        # Compute the loads based on the xyz forces and the computed moments
        loadsA = outputs['loadsA']
        loadsA[:] = 0. # apparently sometimes these don't get zero'd out so we'll do it explicitly
        sec_forces_sum = np.sum(sec_forces, axis=0)
        loadsA[:-1, :] = 0.5 * sec_forces_sum[:, :]
        loadsA[ 1:, :] = loadsA[ 1:, :] + 0.5 * sec_forces_sum[:, :]
        
        loadsB = outputs['loadsB']
        loadsB[:] = 0. # apparently sometimes these don't get zero'd out so we'll do it explicitly
        loadsB[:-1, :] = 0.5 * moment
        loadsB[ 1:, :] = loadsB[ 1:, :] + 0.5 * moment

        outputs['loadsA'] = loadsA
        outputs['loadsB'] = loadsB
        
        outputs['loads'][:,:3] = loadsA # everything on the first 3 columns  
        outputs['loads'][:,3:] = loadsB # everything on the last 3 columns
        
    def compute_partials(self, inputs, J):
        mesh = inputs['def_mesh']
        sec_forces = inputs['sec_forces']
        ny = self.ny
        nx = self.nx
        
        # -------------------------------dloadsB__dsec_forces--------------------------------------
        # dmoment__dsec_forces
        # Compute the aerodynamic centers at the quarter-chord point of each panel
        a_pts = 0.5 * (1-self.w1) * mesh[:-1, :-1, :] + \
                0.5 *   self.w1   * mesh[1:, :-1, :] + \
                0.5 * (1-self.w1) * mesh[:-1,  1:, :] + \
                0.5 *   self.w1   * mesh[1:,  1:, :]
        # Compute the structural midpoints based on the fem_origin location
        s_pts = 0.5 * (1-self.w2) * mesh[0, :-1, :] + \
                0.5 *   self.w2   * mesh[-1, :-1, :] + \
                0.5 * (1-self.w2) * mesh[0,  1:, :] + \
                0.5 *   self.w2   * mesh[-1,  1:, :]
        # Find the moment arm between the aerodynamic centers of each panel
        # and the FEM elements
        diff = a_pts - s_pts
        rows = np.zeros((9*(ny-1)*(nx-1)))
        cols = np.zeros((9*(ny-1)*(nx-1)))
        data = np.zeros((9*(ny-1)*(nx-1)))
        tmp_row = np.array([0,0,0,1,1,1,2,2,2])
        tmp_col = np.array([0,1,2,0,1,2,0,1,2])
        # setup rows and cols
        for nxitr in range(nx-1):
            for nyitr in range(ny-1):
                index = 9*(nyitr+(ny-1)*nxitr)     
                rows[index:index+9] = tmp_row+nyitr*3.    
                cols[index:index+9] = tmp_col+nyitr*3.+3*(ny-1)*nxitr 
                data[index:index+9] = _skew(diff[nxitr,nyitr,:]).flatten()
        dmoment__dsec_forces = coo_matrix((data, (rows,cols)), shape=(3*(ny-1), 3*(ny-1)*(nx-1)))
        
        dloadsB__dsec_forces = self.dloadsB__dmoment * dmoment__dsec_forces
        dloadsB__dsec_forces = dloadsB__dsec_forces.tocoo()
        J['loadsB','sec_forces'] = dloadsB__dsec_forces.data
        
        # --------------------------------dloadsB__ddef_mesh-------------------------------------
        # dmoment__ddiff must to be calulated each time b/c it's a function of sec_forces
        rows = np.zeros((9*(ny-1)*(nx-1)))
        cols = np.zeros((9*(ny-1)*(nx-1)))
        data = np.zeros((9*(ny-1)*(nx-1)))
        tmp_row = np.array([0,0,0,1,1,1,2,2,2])
        tmp_col = np.array([0,1,2,0,1,2,0,1,2])
        # setup rows and cols
        for nxitr in range(nx-1):
            for nyitr in range(ny-1):
                index = 9*(nyitr+(ny-1)*nxitr)     
                rows[index:index+9] = tmp_row+nyitr*3.    
                cols[index:index+9] = tmp_col+nyitr*3.+3*(ny-1)*nxitr 
                data[index:index+9] = -_skew(sec_forces[nxitr,nyitr,:]).flatten()
        dmoment__ddiff = coo_matrix((data, (rows,cols)), shape=(3*(ny-1), 3*(ny-1)*(nx-1)))
    
    
        # multiply the matrix with those previous calculated in setup
        dloadsB__ddef_mesh = self.dloadsB__dmoment * dmoment__ddiff * self.ddiff__ddef_mesh
        dloadsB__ddef_mesh = dloadsB__ddef_mesh.tocoo()
        J['loadsB','def_mesh'] = dloadsB__ddef_mesh.data 
        
        dloads__ddef_mesh = self.dloads__dloadsB * dloadsB__ddef_mesh
        dloads__ddef_mesh = dloads__ddef_mesh.tocoo()
        J['loads','def_mesh'] = dloads__ddef_mesh.data
        
        dloads__dsec_forces = self.dloads__dloadsA__dsec_forces + self.dloads__dloadsB * dloadsB__dsec_forces
        dloads__dsec_forces = dloads__dsec_forces.tocoo()
        J['loads','sec_forces'] = dloads__dsec_forces.data
