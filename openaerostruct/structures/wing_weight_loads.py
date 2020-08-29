import numpy as np

from scipy.sparse import coo_matrix, diags

import openmdao.api as om
from openaerostruct.structures.utils import norm
from openaerostruct.utils.constants import grav_constant


class StructureWeightLoads(om.ExplicitComponent):
    """
    Compute the nodal loads from the weight of the wing structure to be applied to the wing
    structure.

    Parameters
    ----------
    element_mass[ny-1] : numpy array
        Weight for each wing-structure segment.
    nodes[ny, 3] : numpy array
        Flattened array with coordinates for each FEM node.
    load_factor : float
        Load factor for the flight point.

    Returns
    -------
    struct_weight_loads[ny, 6] : numpy array
        Flattened array containing the loads applied on the FEM component,
        computed from the weight of the wing-structure segments.
    element_lengths[ny-1] : numpy array
        Lengths of the FEM finite elements.
    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        self.surface = surface = self.options['surface']
        self.ny = surface['mesh'].shape[1]

        self.add_input('element_mass', val=np.zeros((self.ny-1)), units='kg')
        self.add_input('nodes', val=np.zeros((self.ny, 3), dtype=complex), units='m')
        self.add_input('load_factor', val=1.0)

        self.add_output('struct_weight_loads', val=np.zeros((self.ny, 6)), units='N')
        self.add_output('element_lengths', val=np.zeros(self.ny-1), units='m')

        nym1 = self.ny-1

        #dloads__dzf
        rows = np.zeros(2*nym1)
        cols = np.zeros(2*nym1)
        data = -np.ones(2*nym1)
        counter = 0
        for i in range(nym1):
            for j in range(2):
                rows[counter] = 6*(i+j) + 2
                cols[counter] = i
                counter += 1

        self.dswl__dzf = coo_matrix((data, (rows, cols)), shape=(self.ny*6,nym1))


        rows = np.zeros(3*self.ny)
        cols = np.zeros(3*self.ny)
        for j in range(3):
            for i in range(self.ny):
                idx = i+self.ny*j
                rows[idx] = 6*i+j+2

        self.dswl__dlf_row = rows
        self.dswl__dlf_col = cols

        self.declare_partials('struct_weight_loads', 'load_factor', rows=rows, cols=cols)
        # self.declare_partials('struct_weight_loads', 'load_factor')


        # dstruct_weight_loads__dnodes (this one is super complicated)
        #     we can pre-compute several of the intermediate partial derivative matrices
        #     others we can pre-compute the rows/cols so we just have to set the data
        #     the for loops for rows/cols are a bit expensive (cause python), but the data computations are fast
        #     so we will pay the one-time setup cost to pre-compute the rows/cols
        rows_del = np.zeros(2*nym1)
        cols_d0 = np.zeros(2*nym1)
        cols_d1 = np.zeros(2*nym1)
        data = np.zeros(2*nym1)

        counter = 0
        for i in range(nym1):
            for j in range(2):
                rows_del[counter] = i
                cols_d0[counter] = 3*(i+j)
                cols_d1[counter] = 3*(i+j)+1
                data[counter] = 2*j-1 #this gives 1 or -1 alternating based on j
                counter += 1

        shape = (nym1,3*self.ny)
        self.ddel0__dnodes = coo_matrix((data, (rows_del, cols_d0)), shape=shape)
        self.ddel1__dnodes = coo_matrix((data, (rows_del, cols_d1)), shape=shape)

        # del__dnodes
        rows_el = np.zeros(6*nym1)
        cols_el = np.zeros(6*nym1)
        data = np.ones(6*nym1) # del__dnodes matrix has only ones in it

        counter = 0
        for i in range(nym1):
            for j in range(6):
                rows_el[counter] = i
                cols_el[counter] = 3*i+j
                counter += 1
        self.del__dnodes = coo_matrix((data, (rows_el, cols_el)), shape=shape)

        nym1_rand = np.random.rand(nym1)

        dzm_dnodes_pattern = ((diags(nym1_rand*nym1_rand)*self.ddel0__dnodes+diags(nym1_rand*nym1_rand)*self.ddel1__dnodes)).tocsr()


        dbm3_dnodes = diags(nym1_rand)*dzm_dnodes_pattern \
                    + diags(nym1_rand)*self.ddel1__dnodes \
                    - diags(nym1_rand)*self.del__dnodes


        dbm4_dnodes = diags(nym1_rand)*dzm_dnodes_pattern \
                    + diags(nym1_rand)*self.ddel0__dnodes \
                    - diags(nym1_rand)*self.del__dnodes

        self.dbm3_dnodes_pattern = (diags(nym1_rand)*dzm_dnodes_pattern \
                                  + diags(nym1_rand)*self.ddel1__dnodes\
                                  - diags(nym1_rand)*self.del__dnodes)
        self.dbm4_dnodes_pattern = (diags(nym1_rand)*dzm_dnodes_pattern \
                                  + diags(nym1_rand)*self.ddel0__dnodes\
                                  - diags(nym1_rand)*self.del__dnodes)
        #dswl__dbm
        rows_1 = np.zeros(2*nym1)
        rows_0 = np.zeros(2*nym1)
        cols = np.zeros(2*nym1)
        data = np.random.rand(2*nym1)
        counter = 0
        for i in range(nym1):
            for j in range(2):
                rows_1[counter] = 6*(i+j)+3
                rows_0[counter] = 6*(i+j)+4
                if counter%2==0: # even
                    data[counter] = 1
                else: # odd
                    data[counter] = -1
                cols[counter] = i
                counter += 1

        self.dswl__dbm3 = coo_matrix((data, (rows_1, cols)), shape=(6*self.ny,nym1)).tocsr()
        self.dswl__dbm4 = coo_matrix((data, (rows_0, cols)), shape=(6*self.ny,nym1)).tocsr()

        self.dswl__dnodes_pattern = (self.dswl__dbm3*dbm3_dnodes+self.dswl__dbm4*dbm4_dnodes).tocoo()

        self.dswl__dew_pattern = (self.dswl__dbm4 + self.dswl__dbm3 + self.dswl__dzf).tocoo()

        self.declare_partials('struct_weight_loads', 'nodes',
                              rows=self.dswl__dnodes_pattern.row,
                              cols=self.dswl__dnodes_pattern.col)

        self.declare_partials('struct_weight_loads', 'element_mass',
                              rows=self.dswl__dew_pattern.row,
                              cols=self.dswl__dew_pattern.col)
        self.set_check_partial_options(wrt='*', method='cs')

    def compute(self, inputs, outputs):

        struct_weights = inputs['element_mass'] * inputs['load_factor'] * grav_constant
        nodes = inputs['nodes']

        element_lengths = norm(nodes[1:, :] - nodes[:-1, :], axis=1)

        # And we also need the deltas between consecutive nodes
        deltas = nodes[1:, :] - nodes[:-1, :]
        # save these slices cause I use them a lot
        del0 = deltas[: , 0]
        del1 = deltas[: , 1]

        # Assume weight coincides with the elastic axis
        z_forces_for_each = struct_weights / 2.
        z_moments_for_each = struct_weights / 12. \
                            * (del0**2 + del1**2)**0.5

        loads = outputs['struct_weight_loads']
        loads *= 0 # need to zero it out, since we're accumulating onto it
        # Why doesn't this trigger when running ../../tests/test_aerostruct_wingbox_+weight_analysis.py???
        if self.under_complex_step:
            loads = np.zeros((self.ny, 6), dtype=complex)

        # Loads in z-direction
        loads[:-1, 2] += -z_forces_for_each
        loads[1:, 2] += -z_forces_for_each

        # Bending moments for consistency
        bm3 = z_moments_for_each * del1 / element_lengths
        loads[:-1, 3] += -bm3
        loads[1:, 3] += bm3

        bm4 = z_moments_for_each * del0 / element_lengths
        loads[:-1, 4] += -bm4
        loads[1:, 4] += bm4
        outputs['struct_weight_loads'] = loads

    def compute_partials(self, inputs, J):

        struct_weights = inputs['element_mass'] * inputs['load_factor'] * grav_constant
        nodes = inputs['nodes']

        element_lengths = norm(nodes[1:, :] - nodes[:-1, :], axis=1)

        # And we also need the deltas between consecutive nodes
        deltas = nodes[1:, :] - nodes[:-1, :]
        # save these slices cause I use them alot
        del0 = deltas[: , 0]
        del1 = deltas[: , 1]

        # Assume weight coincides with the elastic axis
        z_moments_for_each = struct_weights / 12. \
                            * (del0**2 + del1**2)**0.5

        dzf__dew = .5*inputs['load_factor'][0]
        dzf__dlf = inputs['element_mass']/2.

        dzm__dew = diags((del0**2 + del1**2)**0.5/12*inputs['load_factor'])
        dzm__dlf =  (del0**2 + del1**2)**.5/12. * inputs['element_mass']

        dbm3__dzm = diags(del1/element_lengths)
        dbm4__dzm = diags(del0/element_lengths)

        # need to convert to lil to re-order the data to match the original row/col indexing from setup
        dswl__dlf = (-self.dswl__dbm3*dbm3__dzm*coo_matrix(dzm__dlf).T +\
                    -self.dswl__dbm4*dbm4__dzm*coo_matrix(dzm__dlf).T +\
                    self.dswl__dzf*coo_matrix(dzf__dlf).T).tolil()

        J['struct_weight_loads', 'load_factor'] = dswl__dlf[self.dswl__dlf_row, self.dswl__dlf_col].toarray().flatten() * grav_constant

        dswl__dew = (-self.dswl__dbm4 * dbm4__dzm * dzm__dew +\
                    -self.dswl__dbm3 * dbm3__dzm * dzm__dew +\
                    self.dswl__dzf  * dzf__dew).tolil()

        data = dswl__dew[self.dswl__dew_pattern.row, self.dswl__dew_pattern.col].toarray().flatten()
        J['struct_weight_loads', 'element_mass'] = data * grav_constant


        # dstruct_weight_loads__dnodes (this one is super complicated)

        # del__dnodes
        # note: del__dnodes matrix already created in setup, just need to set data
        raw_data = diags(1/element_lengths)*(nodes[1:, :] - nodes[:-1, :])
        data = np.hstack((-raw_data,raw_data)).flatten()
        self.del__dnodes.data = data

        dzm_dnodes = diags(struct_weights/12*(del0**2 + del1**2)**-.5)*\
                     (diags(del0)*self.ddel0__dnodes  + diags(del1)*self.ddel1__dnodes)

        dbm3_dnodes = diags(del1/element_lengths)*dzm_dnodes \
                    + diags(z_moments_for_each/element_lengths)*self.ddel1__dnodes \
                    - diags(z_moments_for_each*del1/element_lengths**2)*self.del__dnodes


        dbm4_dnodes = diags(del0/element_lengths)*dzm_dnodes \
                    + diags(z_moments_for_each/element_lengths)*self.ddel0__dnodes \
                    - diags(z_moments_for_each*del0/element_lengths**2)*self.del__dnodes


        # this is kind of dumb, but I need lil cause I have to re-index to preserve order
        #     the coo column ordering doesn't seem to be deterministic
        #     so I use the original row/col from the pattern as index arrays to
        #     pull the data out in the correct order
        dswl__dnodes= (self.dswl__dbm3*dbm3_dnodes+self.dswl__dbm4*dbm4_dnodes).tolil()
        data = dswl__dnodes[self.dswl__dnodes_pattern.row, self.dswl__dnodes_pattern.col].toarray().flatten()
        J['struct_weight_loads', 'nodes'] = -data
