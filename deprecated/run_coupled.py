# Main python script to test OpenAeroStruct coupled system components

from __future__ import print_function
import coupled
import numpy
import warnings
import time

# to disable openmdao warnings which will create an error in Matlab
warnings.filterwarnings('ignore')
numpy.set_printoptions(precision=8)

def test_timing(num_inboard=3, num_outboard=4, n=100):
    print('n=',n)
    print('Run coupled.setup()...')
    tic = time.clock()
    for i in xrange(n):
        def_mesh, params = coupled.setup(num_inboard, num_outboard)
    toc11 = time.clock() - tic
    print('Time per eval: {0} s'.format(toc11 / n))

    print('\nRun coupled.aero()...')
    tic = time.clock()
    for i in xrange(n):
        loads = coupled.aero(def_mesh, params)
    toc21 = time.clock() - tic
    print('Time per eval: {0} s'.format(toc21 / n))

    print('\nRun coupled.struct()...')
    tic = time.clock()
    for i in xrange(n):
        def_mesh = coupled.struct(loads, params)
    toc31 = time.clock() - tic
    print('Time per eval: {0} s'.format(toc31 / n))


def time_iterations(num_inboard=3, num_outboard=4, n=100):
    print('\n...Time iterations... ')
    print('Run coupled loop ...')
    def_mesh, params = coupled.setup(num_inboard, num_outboard)
    tic = time.clock()
    for i in xrange(n):
        loads = coupled.aero(def_mesh, params)
        def_mesh = coupled.struct(loads, params)
    toc1 = time.clock() - tic
    print('Time per iteration: {0} s  '.format(toc1/n))


def main_coupled(num_inboard=2, num_outboard=3, check=False):

    print('\nRun coupled.setup()...')
    def_mesh, params = coupled.setup(num_inboard, num_outboard, check)
    print('def_mesh...  def_mesh.shape =', def_mesh.shape)
    print(def_mesh)
    def_mesh_aerostruct = def_mesh


    print('\nRun coupled.aero()...')
    loads = coupled.aero(def_mesh, params)
    print('loads matrix... loads.shape =', loads.shape)
    print(loads)

    print('\nRun struct()...')
    def_mesh = coupled.struct(loads, params)
    print('def_mesh...  def_mesh.shape =', def_mesh.shape)
    print(def_mesh)


if __name__ == '__main__':

    try:
        import lib
        fortran_flag = True
    except:
        fortran_flag = False
    print('Use Fortran: {0}'.format(fortran_flag))

    npts = [3, 5]
    n = 100
    n_inboard = npts[0]
    n_outboard = npts[1]

    main_coupled(n_inboard, n_outboard)
    # test_timing(n_inboard, n_outboard, n)
    # time_iterations(n_inboard, n_outboard, n)
