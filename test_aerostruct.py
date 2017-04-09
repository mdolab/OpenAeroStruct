# Main python script to test OpenAeroStruct coupled system components

from __future__ import print_function
# import coupled
import numpy
import aerostruct
# from cytn import aerostruct_cython
# import warnings
import sys
import time
import cProfile

# to disable openmdao warnings which will create an error in Matlab
numpy.set_printoptions(precision=8)

def test_timing(num_inboard=3, num_outboard=4, n=100):
    print('n=',n)
    print('Run coupled.setup()...')
    tic = time.clock()
    for i in xrange(n):
        def_mesh_coupled, params_coupled = coupled.setup(num_inboard, num_outboard)
    toc11 = time.clock() - tic
    print('Time per eval: {0} s'.format(toc11 / n))
    # toc11 = 0.0081286546 * n

    print('Run aerostruct.setup()...')
    tic = time.clock()
    for i in xrange(n):
        def_mesh_aerostruct, params_aerostruct = aerostruct.setup(num_inboard, num_outboard)
    toc12 = time.clock() - tic
    print('Time per eval: {0} s  --> {1}x faster than coupled'.format(toc12 / n, toc11/toc12))
    # toc12 = 0.00325556601613 * n

    print('Run aerostruct_cython.setup()...')
    tic = time.clock()
    for i in xrange(n):
        def_mesh_cython, params_cython = aerostruct_cython.setup(num_inboard, num_outboard)
    toc13 = time.clock() - tic
    print('Time per eval: {0} s  --> {1}x faster than coupled'.format(toc13 / n, toc11/toc13))
    print('                      --> {0}x faster than aerostruct'.format(toc12/toc13))
    #

    print('\nRun coupled.aero()...')
    tic = time.clock()
    for i in xrange(n):
        loads_coupled = coupled.aero(def_mesh_coupled, params_coupled)
    toc21 = time.clock() - tic
    print('Time per eval: {0} s'.format(toc21 / n))
    print('Run aerostruct.aerodynamics()...')
    tic = time.clock()
    for i in xrange(n):
        loads_aerostruct = aerostruct.aerodynamics(def_mesh_aerostruct, params_aerostruct)
    toc22 = time.clock() - tic
    print('Time per eval: {0} s  --> {1}x faster than coupled'.format(toc22 / n, toc21/toc22))
    print('Run aerostruct_cython.aerodynamics()...')
    tic = time.clock()
    for i in xrange(n):
        loads_cython = aerostruct_cython.aerodynamics(def_mesh_cython, params_cython)
    toc23 = time.clock() - tic
    print('Time per eval: {0} s  --> {1}x faster than coupled'.format(toc23 / n, toc21/toc23))
    print('                      --> {0}x faster than aerostruct'.format(toc22/toc23))

    #
    # print('\nRun coupled.struct()...')
    # tic = time.clock()
    # for i in xrange(n):
    #     def_mesh2 = coupled.struct(loads_coupled, params)
    # toc31 = time.clock() - tic
    # print('Time per eval: {0} s'.format(toc31 / n))
    # print('Run aerostruct.structures()...')
    # tic = time.clock()
    # for i in xrange(n):
    #     def_mesh2 = aerostruct.structures(loads_aerostruct, params)
    # toc32 = time.clock() - tic
    # print('Time per eval: {0} s  --> {1}x faster than coupled'.format(toc32 / n, toc31/toc32))
    # print('Run aerostruct_cython.structures()...')
    # tic = time.clock()
    # for i in xrange(n):
    #     def_mesh2 = aerostruct_cython.structures(loads_cython, params)
    # toc33 = time.clock() - tic
    # print('Time per eval: {0} s  --> {1}x faster than coupled'.format(toc33 / n, toc31/toc33))


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

    print('Run aerostruct loop ...')
    def_mesh, params = aerostruct.setup(num_inboard, num_outboard)
    tic = time.clock()
    for i in xrange(n):
        loads = aerostruct.aerodynamics(def_mesh, params)
        def_mesh = aerostruct.structures(loads, params)
    toc2 = time.clock() - tic
    print('Time per iteration: {0} s  --> {1}x faster'.format(toc2/n, toc1/toc2))


def test_accuracy(num_inboard=3, num_outboard=4):

    print('\nRun coupled.setup()...')
    def_mesh, params = coupled.setup(num_inboard, num_outboard)
    print('def_mesh...  def_mesh.shape =', def_mesh.shape)
    # print(def_mesh)

    print('\nRun aerostruct.setup()...')
    def_mesh, params = aerostruct.setup(num_inboard, num_outboard)
    print('def_mesh...  def_mesh.shape =', def_mesh.shape)
    # print(def_mesh)

    print('\nRun coupled.aero()...')
    loads = coupled.aero(def_mesh, params)
    print('loads matrix... loads.shape =', loads.shape)
    # print(loads)
    loads_coupled = loads

    print('\nRun aerostruct.aerodynamics()...')
    loads = aerostruct.aerodynamics(def_mesh, params)
    print('loads matrix... loads.shape =', loads.shape)
    # print(loads)
    loads_aerostruct = loads

    print('\nloads Difference:')
    # print(loads_coupled - loads_aerostruct)

    print('\nRun coupled.struct()...')
    def_mesh = coupled.struct(loads_coupled, params)
    print('def_mesh...  def_mesh.shape =',def_mesh.shape)
    # print(def_mesh)
    def_mesh_coupled = def_mesh

    print('\nRun aerostruct.structures()...')
    def_mesh = aerostruct.structures(loads_aerostruct, params)
    print('def_mesh...  def_mesh.shape =',def_mesh.shape)
    # print(def_mesh)
    def_mesh_aerostruct = def_mesh

    print('\ndef_mesh Difference:')
    # print(def_mesh_coupled - def_mesh_aerostruct)


def timings_aerodynamics(num_inboard=2, num_outboard=3,n=100):
    print('n=',n)
    print('\nRun aerostruct.setup()...')
    def_mesh_cython, params_cython = aerostruct_cython.setup(num_inboard, num_outboard)
    # print('def_mesh...  def_mesh.shape =', def_mesh.shape)
    # print(def_mesh)

    toc21 = 0.0156890999*n  # coupled.aero()
    toc22 = 0.0016663263704*n # aerostruct.aerodynamics()

    print('Run aerostruct_cython.aerodynamics()...')
    tic = time.clock()
    for i in xrange(n):
        loads_cython = aerostruct_cython.aerodynamics(def_mesh_cython, params_cython)
    toc23 = time.clock() - tic
    print('Time per eval: "{:8.8f}". s  --> "{:8.5}"x faster than coupled'.format(toc23 / n, toc21/toc23))
    print('                                --> "{:8.5f}"x faster than aerostruct'.format(toc22/toc23))

    # print('\nRun aerostruct.structures()...')
    # def_mesh = aerostruct.structures(loads_aerostruct, params)
    # print('def_mesh...  def_mesh.shape =',def_mesh.shape)
    # # print(def_mesh)
    # def_mesh_aerostruct = def_mesh



def main_aerostruct(num_x=2, num_y=7):
    print('\nRun aerostruct.setup()...')
    def_mesh, params = aerostruct.setup(num_x, num_y)
    # print('def_mesh...  def_mesh.shape =', def_mesh.shape)
    print(def_mesh)

    print('\nRun aerostruct.aerodynamics()...')
    loads = aerostruct.aerodynamics(def_mesh, params)
    print('loads matrix... loads.shape =', loads.shape)
    print(loads)
    loads_aerostruct = loads
    #
    print('\nRun aerostruct.structures()...')
    def_mesh = aerostruct.structures(loads_aerostruct, params)
    print('def_mesh...  def_mesh.shape =',def_mesh.shape)
    print(def_mesh)
    def_mesh_aerostruct = def_mesh


def profile_aerostruct(num_x=2, num_y=7, n=100):
    def_mesh, params = aerostruct.setup()
    for i in range(n):
        # if i % 10 == 0:
        #     print(i)
        loads = aerostruct.aerodynamics(def_mesh, params)
        loads_aerostruct = loads
        def_mesh = aerostruct.structures(loads_aerostruct, params)
        def_mesh_aerostruct = def_mesh


def main_coupled(num_inboard=2, num_outboard=3, check=False):

    print('\nRun coupled.setup()...')
    def_mesh, params = coupled.setup(num_inboard, num_outboard, check)
    print('def_mesh...  def_mesh.shape =', def_mesh.shape)
    print(def_mesh)

    print('\nRun coupled.aero()...')
    loads = coupled.aero(def_mesh, params)
    print('loads matrix... loads.shape =', loads.shape)
    print(loads)

    print('\nRun struct()...')
    def_mesh = coupled.struct(loads, params)
    print('def_mesh...  def_mesh.shape =', def_mesh.shape)
    print(def_mesh)


def main_cython(num_inboard=2, num_outboard=3):

    print('\nRun aerostruct_cython.setup()...')
    def_mesh, params = aerostruct_cython.setup(num_inboard, num_outboard)
    # print('def_mesh...  def_mesh.shape =', def_mesh.shape)
    # print(def_mesh)

    print('\nRun aerostruct.aerodynamics()...')
    loads = aerostruct_cython.aerodynamics(def_mesh, params)
    # print('loads matrix... loads.shape =', loads.shape)
    # print(loads)
    # loads__cython = loads
    #
    # print('\nRun aerostruct.structures()...')
    # def_mesh = aerostruct.structures(loads_aerostruct, params)
    # print('def_mesh...  def_mesh.shape =',def_mesh.shape)
    # # print(def_mesh)
    # def_mesh_aerostruct = def_mesh


if __name__ == '__main__':

    try:
        import lib
        fortran_flag = True
    except:
        fortran_flag = False
    print('Use Fortran: {0}'.format(fortran_flag))

    npts = [3, 5]
    n = 1000
    n_inboard = npts[0]
    n_outboard = npts[1]

    # main_coupled(n_inboard, n_outboard)
    cProfile.run('profile_aerostruct(n=n)','aerostruct.prof')
    # main_aerostruct()
    # main_cython(n_inboard, n_outboard)
    # test_accuracy(n_inboard, n_outboard)
    # test_timing(n_inboard, n_outboard, n)
    # timings_aerodynamics(n_inboard, n_outboard, n)
    # time_iterations(n_inboard, n_outboard, n)
