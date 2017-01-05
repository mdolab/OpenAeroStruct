# Main python script to test OpenAeroStruct coupled system components

from __future__ import print_function
import pyximport; pyximport.install()
from OpenAeroStruct import coupled
import aerostruct
import aerostruct_cython
import numpy

import warnings
import sys
import time

# to disable openmdao warnings which will create an error in Matlab
warnings.filterwarnings("ignore")
numpy.set_printoptions(precision=8)


def test_timing(num_inboard=3, num_outboard=4, n=100):
    print('\nRun coupled.setup()...')
    tic = time.clock()
    for i in xrange(n):
        def_mesh, params = coupled.setup(num_inboard, num_outboard)
    toc1 = time.clock() - tic
    print('Time per eval: {0} s'.format(toc1 / n))

    print('\nRun aerostruct.setup()...')
    tic = time.clock()
    for i in xrange(n):
        def_mesh, params = aerostruct.setup(num_inboard, num_outboard)
    toc2 = time.clock() - tic
    print('Time per eval: {0} s  --> {1}x faster'.format(toc2 / n, toc1/toc2))

    print('\nRun coupled.aero()...')
    tic = time.clock()
    for i in xrange(n):
        loads_coupled = coupled.aero(def_mesh, params)
    toc3 = time.clock() - tic
    print('Time per eval: {0} s'.format(toc3 / n))

    print('\nRun aerostruct.aerodynamics()...')
    tic = time.clock()
    for i in xrange(n):
        loads_aerostruct = aerostruct.aerodynamics(def_mesh, params)
    toc4 = time.clock() - tic
    print('Time per eval: {0} s  --> {1}x faster'.format(toc4 / n, toc3/toc4))

    print('\nRun coupled.struct()...')
    tic = time.clock()
    for i in xrange(n):
        def_mesh2 = coupled.struct(loads_coupled, params)
    toc5 = time.clock() - tic
    print('Time per eval: {0} s'.format(toc5 / n))

    print('\nRun aerostruct.structures()...')
    tic = time.clock()
    for i in xrange(n):
        def_mesh2 = aerostruct.structures(loads_aerostruct, params)
    toc6 = time.clock() - tic
    print('Time per eval: {0} s  --> {1}x faster'.format(toc6 / n, toc5/toc6))


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
    n = 50  # number of times to run each function

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
    print(loads_coupled - loads_aerostruct)

    print('\nRun coupled.struct()...')
    def_mesh = coupled.struct(loads_coupled, params)
    print('def_mesh...  def_mesh.shape =',def_mesh.shape)
    print(def_mesh)
    def_mesh_coupled = def_mesh

    print('\nRun aerostruct.structures()...')
    def_mesh = aerostruct.structures(loads_aerostruct, params)
    print('def_mesh...  def_mesh.shape =',def_mesh.shape)
    print(def_mesh)
    def_mesh_aerostruct = def_mesh

    print('\ndef_mesh Difference:')
    print(def_mesh_coupled - def_mesh_aerostruct)


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

if __name__ == '__main__':

    try:
        from .. import lib
        fortran_flag = True
    except:
        fortran_flag = False
    print('Use Fortran: {0}'.format(fortran_flag))

    npts = [3, 5]
    n = 300
    n_inboard = npts[0]
    n_outboard = npts[1]

    # main_coupled(n_inboard, n_outboard)
    # test_accuracy(n_inboard, n_outboard)
    # test_timing(n_inboard, n_outboard, n)
    time_iterations(n_inboard, n_outboard, n)
