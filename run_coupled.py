# Main python script to test OpenAeroStruct coupled system components

from __future__ import print_function
import coupled
import aerostruct
import warnings
import sys

# to disable openmdao warnings which will create an error in Matlab
warnings.filterwarnings("ignore")


def main(num_inboard=2, num_outboard=3, check=False):

    print('\nRun coupled.setup()...')
    def_mesh, params = coupled.setup(num_inboard, num_outboard, check)
    print('def_mesh...  def_mesh.shape =', def_mesh.shape)
    print(def_mesh)

    print('\nRun aerostruct.setup()...')
    def_mesh, params = aerostruct.setup(num_inboard, num_outboard)
    print('def_mesh...  def_mesh.shape =', def_mesh.shape)
    print(def_mesh)

    print('\nRun coupled.aero()...')
    loads = coupled.aero(def_mesh, params)
    print('loads matrix... loads.shape =', loads.shape)
    print(loads)

    print('\nRun aerostruct.aero()...')
    loads = aerostruct.aero(def_mesh, params)
    print('loads matrix... loads.shape =', loads.shape)
    print(loads)

    # print('\nRun struct()...')
    # def_mesh = coupled.struct(loads, params)
    # print('def_mesh...  def_mesh.shape =',def_mesh.shape)
    # print(def_mesh)

if __name__ == '__main__':
    main()
