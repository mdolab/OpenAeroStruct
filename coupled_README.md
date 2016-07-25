# OpenAeroStruct

OpenAeroStruct is a lightweight Python tool to perform aerostructural optimization of lifting surfaces using OpenMDAO.
It uses a vortex lattice method (VLM) expanded from Phillips' Modern Adaption of Prandtl's Classic Lifting Line Theory (http://arc.aiaa.org/doi/pdfplus/10.2514/2.2649) for the aerodynamics analysis and a spatial beam model with 6-DOF per element for the structural analysis.
The `coupled.py` module isolates the aerodynamic and structure analysis into a coupled system comprising two "black boxes". 
The coupled system modules can be called from Matlab using the `.m` files in the repository. 

![Optimized CRM-type wing with 30 panels](/example.png?raw=true "Example Optimization Result and Visualization")

## Installation and Configuration

To use OpenAeroStruct, you must first install the following software and dependencies:
- Python >=2.7.9 or >=3.4.3 
- Numpy >=1.9.2 
- Scipy >=0.15.1 
- OpenMDAO >=1.7.0 
- Matlab >=2014b for using the Matlab wrapper functions 

Python, Numpy, and Scipy can be easily installed together using Anaconda, which can be downloaded here: https://www.continuum.io/downloads

To use OpenAeroStruct, you must first install OpenMDAO 1.7.0 by following the instructions here: https://github.com/openmdao/openmdao. If you are unfamiliar with OpenMDAO and wish to modify the internals of OpenAeroStruct, you should examine the OpenMDAO documentation at http://openmdao.readthedocs.io/en/1.7.0/. The tutorials provided with OpenMDAO, especially The Sellar Problem, are helpful to understand the basics of using OpenMDAO to solve an optimization problem.

By default, the Python package manager `pip` comes installed with Anaconda. OpenMDAO can be easily installed by opening the Anaconda prompt and entering

    pip install openmdao 

Next, clone the OpenAeroStruct repository from GitHub to get the required files:

    git clone https://github.com/johnjasa/OpenAeroStruct.git

Lastly, from within the OpenAeroStruct folder, make the Fortran files:

    make

Note that the code will run without compiling the Fortran library, but it will run significantly faster when using Fortran.

## Usage

`run_vlm.py` is for aero-only analysis and optimization. It can use a single lifting surface or multiple separate lifting surfaces.

`run_spatialbeam.py` is for structural-only analysis and optimization. It can use a single structural component or multiple structural components, where each component represents a spar within a lifting surface.

`run_aerostruct.py` performs aerostructural analysis and optimization.


For each case, you can view the optimization results using `plot_all.py`. Examine its docstring for keyword information.

An example workflow would be:

    python run_aerostruct.py 1
    python plot_all.py as

The keywords used for each file are explained in their respective docstrings at the top of the file.

## Calling OpenAeroStruct coupled system modules from Matlab

### Configuration

Once OpenMDAO and its dependencies have been installed, you can run the aerodynamics and structure modules from Matlab using wrapper functions. You need Matlab version 2014b or greater and it must be of the same architecture (either 32-bit or 64-bit) as Python installed on your system. Run `pyversion` from the Matlab console to confirm if your Matlab/Python is configured correctly. It should automatically detect the Python executable file. [Link to Matlab documentation on `pyversion`](http://www.mathworks.com/help/matlab/ref/pyversion.html)

```
>> pyversion

       version: '2.7'
    executable: 'C:\Users\samfriedman\Anaconda2_64\python.EXE'
       library: 'C:\Users\samfriedman\Anaconda2_64\python27.dll'
          home: 'C:\Users\samfriedman\Anaconda2_64'
      isloaded: 0
```

If the Python executable isn't specified, or if you need to use a non-default Python version, then call `pyversion` followed by the full path to the executable file.

    >> pyversion 'C:\Python33\python.exe'

or

    >> pyversion '/usr/bin/python'

Matlab loads the Python interpreter when a valid Python command is entered. This action sets the `pyversion` output variable `isloaded` to 1. The path to the Python executable can only be changed when Python isn't loaded in Matlab. To change the path after Python is loaded you must restart Matlab.

Here is an [example Python command](http://www.mathworks.com/help/matlab/matlab_external/call-user-defined-custom-module.html) to use for testing:

```
>> N = py.list({'Jones','Johnson','James'})

N = 

  Python list with no properties.

    ['Jones', 'Johnson', 'James']
```

### Using the Matlab wrappers for the coupled module

Refer to the Matlab script `run_coupled.m` for an example on solving the aerostructural coupled system in Matlab.

To write your own script or function, you must first call `coupled_setup.m` to create the initial wing mesh and the Python dict object containing system parameters needed for the aerodynamics and structures modules. You must specify the number of spanwise inboard and outboard points on the airplane wing as integer arguments to the function. The points will be mirrored to the other side of the wing to produce a full mesh of N = 2(n_inb+n_outb)-3 spanwise points. Each spanwise point has two chordwise points, one of the leading edge and the other on the trailing edge of the wing. The initial mesh sets all points to an elevation of 0.    

To produce the array of loads on the wing, pass the wing mesh and the params dict to the `coupled_aero.m` function.

```
loads = coupled_aero(mesh, params);
```

The `loads` array is an (N, 6) matrix of the force and moment loads on the wing applied at the control points. The control points are along each spanwise pair of points at the 0.35*chord. Columns 1-3 are the x, y, and z components of the force vectors, while columns 4-6 are the x, y, and z components of the moment vectors.

To find the new mesh based upon the aerodynamic loads, call the `coupled_struct.m` function.

```
mesh = coupled_struct(loads, params);
```

This function returns the wing mesh after calculating the displacements caused by the aerodynamic loads.

Other Matlab files of interest:
- `draw_crm.m` makes a plot of one side of the CRM base wing.
- `mat2np.m` converts a Matlab array to a Numpy ndarry.
- `np2mat.m` converts a Numpy ndarry to a Matlab array.
- `coupled_plotdata.m` plots the wings mesh, force vectors, and moment vectors in 3D. This is still a work in progress.

For a full example of solving the aerostructural coupled system by iterating until both arrays convergence, see the `run_coupled.m` Matlab script.

## Known Issues

* Aerostructural optimization sometimes fails to converge for certain geometries. The example provided in `run_aerostruct.py` should converge.
* Aerostructural optimization using multiple lifting surfaces does not converge.
* The residual of the structural system solution for very large problems is sometimes too large and prevents convergence of the optimization problem.
* Internal documentation is lacking.
* `plot_all.py` does not correctly display multiple structural components, but does work for multiple surfaces if using only aerodynamic optimization.
* Multiple surface optimization for structures is not optimally coded.
