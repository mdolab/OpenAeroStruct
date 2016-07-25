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

## Calling OpenAeroStruct coupled system modules from Matlab

### Configuration

Once OpenMDAO and its dependencies have been installed, you can run the aerodynamics and structure modules from Matlab. You need Matlab version 2014b or greater and it must be of the same architecture (either 32-bit or 64-bit) as Python installed on your system. Run `pyversion` from the Matlab console to confirm if your Matlab/Python is configured correctly. It should automatically detect the Python executable file.

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

There are three wrapper functions and two utility functions that are required to call the coupled module from Matlab. Matlab Here is a working version of the isolated coupled functions, `coupled.setup()`, `coupled.aero()`, and `coupled.struct()`.
There are also Matlab files to call these files from Matlab.  

There are two scripts that can be used to test these new functions:
- `run_coupled.py` : this sets up the coupled system problem and performs one iteration.
It creates the `def_mesh` and a dict `param` containing system parameters. Then it gets the `loads` array from the `coupled.aero(def_mesh, params)` function. Then it gets the new mesh from `def_mesh = coupled.struct(loads, params)` function. 
- `run_coupled.m` : This does the same thing as the Python version, except that it iterates through the coupled loop until the two arrays converge. 

Other Matlab files of interest:
- `draw_crm.m` makes a plot of one side of the CRM base wing.
- `mat2np.m` converts a Matlab array to a Numpy ndarry.
- `np2mat.m` converts a Numpy ndarry to a Matlab array.
- `coupled_plotdata.m` plots the wings mesh, force vectors, and moment vectors in 3D. This one is still a work in progress.


## Usage

`run_vlm.py` is for aero-only analysis and optimization. It can use a single lifting surface or multiple separate lifting surfaces.

`run_spatialbeam.py` is for structural-only analysis and optimization. It can use a single structural component or multiple structural components, where each component represents a spar within a lifting surface.

`run_aerostruct.py` performs aerostructural analysis and optimization.


For each case, you can view the optimization results using `plot_all.py`. Examine its docstring for keyword information.

An example workflow would be:

    python run_aerostruct.py 1
    python plot_all.py as

The keywords used for each file are explained in their respective docstrings at the top of the file.

## Known Issues

* Aerostructural optimization sometimes fails to converge for certain geometries. The example provided in `run_aerostruct.py` should converge.
* Aerostructural optimization using multiple lifting surfaces does not converge.
* The residual of the structural system solution for very large problems is sometimes too large and prevents convergence of the optimization problem.
* Internal documentation is lacking.
* `plot_all.py` does not correctly display multiple structural components, but does work for multiple surfaces if using only aerodynamic optimization.
* Multiple surface optimization for structures is not optimally coded.
