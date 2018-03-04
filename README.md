[![Build Status](https://travis-ci.org/mdolab/OpenAeroStruct.svg?branch=master)](https://travis-ci.org/mdolab/OpenAeroStruct)
[![Coverage Status](https://coveralls.io/repos/github/mdolab/OpenAeroStruct/badge.svg?branch=master)](https://coveralls.io/github/mdolab/OpenAeroStruct?branch=master)
[![Documentation Status](https://readthedocs.org/projects/openaerostruct/badge/?version=latest)](http://openaerostruct.readthedocs.io/en/latest/?badge=latest)
# OpenAeroStruct

OpenAeroStruct is a lightweight Python tool that performs aerostructural optimization of lifting surfaces using OpenMDAO. It uses a vortex lattice method (VLM) for the aerodynamics analysis and a spatial beam model with 6-DOF per element for the structural analysis.

Documentation is available [here](http://openaerostruct.readthedocs.io/en/latest/).

Please see the [SMO journal paper](https://link.springer.com/article/10.1007%2Fs00158-018-1912-8) for more information and please cite this article if you use OpenAeroStruct in your research. Here's an open-access copy of the paper: http://rdcu.be/Gtl1

![Optimized CRM-type wing with 30 panels](/example.png?raw=true "Example Optimization Result and Visualization")

## Installation

To use OpenAeroStruct, you must first install OpenMDAO 1.7.3 by following the instructions here: https://github.com/openmdao/openmdao. If you are unfamiliar with OpenMDAO and wish to modify the internals of OpenAeroStruct, you should examine the OpenMDAO documentation at http://openmdao.readthedocs.io/en/1.7.3/. The tutorials provided with OpenMDAO, especially the Sellar Problem, are helpful to understand the basics of using OpenMDAO to solve an optimization problem. Note that OpenMDAO 1.7.3 is the most recent version that has been tested and confirmed working with OpenAeroStruct.

Next, clone this repository:

    git clone https://github.com/mdolab/OpenAeroStruct.git

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
    python plot_all.py aerostruct.db

The first command performs aerostructural optimization and the second visualizes the optimization history.

The keywords used for each file are explained in their respective docstrings at the top of the file.

If you wish to examine the code in more depth, see `run_classes.py` and the methods it calls. These methods interface directly with OpenMDAO.

## Notes

This current version of this repository has grown past the previous Matlab implementation. If you are looking for a Matlab-capable version, please see https://github.com/samtx/OpenAeroStruct for the latest version.

## Known Issues

* The increase in accuracy of results when using a cosine-spaced mesh is not as large as it should be.
* Aerostructural optimization sometimes fails to converge for certain geometries. The example provided in `run_aerostruct.py` should converge. The structural and aerodynamic values must make sense together, e.g. the beam thickness and radius must be able to support the aerodynamic loads.
