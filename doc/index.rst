.. OpenAeroStruct documentation master file

OpenAeroStruct Documentation
==========================================

OpenAeroStruct is a lightweight tool that performs aerostructural optimization using OpenMDAO.
It couples a vortex-lattice method (VLM) and a 6 degrees of freedom 3-dimensional spatial beam model to simulate aerodynamic and structural analyses using lifting surfaces.
These simulations are wrapped with an optimizer using NASA's OpenMDAO framework.
The analysis and optimization results can be visualized using included tools, producing figures such as this:

.. image:: ../example.png

Installation
-----------------

To use OpenAeroStruct, you must first install OpenMDAO 1.7.3 by following the instructions here: https://github.com/openmdao/openmdao1. If you are unfamiliar with OpenMDAO and wish to modify the internals of OpenAeroStruct, you should examine the OpenMDAO documentation at http://openmdao.readthedocs.io/en/1.7.3/. The tutorials provided with OpenMDAO, especially the Sellar Problem, are helpful to understand the basics of using OpenMDAO to solve an optimization problem. Note that OpenMDAO 1.7.3 is the most recent version that has been tested and confirmed working with OpenAeroStruct.

Next, clone this repository:

.. code-block:: bash

    git clone https://github.com/mdolab/OpenAeroStruct.git

If you are using MacOS, move the `config-macOS.mk` file to the `config/` folder and rename it to `config.mk` to replace the existing `config.mk` file. The default `config.mk` file works for Ubuntu users. Unfortunately, there is no documented Windows support for the Fortran code.

If you're on Mac or Ubuntu, run this command from within the OpenAeroStruct folder to make the Fortran files:

.. code-block:: bash

    make

Note that the code will run without compiling the Fortran library, but it will run significantly faster when using Fortran.

We include a script that runs a variety of aerodynamic, structural, and aerostructural tests. To check your installation, run:

.. code-block:: bash

    python test_suite.py

Usage
-----

`run_vlm.py` is for aero-only analysis and optimization. It can use a single lifting surface or multiple separate lifting surfaces.

`run_spatialbeam.py` is for structural-only analysis and optimization. It can use a single structural component or multiple structural components, where each component represents a spar within a lifting surface.

`run_aerostruct.py` performs aerostructural analysis and optimization.


For each case, you can view the optimization results using `plot_all.py`. Examine its docstring for keyword information.

An example workflow would be:

.. code-block:: bash

    python run_aerostruct.py 1
    python plot_all.py aerostruct.db

The first command performs aerostructural optimization and the second visualizes the optimization history.

The keywords used for each file are explained in their respective docstrings at the top of the file.

If you wish to examine the code in more depth, see `run_classes.py` and the methods it calls. These methods interface directly with OpenMDAO.

Notes
-----

A tutorial on how to implement OpenAeroStruct in :ref:`Matlab <Matlab>` is available in the
Tutorials menu.

Known Issues
------------

* The increase in accuracy of results when using a cosine-spaced mesh is not as large as it should be.
* Aerostructural optimization sometimes fails to converge for certain geometries. The example provided in `run_aerostruct.py` should converge. The structural and aerodynamic values must make sense together, e.g. the beam thickness and radius must be able to support the aerodynamic loads.



Tutorials and Indices
==================

.. toctree::
   :maxdepth: 2

   tutorials.rst

Check out the module index below to see the internal methods within each file and how they're used.

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
