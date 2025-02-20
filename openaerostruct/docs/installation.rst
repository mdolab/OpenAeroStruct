.. _Installation:

Installation
============

To use OpenAeroStruct, you must first install Python 3.

The easiest way to get started is to install OpenAeroStruct from PyPI :

.. code-block:: bash

    pip install openaerostruct

If you'd like easier access to the examples and source code, you can also install OpenAeroStruct by cloning the OpenAeroStruct repository:

.. code-block:: bash

    git clone https://github.com/mdolab/OpenAeroStruct.git

Then from within the OpenAeroStruct folder, pip install the package:

.. code-block:: bash

	cd openaerostruct
	pip install -e .

Both methods will automatically install the dependencies: numpy, scipy, matplotlib, and OpenMDAO.

The oldest and latest versions of the dependencies that we test regularly are the following (other versions may work, but no guarantees):

.. list-table::
    :header-rows: 1

    * - Dependency
      - oldest
      - latest
    * - Python
      - 3.9
      - 3.11
    * - NumPy
      - 1.21
      - 1.26
    * - SciPy
      - 1.7
      - 1.15
    * - OpenMDAO
      - 3.35
      - latest
    * - pyGeo (optional)
      - 1.15.0
      - latest
    * - OpenVSP (optional)
      - 3.33
      - 3.33
    * - MPhys (optional)
      - 2.0
      - latest

Numpy 2.0 or later should also work, but we currently do not run tests with Numpy 2.

If you are unfamiliar with OpenMDAO and wish to modify the internals of OpenAeroStruct, you should examine the OpenMDAO documentation at https://openmdao.org/newdocs/versions/latest/main.html. The tutorials provided with OpenMDAO are helpful to understand the basics of using OpenMDAO to solve an optimization problem.

Advanced Options
~~~~~~~~~~~~~~~~

To run the tests on your machine, use the [test] option. This will install the pytest package.

.. code-block:: bash

    pip install -e .[test]

Then run the tests from the OpenAeroStruct root directory by calling:

.. code-block:: bash

    testflo -v .

To install the dependencies to build the documentation locally, run:

.. code-block:: bash

    pip install -e .[docs]
