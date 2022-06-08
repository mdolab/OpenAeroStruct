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

The latest version of OpenAeroStruct supports OpenMDAO version 3.2.0 to 3.16.0; 3.16.0 is recommended.
If you wish to install a specific version of OpenMDAO, follow the instructions at https://github.com/OpenMDAO/OpenMDAO/.

If you are unfamiliar with OpenMDAO and wish to modify the internals of OpenAeroStruct, you should examine the OpenMDAO documentation at http://openmdao.org/twodocs/versions/latest/index.html. The tutorials provided with OpenMDAO are helpful to understand the basics of using OpenMDAO to solve an optimization problem.

Advanced Options
~~~~~~~~~~~~~~~~

To run the tests on your machine, use the [test] option. This will install the `testflo <https://github.com/OpenMDAO/testflo>`_ package.

.. code-block:: bash

    pip install -e .[test]


To install the dependencies to build the documentation locally, run:

.. code-block:: bash

    pip install -e .[docs]

The documentation build requires OpenMDAO 3.9.2 or older.