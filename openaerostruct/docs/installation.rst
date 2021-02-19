.. _Installation:

Installation
============

To use OpenAeroStruct, you must first install OpenMDAO (versions 3.2.0 to 3.6.0 work; 3.6.0 is recommended) by following the instructions at https://github.com/OpenMDAO/OpenMDAO/.
If you are unfamiliar with OpenMDAO and wish to modify the internals of OpenAeroStruct, you should examine the OpenMDAO documentation at http://openmdao.org/twodocs/versions/latest/index.html. The tutorials provided with OpenMDAO are helpful to understand the basics of using OpenMDAO to solve an optimization problem.

Next, clone the OpenAeroStruct repository:

.. code-block:: bash

    git clone https://github.com/mdolab/OpenAeroStruct.git

Then from within the OpenAeroStruct folder, pip install the package:

.. code-block:: bash

    pip install -e .
