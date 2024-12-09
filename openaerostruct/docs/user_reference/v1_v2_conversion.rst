.. _V1_V2_Conversion:

OpenAeroStruct v1 to v2 Conversion
==================================

There are quite a few differences between OpenAeroStruct v1 and v2, though the underlying analyses remain largely the same.
In this document, we'll example code in each subsection using v1 on the top and the corresponding code to perform the same computations in v2 on the bottom.
This can serve as a Rosetta Stone to translate older v1 scripts to run in v2.

In general, much more of the code is exposed to the user in v2.
This was done to make the settings more explicit and to make it clear what values the user needs to supply.
For example, in v1 there was a default surface dictionary that contained some nominal values, but in v2 there are no defaults for any values in the surface dictionaries so the user must supply everything necessary.

We'll first go through individual sets of commands then present a full example script that performs the same analyses in both versions.

Instantiate the problem
-----------------------

v1 Script
^^^^^^^^^

.. code-block:: python

  from __future__ import division, print_function
  from time import time
  import numpy as np

  from OpenAeroStruct import OASProblem


  # Set problem type
  prob_dict = {'type' : 'aero',
                'optimize' : True}

  # Instantiate problem and add default surface
  OAS_prob = OASProblem(prob_dict)

v2 Script
^^^^^^^^^

.. literalinclude:: /../../tests/integration_tests/test_v1_aero_opt.py
  :start-after: checkpoint 0
  :end-before: checkpoint 1
  :dedent: 8

Create the surface and add it to the problem
--------------------------------------------

v1 Script
^^^^^^^^^

.. code-block:: python

  # Create a dictionary to store options about the surface
  surf_dict = {'num_y' : 5,
                'num_x' : 3,
                'wing_type' : 'rect',
                'symmetry' : True,
                'num_twist_cp' : 2}

  # Add the specified wing surface to the problem
  OAS_prob.add_surface(surf_dict)

v2 Script
^^^^^^^^^

.. literalinclude:: /../../tests/integration_tests/test_v1_aero_opt.py
  :start-after: checkpoint 2
  :end-before: checkpoint 3
  :dedent: 8


Set up the problem, add design variables, and run the optimization
------------------------------------------------------------------

v1 Script
^^^^^^^^^

.. code-block:: python

  # Setup problem and add design variables, constraint, and objective
  OAS_prob.add_desvar('wing.twist_cp', lower=-10., upper=15.)
  OAS_prob.add_desvar('wing.sweep', lower=10., upper=30.)
  OAS_prob.add_desvar('wing.dihedral', lower=-10., upper=20.)
  OAS_prob.add_constraint('wing_perf.CL', equals=0.5)
  OAS_prob.add_objective('wing_perf.CD', scaler=1e4)
  OAS_prob.setup()

  # Actually run the problem
  OAS_prob.run()

v2 Script
^^^^^^^^^

.. literalinclude:: /../../tests/integration_tests/test_v1_aero_opt.py
  :start-after: checkpoint 4
  :end-before: checkpoint 5
  :dedent: 8

Full run scripts
----------------

v1 Script
^^^^^^^^^

.. code-block:: python

  from __future__ import division, print_function
  from time import time
  import numpy as np

  from OpenAeroStruct import OASProblem


  # Set problem type
  prob_dict = {'type' : 'aero',
                'optimize' : True}

  # Instantiate problem and add default surface
  OAS_prob = OASProblem(prob_dict)

  # Create a dictionary to store options about the surface
  surf_dict = {'num_y' : 5,
                'num_x' : 3,
                'wing_type' : 'rect',
                'symmetry' : True,
                'num_twist_cp' : 2}

  # Add the specified wing surface to the problem
  OAS_prob.add_surface(surf_dict)

  # Setup problem and add design variables, constraint, and objective
  OAS_prob.add_desvar('wing.twist_cp', lower=-10., upper=15.)
  OAS_prob.add_desvar('wing.sweep', lower=10., upper=30.)
  OAS_prob.add_desvar('wing.dihedral', lower=-10., upper=20.)
  OAS_prob.add_constraint('wing_perf.CL', equals=0.5)
  OAS_prob.add_objective('wing_perf.CD', scaler=1e4)
  OAS_prob.setup()

  # Actually run the problem
  OAS_prob.run()

v2 Script
^^^^^^^^^

.. literalinclude:: /../../tests/integration_tests/test_v1_aero_opt.py
  :start-after: checkpoint 0
  :end-before: checkpoint 5
  :dedent: 8