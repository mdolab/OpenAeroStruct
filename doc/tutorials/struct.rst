.. _Struct:

=========
Structural Optimization
=========

The following Python script performs structural optimization to minimize weight while varying thickness subject to a stress failure constraint.

.. code-block:: python

  from __future__ import division, print_function
  from OpenAeroStruct import OASProblem

  # Set problem type
  prob_dict = {'type' : 'struct',
               'optimize' : True}

  # Instantiate problem and add default surface
  OAS_prob = OASProblem(prob_dict)
  OAS_prob.add_surface({'name' : 'wing',
                        'num_y' : 11,
                        'symmetry' : True})

  # Setup problem and add design variables, constraint, and objective
  OAS_prob.setup()
  OAS_prob.add_desvar('wing.thickness_cp', lower=0.001, upper=0.25, scaler=1e2)
  OAS_prob.add_constraint('wing.failure', upper=0.)
  OAS_prob.add_objective('wing.weight', scaler=1e-3)

  # Actually run the problem
  OAS_prob.run()
  print("\nWing weight:", OAS_prob.prob['wing.weight'])

Which should output:

.. code-block:: console

  Wing weight: 547.375784603

We will now go through each block of code to explain what is going on within OpenAeroStruct.

.. code-block:: python

  from __future__ import print_function
  from OpenAeroStruct import OASProblem

We import the OASProblem class from OpenAeroStruct, which is how we access the methods within OpenAeroStruct.
Additionally, we import print_function to ensure compatibility between Python 2 and Python 3.

.. code-block:: python

  # Set problem type
  prob_dict = {'type' : 'struct',
               'optimize' : True}

We then create a dictionary containing options for the problem we want to solve.
We define our problem as struct-only and that we want to perform optimization.
Please see :meth:`OASProblem.get_default_prob_dict` within :doc:`../source/run_classes` to see the defaults for the problem options dictionary.

.. code-block:: python

  # Instantiate problem and add default surface
  OAS_prob = OASProblem(prob_dict)
  OAS_prob.add_surface({'name' : 'wing',
                        'num_y' : 11,
                        'symmetry' : True})

Next, we add a single lifting surface to the problem.
Even though this is a structures-only problem, we add a lifting surface to define the structure in a manner consistent with the aerostructural case.
This means that we will create a tubular spar based on the lifting surface's span with its element radii set from the lifting surface's chord.

In this case, we provide a name and tell OpenAeroStruct to explicitly model only one half of the surface and compute the effects from the other half of the surface.
This is less computationally expensive than modeling the entire surface.

We then provide the number of spanwise ('num_y') mesh points to use for the surface.
Note that in the aerodynamic case, we specified 'num_x', but here we do not.
'num_x' has no bearing on the structural analysis because there is only one beam regardless of the number of chordwise panels.
These numbers correspond to the entire surface even though we are using symmetric effects.

.. code-block:: python

  # Setup problem and add design variables, constraint, and objective
  OAS_prob.setup()
  OAS_prob.add_desvar('wing.thickness_cp', lower=0.001, upper=0.25, scaler=1e2)
  OAS_prob.add_constraint('wing.failure', upper=0.)
  OAS_prob.add_objective('wing.weight', scaler=1e-3)

First we set up the problem using OASProblem's built-in method and add optimization parameters.
We set our design variables as the b-spline control points for the thickness distribution with bounds between 0.001 and 0.25 meters.
We then set the constraint to not allow the KS aggregated stress measures to fail while we minimize structural weight.

Note that the objective and thickness control points have a scaler value which internally multiplies the values that the optimizer sees.
This is necessary because the optimization problem is better conditioned if the design variables, constraints, and objective are on the same order of magnitude.
The correct scaling parameters are difficult to know before examining the possible design space, so some experimentation may be necessary to find the best scalers.

.. code-block:: python

  # Actually run the problem
  OAS_prob.run()
  print("\nWing weight:", OAS_prob.prob['wing.weight'])

Lastly, we actually run the optimization and print the resulting minimized weight.

We can then visualize the results by running

.. code-block:: bash

  python plot_all.py struct.db
