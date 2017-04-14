.. _Aero:

Aerodynamic Optimization
========================

The following Python script performs aerodynamic optimization to minimize drag while varying twist subject to a lift constraint.

.. code-block:: python

  from __future__ import print_function
  from OpenAeroStruct import OASProblem

  # Set problem type and instantiate problem
  prob_dict = {'type' : 'aero',
               'optimize' : True}
  OAS_prob = OASProblem(prob_dict)

  # Add lifting surface
  surf_dict = {'name' : 'wing',
               'symmetry' : True,
               'num_y' : 11,
               'num_x' : 3}
  OAS_prob.add_surface(surf_dict)

  # Add design variables, constraint, and objective and setup problem
  OAS_prob.add_desvar('wing.twist_cp', lower=-10., upper=15.)
  OAS_prob.add_constraint('wing_perf.CL', equals=0.5)
  OAS_prob.add_objective('wing_perf.CD', scaler=1e4)
  OAS_prob.setup()

  # Actually run the problem
  OAS_prob.run()

  print("\nWing CL:", OAS_prob.prob['wing_perf.CL'])
  print("Wing CD:", OAS_prob.prob['wing_perf.CD'])

Which should output the optimization results and then these lines:

.. code-block:: console

  Wing CL: 0.5
  Wing CD: 0.00702481143794

We will now go through each block of code to explain what is going on within OpenAeroStruct.

.. code-block:: python

  from __future__ import print_function
  from OpenAeroStruct import OASProblem

We import the OASProblem class from OpenAeroStruct, which is how we access the methods within OpenAeroStruct.
Additionally, we import print_function to ensure compatibility between Python 2 and Python 3.

.. code-block:: python

  # Set problem type and instantiate problem
  prob_dict = {'type' : 'aero',
               'optimize' : True}
  OAS_prob = OASProblem(prob_dict)

We then create a dictionary containing options for the problem we want to solve.
We define our problem as aero-only and specify that we want to perform optimization.
Please see :meth:`OASProblem.get_default_prob_dict` within the :doc:`../source/run_classes` to see the defaults for the problem options dictionary.

.. code-block:: python

  # Add lifting surface
  surf_dict = {'name' : 'wing',
               'symmetry' : True,
               'num_y' : 11,
               'num_x' : 3}
  OAS_prob.add_surface(surf_dict)

Next, we add a single lifting surface to the problem.
In this case, we provide a name and tell OpenAeroStruct to explicitly model only one half of the surface and compute the effects from the other half of the surface.
This is computationally cheaper than modeling the entire surface.

We then provide the number of spanwise (num_y) and chordwise (num_x) mesh points to use for the surface.
These numbers correspond to the entire surface even though we are using symmetric effects.
So, this wing has 10 spanwise panels and 2 chordwise panels, but we only model 5 spanwise panels and 2 chordwise panels, as shown below.

.. image:: aero_sample.png

.. code-block:: python

  # Add design variables, constraint, and objective and setup problem
  OAS_prob.add_desvar('wing.twist_cp', lower=-10., upper=15.)
  OAS_prob.add_constraint('wing_perf.CL', equals=0.5)
  OAS_prob.add_objective('wing_perf.CD', scaler=1e4)
  OAS_prob.setup()

First we set up the problem using OASProblem's built-in method and add optimization parameters.
We set our design variables as the b-spline control points for the twist distribution with bounds at -10 and 15 degrees.
We then set the constraint to keep CL = 0.5 and the objective to minimize CD.

Note that the objective has a scaler value which internally multiplies the values that the optimizer sees.
This is necessary because the optimization problem is better conditioned if the design variables, constraints, and objective are on the same order of magnitude.
The correct scaling parameters are difficult to know before examining the possible design space, so some experimentation may be necessary to find the best scalers.

.. code-block:: python

  # Actually run the problem
  OAS_prob.run()

  print("\nWing CL:", OAS_prob.prob['wing_perf.CL'])
  print("Wing CD:", OAS_prob.prob['wing_perf.CD'])

Lastly, we actually run the optimization and print the resulting CL and CD.

We can then visualize the results by running

.. code-block:: bash

  python plot_all.py aero.db
