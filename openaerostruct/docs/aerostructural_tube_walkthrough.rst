.. _Aerostructural_Walkthrough:

Aerostructural Optimization with Tubular Spar
=============================================

With aerodynamic- and structural-only analyses done, we now examine an aerostructural design problem.
The construction of the problem follows the same logic as outlined in :ref:`Aerodynamic_Optimization_Walkthrough`, though with some added details.
For example, we use an `AerostructPoint` group instead of an `AeroGroup` because it contains the additional components needed for aerostructural optimization.
Additionally, we have more variable connections due to the more complex problem formulation.

For more details about ``mesh_dict`` and ``surface`` in the following script, see :ref:`Mesh and Surface Dict`.

.. literalinclude:: /../tests/test_aerostruct.py
   :start-after: checkpoint 0
   :end-before: checkpoint 1
   :dedent: 8

.. code-block::

  [0.]

.. literalinclude:: /../tests/test_aerostruct.py
   :start-after: checkpoint 2
   :end-before: checkpoint 3
   :dedent: 8

.. code-block::

  Optimization terminated successfully    (Exit mode 0)
            Current function value: [0.97696333]
            Iterations: 42
            Function evaluations: 47
            Gradient evaluations: 42
  Optimization Complete
  -----------------------------------

.. code-block:: python

  print(print(prob["AS_point_0.fuelburn"][0]))

.. code-block::

  97696.3325251465
