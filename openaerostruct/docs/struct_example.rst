.. _Structural_Optimization_Example:

Structural Optimization Example
===============================

OpenAeroStruct can also handle structural-only optimization problems.
Here we prescribe a load on the spar and allow the optimizer to vary the structural thickness to minimize weight subject to failure constraints.
Although doing structural-only optimizations is relatively rare, this is a building block towards aerostructural optimization.

.. literalinclude:: /../tests/test_struct.py
  :start-after: checkpoint 0
  :end-before: checkpoint 1
  :dedent: 8

.. code-block::

  Optimization terminated successfully    (Exit mode 0)
            Current function value: [0.71088468]
            Iterations: 14
            Function evaluations: 22
            Gradient evaluations: 14
  Optimization Complete
  -----------------------------------

.. code-block:: python

  print(prob["wing.structural_mass"][0])

.. code-block::

  71088.46823994313