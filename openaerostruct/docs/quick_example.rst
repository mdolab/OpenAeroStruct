.. _Quick_Example:

Quick Example
=============

Here is an example run script to perform aerodynamic optimization.
We'll go into more detail later about how to set up a model, define the optimization problem, and postprocess results.

.. literalinclude:: /../tests/test_aero.py
    :start-after: checkpoint 0
    :end-before: checkpoint 1
    :dedent: 8

.. code-block::

    Optimization terminated successfully    (Exit mode 0)
            Current function value: [333.89699872]
            Iterations: 18
            Function evaluations: 18
            Gradient evaluations: 18
    Optimization Complete
    -----------------------------------

.. code-block:: python

    print(prob["aero_point_0.wing_perf.CD"][0])

.. code-block::

    0.0333896998716508

.. code-block:: python

    print(prob["aero_point_0.wing_perf.CL"][0])

.. code-block::

    0.49999999999999667

.. code-block:: python

    print(prob["aero_point_0.CM"][1])

.. code-block::

    -1.788555037237238

