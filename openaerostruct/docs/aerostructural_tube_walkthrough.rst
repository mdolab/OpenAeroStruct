.. _Aerostructural_Walkthrough:

Aerostructural Optimization with Tubular Spar
=============================================

With aerodynamic- and structural-only analyses done, we now examine an aerostructural design problem.
The construction of the problem follows the same logic as outlined in :ref:`Aerodynamic_Optimization_Walkthrough`, though with some added details.
For example, we use an `AerostructPoint` group instead of an `AeroGroup` because it contains the additional components needed for aerostructural optimization.
Additionally, we have more variable connections due to the more complex problem formulation.

For more details about ``mesh_dict`` and ``surface`` in the following script, see :ref:`Mesh and Surface Dict`.

.. literalinclude:: /../../tests/integration_tests/test_aerostruct.py
    :start-after: docs checkpoint 0
    :end-before: docs checkpoint 1
    :dedent: 8
