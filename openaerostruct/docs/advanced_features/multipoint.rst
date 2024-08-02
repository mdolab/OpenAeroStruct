.. _Multipoint Optimization:

Multipoint Optimization
=======================

To simulate multiple flight conditions in a single analysis or optimization, you can add multiple `AeroPoint` or `AerostructPoint` groups to the problem.
This allows you to analyze the performance of the aircraft at multiple flight conditions simultaneously, such as at different cruise and maneuver conditions.


Aerodynamic Optimization Example
--------------------------------
We optimize the aircraft at two cruise flight conditions below.

.. literalinclude:: /../../tests/integration_tests/test_multipoint_aero.py
    :start-after: docs checkpoint 0
    :end-before: docs checkpoint 1
    :dedent: 8


Aerostructural Optimization Example (Q400)
------------------------------------------

This is an additional example of a multipoint aerostructural optimization with the wingbox model using a wing based on the Bombardier Q400.
Here we also create a custom mesh instead of using one provided by OpenAeroStruct.
Make sure you go through the :ref:`Aerostructural_with_Wingbox_Walkthrough` before trying to understand this example.

.. embed-code::
    advanced_features/scripts/wingbox_mpt_Q400_example.py

The following shows a visualization of the results.
As can be seen, there is plenty of room for improvement.
A finer mesh and a tighter optimization tolerance should be used.

.. image:: /advanced_features/figs/wingbox_Q400.png