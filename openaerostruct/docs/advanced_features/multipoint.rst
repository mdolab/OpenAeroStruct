.. _Multipoint Optimization:

Multipoint Optimization
=======================

To simulate multiple flight conditions in a single analysis or optimization, you can add multiple `AeroPoint` or `AerostructPoint` groups to the problem.
This allows you to analyze the performance of the aircraft at multiple flight conditions simultaneously, such as at different cruise and maneuver conditions.


Aerodynamic Optimization Example
--------------------------------
We optimize the aircraft at two cruise flight conditions below.

.. literalinclude:: /../tests/test_multipoint_aero.py
    :start-after: checkpoint 0
    :end-before: checkpoint 1
    :dedent: 8

.. code-block:: python

    print(prob["aero_point_0.wing_perf.CL"][0])
    print(prob["aero_point_0.wing_perf.CD"][0])
    print(prob["aero_point_1.wing_perf.CL"][0])
    print(prob["aero_point_1.wing_perf.CD"][0])


Aerostructural Optimization Example (Q400)
------------------------------------------

This is an additional example of a multipoint aerostructural optimization with the wingbox model using a wing based on the Bombardier Q400.
Here we also create a custom mesh instead of using one provided by OpenAeroStruct.
Make sure you go through the :ref:`Aerostructural_with_Wingbox_Walkthrough` before trying to understand this example.

.. literalinclude:: /advanced_features/scripts/wingbox_mpt_Q400_example.py

The following shows a visualization of the results.
As can be seen, there is plenty of room for improvement.
A finer mesh and a tighter optimization tolerance should be used.

.. image:: /advanced_features/figs/wingbox_Q400.png