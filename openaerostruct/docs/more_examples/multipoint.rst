.. _Multipoint:

Multipoint
==========

To simulate multiple flight conditions in a single analysis or optimization, you can add multiple `AeroPoint` or `AerostructPoint` groups to the problem.
This allows you to analyze the performance of the aircraft at multiple flight conditions simultaneously, such as at different cruise and maneuver conditions.
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
