.. _Multipoint:

Multipoint
==========

To simulate multiple flight conditions in a single analysis or optimization, you can add multiple `AeroPoint` or `AerostructPoint` groups to the problem.
This allows you to analyze the performance of the aircraft at multiple flight conditions simultaneously, such as at different cruise and maneuver conditions.
We optimize the aircraft at two cruise flight conditions below.

.. embed-code::
    openaerostruct.tests.test_multipoint_aero.Test.test
