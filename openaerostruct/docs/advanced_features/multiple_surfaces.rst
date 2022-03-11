.. _Multiple_Lifting_Surfaces:

Multiple Lifting Surfaces
=========================

It's easily possible to simulate multiple lifting surfaces simultaneously in OpenAeroStruct.
The most straightforward example is a wing and a tail for a conventional airplane, as shown below, though OpenAeroStruct can handle any arbitrary collection of lifting surfaces.

.. literalinclude:: /../tests/test_multiple_aero_analysis.py
    :start-after: checkpoint 0
    :end-before: checkpoint 1
    :dedent: 8

.. code-block:: python

    print(prob["aero_point_0.wing_perf.CD"][0])
    print(prob["aero_point_0.wing_perf.CL"][0])
    print(prob["aero_point_0.CM"][1])
