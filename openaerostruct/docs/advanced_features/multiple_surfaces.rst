.. _Multiple_Lifting_Surfaces:

Multiple Lifting Surfaces
=========================

It's easily possible to simulate multiple lifting surfaces simultaneously in OpenAeroStruct.
The most straightforward example is a wing and a tail for a conventional airplane, as shown below, though OpenAeroStruct can handle any arbitrary collection of lifting surfaces.

.. literalinclude:: /../../tests/integration_tests/test_multiple_aero_analysis.py
    :start-after: docs checkpoint 0
    :end-before: docs checkpoint 1
    :dedent: 8
