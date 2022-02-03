.. _Ground Effect:

Ground Effect
=============

For certain flight conditions or types of aircraft, incorporating ground effect may be important.
Mathematically, the ground effect is simply an extra boundary condition imposed such that the velocities normal to the ground plane are zero.
In vortex lattice methods, this is accomplished by mirroring a second copy of the mesh across the ground plane such that the vortex induced velocities cancel out at the ground.
Some VLM solvers, such as AVL, model ground effect by mirroring across an x-y plane.
This is simple to implement but is not strictly correct because of the influence of the angle of attack.

In OpenAeroStruct, the ground effect mirroring plane is parallel to the freestream (influenced by angle of attack).
This means that configurations deep in ground effect (small altitude compared to wingspan) or at higher angles of attack will obtain the correct ground effect correction.

.. image:: groundplane.svg
    :width: 600

To enable ground effect, add a :code:`groundplane: True` attribute to your aerosurfaces, like so:

.. literalinclude:: /../tests/test_aero_ground_effect.py
    :start-after: checkpoint 0
    :end-before: checkpoint 1
    :dedent: 8

.. code-block::

    Optimization terminated successfully    (Exit mode 0)
            Current function value: [333.89684173]
            Iterations: 18
            Function evaluations: 18
            Gradient evaluations: 18
    Optimization Complete
    -----------------------------------

.. code-block:: python

    print(prob["aero_point_0.wing_perf.CD"][0])

.. code-block::

    0.03338968417307574

.. code-block:: python

    print(prob["aero_point_0.CM"][1]

.. code-block::

    0.4999999999999967

.. code-block:: python

    prob["height_agl"] = 10.0
    prob.run_driver()

.. code-block::

    Optimization terminated successfully    (Exit mode 0)
            Current function value: [291.45613949]
            Iterations: 19
            Function evaluations: 19
            Gradient evaluations: 19
    Optimization Complete
    -----------------------------------

.. code-block:: python

    print(prob["aero_point_0.wing_perf.CD"][0])

.. code-block::

    0.029145613948518823

.. code-block:: python

    print(prob["aero_point_0.wing_perf.CL"][0])

.. code-block::

    0.4999999999983347

.. code-block:: python

    print(prob["aero_point_0.CM"][1])

.. code-block::

    -1.7719184423417516

.. code-block:: python

    totals = prob.check_totals(
        of=["aero_point_0.wing_perf.CD", "aero_point_0.wing_perf.CL"],
        wrt=["wing.twist_cp", "height_agl"],
        compact_print=True,
        out_stream=None,
    )
    assert_check_totals(totals, atol=1e-2, rtol=1e-5)

If groundplane is turned on for an AeroPoint or AeroStructPoint, a new input will be created (height_agl) which represents the distance from the origin (in airplane coordinates) to the ground plane.
The default value, 8000 meters, produces essentially zero ground effect.

Note that symmetry must be turned on for the ground effect correction to be used.
Also, crosswind (beta) may not be used when ground effect is turned on.
Finally, take care when defining geometry and run cases that your baseline mesh does not end up below the ground plane.
This can occur for wings with long chord, anhedral, shear, tail surfaces located far behind the wing, high angles of attack, or some combination.

The following plots (generated using the :code:`examples/drag_polar_ground_effect.py` file) illustrate the effect of the ground plane on a rectangular wing with aspect ratio 12.
As the wing approaches the ground, induced drag is significantly reduced compared to the free-flight induced drag.
These results are consistent with published values in the literature, for example "Lifting-Line Predictions for Induced Drag and Lift in Ground Effect" by Phillips and Hunsaker.

.. image:: ground_effect_correction.png
    :width: 600

.. image:: ground_effect_polars.png
    :width: 600


