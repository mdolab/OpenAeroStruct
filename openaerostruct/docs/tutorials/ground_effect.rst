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

.. embed-code::
    openaerostruct.tests.test_aero_ground_effect.Test.test
    :layout: interleave

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


