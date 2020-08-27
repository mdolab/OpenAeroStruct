.. _Custom_Mesh:

User-Provided Mesh Example
==========================

Here is an example script with a custom mesh provided as an array of coordinates.
This should help you understand how meshes are defined in OpenAeroStruct and how to create them for your own custom planform shapes.
This is an alternative to the helper-function approach described in :ref:`Aerodynamic_Optimization_Walkthrough`.

The following shows the portion of the example script in which the user provides the coordinates for the mesh.
This example is for a wing with a kink and two distinct trapezoidal segments.

.. literalinclude:: two_part_wing_custom_mesh.py
   :start-after: checkpoint 0
   :end-before: checkpoint 1

The following shows a visualization of the mesh.

.. image:: two_part_mesh.png

The complete script for the optimization is as follows.
Make sure you go through the :ref:`Aerostructural_with_Wingbox_Walkthrough` before trying to understand this setup.

.. embed-code::
    two_part_wing_custom_mesh.py
    :layout: interleave

There is plenty of room for improvement.
A finer mesh and a tighter optimization tolerance should be used.