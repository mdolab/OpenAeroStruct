.. _Openvsp_Mesh:

OpenVSP-generated Mesh Example
==============================

OpenAeroStruct also has the ability to generate VLM meshes from an OpenVSP ``.vsp3`` file.
In order to use this feature, users must have the OpenVSP Python API installed in their Python environment.
For instructions on how to install OpenVSP's Python API see `here <https://github.com/tadkollar/RevHack2020/blob/master/problems/oas_stability_derivs/openmdao_with_vsp.md>`_.

Here is an example script with an aero surface mesh generated for a Boeing 777 aircraft defined in an OpenVSP model.
This should help you understand how meshes are defined in OpenAeroStruct and how to create them for your own custom planform shapes.
This is yet another alternative to the mesh definition approaches described in :ref:`Aerodynamic_Optimization_Walkthrough` and :ref:`Custom_Mesh`.

The OpenVSP model used in this example is shown below:

.. image:: /advanced_features/figs/vsp777.png

The following shows the portion of the example script in which the user provides the OpenVSP model for the mesh.
In this example we'll only be importing the wing, horizontal, and vertical tail from the OpenVSP model shown above.

.. literalinclude:: /advanced_features/scripts/run_vsp_777.py
   :start-after: checkpoint 0
   :end-before: checkpoint 1

The following shows a visualization of the mesh.

.. image:: /advanced_features/figs/oas_vsp_mesh.png

The complete script for the aerodynamic analysis is as follows.

.. literalinclude:: /advanced_features/scripts/run_vsp_777.py

The user may have noticed that the VLM mesh density is never explicitly defined in the script.
The discretization of the VLM is based on the tessellation refinement define in OpenVSP model.
To get an idea of what the discretization will look like, users can turn on the wireframe model for each surface in OpenVSP.
The spanwise and chordwise lines will roughly correspond to panels in the exported VLM mesh.
More or fewer VLM panels can be added in the chordwise and spanwise direction for each surface by toggling the following OpenVSP parameters.

.. image:: /advanced_features/figs/vsp_chordwise.png
.. image:: /advanced_features/figs/vsp_spanwise.png
