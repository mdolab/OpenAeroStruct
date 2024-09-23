.. _Custom_Mesh:

Multi-section Mesh Example
==========================

These simple example scipts demonstrate the usage of the multi-section wing geometry features in OpenAeroStruct.
We first start with the induced drag minimization of a simple two-section symmetrical wing.

Let's start with the necessary imports.

.. literalinclude:: /advanced_features/scripts/basic_2_sec.py
   :start-after: checkpoint 0
   :end-before: checkpoint 1


Then we are ready to discuss the multi-section parameterization and multi-section surface dictionary.

.. literalinclude:: /advanced_features/scripts/basic_2_sec.py
   :start-after: checkpoint 1
   :end-before: checkpoint 2

Next, we setup the flow variables.

.. literalinclude:: /advanced_features/scripts/basic_2_sec.py
   :start-after: checkpoint 2
   :end-before: checkpoint 3


Enforcing C0 continuity is important part of a multi-section wing design optimziation problem.
This section describes one of the two ways we can do that in OpenAeroStruct.

.. literalinclude:: /advanced_features/scripts/basic_2_sec.py
   :start-after: checkpoint 3
   :end-before: checkpoint 4


We can now create the aerodynamic analysis group.

.. literalinclude:: /advanced_features/scripts/basic_2_sec.py
   :start-after: checkpoint 4
   :end-before: checkpoint 5

Connecting the geometry and analysis groups requires care when using the multi-section parameterization.

.. literalinclude:: /advanced_features/scripts/basic_2_sec.py
   :start-after: checkpoint 5
   :end-before: checkpoint 6

We can now setup our optimization problem and run it.

.. literalinclude:: /advanced_features/scripts/basic_2_sec.py
   :start-after: checkpoint 6
   :end-before: checkpoint 7

We then finish by plotting the result.

.. literalinclude:: /advanced_features/scripts/basic_2_sec.py
   :start-after: checkpoint 7
   :end-before: checkpoint 8





The following shows a visualization of the resulting mesh.

.. image:: /advanced_features/figs/multi_section_2_sym.png
