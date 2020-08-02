.. _Aerostructural_with_wingbox_(Q400):

Aerostructural with wingbox (Q400)
==================================

This is an additional example of a multipoint aerostructural optimization with the wingbox model using a wing based on the Bombardier Q400.
Here we also create a custom mesh instead of using one provided by OpenAeroStruct.
Make sure you go through the :ref:`Aerostructural_with_Wingbox_Walkthrough` before trying to understand this example.

.. embed-code::
    openaerostruct.docs.wingbox_mpt_Q400_example
    :layout: interleave

The following shows a visualization of the results.
As can be seen, there is plenty of room for improvement.
A finer mesh and a tighter optimization tolerance should be used.

.. image:: ../wingbox_Q400.png
