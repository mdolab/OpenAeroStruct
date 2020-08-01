.. OpenAeroStruct documentation master file

OpenAeroStruct Documentation
============================

OpenAeroStruct is a lightweight tool that performs aerostructural optimization using OpenMDAO.
It couples a vortex-lattice method (VLM) and a 6 degrees of freedom (per node) 3-dimensional spatial beam model to simulate aerodynamic and structural analyses using lifting surfaces.
These simulations are wrapped with an optimizer using NASA's OpenMDAO framework.
The analysis and optimization results can be visualized using included tools, producing figures such as this:

.. figure:: example.png
   :align: center
   :width: 100%
   :alt: sample visualization of aerostructural system

Walkthroughs and Examples
=========================

These first few doc pages go into detail about how to set up and run a problem in OpenAeroStruct.
Please review these at a minimum to understand how aerodynamic, structural, and aerostructural problems are constructed.

.. toctree::
   :maxdepth: 2

   installation.rst
   quick_example.rst
   aero_walkthrough.rst
   struct_example.rst
   aerostructural_walkthrough.rst
   aerostructural_wingbox_walkthrough.rst
   custom_mesh_example.rst

Once you have reviewed and understand these examples, you can move on to some more advanced examples below.

.. toctree::
   :maxdepth: 2

   more_examples.rst

Other Useful Docs
=================

.. toctree::
   :maxdepth: 2

   specialty_topics.rst
   v1_v2_conversion.rst
   how_to_contribute.rst

Source Docs
===========

.. toctree::
   :maxdepth: 1

   _srcdocs/index.rst

Notes
=====

This current version of this repository has grown past the previous Matlab implementation. If you are looking for a Matlab-capable version, please see https://github.com/samtx/OpenAeroStruct for the latest version.
