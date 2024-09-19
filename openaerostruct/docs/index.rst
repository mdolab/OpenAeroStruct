.. OpenAeroStruct documentation master file

OpenAeroStruct Documentation
============================

OpenAeroStruct is a lightweight tool that performs aerostructural optimization using OpenMDAO.
It couples a vortex-lattice method (VLM) and a 6 degrees of freedom (per node) 3-dimensional spatial beam model to simulate aerodynamic and structural analyses using lifting surfaces.
These simulations are wrapped with an optimizer using NASA's OpenMDAO framework.
The analysis and optimization results can be visualized using included tools, producing figures such as this:

.. figure:: /figures/example.png
   :align: center
   :width: 100%
   :alt: sample visualization of aerostructural system

   Aerostructural optimization of the Common Research Model (CRM) wing.

.. figure:: /figures/aerostruct_xdsm.png
   :align: center
   :width: 70%
   :alt: sample XDSM of aerostructural system

   The eXtended Design Structure Matrix (XDSM) of the aerostructural system.

Walkthroughs and Examples
=========================

These first few doc pages go into detail about how to set up and run a problem in OpenAeroStruct.
Please review these at a minimum to understand how aerodynamic, structural, and aerostructural problems are constructed.

.. toctree::
   :maxdepth: 1

   installation.rst
   quick_example.rst
   aero_walkthrough.rst
   struct_example.rst
   aerostructural_index.rst


Advanced Features
=================
Once you have reviewed and understand basic walkthroughs, you can move on to some more advanced features below.

.. toctree::
   :maxdepth: 2

   advanced_features.rst

User Reference
==============
Other reference guide can be found below.

.. toctree::
   :maxdepth: 2

   user_reference.rst

How to Contribute
=================

.. toctree::
   :maxdepth: 1

   how_to_contribute.rst
   
Source Docs
===========

.. toctree::
   :maxdepth: 1

   _srcdocs/index.rst

Notes
=====

This current version of this repository has grown past the previous Matlab implementation. If you are looking for a Matlab-capable version, please see https://github.com/samtx/OpenAeroStruct for the latest version.
