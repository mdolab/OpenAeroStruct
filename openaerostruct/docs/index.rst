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

---------------
Please Cite Us!
---------------

If you use OpenAeroStruct, please cite the `following paper <https://www.researchgate.net/publication/322991521_Open-source_coupled_aerostructural_optimization_using_Python>__`:

John P. Jasa, John T. Hwang, and Joaquim R. R. A. Martins. "Open-source coupled aerostructural optimization using Python." Structural and Multidisciplinary Optimization 57.4 (2018): 1815-1827. DOI: 10.1007/s00158-018-1912-8

.. code-block:: bibtex

   @article{Jasa2018a,
	   Author = {John P. Jasa and John T. Hwang and Joaquim R. R. A. Martins},
	   Doi = {10.1007/s00158-018-1912-8},
	   Journal = {Structural and Multidisciplinary Optimization},
	   Month = {April},
	   Number = {4},
	   Pages = {1815--1827},
	   Publisher = {Springer},
	   Title = {Open-source coupled aerostructural optimization using {Python}},
	   Volume = {57},
	   Year = {2018}
   }


If you use the wingbox structural model, fuel-weight inertial loads, or structural-weight inertial loads, please also cite the `following paper <https://www.researchgate.net/publication/327654423_Low-Fidelity_Aerostructural_Optimization_of_Aircraft_Wings_with_a_Simplified_Wingbox_Model_Using_OpenAeroStruct>__`:

Shamsheer S. Chauhan and Joaquim R. R. A. Martins, “Low-Fidelity Aerostructural Optimization of Aircraft Wings with a Simplified Wingbox Model Using OpenAeroStruct,” Proceedings of the 6th International Conference on Engineering Optimization, EngOpt 2018, Springer, Lisbon, Portugal, September 2018, pp. 418–431. doi:10.1007/978-3-319-97773-7 38

.. code-block:: bibtex

   @inproceedings{Chauhan2018b,
	   Author = {Shamsheer S. Chauhan and Joaquim R. R. A. Martins},
	   Address = {Lisbon, Portugal},
	   Booktitle = {Proceedings of the 6th International Conference on Engineering Optimization, EngOpt 2018},
	   Doi = {10.1007/978-3-319-97773-7_38},
	   Pages = {418-431},
	   Publisher = {Springer},
	   Title = {Low-Fidelity Aerostructural Optimization of Aircraft Wings with a Simplified Wingbox Model Using {OpenAeroStruct}},
	   Year = {2018}
   }


If using point-mass loads or thrust loads, please cite the following `paper (https://www.researchgate.net/publication/333806174_How_Certain_Physical_Considerations_Impact_Aerostructural_Wing_Optimization)__`.

John P. Jasa, Shamsheer S. Chauhan, Justin S. Gray, and Joaquim R. R. A. Martins, “How Certain Physical Considerations Impact Aerostructural Wing Optimization,” AIAA/ISSMO Multidisciplinary Analysis and Optimization Conference, Dallas, TX, 2019. doi:10.2514/6.2019-3242

.. code-block:: bibtex

   @inproceedings{Jasa2019c,
	   Author = {John P. Jasa and Shamsheer S. Chauhan and Justin S. Gray and Joaquim R. R. A. Martins},
	   Address = {Dallas, TX},
	   Booktitle = {AIAA/ISSMO Multidisciplinary Analysis and Optimization Conference},
	   Doi = {10.2514/6.2019-3242},
	   Title = {How Certain Physical Considerations Impact Aerostructural Wing Optimization},
	   Month = {June},
	   Year = {2019}}
   }


Notes
=====

This current version of this repository has grown past the previous Matlab implementation. If you are looking for a Matlab-capable version, please see https://github.com/samtx/OpenAeroStruct for the latest version.
