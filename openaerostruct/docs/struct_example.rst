.. _Structural_Optimization_Example:

Structural Optimization
=======================

OpenAeroStruct can also handle structural-only optimization problems.
Here we prescribe a load on the spar and allow the optimizer to vary the structural thickness to minimize weight subject to failure constraints.
Although doing structural-only optimizations is relatively rare, this is a building block towards aerostructural optimization.

For more details about ``mesh_dict`` and ``surf_dict`` in the following script, see :ref:`Mesh and Surface Dict`.

.. embed-code::
  openaerostruct.tests.test_struct.Test.test
  :layout: interleave