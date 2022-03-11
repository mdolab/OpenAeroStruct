.. _Aerostruct_ffd:

Aerostructural FFD
==================

In additional to OpenAeroStruct's internal geometry manipulation group, you can also use `pyGeo` to perform free-form deformation (FFD) manipulation on the mesh.
This allows for more general design shape changes and helps sync up geometry changes between meshes from different levels of fidelity.

.. warning::
  This example requires `pyGeo`, an external code developed by the MDO Lab. You can install it from `here <https://github.com/mdolab/pygeo>`_.

.. literalinclude:: /../tests/test_aerostruct_ffd.py
  :start-after: checkpoint 0
  :end-before: checkpoint 1
  :dedent: 8

.. code-block:: python

  print(prob["AS_point_0.fuelburn"][0])
