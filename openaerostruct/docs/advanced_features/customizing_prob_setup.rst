.. _Customizing_Prob_Setup:

Customizing Problem Setup
=========================

Specifying a Reference Area
---------------------------

In order to use a specific reference area for calculating coefficients, rather than the default approach of area-weighted averaging of each lifting surface, a few modifications must be done. First, an option need to be passed to either ``AeroPoint`` or ``AeroStructPoint`` when it's intialized:

.. code-block:: python

  aero_group = AeroPoint(surfaces=surfaces, user_specified_Sref=True)

or in the aerostructural case,

.. code-block:: python

  AS_point = AerostructPoint(surfaces=surfaces,user_specified_Sref=True)

Next, a new independent variable called ``S_ref_total`` needs to be created in the run script. Here it is set to the value of ``areaRef``. Then, for each analysis point, this variable needs to be connected to each analysis point for performance computation.

.. code-block:: python

  indep_var_comp.add_output('S_ref_total', val=areaRef, units='m**2')
  prob.model.connect('S_ref_total', point_name + '.S_ref_total')