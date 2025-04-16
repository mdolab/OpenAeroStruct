.. _Custom_Mesh:

Multi-section Surfaces
==========================

OpenAeroStruct features the ability to specify surfaces as a series of sequentially connected sections. 
Rather than controling the geometry of the surface as a whole, the optimizer can control the geometric parameters of each section individually using any of the geometric transformations available in OpenAeroStruct. 
This feature was developed for the modular wing morphing applications but can be useful in situations where the user wishes to optimize a particular wing section while leaving the others fixed. 

*Please note that aerostructural optimization with multi-section wings is currently not supported*

This example script demonstrate the usage of the multi-section wing geometry features in OpenAeroStruct.
We first start with the induced drag minimization of a simple two-section symmetrical wing.

Let's start with the necessary imports.

.. literalinclude:: /advanced_features/scripts/basic_2_sec_construction.py
   :start-after: checkpoint 0
   :end-before: checkpoint 1


Then we are ready to discuss the multi-section parameterization and multi-section surface dictionary.

.. literalinclude:: /advanced_features/scripts/basic_2_sec_construction.py
   :start-after: checkpoint 1
   :end-before: checkpoint 2

Next, we setup the flow variables.

.. literalinclude:: /advanced_features/scripts/basic_2_sec_construction.py
   :start-after: checkpoint 2
   :end-before: checkpoint 3


Giving the optimizer control of each wing section without any measure to constrain them can lead to undesirable geometries, such as wing sections separating from each other or large kinks appearing along the span, that cause numerical deficiencies in the optimization. 
This concern is addressed by enforcing C0 continuity along the span-wise junctions between the sections. 
C0 continuity can be enforced for any geometric parameter OAS controls for a given section. 
For example, the tip chord or twist of a given section should match the root chord or twist of the subsequent section. 
There are two approaches for enforcing C0 continuity between sections: a construction-based approach and a constraint-based approach. 

The constaint-based approach involves explicitly setting a position constraint that joins the surface points at section junctions to a specified tolerance. 
This constraint is useful if the user desires a high degree of customization in how section are attached to each other.
The number of linear constraints can quickly grow for problems with many sections.
By reducing the degree of customization in section attachment, it is possible to eliminate these additional constraints and maintain C0 continuity using the construction-based approach.
The remainder of this section describes the construction-based approach while the constraint-based approach is described in the subsequent section.


Continuity by construction is enforced by assigning the B-spline control points located at section edges to the same independent variable controlled by the optimizer. 
Enforcing C0 continuity by construction applies to geometric parameters that employ B-spline parametrization in OAS. 
These include chord distribution, twist distribution, shear distribution in all three directions, and thickness and radius distributions for structural spars. 


.. literalinclude:: /advanced_features/scripts/basic_2_sec_construction.py
   :start-after: checkpoint 3
   :end-before: checkpoint 4

Geometric transformations parameterized by a single value, such as sweep, span, and dihedral, can still be used in the construction approach.
The multi-section unification component built into the multi-section geometry group automatically shifts sections outboard/inboard to maintain leading edge coincidence between sections.
This shifting is able to keep the sections together during span, sweep, and dihedral changes at the leading edge but enforcing full C0 continuity may still require the construction based approach described above.

.. note::
  Mesh shifting cannot compensate for taper transformations. 
  Instead, the constraint-based approach should be used with the scaler taper variable.
  Alternatively, the effect of linear taper can be emulated by using two chord B-spline control points.

.. note::
   Mesh shifting could interfere with using the constraint-based approach to maintain C0 continuity during span, sweep, and dihedral transformations.
   It is advised to disable mesh shifting by passing the following option into the multi-section geometery group.

   .. code-block:: python
      
      shift_uni_mesh = False


We can now continue by creating the aerodynamic analysis group.

.. literalinclude:: /advanced_features/scripts/basic_2_sec_construction.py
   :start-after: checkpoint 4
   :end-before: checkpoint 5

Connecting the geometry and analysis groups requires care when using the multi-section parameterization.
While the multi-section geometry group is effectively a drop-in replacement for the standard geometery group for multi-section wings, it does require that the unified mesh component be connected to the AeroPoint analysis.

.. literalinclude:: /advanced_features/scripts/basic_2_sec_construction.py
   :start-after: checkpoint 5
   :end-before: checkpoint 6

.. note::
   When using a thickness to chord ratio B-spline to account for viscous effects the user should be careful to connect the unified thickess to chord ratio B-spline automatically generated by the multi-section geometery group.

   .. literalinclude:: /advanced_features/scripts/basic_2_sec_visc.py
      :start-after: checkpoint 0
      :end-before: checkpoint 1

We can now setup our optimization problem and run it.

.. literalinclude:: /advanced_features/scripts/basic_2_sec_construction.py
   :start-after: checkpoint 6
   :end-before: checkpoint 7

We then finish by plotting the result.

.. literalinclude:: /advanced_features/scripts/basic_2_sec_construction.py
   :start-after: checkpoint 7
   :end-before: checkpoint 8


The following shows the resulting optimized mesh. 
The result is identical regardless of if the constraint-based or construction-based joining approahces are used.

.. image:: /advanced_features/figs/multi_section_2_sym.png


Constraint-based Approach
==========================

Setting up the multi-section goemetry group with the constrint-based approach for section joining requires a slightly different inital setup.
The constraint can be enforced at the leading and/or trailing edge of each segment junction in any combination of x, y, or z directions.
A fully differentiated implementation that facilitates setting this constraint is incorporated into OAS. This approach is robust but introduces at least two linear constraints per segment junction. 
The steps required to setup the multi-section geometry group by constraint feature differences from the construction-based approach starting at the geometry group setup.


.. literalinclude:: /advanced_features/scripts/basic_2_sec.py
   :start-after: checkpoint 3
   :end-before: checkpoint 4



