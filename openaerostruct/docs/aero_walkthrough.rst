.. _Aerodynamic_Optimization_Walkthrough:

Aerodynamic Optimization Walkthrough
====================================

This page documents in detail how we perform aerodynamic optimization in OpenAeroStruct.
You should read this page in its entirety to get a feel for how the model is set up and used.

OpenMDAO basics
---------------

OpenAeroStruct is a tool written using the `OpenMDAO <http://openmdao.org/>`_ framework.
OpenMDAO is an open-source high-performance computing platform for systems analysis and multidisciplinary optimization, written in Python.
Visit OpenMDAO's `documentation <http://openmdao.org/twodocs/versions/latest/index.html>`_ for information on its capabilities and API.

Here is an extremely quick rundown of the basic terminology within OpenMDAO:

- Component: smallest unit of computational work
- Group: container used to build up model hierarchies
- Driver: controls iterative model execution process (e.g., Optimizer, DOE)
- Problem: top level class that contains everything and provide the model execution API

.. figure:: problem_diagram.png
   :align: center
   :width: 50%
   :alt: diagram of the problem structure


Intro to the structure of the OpenAeroStruct problem
----------------------------------------------------

In general, you'll follow these four steps to set up and run a problem in OpenAeroStruct:

1. Define your lifting surface
2. Initialize your problem and add flow conditions
3. Add your design variables, constraints, and objective
4. Set up and run the optimization problem

We'll now investigate these steps individually using an aerodynamic optimization case as an example.
The full run script for this case is in :ref:`Quick_Example`.

1. Define your lifting surface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In the simplest case, you can add a single lifting surface that represents the wing of an aircraft.
To define a lifting surface, you first need to produce a computational mesh that represents that surface.

OpenAeroStruct contains a helper function to create these meshes or you can create your own array of points through another method.
To use OpenAeroStruct's helper function, you need to give it the number of spanwise points, `num_y`, the number of chordwise points, `num_x`, and the type of baseline planform to start with, `wing_type`.
There are two main options for `wing_type`: `rect`, which is a rectangular wing, and `CRM`, which is based on the Common Research Model wing.
This example uses the `CRM` option.
See the `run_scaneagle.py` example in the `examples/` directory for an example that uses the `rect` option.
The `rect` option is convenient for simple rectangular and trapezoidal planforms.
Alternatively, if you want to create your own custom mesh, see :ref:`Custom_Mesh`.

In the code block shown below, we call the helper function to define a mesh and get a starting twist distribution.

.. embed-code::
    aero_walkthrough/part_1.py

There are many options for each surface, and they are loosely organized into the following categories:

- Wing definition (mesh, wing position, symmetry option, etc)
- Geometric variable definitions (span, dihedral, sweep, twist, etc)
- Aerodynamic performance (CL and CD at angle of attack=0)
- Airfoil properties (turbulence transition point, t/c, location of max t)
- Structural properties (E, G, yield stress, location of spar, etc)
- Options for constraints (KS aggregation, monotonic design variables)

.. embed-code::
    aero_walkthrough/part_2.py

2. Initialize your problem and add problem conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Next, we need to initialize the OpenMDAO problem, add flow conditions, and add the groups that contain the analysis components.
In a more complex model, these flow conditions might come from a different OpenMDAO component, but here we hook them up into the model using an independent variable component, or `IndepVarComp`.
Set the values for these parameters that you want to use here.
We then add this component to the OpenMDAO model.

.. embed-code::
    aero_walkthrough/part_3.py

We now need to provide the geometry and analysis groups to the OpenMDAO problem.

We instantiate and add the `Geometry` group, which computes the new mesh shape based on the current values of the design parameters.
In an optimization context, the driver will change these values, and the geometry group computes the new mesh to use in the analysis components.

We then add an `AeroPoint` group, which contains the analysis components to compute the aerodynamic performance of the lifting surface.
Additionally, we promote the flow condition variables from the group up to the model level.
This means that the values in our IndepVarComp can pass data into this `AeroPoint` group, which is how the aerodynamic analysis knows which flow conditions to use.

We need to connect some of the variables from the `Geometry` group into the `AeroPoint` group.
These connections allow information about the mesh to flow through the model correctly.

.. embed-code::
    aero_walkthrough/part_4.py

3. Add your design variables, constraints, and objective
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. note::
  This step is only necessary when performing an optimization and is not needed for only analysis.

With the problem and surfaces defined, we can now add a description of the
optimization problem.

You can use a certain set of parameters as design variables, including:

- angle of attack
- taper
- span
- dihedral
- sweep
- chord distribution
- twist distribution
- shear deformation in x direction
- shear deformation in y direction
- structural spar radius distribution
- structural spar thickness distribution

For the constraints and objective, you can use any outputted variable.
Common constraints include:

- structural failure
- CL = fixed value
- monotonic constraint on spanwise variable (e.g. chord can only get decrease as you go outboard)

Common objectives include:

- weight
- fuel burn
- CL
- CD

We also tell the OpenMDAO problem to record information about each optimization iteration.
This will allow us to visualize the history during and after the optimization.

.. embed-code::
    aero_walkthrough/part_5.py

4. Set up and run the optimization problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With the problem defined, we can now actually run the optimization.
If you only wanted to perform analysis, not optimization, you could use `prob.run_model()` instead of `prob.run_driver()` in the code below.

The code below find the lowest `CD` value while providing a certain amount of lift by constraining `CL`.

.. embed-code::
    aero_walkthrough/part_6.py

.. embed-code::
    openaerostruct.tests.test_aero.Test.test
    :layout: output



Investigation of the problem structure -- N2 diagram
----------------------------------------------------

We'll now take a moment to explain the organization of the aerodynamic model.

.. embed-n2::
    aero_walkthrough/generate_n2.py

Mouse over components and parameters to see the data-passing connections between them.
You can expand this view, click on boxes to zoom in, or right-click to collapse boxes.
This shows the layout of the components within the OpenAeroStruct model.
There's also a help button (the ? mark) on the far right of the top toolbar with information about more features.

To create this diagram for any OpenMDAO problem, add these two lines after you call `prob.setup()`:

.. code-block:: python

  from openmdao.api import n2
  n2(prob)

Use any web browser to open the `.html` file and you can examine your problem layout.
This diagram shows groups in dark blue, components in light blue, as organized by your actual problem hierarchy.
Parameters (inputs and outputs) are shown on the diagonal, with off-diagonal terms representing where an output from a component is passed as input to another component.
Any red parameters shown are unconnected inputs, which means that those parameters are not receiving any data from another component.
In general this is bad -- it suggests that the model is not set up properly -- but if you want to use the default value of those inputs, then it's OK.
For example, we currently have unconnected inputs in the geometry group within OpenAeroStruct, but this is allowable because the default values used in that group are correct and do not modify the mesh.

How to visualize results
------------------------

You can visualize the lifting surface and structural spar using:

.. code-block:: console

  plot_wing aero.db

Here you'll use `aero.db` or the filename for where you saved the problem data.
This will produce a window where you can see how the lifting surface and design variables change with each iteration, as shown below.
You can monitor the results from your optimization as it progresses by checking the `Automatically refresh` button.

.. image:: aero.png
