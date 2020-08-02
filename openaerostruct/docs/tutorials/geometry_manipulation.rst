.. _Geometry_Creation_and_Manipulation:

==================================
Geometry Creation and Manipulation
==================================

Geometry Creation
-----------------

OpenAeroStruct contains two types of default surfaces: a simple rectangular lifting surface, and a wing modeled after a B777-type aircraft called the Common Research Model (CRM).
See the `run_scaneagle.py` example in the `examples/` directory for an example that uses the `rect` option.
The `rect` option is convenient for simple rectangular and trapezoidal planforms.
See :ref:`Aerodynamic_Optimization_Walkthrough` and :ref:`Aerostructural_with_Wingbox_Walkthrough` for examples that use the CRM-based planforms.

Alternatively, if you want to create your own custom mesh, see :ref:`Custom_Mesh`.

Explanation of design variables
-------------------------------

We currently have eight design variables, which we will detail below. Sweep, taper, and dihedral each have one design variable defined for a single lifting surface, while twist, chord, xshear, zshear, thickness, and radius are arrays that contain values for specific spanwise locations along the surface (there are additional design variables when using the wingbox structural model; see :ref:`Aerostructural_with_Wingbox_Walkthrough` for more).
Along with being used as design variables during optimization, the planform design variables can also be used to manipulate the default surfaces already provided in OpenAeroStruct (i.e., `rect` and `CRM`) to get a desired initial geometry (see the `run_scaneagle.py` example in the `examples/` directory).

Instead of directly controlling the nodal mesh points for the design variables defined as arrays along the span, we vary b-spline control points which influence a b-spline interpolation.
For example, in the figure below, we control the green points as our design variables, which are the b-spline knots.
The blue curve is interpolated from the green points and the blue curve is what would modify the aerodynamic mesh.
In this way we can choose the number of design variables independently of the fidelity of the aerodynamic mesh.

.. image:: bspline.svg

In each of the images below, the left wing is the initial rectangular wing and the right wing is the perturbed wing based on the specific design variable.

Angle of attack
~~~~~~~~~~~~~~~
Angle of attack, or alpha, does not change the geometry of the lifting surface at all.
Instead, it changes the incidence angle of the incoming flow.
This affects the right-hand-side of the linear system that we solve to obtain the aerodynamic circulations.

Taper
~~~~~

The taper design variable is the taper ratio of the wing which linearly varies the chord of the wing.
A value of 1.0 corresponds to a straight rectangular wing, whereas a value less than 1 corresponds to a tapered wing.
Values greater than 1 are possible.

.. image:: taper.svg

Sweep
~~~~~

The sweep design variable performs shearing sweep in which the planform area and the y-coordinate extents remain constant.
Positive angles sweep back.

.. image:: sweep.svg

Dihedral
~~~~~~~~

Positive dihedral rotates the wing such that the tip is higher than the root.
As with taper and sweep, this linearly varies across the span.

.. image:: dihedral.svg

Chord
~~~~~

Whereas taper ratio can only vary the chord linearly across the span, the chord design variable allows for arbitrary chord distributions, as shown below.

.. image:: chord.svg

X Shear
~~~~~~~

This design variable changes the x-coordinate at each certain spanwise location, as shown below.
It can be any arbitrary distribution.
This is a more general form of the sweep variable.

.. image:: xshear.svg

Z Shear
~~~~~~~

This design variable changes the z-coordinate at each certain spanwise location, as shown below.
It can be any arbitrary distribution.
This is a more general form of the dihedral variable.

.. image:: zshear.svg

Twist
~~~~~

Below we show a wing with linear twist variation along the span.
OpenAeroStruct is capable of arbitrary twist distributions.

.. image:: twist.svg

Thickness
~~~~~~~~~

Control the thickness of the structural spar; can have any distribution.
The thickness is added internally to the tubular spar, so we must impose a non-intersection
constraint when optimizing that limits the thickness so it does not go past the physical boundary
of a solid cylinder.
Here we can't reliably see thickness changes, so the color of the spar corresponds to thickness.

.. image:: thickness.svg

Radius
~~~~~~

Control the radius of the structural spar; can have any distribution.
With an aerostructural case, it would make physical sense to have some limit on the radius such
that the spar is not larger than the thickness of the airfoil.
You can set this manually when you set the design variable or you can use the experimental
`SparWithinWing` component.

.. image:: radius.svg

Multiple lifting surfaces
-------------------------

So far we have only discussed cases with a single lifting surface, though OpenAeroStruct can handle multiple surfaces.
For example, you could have a case with a main wing surface and a tail surface as shown below.

.. image:: wing_and_tail.png

Most components operate only on one lifting surface without regard for the others in the problem.
Only two components need to have information from all lifting surfaces -- `AssembleAIC` and `VLMForces`.
`AssembleAIC` considers all lifting surfaces when it constructs the aerodynamic influence coefficient (AIC) matrix.

Utility Scripts
---------------
A few useful scripts can be found in geometry/utils.py, such as
:py:meth:`writing the mesh to a Tecplot file <openaerostruct.geometry.utils.writeMesh>`
and
:py:meth:`mirroring half-meshes to obtain the full mesh <openaerostruct.geometry.utils.getFullMesh>`.
