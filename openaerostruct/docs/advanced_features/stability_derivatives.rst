.. _Stability_derivatives:

Stability Derivatives
=====================

OpenAeroStruct can compute stability derivatives and use them for design optimization, but the implementation is a bit tricky.
This example shows how to compute the longitudinal stability derivatives and static margins for an aerodynamic lifting surface, but this example can be extended to aerostructural analysis and design optimization.

A challenge associated with stability derivatives in OpenAeroStruct is that for gradient-based optimization, we need derivatives of stability derivatives w.r.t. design variables (e.g., derivatives of CL_alpha w.r.t. wing sweep).
The derivatives of stability derivatives are second-order derivatives of OpenAeroStruct's outputs (e.g., CL), which cannot be computed analytically in OpenAeroStruct.
In fact, this is an OpenMDAO's limitation because the background theory of OpenMDAO is only applicable to first-order derivatives.

To overcome this limitation, we compute the stability derivatives using finite difference, and then analytically differentiate the finite difference component using OpenMDAO.
This can be implemented in multiple ways, but in this example, we use an approach similar to multipoint optimization, where we create two `AeroPoint` instances---one for the specified flight condition, and the other for the perturbed angle of attack for finite difference derivatives.
The script is similar to the :ref:`Multipoint Optimization` example, but in this example, the flight conditions for two AeroPoints are not independent, and we add a few `ExecComps` to compute perturbed angle of attack, stability derivatives, and static margin.

The following script computes CL_alpha, CM_alpha, and static margin for a swept wing.
You can play around with different sweep angles and see how the static margin changes as a function of the sweep angle.
Note that this example does not solve the trim (CM=0).


.. embed-code::
   ../examples/stability_derivatives.py
