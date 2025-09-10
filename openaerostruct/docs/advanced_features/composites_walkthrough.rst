.. _Composites Walkthrough:

Composite material
==================

This page will walk you through the composites model in OpenAeroStruct.
The composites module allows you to define composite material properties and laminate layups for the wing structure.

Model setup
-----------

First, define the following parameters in the `surface` dictionary:

- ``useComposite`` is a boolean variable that is set to True to enable the composites model criteria.
- ``safety_factor`` is the factor of safety used for determining the Tsai-Wu based failure
- ``ply_angles`` is a list of the angles of the plies with respect to the x-axis.
- ``ply_fractions`` is a list of the ply fractions of the plies. (Should be the same length as ``ply_angles``, with the sum of the fractions equal to 1).
- ``E1`` is the modulus of elasticity in the fiber direction.
- ``E2`` is the modulus of elasticity in the transverse direction.
- ``G12`` is the shear modulus.
- ``nu12`` is the Poisson's ratio.
- ``sigma_t1`` is the tensile strength in the fiber direction.
- ``sigma_c1`` is the compressive strength in the fiber direction.
- ``sigma_t2`` is the tensile strength in the transverse direction.
- ``sigma_c2`` is the compressive strength in the transverse direction.
- ``sigma_12max`` is the maximum shear strength.

.. note::
    The composites failure model doesn't use the ``strength_factor_for_upper_skin`` option from the surface dictionary.
    If you want to apply a knockdown factor on the compressive strength to account for buckling, you should scale down the values of ``sigma_c1`` and ``sigma_c2``.

Next, call a utility function to compute the effective E and G for the composite material.
The following function adds ``E`` and ``G`` to ``surface``.

.. code-block:: python

    from openaerostruct.structures.utils import compute_composite_stiffness  # noqa: E402
    compute_composite_stiffness(surf_dict)

The rest of the model setup is the same as the original metallic problem.
OpenAeroStruct will compute the failure metric based on the Tsai-Wu failure criteria instead of the von Mises failure criteria when we set ``useComposite`` to True.
But this is done automatically within ``AerostructPoint``, so you don't need to do anything special.

Theory
------

Approximation of the Moduli of Elasticity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Currently, the moduli of elasticity of the entire FEM spatial beam model are assumed to be isotropic
in 2D plane so as to not change the entire model and is left for the future works.
The values of the moduli of elasticity are found using the following procedure.

The unidirectional ply properties are used to find the stiffness matrix of the plies:

.. math::

    Q = \begin{bmatrix}
    \frac{E_1}{1-\nu_{12}\nu_{21}} & \frac{\nu_{21}E_2}{1-\nu_{12}\nu_{21}} & 0 \\
    \frac{\nu_{12}E_1}{1-\nu_{12}\nu_{21}} & \frac{E_2}{1-\nu_{12}\nu_{21}} & 0 \\
    0 & 0 & G_{12}
    \end{bmatrix}

where :math:`E_1` and :math:`E_2` are the moduli of elasticity in the fiber direction and transverse direction, respectively,
:math:`\nu_{12}` and :math:`\nu_{21}` are the Poisson's ratios, and :math:`G_{12}` is the shear modulus.

The transformed reduced stiffness matrix is found using the following equation:

.. math::

    \bar{Q} = T^T Q T

where :math:`T` is the transformation matrix. The transformation matrix is found using the following equation:

.. math::

    T = \begin{bmatrix}
    \cos \theta & \sin \theta & 0 \\
    -\sin \theta & \cos \theta & 0 \\
    0 & 0 & 1
    \end{bmatrix}

where :math:`\theta` is the angle of the ply with respect to the x-axis.

The effective reduced stiffness matrix for the laminate is found using the weighted sum of the reduced stiffness matrices of the plies,
using their respective ply fraction constituition:

.. math::

    Q_{eff} = \sum_{i=1}^{n} f_i Q_{bar_i}

where :math:`f_i` is the ply fraction of the :math:`i^{th}` ply.

The effective compliance matrix is found using the following equation:

.. math::

    S_{eff} = Q_{eff}^{-1}

The effective laminate properties are found using the following equations:

.. math::
    E = \frac{1}{S_{eff_{11}}}\\
    G = \frac{1}{S_{eff_{66}}}

These moduli of elasticility values are hence used to determine the stiffness matrix of the entire FEM spatial beam model.

Tsai-Wu Failure Criteria
~~~~~~~~~~~~~~~~~~~~~~~~~

Thereafter, at the 4 critical points in the wingbox (mentioned in the aerostruct-wingbox walkthrough),
the strains are calculated for each of the constituent plies by transforming the strains at the critical points to the laminate coordinate system. This is done using the following equation:

.. math::

    \begin{pmatrix}
    \epsilon_1 \\
    \epsilon_2 \\
    \gamma_{12}
    \end{pmatrix}
    =
    [T]
    \begin{pmatrix}
    \epsilon_x \\
    \epsilon_y \\
    \gamma_{xy}
    \end{pmatrix}

The strains are then used to calculate the stresses in the laminate using the following equation:

.. math::

    \begin{pmatrix}
    \sigma_1 \\
    \sigma_2 \\
    \tau_{12}
    \end{pmatrix}
    =
    [Q]
    \begin{pmatrix}
    \epsilon_1 \\
    \epsilon_2 \\
    \gamma_{12}
    \end{pmatrix}

These local axial and shear stresses are then utilized to calculate the value of the **Strength Ratios**, where the coefficients are defined by:

.. math::

    F_{11} = \frac{1}{S_L^{(+)} S_L^{(-)}} \quad \text{and} \quad F_1 = \frac{1}{S_L^{(+)}} - \frac{1}{S_L^{(-)}}

.. math::

    F_{22} = \frac{1}{S_T^{(+)} S_T^{(-)}} \quad \text{and} \quad F_2 = \frac{1}{S_T^{(+)}} - \frac{1}{S_T^{(-)}}

.. math::

    F_{66} = \frac{1}{2 S_{LT}^{2}}

where :math:`S_L^{(+)} \text{and} S_L^{(-)}` are the longitudinal strengths in tension and compression respectively,
:math:`S_T^{(+)} \text{and} S_T^{(-)}` are the transverse strengths in tension and compression respectively and
:math:`S_{LT}^{(+)}` is the shear strength of a ply. The strength ratios are then used to calculate the Tsai-Wu based failure criteria for each ply.
The Tsai-Wu failure criteria is given by:

.. math::

    F_1 \sigma_1 + F_2 \sigma_2 + F_{11} \sigma_1^2 + F_{22} \sigma_2^2 + F_{66} \tau_{12}^2 = 1

In order to implement the safety factor in the Tsai-Wu failure criteria, the equation is re-written as:

.. math::
    a &= F_1 \sigma_1 + F_2 \sigma_2 \\
    b &= F_{11} \sigma_1^2 + F_{22} \sigma_2^2 + F_{12} \sigma_1 \sigma_2

We hence calculate the **Strength Ratios** using the formula:

.. math::

    SR = \frac{1}{2} (a + \sqrt{a^2 + 4 b})

The strength ratio values hence calculated for each ply (determined by the length of ``ply_angles``) at each critical point (4 total),
(hence 4 x ``numplies`` strength ratio values for each beam element) for all beam elements are aggregated using a **KS Aggregate** function:

.. math::

    \hat{g}_{KS}(\rho, g) = \max_j g_j + \frac{1}{\rho} \ln \left( \sum_{j=1}^{n_g} \exp \left( \rho (g_j - \max_j g_j) \right) \right)


where :math:`g` is :math:`\left( \frac{SR}{SR_{\text{lim}}} - 1 \right)` value for each ply and :math:`SR_{\text{lim}}` is defined as:

.. math::

    SR_{\text{lim}} = \frac{1}{FOS}


The failure is determined by the value of :math:`\hat{g}_{KS}(\rho, g)` exceeding 0.

Results
-------

The effect of using composites can be seen in the following figure. A Pareto-optimal front is generated for the wingbox model using Isotropic (Aluminum) and Orthotropic (Carbon Fiber Reinforced Polymer) materials.

.. image:: /advanced_features/figs/compositeModelPareto.png
   :width: 600
   :align: center


Complete script
---------------

.. embed-code::
  openaerostruct.examples.run_aerostruct_uCRM_composite_wingbox
