.. _Mesh and Surface Dict:

Mesh and Surface Dictionaries
=============================

Mesh Dict
---------

Here is a list of the keys and default values of the ``mesh_dict``, which is used to generate a mesh, e.g. ``mesh = generate_mesh(mesh_dict)``.

.. list-table:: Mesh definition
    :widths: 20 20 5 55
    :header-rows: 1

    * - Key
      - Default value
      - Units
      - Description
    * - num_x
      - 3
      - 
      - Number of chordwise vertices.
    * - num_y
      - 5
      - 
      - Number of spanwise vertices for the entire wing. When ``symmetry = True``, the number of vertices for a half wing will be ``(num_y + 1) / 2``.
    * - span
      - 10.0
      - m
      - Full wingspan, even for symmetric cases. 
    * - root_chord
      - 1.0
      - m
      - Root chord length.
    * - span_cos_spacing
      - 0
      - 
      - Spanwise cosine spacing. 0 for uniform spanwise panels, 1 for cosine-spaced panels.
    * - chord_cos_spacing
      - 0
      - 
      - Chordwise cosine spacing.
    * - wing_type
      - "rect"
      - 
      - Initial shape of the wing. ``"rect"`` or ``"CRM"``.
    * - symmetry
      - True
      - 
      - If true, OAS models half of the wing reflected across the plane ``y = 0``.
    * - offset
      - np.array([0, 0, 0])
      - m
      - Coordinates to offset the surface from its default location.
    * - num_twist_cp
      - 2
      - 
      - Number of twist control points. Only relevant when ``wing_type = "CRM"``.
 

Surface Dict
------------
Here is a non-exclusive list of the surface dict keys.
The surface dict will be provided to Groups, including ``Geometry``, ``AeroPoint``, and ``AerostructGeometry``.

.. list-table:: Surface definition
    :widths: 20 20 5 55
    :header-rows: 1

    * - Key
      - Example value
      - Units
      - Description
    * - name
      - "wing"
      - 
      - Name of the surface.
    * - symmetry
      - True or False
      - 
      - If true, OAS models half of the wing reflected across the plane ``y = 0``.
    * - S_ref_type
      - "wetted" or "projected"
      - 
      - How we compute the wing reference area.
    * - mesh
      - 3D ndarray or "gen-mesh"
      - m
      - ``x, y, z`` coordinates of the mesh vertices, can be created by ``generate_mesh``. Can also specify "gen-mesh" to generate a mesh based on taper, sweep, and span variables given in the surface dictionary.
    * - span
      - 10.0
      - m
      - Wing span.
    * - taper
      - 0.5
      - 
      - Wing taper ratio.
    * - sweep
      - 10.0
      - deg
      - Wing sweep angle.
    * - dihedral
      - 5.0
      - deg
      - Wing dihedral.
    * - twist_cp
      - np.array([0, 5])
      - deg
      - B-spline control points for twist distribution. Array convention is ``[wing tip, ..., root]`` in symmetry cases, and ``[tip, ..., root, ... tip]`` when ``symmetry = False``.
    * - chord_cp
      - np.array([0.1, 5])
      - 
      - B-spline control points for chord distribution. This is a chord scaler applied to the initial mesh, not the chord value [m] itself. Array convention is the same as ``twist_cp``.
    * - xshear_cp
      - np.array([0.1, 0.2])
      - m
      - B-spline control points for the x-wise shear deformation of the wing.
    * - yshear_cp
      - np.array([0.1, 0.2])
      - m
      - B-spline control points for the y-wise shear deformation of the wing.
    * - zshear_cp
      - np.array([0.1, 0.2])
      - m
      - B-spline control points for the z-wise shear deformation of the wing.
    * - ref_axis_pos
      - 0.25
      - 
      - Position of reference axis along the chord about which to apply twist, chord, taper, and span geometry transformations. 1 is the trailing edge, 0 is the leading edge.

.. list-table:: Multi-section Surface definition
    :widths: 20 20 5 55
    :header-rows: 1

    * - Key
      - Example value
      - Units
      - Description
    * - is_multi_section
      - True
      - 
      - This key must be present and set to True for the AeroPoint to correctly interpret this surface as multi-section.
    * - num_sections
      - 2
      - 
      - The number of sections in the multi-section surface.
    * - sec_name
      - ["sec0","sec1","sec2"]
      - 
      - Names of the individual sections. Each section must be named and the list length must match the specified number of sections.
    * - meshes
      - "gen-meshes" or [mesh1,mesh2,...]
      - 
      - Supply a list of meshes for each section or "gen-meshes" for automatic mesh generation
    * - root_chord
      - 1.0
      - m
      - Root chord length of the section indicated as "root section"(required if using the built-in mesh generator)
    * - span
      - [10.0,10.0]
      - m
      - Wing span for each section. The list length must match the specified number of sections.
    * - ny
      - [21,21]
      - 
      - Number of spanwise points for each section. The list length must match the specified number of sections. (required if using the built-in multi-section mesh generator)
    * - nx
      - 10
      - 
      - Number of chordwise points. Same for all sections.(required if using the built-in multi-section mesh generator)
    * - bpanels
      - [10,10]
      - 
      - Number of spanwise panels for each section. The list length must match the specified number of sections. An alternative to specifying nx.
    * - cpanels
      - [10,10]
      - 
      - Number of chordwise panels for each section. An alternative to specifying ny.
    * - root_section
      - 1
      - 
      - Root chord length of the section indicated as "root section"(required if using the built-in mesh generator)

.. list-table:: Aerodynamics definitions
    :widths: 20 20 5 55
    :header-rows: 1

    * - Key
      - Example value
      - Units
      - Description
    * - CL0
      - 0.0
      - 
      - Lift coefficient of the surface at 0 angle of attack.
    * - CD0
      - 0.015
      - 
      - Drag coefficient of the surface at 0 angle of attack.
    * - with_viscous
      - True or False
      - 
      - If true, compute viscous drag
    * - with_wave
      - True or False
      - 
      - If true, compute wage drag
    * - groundplane
      - True or False
      - 
      - If true, compute ground effect.
    * - k_lam
      - 0.05
      - 
      - Airfoil property for viscous drag calculation. Percentage of chord with lanimar flow.
    * - t_over_c_cp
      - np.array([0.12, 0.12])
      - 
      - B-spline control points for airfoil thickness-over-chord ratio
    * - c_max_t
      - 0.303
      - 
      - Chordwise nondimensional location of the maximum airfoil thickness.

.. list-table:: Structure definitions
    :widths: 20 20 5 55
    :header-rows: 1

    * - Key
      - Example value
      - Units
      - Description
    * - fem_model_type
      - "tube" or "wingbox"
      - 
      - Structure model.
    * - E
      - 73.1e9
      - Pa
      - Young's modulus
    * - G
      - 27.5e9
      - Pa
      - Shear modulus
    * - yield
      - 420.0e6 / 1.5
      - Pa
      - Allowable yield stress including the safety factor.
    * - mrho
      - 2.78e3
      - kg/m^3
      - Material density
    * - fem_origin
      - 0.35
      - 
      - Normalized chordwise location of the spar
    * - wing_weight_ratio
      - 2.0
      - 
      - Ratio of the total wing weight (including non-structural components) to the wing structural weight.
    * - exact_failure_constraint
      - True or False
      - 
      - If False, we use KS function to aggregate the stress constraint.
    * - struct_weight_relief
      - True or False
      - 
      - Set True to add the weight of the structure to the loads on the structure.
    * - distributed_fuel_weight
      - True or False
      - 
      - Set True to distribute the fuel weight across the entire wing.
    * - fuel_density
      - 803.0
      - kg/m^3
      - Fuel density only needed if the fuel-in-wing volume constraint is used)
    * - Wf_reserve
      - 15000.0
      - kg
      - Reserve fuel mass
    * - n_point_masses
      - 1
      - 
      - Number of point masses in the system (for example, engine)


.. list-table:: Structure parameterization for tubular spar
    :widths: 20 20 5 55
    :header-rows: 1

    * - Key
      - Example value
      - Units
      - Description
    * - thickness_cp
      - np.array([0.01, 0.02])
      - m
      - B-spline control point of the tube thickness distribution.
    * - radius_cp
      - np.array([0.1, 0.2])
      - m
      - B-spline control point of the tube radius distribution.

.. list-table:: Structure parameterization for wingbox
    :widths: 20 20 5 55
    :header-rows: 1

    * - Key
      - Example value
      - Units
      - Description
    * - spar_thickness_cp
      - np.array([0.004, 0.01])
      - m
      - Control point of spar thickness distribution.
    * - skin_thickness_cp
      - np.array([0.005, 0.02])
      - m
      - Control point of skin thickness distribution.
    * - original_wingbox _airfoil_t_over_c
      - 0.12
      - 
      - Thickness-over-chord ratio of airfoil provided for the wingbox cross-section.
    * - strength_factor _for_upper_skin
      - 1.0
      - 
      - A factor to adjust the yield strength of the upper skin relative to the lower skin.
    * - data_x_upper
      - 1D ndarray
      - 
      - ``x`` coordinates of the wingbox cross-section's upper surface for an airfoil with the chord scaled to 1.
    * - data_y_upper
      - 1D ndarray
      - 
      - ``y`` coordinates of the wingbox cross-section's upper surface
    * - data_x_lower
      - 1D ndarray
      - 
      - ``x`` coordinates of the wingbox cross-section's lower surface
    * - data_y_lower
      - 1D ndarray
      - 
      - ``y`` coordinates of the wingbox cross-section's lower surface

.. list-table:: FFD parameters
    :widths: 20 20 5 55
    :header-rows: 1

    * - Key
      - Example value
      - Units
      - Description
    * - mx
      - 2
      - 
      - Number of the FFD control points in the x direction.
    * - my
      - 2
      - 
      - Number of the FFD control points in the y direction.
..
  TODO: list default values (if any), and whethre each key is required or optional.