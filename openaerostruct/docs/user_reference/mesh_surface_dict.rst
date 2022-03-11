.. _Mesh and Surface Dict:

Mesh and Surface Dictionaries
=============================

Mesh Dict
---------

Here is a list of the keys and default values of the ``mesh-dict`` used to generate a mesh by e.g. `mesh = generate_mesh(mesh_dict)`.

.. list-table::
    :widths: 20 20 60
    :header-rows: 1

    * - Key
      - Default value
      - Description
    * - ``num_x``
      - 3
      - Number of chordwise vertices. Needs to be 2 or an odd number.
    * - ``num_y``
      - 5
      - Number of spanwise vertices for the entire wing. When `symmetry = True`, the numbe of vertices for a half wing will be ``(num_y + 1) / 2``.
    * - ``span``
      - 10.0 (m)
      - Full wingspan, even for symmetric cases. 
    * - ``root_chord``
      - 1.0 (m)
      - Root chord length.
    * - ``span_cos_spacing``
      - 0
      - Spanwise cosine spacing. 0 for uniform spanwise panels, 1 for cosine-spaced panels.
    * - ``chord_cos_spacing``
      - 0
      - Chordwise cosine spacing.
    * - ``wing_type``
      - "rect"
      - Initial shape of the wing. ``"rect"`` or ``"CRM"``.
    * - ``symmetry``
      - True
      - If true, OAS models half of the wing reflected across the plane ``y = 0``.
    * - ``offset``
      - ``np.array([0, 0, 0])``
      - Coordinates to offset the surface from its default location.
 

Surface Dict
------------
TODO