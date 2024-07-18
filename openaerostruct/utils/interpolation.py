def get_normalized_span_coords(surface, mid_panel=False):
    """Get the normalised coordinates used for interpolating values along the wingspan

    These normalized coordinates range from 0 at the tip of the wing to 1 at the root.

    Parameters
    ----------
    surface : OpenAeroStruct surface dictionary
        Surface to generate coordinates for
    mid_panel : bool, optional
        Whether the normalized coordinate should be of the panel midpoints rather than the mesh nodes, by default False

    Returns
    -------
    np.array
        Normalized coordinate values
    """
    spanwise_coord = surface["mesh"][0, :, 1]
    span_range = spanwise_coord[-1] - spanwise_coord[0]
    span_offset = spanwise_coord[0]
    if mid_panel:
        x_real = (spanwise_coord[:-1] + spanwise_coord[1:]) / 2
    else:
        x_real = spanwise_coord
    x_norm = (x_real - span_offset) / span_range
    return x_norm
