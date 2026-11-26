"""
interpolation.py

Interpolation of element fields onto the mesh nodes and back.

This implements the smoothing used by Scanco Medical's fe_post when its
interpol_option is set to 1, so that quantities reported by n88postfaim and
n88pistoia can be compared directly against results produced by that workflow.

Copyright (c) 2010-2025, Numerics88 Solutions.
http://www.numerics88.com/
See LICENSE for details.
"""

from __future__ import division
import numpy


def _node_column (node_numbers, j, one_based):
    """Return local node column j as a 0-based index array.

    Read one column at a time so that a netCDF variable can be passed in
    directly: materialising the whole (n_elements, 8) table costs about 750 MB
    on a 23M element model, which neither postfaim nor pistoia otherwise needs.
    """
    col = numpy.asarray (node_numbers[:, j])
    if one_based:
        col = col - 1
    return col


def smooth_element_field (values, node_numbers, number_of_nodes,
                          one_based=False, nodes_per_element=None):
    """Smooth a per-element field by averaging through the mesh nodes.

    Element values are extrapolated to the element nodes, averaged between
    neighbouring elements, and interpolated back to the element centres:

        node value    = mean of the values of the elements touching that node
        element value = mean of its own node values

    Only elements that exist contribute to a node average, so free surfaces are
    handled correctly without any special case - which matters, because
    surfaces are exactly where extreme values of quantities such as the
    energy-equivalent strain tend to concentrate.

    On a regular voxel mesh, in the interior, this is exactly the separable
    [1,2,1]/4 binomial kernel: weight 8/64 for the element itself, 4/64 for
    each of its 6 face neighbours, 2/64 for each of its 12 edge neighbours and
    1/64 for each of its 8 corner neighbours. A uniform field is reproduced
    exactly, including at surfaces, so the operation cannot introduce a bias.

    The operation is linear, and is applied to whatever quantity is being
    reported: for a derived scalar such as von Mises stress or the
    energy-equivalent strain, smooth the derived scalar rather than the fields
    it was derived from, since the derivation itself is not linear.

    Parameters
    ----------
    values : array of shape (n_elements,) or (n_elements, n_components)
        The element field. A MaskedArray may be passed, in which case masked
        elements are excluded from the node averages entirely - treated as
        absent, exactly like the void outside a free surface - and the same
        mask is applied to the result.
    node_numbers : array of shape (n_elements, nodes_per_element)
        Element connectivity. May be a netCDF variable supporting [:, j]
        slicing. Pass one_based=True to hand in an n88model NodeNumbers
        variable unmodified.
    number_of_nodes : int
        Number of nodes in the mesh.
    one_based : bool
        Subtract 1 from node indices as they are read.
    nodes_per_element : int
        Override, for the case that node_numbers has no usable shape attribute.

    Returns
    -------
    array of the same shape as values, dtype float64. A MaskedArray if values
    was one.
    """
    mask = numpy.ma.getmask (values)
    v = numpy.ma.getdata (values).astype (numpy.float64, copy=True)

    one_dimensional = (v.ndim == 1)
    if one_dimensional:
        v = v.reshape (-1, 1)
    if mask is numpy.ma.nomask:
        component_mask = None
    elif one_dimensional:
        component_mask = mask.reshape (-1, 1)
    else:
        component_mask = mask

    npe = nodes_per_element or node_numbers.shape[1]
    n_elements, n_components = v.shape

    result = numpy.empty_like (v)
    for c in range (n_components):
        column = v[:, c]
        if component_mask is None:
            weights = None
        else:
            # Masked elements are absent, not zero.
            weights = numpy.logical_not (component_mask[:, c]).astype (numpy.float64)
            column = column * weights

        # Element -> node. Accumulate one local node index at a time so that no
        # (n_elements * nodes_per_element) temporary is ever built.
        accumulated = numpy.zeros (number_of_nodes, dtype=numpy.float64)
        count = numpy.zeros (number_of_nodes, dtype=numpy.float64)
        for j in range (npe):
            col = _node_column (node_numbers, j, one_based)
            accumulated += numpy.bincount (col, weights=column,
                                           minlength=number_of_nodes)
            count += numpy.bincount (col, weights=weights,
                                     minlength=number_of_nodes)
            del col
        numpy.divide (accumulated, numpy.maximum (count, 1.0), out=accumulated)

        # Node -> element centre.
        smoothed = numpy.zeros (n_elements, dtype=numpy.float64)
        for j in range (npe):
            col = _node_column (node_numbers, j, one_based)
            smoothed += accumulated[col]
            del col
        smoothed /= npe
        result[:, c] = smoothed

    if one_dimensional:
        result = result[:, 0]
    if mask is not numpy.ma.nomask:
        result = numpy.ma.MaskedArray (result, mask)
    return result
