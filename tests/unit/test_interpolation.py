"""
test_interpolation.py

Unit tests for n88tools.interpolation.smooth_element_field.

Copyright (c) 2010-2025, Numerics88 Solutions.
http://www.numerics88.com/
See LICENSE for details.
"""

from __future__ import division
import unittest
import numpy

from n88tools.interpolation import smooth_element_field


def regular_grid (nx, ny, nz):
    """Connectivity of a full nx*ny*nz hexahedral voxel grid.

    Returns (node_numbers, number_of_nodes), 0-based, with the local node
    ordering used by n88model files.
    """
    nnx, nny = nx + 1, ny + 1
    def node (i, j, k):
        return (k*nny + j)*nnx + i
    connectivity = numpy.empty ((nx*ny*nz, 8), dtype=int)
    e = 0
    for k in range (nz):
        for j in range (ny):
            for i in range (nx):
                connectivity[e] = (node(i,j,k),     node(i+1,j,k),
                                   node(i+1,j+1,k), node(i,j+1,k),
                                   node(i,j,k+1),   node(i+1,j,k+1),
                                   node(i+1,j+1,k+1), node(i,j+1,k+1))
                e += 1
    return connectivity, nnx*nny*(nz + 1)


def binomial_reference (values, nx, ny, nz):
    """The separable [1,2,1]/4 kernel applied in each direction.

    This is what the node averaging reduces to in the interior of a regular
    grid, and is computed here by a completely independent route.
    """
    kernel = numpy.array ([1.0, 2.0, 1.0]) / 4
    field = values.reshape (nz, ny, nx).copy()
    for axis in (0, 1, 2):
        field = numpy.apply_along_axis (
            lambda line: numpy.convolve (line, kernel, mode="same"), axis, field)
    return field


class TestSmoothElementField (unittest.TestCase):

    def setUp (self):
        self.nx = self.ny = self.nz = 7
        self.connectivity, self.number_of_nodes = regular_grid (self.nx, self.ny, self.nz)
        self.number_of_elements = self.nx * self.ny * self.nz
        self.values = numpy.random.RandomState(0).normal (size=self.number_of_elements)

    def smooth (self, values, **kwargs):
        return smooth_element_field (values, self.connectivity,
                                     self.number_of_nodes, **kwargs)

    def test_uniform_field_is_preserved_exactly (self):
        # Including at free surfaces, where a naive fixed-stencil kernel would
        # pull values toward the value assumed for the void outside.
        uniform = numpy.full (self.number_of_elements, 3.25)
        numpy.testing.assert_array_equal (self.smooth (uniform), uniform)

    def test_interior_equals_binomial_kernel (self):
        smoothed = self.smooth (self.values).reshape (self.nz, self.ny, self.nx)
        reference = binomial_reference (self.values, self.nx, self.ny, self.nz)
        interior = (slice(1,-1),) * 3
        numpy.testing.assert_allclose (smoothed[interior], reference[interior],
                                       rtol=0, atol=1e-12)

    def test_range_is_not_extended (self):
        # An average of averages cannot leave the original range.
        smoothed = self.smooth (self.values)
        self.assertGreaterEqual (smoothed.min(), self.values.min())
        self.assertLessEqual (smoothed.max(), self.values.max())

    def test_one_based_node_numbers (self):
        numpy.testing.assert_array_equal (
            smooth_element_field (self.values, self.connectivity + 1,
                                  self.number_of_nodes, one_based=True),
            self.smooth (self.values))

    def test_nodes_per_element_override (self):
        class NoShape (object):
            def __init__ (self, array):
                self.array = array
            def __getitem__ (self, key):
                return self.array[key]
        numpy.testing.assert_array_equal (
            smooth_element_field (self.values, NoShape (self.connectivity),
                                  self.number_of_nodes, nodes_per_element=8),
            self.smooth (self.values))

    def test_components_are_independent (self):
        multi = numpy.random.RandomState(1).normal (size=(self.number_of_elements, 6))
        smoothed = self.smooth (multi)
        self.assertEqual (smoothed.shape, multi.shape)
        for c in range (6):
            numpy.testing.assert_allclose (smoothed[:,c], self.smooth (multi[:,c]))

    def test_integer_input_is_promoted (self):
        smoothed = self.smooth (numpy.arange (self.number_of_elements))
        self.assertEqual (smoothed.dtype, numpy.float64)

    def test_input_is_not_modified (self):
        original = self.values.copy()
        self.smooth (self.values)
        numpy.testing.assert_array_equal (self.values, original)

    def test_masked_elements_are_absent_not_zero (self):
        # Smoothing a masked field must give the same answer as deleting those
        # elements from the mesh: a masked element is void, not a zero value.
        mask = numpy.zeros (self.number_of_elements, dtype=bool)
        mask[numpy.random.RandomState(2).choice (
            self.number_of_elements, self.number_of_elements//3, replace=False)] = True
        kept = numpy.logical_not (mask)
        smoothed = self.smooth (numpy.ma.MaskedArray (self.values, mask))
        submesh = smooth_element_field (self.values[kept], self.connectivity[kept],
                                        self.number_of_nodes)
        numpy.testing.assert_array_equal (numpy.ma.getmask (smoothed), mask)
        numpy.testing.assert_allclose (numpy.ma.getdata (smoothed)[kept], submesh)

    def test_uniform_field_with_mask_is_preserved_exactly (self):
        mask = numpy.zeros (self.number_of_elements, dtype=bool)
        mask[::3] = True
        uniform = numpy.ma.MaskedArray (
            numpy.full (self.number_of_elements, 3.25), mask)
        smoothed = self.smooth (uniform)
        numpy.testing.assert_array_equal (
            numpy.ma.getdata (smoothed)[numpy.logical_not (mask)], 3.25)

    def test_disconnected_element_is_unchanged (self):
        # A single element touching nothing else is its own only neighbour.
        connectivity = numpy.arange (8).reshape (1, 8)
        numpy.testing.assert_allclose (
            smooth_element_field (numpy.array ([5.0]), connectivity, 8),
            numpy.array ([5.0]))


if __name__ == "__main__":
    unittest.main()
