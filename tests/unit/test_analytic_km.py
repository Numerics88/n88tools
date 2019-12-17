from __future__ import division
import unittest
from n88tools import analytic_km
import numpy as np

class TestAnalyticKM(unittest.TestCase):

    def test_instantiate_Hexahedron(self):
        '''Can instantiate Hexahedron'''
        a = np.ones((3))
        _ = analytic_km.Hexahedron(a)

    def test_instantiate_IsotropicHexahedron(self):
        '''Can instantiate IsotropicHexahedron'''
        E = 10
        nu = 0.3
        a = np.ones((3))
        _ = analytic_km.IsotropicHexahedron(E, nu, a)

    def test_instantiate_OrthotropicHexahedron(self):
        '''Can instantiate OrthotropicHexahedron'''
        E = 10 * np.ones((3))
        nu = 0.3 * np.ones((3))
        G = 5 * np.ones((3))
        a = np.ones((3))
        _ = analytic_km.OrthotropicHexahedron(E, nu, G, a)

    def test_instantiate_AnisotropicHexahedron(self):
        '''Can instantiate AnisotropicHexahedron'''
        D = np.ones((6,6))
        a = np.ones((3))
        _ = analytic_km.AnisotropicHexahedron(D, a)


if __name__ == '__main__':
    unittest.main()
