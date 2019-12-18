from __future__ import division
import unittest
from n88tools import transformations
import numpy as np
import math

class TestTransformations(unittest.TestCase):

    def test_rotationZ_phi_zero(self):
        '''Zero phi returns identity for rotationZ'''
        M = transformations.rotationZ(0.0)
        self.assertTrue(
            np.allclose(M, np.eye(3)),
            'rotationZ(0) does not return identity'
        )

    def test_rotationY_phi_zero(self):
        '''Zero phi returns identity for rotationY'''
        M = transformations.rotationY(0.0)
        self.assertTrue(
            np.allclose(M, np.eye(3)),
            'rotationY(0) does not return identity'
        )

    def test_rotationX_phi_zero(self):
        '''Zero phi returns identity for rotationX'''
        M = transformations.rotationX(0.0)
        self.assertTrue(
            np.allclose(M, np.eye(3)),
            'rotationX(0) does not return identity'
        )

    def test_rotationZ_phi_180(self):
        '''180 phi returns negative identity for rotationZ'''
        M = transformations.rotationZ(math.radians(180.0))
        T = np.eye(3)
        T[0,0]*=-1
        T[1,1]*=-1

        self.assertTrue(np.allclose(M,T))

    def test_rotationY_phi_180(self):
        '''180 phi returns negative identity for rotationY'''
        M = transformations.rotationY(math.radians(180.0))
        T = np.eye(3)
        T[0,0]*=-1
        T[2,2]*=-1

        self.assertTrue(np.allclose(M,T))

    def test_rotationX_phi_180(self):
        '''180 phi returns negative identity for rotationX'''
        M = transformations.rotationX(math.radians(180.0))
        T = np.eye(3)
        T[1,1]*=-1
        T[2,2]*=-1

        self.assertTrue(np.allclose(M,T))

    def test_rotationZ(self):
        '''rotationZ returns correct for phi'''
        phi=18.967
        M = transformations.rotationZ(math.radians(phi))
        T = np.array([
            [0.9457060, -0.3250235,  0.0000000],
            [0.3250235,  0.9457060,  0.0000000],
            [0.0000000,  0.0000000,  1.0000000]
        ])

        self.assertTrue(np.allclose(M,T))

    def test_rotationY(self):
        '''rotationY returns correct for phi'''
        phi=18.967
        M = transformations.rotationY(math.radians(phi))
        T = np.array([
            [0.9457060,  0.0000000,  0.3250235],
            [0.0000000,  1.0000000,  0.0000000],
            [-0.3250235,  0.0000000,  0.9457060]
        ])

        self.assertTrue(np.allclose(M,T))

    def test_rotationX(self):
        '''rotationX returns correct for phi'''
        phi=18.967
        M = transformations.rotationX(math.radians(phi))
        T = np.array([
            [1.0000000,  0.0000000,  0.0000000],
            [0.0000000,  0.9457060, -0.3250235],
            [0.0000000,  0.3250235,  0.9457060]
        ])

        self.assertTrue(np.allclose(M,T))

if __name__ == '__main__':
    unittest.main()
