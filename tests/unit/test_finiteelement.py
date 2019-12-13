from __future__ import division
import unittest
from n88tools import finiteelement
import numpy as np

class TestFiniteElement(unittest.TestCase):

    def test_stress_strain_isotropic(self):
        '''Create an isotropic stiffness tensor from (E,nu)'''
        E = 1.0
        nu = 0.3
        D = finiteelement.stress_strain_isotropic(E, nu)

        M = np.array([
            [ 1.34615385,  0.57692308,  0.57692308,  0.,          0.,          0.        ],
            [ 0.57692308,  1.34615385,  0.57692308,  0.,          0.,          0.        ],
            [ 0.57692308,  0.57692308,  1.34615385,  0.,          0.,          0.        ],
            [ 0.,          0.,          0.,          0.38461538,  0.,          0.        ],
            [ 0.,          0.,          0.,          0.,          0.38461538,  0.        ],
            [ 0.,          0.,          0.,          0.,          0.,          0.38461538]
        ])

        self.assertTrue(np.allclose(D, M))

    def test_stress_strain_orthotropic(self):
        '''Create an isotropic stiffness tensor from (E,nu,G)'''
        E = np.array([1, 2, 3])
        nu = np.array([0.2, 0.3, 0.4])
        G = np.array([10, 20, 30])
        D = finiteelement.stress_strain_orthotropic(E, nu, G)

        M = np.array([
            [  1.73431734,   1.58671587,   0.99630996,   0.,           0.,           0.        ],
            [  1.58671587,   3.57933579,   1.54981550,   0.,           0.,           0.        ],
            [  0.99630996,   1.54981550,   3.76383764,   0.,           0.,           0.        ],
            [  0.,           0.,           0.,          10.,           0.,           0.        ],
            [  0.,           0.,           0.,           0.,          20.,           0.        ],
            [  0.,           0.,           0.,           0.,           0.,          30.        ]
        ])

        self.assertTrue(np.allclose(D, M))

if __name__ == '__main__':
    unittest.main()