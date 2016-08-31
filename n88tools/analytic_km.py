"""
analytic_km

Classes for generating single-element stiffness matrices.

Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
http://www.numerics88.com/
See LICENSE for details.
"""


from __future__ import division

import numpy
from numpy.core import *


class Hexahedron:

    I1 = array(( 1/3,-1/2, 1/6, 1/2,-1/2, 1,-1/2, -1,
                 1/6,-1/2, 1/3, 1/2, 1/2,-1, 1/2,  1))
    I1.resize((2,2,2,2))
    H = array(((0,0,0),
               (0,0,1),
               (1,0,1),
               (1,0,0),
               (0,1,0),
               (0,1,1),
               (1,1,1),
               (1,1,0)))
    i,m,j,n = numpy.mgrid[0:8,0:3,0:8,0:3]
    I_unscaled = (I1[H[i,0], 1*(m==0), H[j,0], 1*(n==0)] *
                  I1[H[i,1], 1*(m==1), H[j,1], 1*(n==1)] *
                  I1[H[i,2], 1*(m==2), H[j,2], 1*(n==2)])

    def __init__(self, a):
        self.I = self.I_unscaled * prod(a) / (a[self.m]*a[self.n])


class IsotropicHexahedron(Hexahedron):

    def __init__(self, E, nu, a):
        Hexahedron.__init__(self, a)

        # Shorten notation
        i,m,j,n = self.i,self.m,self.j,self.n
        I = self.I

        # Precalculate some values
        m_1 = (m+1)%3
        m_2 = (m+2)%3
        n_1 = (n+1)%3
        n_2 = (n+2)%3
        m_n = (m-n)%3
        c1 = E*(1-nu)/((1+nu)*(1-2*nu))
        c2 = 0.5*E/(1+nu)
        c3 = E*nu/((1+nu)*(1-2*nu))

        self.k = ((m_n==0) * (c1 * I[i,m,j,n] +
                              c2 * (I[i,m_1,j,n_1] + I[i,m_2,j,n_2])) +
                  (m_n!=0) * (c3*I[i,m,j,n] + c2*I[i,n,j,m]))
        self.k.resize((24,24))


class OrthotropicHexahedron(Hexahedron):

    def __init__(self, E, nu, G, a):
        Hexahedron.__init__(self, a)

        # Change flat independent parameters to tensors
        flat_nu = nu
        nu = zeros((3,3), float)
        nu[1,2] = flat_nu[0]
        nu[2,0] = flat_nu[1]
        nu[0,1] = flat_nu[2]
        nu[2,1] = nu[1,2] * E[2] / E[1]
        nu[0,2] = nu[2,0] * E[0] / E[2]
        nu[1,0] = nu[0,1] * E[1] / E[0]
        flat_G = G
        G = zeros((3,3), float)
        G[1,2] = flat_G[0]
        G[2,0] = flat_G[1]
        G[0,1] = flat_G[2]
        G[2,1] = G[1,2]
        G[0,2] = G[2,0]
        G[1,0] = G[0,1]

        # Shorten notation
        i,m,j,n = self.i,self.m,self.j,self.n
        I = self.I

        # Precalculate some values
        m_1 = (m+1)%3
        m_2 = (m+2)%3
        n_1 = (n+1)%3
        n_2 = (n+2)%3
        m_n = (m-n)%3
        p = (3 - m - n)%3   # The %3 is only to keep it in range for irrelevant indices
        xi = (1 - nu[1,2]*nu[2,1] - nu[2,0]*nu[0,2] - nu[0,1]*nu[1,0]
                - nu[1,2]*nu[2,0]*nu[0,1] - nu[2,1]*nu[0,2]*nu[1,0])
        c1 = E[m]*(1 - nu[m_1,m_2]*nu[m_2,m_1])/xi
        c2 = E[m]*(nu[n,m] + nu[p,m]*nu[n,p])/xi

        self.k = ((m_n==0) * (c1 * I[i,m,j,n] +
                              G[m,m_1] * I[i,m_1,j,n_1] +
                              G[m_2,m] * I[i,m_2,j,n_2]) +
                  (m_n!=0) * (c2*I[i,m,j,n] + G[m,n]*I[i,n,j,m]))
        self.k.resize((24,24))


class AnisotropicHexahedron(Hexahedron):

    def __init__(self, D, a):
        Hexahedron.__init__(self, a)

        # Shorten notation - also note that variables are differently named
        s,i,t,j = self.i,self.m,self.j,self.n
        I = self.I

        # Precalculate some values
        i_1 = (i+1)%3
        i_2 = (i+2)%3
        j_1 = (j+1)%3
        j_2 = (j+2)%3

        self.k = (
          # Term for D1.  i.e. in the upper left of D
            D[i,j] * I[s,i,t,j]
          # Terms for D2.  i.e. in the upper right of D
          + D[i,j_1+3] * I[s,i,t,j_2]
          + D[i,j_2+3] * I[s,i,t,j_1]
          # Terms for D3.  i.e. in the lower left of D
          + D[i_1+3,j] * I[s,i_2,t,j]
          + D[i_2+3,j] * I[s,i_1,t,j]
          # Terms for D4.  i.e. in the lower right of D
          + D[i_1+3,j_1+3] * I[s,i_2,t,j_2]
          + D[i_2+3,j_2+3] * I[s,i_1,t,j_1]
          + D[i_1+3,j_2+3] * I[s,i_2,t,j_1]
          + D[i_2+3,j_1+3] * I[s,i_1,t,j_2]
          )
        self.k.resize((24,24))
