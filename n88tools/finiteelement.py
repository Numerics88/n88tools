"""
finiteelement.py

A number of functions that are useful for finite element calculations.

Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
http://www.numerics88.com/
See LICENSE for details.
"""

from __future__ import division

import sys
import numpy
from numpy.core import *
from numpy import linalg
from math import *


def stress_strain_isotropic (E, nu):
    s = nu/(1-nu)
    t = (1-2*nu)/(2*(1-nu))
    D = array ([[ 1, s, s, 0, 0, 0],
                [ s, 1, s, 0, 0, 0],
                [ s, s, 1, 0, 0, 0],
                [ 0, 0, 0, t, 0, 0],
                [ 0, 0, 0, 0, t, 0],
                [ 0, 0, 0, 0, 0, t]]
               )
    D *= E*(1-nu)/((1+nu)*(1-2*nu))
    return D


def stress_strain_orthotropic (E, nu, G):
    nu12 = nu[0]
    nu20 = nu[1]
    nu01 = nu[2]
    nu21 = (E[2]/E[1])*nu12
    nu02 = (E[0]/E[2])*nu20
    nu10 = (E[1]/E[0])*nu01
    D = zeros((6,6), float)
    D[:3,:3] = array ([[     1/E[0], -nu10/E[1], -nu20/E[2]],
                       [ -nu01/E[0],     1/E[1], -nu21/E[2]],
                       [ -nu02/E[0], -nu12/E[1],     1/E[2]]])
    D[3,3] = 1/G[0]
    D[4,4] = 1/G[1]
    D[5,5] = 1/G[2]
    return numpy.linalg.inv(D)


def generate_local_stiffness(reader):
    """"Generates mat_km, a dictionary of local stiffness matrices indexed
    by material id. The local stiffness matrix is the stiffness matrix for
    a single element.
    
    reader must be an opened N88ModelReader object.
    """
    import analytic_km
    
    mat_D = {}
    mat_km = {}
    a = reader.ElementSize
    for id,name in reader.MaterialTable.iteritems():
        material = reader.MaterialDefinitions[name]
        if material['Type'][-9:] == "Isotropic":
            E = material['E']
            nu = material['nu']
            mat_D[id] = stress_strain_isotropic (E, nu)
            element = analytic_km.IsotropicHexahedron(E, nu, a)
            mat_km[id] = element.k
        elif material['Type'][-11:] == "Orthotropic":
            E = material['E']
            nu = material['nu']
            G = material['G']
            mat_D[id] = stress_strain_orthotropic (E, nu, G)
            element = analytic_km.OrthotropicHexahedron(E, nu, G, a)
            mat_km[id] = element.k
        else:
            print "ERROR: Unhandled material type"
            sys.exit(-1)
    return mat_D, mat_km


def renumber_material_ids (mat_D_in, mat_km_in, g_mat_in):
    """Renumbers the material ids, so that they are sequential 0,1,2...

    mat_D_in and mat_km_in must be dictionaries, indexed by arbitrary
    numerical material ID.
    
    Output is mat_D, mat_km, g_mat, where mat_g and mat_km are 3D arrays,
    and g_mat contains the new sequential material IDs, which can be used
    to index into mat_D and mat_km.
    """
    ntype = len(mat_D_in)
    ndof = 24
    nst = 6
    keys = mat_D_in.keys()
    D_values = mat_D_in.values()
    km_values = mat_km_in.values()
    mat_D = zeros((ntype,nst,nst), float64)
    mat_km = zeros((ntype,ndof,ndof), float64)
    for i in range(ntype):
        mat_D[i] = D_values[i]
        mat_km[i] = km_values[i]
    nels = len(g_mat_in)
    g_mat = zeros(nels, int)
    # There might be a faster way to perform the following operation
    for i in range(ntype):
        g_mat += i*(g_mat_in==keys[i])
    id_key = {}
    for i in range(ntype):
        id_key[keys[i]] = i
    return id_key, mat_D, mat_km, g_mat


def generate_global_stiffness(mat_km, g_mat, g_num, nn, sparse=True, fast=True):
    """Generates the Global Stiffness Matrix K.

    Takes no account of constraints.
    g_num must be 0-indexed.
    """
    nels = g_num.shape[0]
    if sparse:
        # Note that duplicate entries in COO format are summed
        # when converting to CSR or CSC.
        import scipy
        from scipy import sparse
        from scipy.sparse import linalg
        if fast:
            from finiteelementfunctions import assembleIJV
            I,J,V = assembleIJV (mat_km, g_mat, g_num)
        else:
            Ie,Je = numpy.mgrid[0:3,0:3]
            Ie = Ie.flatten()
            Je = Je.flatten()
            I = zeros(nels*24**2, int)
            J = zeros(nels*24**2, int)
            V = zeros(nels*24**2, float64)
            k = 0
            for e in xrange(nels):
                m = g_mat[e]
                for i in xrange(8):
                    iglobal = g_num[e,i]
                    for j in xrange(8):
                        jglobal = g_num[e,j]
                        I[k:k+9] = Ie + 3*iglobal
                        J[k:k+9] = Je + 3*jglobal
                        V[k:k+9] = mat_km[m,3*i:3*(i+1),3*j:3*(j+1)].flatten()
                        k += 9
        K = sparse.coo_matrix((V,(I,J)), shape=(nn*3,nn*3), dtype=float64)
    else:
        # Not sparse
        from numpy import linalg
        K = zeros((nn*3,nn*3), float64)
        for e in xrange(nels):
            m = g_mat[e]
            for i in xrange(8):
                iglobal = g_num[e,i]
                for j in xrange(8):
                    jglobal = g_num[e,j]
                    # Here we add the whole 3x3 node-node matrix for (i,j) into the relevant
                    # part of K
                    K[3*iglobal:3*(1+iglobal),3*jglobal:3*(1+jglobal)] += \
                         mat_km[m,3*i:3*(i+1),3*j:3*(j+1)]
    return K


def generate_force_terms(nodeid, sense, values, nn):
    indices = 3*nodeid + sense
    force_terms = zeros(3*nn, float)
    force_terms[indices] = values
    return force_terms


def eliminate_known_values(A, b, indices, values, sparse=True):
    """Eliminate known values from linear equations.

    In two dimensions, the matrix equation is

        Ax = b
    
    Given a set of number of known values for specific x[i], eliminate the
    corresponding columns of A and recalculate b correspondingly.

    Arguments:
    A -- input matrix, with shape (m,n).
    b -- R.H.S. of equation; a vector of length m.
    indices -- a set of integer indices, in the range [0,m), that are to be
               eliminated.  indices must be unique-valued and sorted.
    values -- a vector of the same length as indices, giving the known value
              of x[i] for each i in indices.

    Returns: (A_reduced, b_reduced, g)
    A_reduced -- A with columns corresponing to indices removed.
    b_reduced -- Recalculated b
    g -- The reverse lookup function.  This is a vector of length
         n-len(indices), with values in the range [0,n), such that
         A_reduced[i,j] = A[i,g[j]] and x_reduced[j] = x[g[j]] 
         
    Note that also works in principle on arbitrary dimensions.  For example, 
    A and x vectors, and b a scalar.
    """
    n = A.shape[-1]
    n_known = len(indices)
    assert(n_known == len(values))
    if n_known == 0:
        return (A, b, arange(n))
    n_reduced = n - n_known
    # Construct reverse lookup function g
    g = zeros(n_reduced, int)
    j = 0
    k = 0
    for i in xrange(n):
        if (k < n_known) and (indices[k] == i):
            k += 1
            # While we're at it, check for consistency
            if (k < n_known) and (indices[k] <= indices[k-1]):
                raise exceptions.ValueError("indices is not unique and sorted")
        else:
            g[j] = i
            j += 1
    assert(j == n_reduced)
    if sparse:
        # Convert to Compressed Sparse Column matrix to efficiently select columns
        Acsc = A.tocsc()
        A_reduced = Acsc[:,g]
        # Note: The * operator for sparse matrices is matrix multiplication
        b_reduced = b - Acsc[:,indices] * values
        return A_reduced, b_reduced, g
    else:  # not sparse            
        A_reduced = A.take(g,axis=-1)
        b_reduced = b - inner(A.take(indices,axis=-1), values)
        return A_reduced, b_reduced, g


def shape_der_hexahedron (a, x):
    """Calculates the derivatives of the shape functions (for box-like hexahedron).
    
    a -- dimensions of element
    x -- coordinates of point inside element to sample the derivatives
    
    Note that to get the equivalent to the S&G shape_der function, you
    need to use a=(2,2,2), and then transpose the output.
    """
    ndim = 3
    nod = 8

    s = 1 - x/a
    sp = -1/a
    t = x/a
    tp = 1/a

    deriv = zeros((nod,ndim),float64)

    deriv[0,0] = sp[0] *  s[1] *  s[2]
    deriv[0,1] =  s[0] * sp[1] *  s[2]
    deriv[0,2] =  s[0] *  s[1] *  sp[2]

    deriv[1,0] = sp[0] *  s[1] *  t[2]
    deriv[1,1] =  s[0] * sp[1] *  t[2]
    deriv[1,2] =  s[0] *  s[1] *  tp[2]
        
    deriv[2,0] = tp[0] *  s[1] *  t[2]
    deriv[2,1] =  t[0] * sp[1] *  t[2]
    deriv[2,2] =  t[0] *  s[1] *  tp[2]
        
    deriv[3,0] = tp[0] *  s[1] *  s[2]
    deriv[3,1] =  t[0] * sp[1] *  s[2]
    deriv[3,2] =  t[0] *  s[1] *  sp[2]

    deriv[4,0] = sp[0] *  t[1] *  s[2]
    deriv[4,1] =  s[0] * tp[1] *  s[2]
    deriv[4,2] =  s[0] *  t[1] *  sp[2]

    deriv[5,0] = sp[0] *  t[1] *  t[2]
    deriv[5,1] =  s[0] * tp[1] *  t[2]
    deriv[5,2] =  s[0] *  t[1] *  tp[2]
        
    deriv[6,0] = tp[0] *  t[1] *  t[2]
    deriv[6,1] =  t[0] * tp[1] *  t[2]
    deriv[6,2] =  t[0] *  t[1] *  tp[2]
        
    deriv[7,0] = tp[0] *  t[1] *  s[2]
    deriv[7,1] =  t[0] * tp[1] *  s[2]
    deriv[7,2] =  t[0] *  t[1] *  sp[2]

    return deriv


def beemat(deriv):
    """Calculate the matrix B, where
    
    strain = B x
    
    for x the displacements of the nodes of a hexahedral element
    with S&G topology. x has dimensions of (nod,nodof) flattened.
    """
    nst = 6
    nod = 8
    nodof = 3
    B = zeros((nst,nod,nodof),float64)
    B[0,:,0] = deriv[:,0]
    B[1,:,1] = deriv[:,1]
    B[2,:,2] = deriv[:,2]
    B[3,:,1] = deriv[:,2]
    B[3,:,2] = deriv[:,1]
    B[4,:,2] = deriv[:,0]
    B[4,:,0] = deriv[:,2]
    B[5,:,0] = deriv[:,1]
    B[5,:,1] = deriv[:,0]
    B.shape = (nst,nod*nodof)
    return B


def invar (stress):
    sigm = (stress[0] + stress[1] + stress[2])/3
    t = sqrt(
           ((stress[0]-stress[1])**2
          + (stress[1]-stress[2])**2
          + (stress[2]-stress[0])**2
          + 6*(stress[3]**2 + stress[4]**2 + stress[5]**2))/3)
    dsbar = sqrt(3/2)*t
    sx = (2*stress[0]-stress[1]-stress[2])/3
    sy = (2*stress[1]-stress[2]-stress[0])/3
    sz = (2*stress[2]-stress[0]-stress[1])/3
    J3 = (sx*sy*sz - sx*stress[3]**2 - sy*stress[4]**2 - sz*stress[5]**2
          + 2*stress[3]*stress[4]*stress[5])
    arg = -3*sqrt(6)*J3/t**3
    if arg >= 1:
        arg = 1
    if arg <= -1:
        arg = -1
    theta = asin(arg)/3
    return sigm, dsbar, theta


def rotate_stress (stress_in, R):
    s = stress_in
    stress = zeros(6, float64)
    stress[:3] = (R[:3,0]**2*s[0]
                + R[:3,1]**2*s[1]
                + R[:3,2]**2*s[2]
                + 2*R[:3,1]*R[:3,2]*s[3]
                + 2*R[:3,2]*R[:3,0]*s[4]
                + 2*R[:3,0]*R[:3,1]*s[5] )
    for i in range(3):
        j = (arange(3)+i)%3
        stress[3+i] = (
            R[j[1],0]*R[j[2],0]*s[0]
          + R[j[1],1]*R[j[2],1]*s[1]
          + R[j[1],2]*R[j[2],2]*s[2]
          + (R[j[1],1]*R[j[2],2] + R[j[2],1]*R[j[1],2])*s[3]
          + (R[j[1],2]*R[j[2],0] + R[j[2],2]*R[j[1],0])*s[4]
          + (R[j[1],0]*R[j[2],1] + R[j[2],0]*R[j[1],1])*s[5] )
    return stress


def principal_stresses (stress):
    M = array ([[stress[0], stress[5], stress[4]],
                [stress[5], stress[1], stress[3]],
                [stress[4], stress[3], stress[2]]])
    return sort(linalg.eigvalsh(M))


def VonMises_yield_function (Y, stress):
    sigma = principal_stresses (stress)
    f = sqrt(((sigma[0] - sigma[1])**2
            + (sigma[1] - sigma[2])**2
            + (sigma[2] - sigma[0])**2)/2) - Y
    return f


def VonMises_dQ (Y, sigm, dsbar, theta):
    dQ_dsigm = 0.0
    dQ_dJ2 = 1.5/dsbar
    dQ_dJ3 = 0.0
    return (dQ_dsigm, dQ_dJ2, dQ_dJ3)


def VonMises_dQ_dsigma (Y, stress):
    dsbar = sqrt(((sigma[0] - sigma[1])**2
                + (sigma[1] - sigma[2])**2
                + (sigma[2] - sigma[0])**2)/2)
    dQ = VonMises_dQ (Y, 0, dsbar, 0)
    M2 = array([[ 2, -1, -1, 0, 0, 0],
                [-1,  2, -1, 0, 0, 0],
                [-1, -1,  2, 0, 0, 0],
                [ 0,  0,  0, 6, 0, 0],
                [ 0,  0,  0, 0, 6, 0],
                [ 0,  0,  0, 0, 0, 6]], float)
    M2 *= 1/3
    return dQ[1]*M2


def VonMises_stable_time_step (E, nu, Y):
    return 4*(1+nu)/(3*E)


def MaximumPrincipalStrain_yield_function (E, nu, epsilon_YT, epsilon_YC, stress):
    sigma = principal_stresses (stress)
    fT = max(sigma[0] - nu*(sigma[1] + sigma[2]),
             sigma[1] - nu*(sigma[2] + sigma[0]),
             sigma[2] - nu*(sigma[0] + sigma[1])) \
         - epsilon_YT * E;
    fC = max(nu*(sigma[1] + sigma[2]) - sigma[0],
             nu*(sigma[2] + sigma[0]) - sigma[1],
             nu*(sigma[0] + sigma[1]) - sigma[2]) \
         - epsilon_YC * E;
    return max (fT, fC)


def MaximumPrincipalStrain_dQ (Y, sigm, dsbar, theta):
    # Uncertain what the best choice for the potential function for
    # maximum principal strain is: for now use von Mises potential function.
    return VonMises_dQ (Y, sigm, dsbar, theta)


def MaximumPrincipalStrain_dQ_dsigma (Y, stress):
    dsbar = sqrt(((sigma[0] - sigma[1])**2
                + (sigma[1] - sigma[2])**2
                + (sigma[2] - sigma[0])**2)/2)
    dQ = MaximumPrincipalStrain_dQ (Y, 0, dsbar, 0)
    M2 = array([[ 2, -1, -1, 0, 0, 0],
                [-1,  2, -1, 0, 0, 0],
                [-1, -1,  2, 0, 0, 0],
                [ 0,  0,  0, 6, 0, 0],
                [ 0,  0,  0, 0, 6, 0],
                [ 0,  0,  0, 0, 0, 6]], float)
    M2 *= 1/3
    return dQ[1]*M2


def MaximumPrincipalStrain_stable_time_step (E, nu, Y):
    # Uncertain what the best choice for the stable time step for
    # maximum principal strain is: for now use von Mises stable time step.
    return VonMises_stable_time_step (E, nu, Y)


def MohrCoulomb_yield_function (c, phi, stress):
    sigma = principal_stresses (stress)
    f = sigma[2] - sigma[0] + (sigma[2] + sigma[0])*sin(phi) - 2*c*cos(phi)
    # Note that the S&G version is 1/2 of this. The yield function really only
    # matters in defining the yield surface f=0, but we will divide by 2 anyway.
    return f/2


def MohrCoulomb_dQ (c, phi, psi, sigm, dsbar, theta):
    dQ_dsigm = sin(psi)
    if sin(theta) > 0.49:
        dQ_dJ2 = (sqrt(3)*0.5-sin(psi)*0.5/sqrt(3))*sqrt(3)*0.5/dsbar
        dQ_dJ3 = 0.
    elif sin(theta) < -0.49:
        dQ_dJ2 = (sqrt(3)*0.5+sin(psi)*0.5/sqrt(3))*sqrt(3)*0.5/dsbar
        dQ_dJ3 = 0.
    else:
        t = sqrt(2/3)*dsbar
        dQ_dJ2 = (cos(theta)/(sqrt(2)*t))*(1 + tan(theta)*tan(3*theta)
                                          + (sin(psi)/sqrt(3))*(tan(3*theta)-tan(theta)))
        dQ_dJ3 = (sqrt(3)*sin(theta) + sin(psi)*cos(theta))/(t**2*cos(3*theta))
    return (dQ_dsigm, dQ_dJ2, dQ_dJ3)


def MohrCoulomb_dQ_dsigma (c, phi, psi, stress):
    sigm, dsbar, theta = invar (stress)
    dQ = MohrCoulomb_dQ (c, phi, psi, sigm, dsbar, theta)
    M1 = array([[1, 1, 1, 0, 0, 0],
                [1, 1, 1, 0, 0, 0],
                [1, 1, 1, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0]], float)
    M1 *= 1/(3*(stress[0] + stress[1] + stress[2]))
    M2 = array([[ 2, -1, -1, 0, 0, 0],
                [-1,  2, -1, 0, 0, 0],
                [-1, -1,  2, 0, 0, 0],
                [ 0,  0,  0, 6, 0, 0],
                [ 0,  0,  0, 0, 6, 0],
                [ 0,  0,  0, 0, 0, 6]], float)
    M2 *= 1/3
    s = array([(3*stress[0] - stress[1] - stress[2])/3,
               (3*stress[1] - stress[2] - stress[0])/3,
               (3*stress[2] - stress[0] - stress[1])/3])
    M3 = zeros((6,6), float)
    M3[:3,:3] = array([[s[0], s[2], s[1]],
                       [s[2], s[1], s[0]],
                       [s[1], s[0], s[2]]])
    M3[:3,3:6] = array([[-2*stress[3],    stress[4],    stress[5]],
                        [   stress[3], -2*stress[4],    stress[5]],
                        [   stress[3],    stress[4], -2*stress[5]]])
    M3[3:6,:3] = transpose(M3[:3,3:6])
    M3[3:6,3:6] = 3*array([[    -s[0], stress[2], stress[1]],
                           [stress[2],     -s[1], stress[0]],
                           [stress[1], stress[0],     -s[2]]])
    M3 *= 1/3
    return dQ[0]*M1 + dQ[1]*M2 + dQ[2]*M3


def MohrCoulomb_stable_time_step (E, nu, c, phi, psi):
    return 4*(1+nu)*(1-2*nu)/(E*(1-2*nu+sin(phi)**2))


def calculate_strain (x, g_g, a):
    nod = 8
    ndim = 3
    nst = 6
    nels = g_g.shape[0]
    # Note that a/2 is coordinates of the center of the element.
    deriv = shape_der_hexahedron(a, a/2)
    B = beemat(deriv)
    x_all_local = x[g_g.flatten()]
    x_all_local.shape = (nels,nod*ndim)
    strain = numpy.tensordot(x_all_local, B, axes=([1],[1]))
    return strain


def calculate_stress (strain, mat_D, g_mat):
    nst = 6
    nels = g_mat.shape[0]
    stress = zeros((nels,nst), float64)
    for iel in xrange(nels):
        stress[iel] = dot(mat_D[g_mat[iel]], strain[iel])
    return stress
    

def calculate_body_load (strain, mat_D, g_mat, g_num, a, nn):
    ndim = 3
    nod = 8
    nels = g_mat.shape[0]
    deriv = shape_der_hexahedron(a, a/2)
    B = beemat(deriv)
    body_load = zeros((nn,ndim), float64)
    for iel in xrange(nels):
        eload = product(a) * dot(transpose(B),
                  dot(mat_D[g_mat[iel]], strain[iel]))
        body_load[g_num[iel,:],:] += eload.reshape((nod,ndim))
    return body_load


def calculate_yield_function (stress, material_table, material_definitions, id_key, g_mat):
    nels = g_mat.size
    f = zeros(nels, float64)
    # Step through each material; treat nonlinear ones
    for id,name in material_table.iteritems():
        material = material_definitions[name]
        if material['Type'] == "VonMisesIsotropic":
            Y = material['Y']
            mask = (g_mat != id_key[id])
            mat_list = numpy.ma.masked_array(arange(nels), mask=mask).compressed()
            for iel in mat_list:
                f[iel] = VonMises_yield_function (Y, stress[iel])
        elif material['Type'] == "MaximumPrincipalStrainIsotropic":
            E = material["E"]
            nu = material["nu"]
            epsilon_YT = material['epsilon_YT']
            epsilon_YC = material['epsilon_YC']
            mask = (g_mat != id_key[id])
            mat_list = numpy.ma.masked_array(arange(nels), mask=mask).compressed()
            for iel in mat_list:
                f[iel] = MaximumPrincipalStrain_yield_function (E, nu, epsilon_YT, epsilon_YC, stress[iel])
        elif material['Type'] == "MohrCoulombIsotropic":
            c = material['c']
            phi = material['phi']
            mask = (g_mat != id_key[id])
            mat_list = numpy.ma.masked_array(arange(nels), mask=mask).compressed()
            for iel in mat_list:
                f[iel] = MohrCoulomb_yield_function (c, phi, stress[iel])
    return f


def calculate_forces (K, x, sparse=True):
    if sparse:
        # Note: The * operator for sparse matrices is matrix multiplication
        b = K * x.flatten()
    else:
        b = dot (K, x.flatten())
    return b
