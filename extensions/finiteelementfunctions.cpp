/*=========================================================================

  Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
  http://www.numerics88.com/
  See LICENSE for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include "finiteelementfunctions.h"
#include "numpy/arrayobject.h"
#include "n88util/array.hpp"


static PyObject *
assembleIJV (PyObject *self, PyObject *args)
  {
  // Note that all variables have to be declared in advance or in scope to
  // allow goto to compile.
  int ok;
  PyObject* mat_km_object = NULL;
  PyArrayObject* mat_km_array = NULL;
  PyObject* g_mat_object = NULL;
  PyArrayObject* g_mat_array = NULL;
  PyObject* g_num_object = NULL;
  PyArrayObject* g_num_array = NULL;
  PyObject * I_array = NULL;
  PyObject * J_array = NULL;
  PyObject * V_array = NULL;
  PyObject* result = NULL;
  size_t nels = 0;
  size_t nod = 0;
  size_t ntype = 0;
  size_t ndof = 0;
  npy_intp dims[3] = {0,0,0};

  ok = PyArg_ParseTuple(args, "OOO", &mat_km_object, &g_mat_object, &g_num_object);
  if (!ok) return NULL;

  mat_km_array = (PyArrayObject *)PyArray_FROM_OTF(
                      mat_km_object, NPY_FLOAT64,  NPY_ARRAY_IN_ARRAY);
  if (mat_km_array == NULL) goto fail;
  if (mat_km_array->nd != 3)
    {
    PyErr_SetString (PyExc_RuntimeError, "mat_km must have dimension 3.");
    goto fail;
    }
  ntype = mat_km_array->dimensions[0];
  ndof = mat_km_array->dimensions[1];
  if (mat_km_array->dimensions[2] != ndof)
    {
    PyErr_SetString (PyExc_RuntimeError, "Second and third dimensions of mat_km must be equal.");
    goto fail;
    }
  if (ndof != 24)
    {
    PyErr_SetString (PyExc_RuntimeError, "Second dimensions of mat_km must have length 24.");
    goto fail;
    }

  g_mat_array = (PyArrayObject *)PyArray_FROM_OTF(
                    g_mat_object, NPY_INT64,  NPY_ARRAY_IN_ARRAY);
  if (g_mat_array == NULL) goto fail;
  if (g_mat_array->nd != 1)
    {
    PyErr_SetString (PyExc_RuntimeError, "g_mat must have dimension 1.");
    goto fail;
    }
  nels = g_mat_array->dimensions[0];

  g_num_array = (PyArrayObject *)PyArray_FROM_OTF(
                   g_num_object, NPY_INT64,  NPY_ARRAY_IN_ARRAY);
  if (g_num_array == NULL) return NULL;
  if (g_num_array->nd != 2)
    {
    PyErr_SetString (PyExc_RuntimeError, "g_num must have dimension 2.");
    goto fail;
    }
  nod = g_num_array->dimensions[1];
  if (nod != 8)
    {
    PyErr_SetString (PyExc_RuntimeError, "g_num second dimension must have length 8.");
    goto fail;      
    }
  if (g_num_array->dimensions[0] != nels)
    {
    PyErr_SetString (PyExc_RuntimeError, "Mismatched sized: g_num and g_mat.");
    goto fail;    
    }

  dims[0] = nels*24*24;
  I_array = PyArray_SimpleNew (1, dims, NPY_INT64);
  if (I_array == NULL) goto fail;
  J_array = PyArray_SimpleNew (1, dims, NPY_INT64);
  if (J_array == NULL) goto fail;
  V_array = PyArray_SimpleNew (1, dims, NPY_FLOAT64);
  if (V_array == NULL) goto fail;
  
  // scope: necessary to allow goto to jump past declarations.
    {
    // Wrap python array data in n88::array objects. This will make it easier
    // to perform operations on the data.
    dims[0] = ntype;
    dims[1] = ndof;
    dims[2] = ndof;
    n88::array<3,npy_float64> mat_km((npy_float64*)(PyArray_DATA(mat_km_array)),
                                     dims[0], dims[1], dims[2]);
    dims[0] = nels;
    n88::array<1,npy_int64> g_mat((npy_int64*)(PyArray_DATA(g_mat_array)), dims[0]);
    dims[0] = nels;
    dims[1] = nod;
    n88::array<2,npy_int64> g_num((npy_int64*)(PyArray_DATA(g_num_array)),
                                   dims[0], dims[1]);
    dims[0] = nels*24*24;
    n88::array<1,npy_int64> I((npy_int64*)(PyArray_DATA(I_array)), dims[0]);
    n88::array<1,npy_int64> J((npy_int64*)(PyArray_DATA(J_array)), dims[0]);
    n88::array<1,npy_float64> V((npy_float64*)(PyArray_DATA(V_array)), dims[0]);

    // Note that PyArray_SimpleNew does not initialize data.
    I.zero();
    J.zero();
    V.zero();

    const npy_int64 Ie[9] = {0, 0, 0, 1, 1, 1, 2, 2, 2};
    const npy_int64 Je[9] = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    
    size_t k = 0;
    for (size_t e=0; e<nels; ++e)
      {
      size_t m = g_mat[e];
      for (size_t i=0; i<8; ++i)
        {
        size_t iglobal = g_num(e,i);
        for (size_t j=0; j<8; ++j)
          {
          size_t jglobal = g_num(e,j);
          for (size_t l=0; l<9; ++l)
            { I[k+l] = Ie[l] + 3*iglobal; }
          for (size_t l=0; l<9; ++l)
            { J[k+l] = Je[l] + 3*jglobal; }
          for (size_t a=0; a<3; ++a)
            for (size_t b=0; b<3; ++b)
              { V[k+3*a+b] = mat_km(m,3*i+a,3*j+b); }
          k += 9;
          }
        }
      }
    }

  result = Py_BuildValue ("(OOO)", I_array, J_array, V_array);
  if (result == NULL) goto fail;

  // I,J,V should be referenced now by 'result'.
  Py_XDECREF (V_array);
  Py_XDECREF (J_array);
  Py_XDECREF (I_array);
  Py_XDECREF (g_num_array);
  Py_XDECREF (g_mat_array);
  Py_XDECREF (mat_km_array);
  return result;

  fail:
    Py_XDECREF (V_array);
    Py_XDECREF (J_array);
    Py_XDECREF (I_array);
    Py_XDECREF (g_num_array);
    Py_XDECREF (g_mat_array);
    Py_XDECREF (mat_km_array);
    return NULL;

  }


static PyMethodDef finiteelementfunctionsMethods[] = 
  {
    {"assembleIJV",  assembleIJV, METH_VARARGS,
     "Some functions to speed up finiteelement module."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
  };


PyMODINIT_FUNC
initfiniteelementfunctions(void)
  {
  (void) Py_InitModule("finiteelementfunctions", finiteelementfunctionsMethods);
  import_array();
  }
