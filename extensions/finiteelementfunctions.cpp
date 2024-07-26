/*=========================================================================

  Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
  http://www.numerics88.com/
  See LICENSE for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/* This removes functions deprecated as of numpy 1.7,
   since such functions will be removed in numpy 2.x */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION


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
  PyArrayObject* I_array = NULL;
  PyArrayObject* J_array = NULL;
  PyArrayObject* V_array = NULL;
  PyObject* result = NULL;
  npy_intp nels = 0;
  npy_intp nod = 0;
  int ntype = 0;
  int ndof = 0;
  npy_intp dims[3] = {0,0,0};

  ok = PyArg_ParseTuple(args, "OOO", &mat_km_object, &g_mat_object, &g_num_object);
  if (!ok) return NULL;

  mat_km_array = (PyArrayObject *)PyArray_FROM_OTF(
                      mat_km_object, NPY_FLOAT64,  NPY_ARRAY_IN_ARRAY);
  if (mat_km_array == NULL) goto fail;
  if (PyArray_NDIM(mat_km_array) != 3)
    {
    PyErr_SetString (PyExc_RuntimeError, "mat_km must have dimension 3.");
    goto fail;
    }
  ntype = PyArray_DIMS(mat_km_array)[0];
  ndof = PyArray_DIMS(mat_km_array)[1];
  if (PyArray_DIMS(mat_km_array)[2] != ndof)
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
  if (PyArray_NDIM(g_mat_array) != 1)
    {
    PyErr_SetString (PyExc_RuntimeError, "g_mat must have dimension 1.");
    goto fail;
    }
  nels = PyArray_DIMS(g_mat_array)[0];

  g_num_array = (PyArrayObject *)PyArray_FROM_OTF(
                   g_num_object, NPY_INT64,  NPY_ARRAY_IN_ARRAY);
  if (g_num_array == NULL) return NULL;
  if (PyArray_NDIM(g_num_array) != 2)
    {
    PyErr_SetString (PyExc_RuntimeError, "g_num must have dimension 2.");
    goto fail;
    }
  nod = PyArray_DIMS(g_num_array)[1];
  if (nod != 8)
    {
    PyErr_SetString (PyExc_RuntimeError, "g_num second dimension must have length 8.");
    goto fail;
    }
  if (PyArray_DIMS(g_num_array)[0] != nels)
    {
    PyErr_SetString (PyExc_RuntimeError, "Mismatched sized: g_num and g_mat.");
    goto fail;
    }

  dims[0] = nels*24*24;
  I_array = (PyArrayObject*)PyArray_SimpleNew (1, dims, NPY_INT64);
  if (I_array == NULL) goto fail;
  J_array = (PyArrayObject*)PyArray_SimpleNew (1, dims, NPY_INT64);
  if (J_array == NULL) goto fail;
  V_array = (PyArrayObject*)PyArray_SimpleNew (1, dims, NPY_FLOAT64);
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
    for (npy_intp e=0; e<nels; ++e)
      {
      size_t m = g_mat[e];
      for (int i=0; i<8; ++i)
        {
        size_t iglobal = g_num(e,i);
        for (int j=0; j<8; ++j)
          {
          size_t jglobal = g_num(e,j);
          for (int l=0; l<9; ++l)
            { I[k+l] = Ie[l] + 3*iglobal; }
          for (int l=0; l<9; ++l)
            { J[k+l] = Je[l] + 3*jglobal; }
          for (int a=0; a<3; ++a)
            for (int b=0; b<3; ++b)
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


#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef finiteelementfunctions =
{
    PyModuleDef_HEAD_INIT,
    "finiteelementfunctions", /* name of module */
    "",          /* module documentation, may be NULL */
    -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    finiteelementfunctionsMethods
};

PyMODINIT_FUNC PyInit_finiteelementfunctions(void)
{
  import_array();
  return PyModule_Create(&finiteelementfunctions);
}

#else

PyMODINIT_FUNC
initfiniteelementfunctions(void)
  {
  (void) Py_InitModule("finiteelementfunctions", finiteelementfunctionsMethods);
  import_array();
  }

#endif
