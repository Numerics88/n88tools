/*=========================================================================

  Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
  http://www.numerics88.com/
  See LICENSE for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include <Python.h>

#ifndef FAIM_finiteelementfunctions_module_h_INCLUDED
#define FAIM_finiteelementfunctions_module_h_INCLUDED

// use C calling conventions - even when using a C++ compiler
#ifdef __cplusplus
extern "C"
{
#endif

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_finiteelementfunctions(void);
#else
PyMODINIT_FUNC initfiniteelementfunctions(void);
#endif

#ifdef __cplusplus
}
#endif

#endif
