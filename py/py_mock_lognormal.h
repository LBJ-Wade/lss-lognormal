#ifndef PY_MOCK_LOGNORMAL_H
#define PY_MOCK_LOGNORMAL_H 1

#include "Python.h"

PyObject* py_mock_lognormal_generate(PyObject* self, PyObject* args);
PyObject* py_mock_lognormal_generate_octant(PyObject* self, PyObject* args);
PyObject* py_random_generate_octant(PyObject* self, PyObject* args);

#endif
