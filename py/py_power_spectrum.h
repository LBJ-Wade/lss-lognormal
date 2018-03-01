#ifndef PY_POWER_SPECTRUM_H
#define PY_POWER_SPECTRUM_H 1

#include "Python.h"

PyObject* py_power_spectrum_alloc(PyObject* self, PyObject* args);
PyObject* py_power_spectrum_call(PyObject* self, PyObject* args);

#endif
