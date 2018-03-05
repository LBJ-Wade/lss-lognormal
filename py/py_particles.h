#ifndef PY_PARTICLES_H
#define PY_PARTICLES_H 1

#include "Python.h"

PyMODINIT_FUNC
py_particles_module_init();

PyObject* py_particles_alloc(PyObject* self, PyObject* args);
PyObject* py_particles_asarray(PyObject* self, PyObject* args);

#endif
