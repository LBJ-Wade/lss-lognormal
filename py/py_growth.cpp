#include "growth.h"
#include "py_growth.h"

PyObject* py_growth_D(PyObject* self, PyObject* args)
{
  // _growth_D(a, omega_m)
  double a, omega_m;
  if(!PyArg_ParseTuple(args, "dd", &a, &omega_m)) {
    return NULL;
  }

  return Py_BuildValue("d", growth_D(a, omega_m));
}

PyObject* py_growth_f(PyObject* self, PyObject* args)
{
  // _growth_f(a, omega_m)
  double a, omega_m;
  if(!PyArg_ParseTuple(args, "dd", &a, &omega_m)) {
    return NULL;
  }

  return Py_BuildValue("d", growth_f(a, omega_m));
}

