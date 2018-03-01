#include <iostream>
#include "power_spectrum.h"
#include "error.h"
#include "py_power_spectrum.h"
#include "py_assert.h"

static void py_power_spectrum_free(PyObject *obj);

using namespace std;

PyObject* py_power_spectrum_alloc(PyObject* self, PyObject* args)
{
  // _power_spectrum_alloc(filename)
  // Create a new power spectrum object
  PyObject* bytes;
  
  if(!PyArg_ParseTuple(args, "O&",
		       PyUnicode_FSConverter, &bytes)) {
    return NULL;
  }

  char* filename;
  Py_ssize_t len;
  PyBytes_AsStringAndSize(bytes, &filename, &len);

  
  PowerSpectrum* ps;

  try {
    ps= new PowerSpectrum(filename);
  }
  catch (ErrorFileNotFound e) {
    PyErr_SetNone(PyExc_FileNotFoundError);
    Py_DECREF(bytes);
    return NULL;
  }
  catch (...){
    PyErr_SetNone(PyExc_OSError);
    Py_DECREF(bytes);
    return NULL;
  }

  Py_DECREF(bytes);
  
  return PyCapsule_New(ps, "_PowerSpectrum", py_power_spectrum_free);
}

void py_power_spectrum_free(PyObject *obj)
{
  // Delete the power spectrum object
  // Called automatically by Python
  PowerSpectrum* const ps=
    (PowerSpectrum*) PyCapsule_GetPointer(obj, "_PowerSpectrum");
  py_assert_void(ps);

  delete ps;
}

PyObject* py_power_spectrum_call(PyObject* self, PyObject* args)
{
  // _power_spectrum_call(_ps, k)
  // return P(k)
  
  PyObject* py_ps;
  double k;
  if(!PyArg_ParseTuple(args, "Od", &py_ps, &k)) {
    return NULL;
  }

  PowerSpectrum* const ps=
    (PowerSpectrum*) PyCapsule_GetPointer(py_ps, "_PowerSpectrum");
  py_assert_ptr(ps);

  return Py_BuildValue("d", ps->P(k));
}
