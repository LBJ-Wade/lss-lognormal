#include "Python.h"

#include "py_power_spectrum.h"
#include "py_growth.h"

using namespace std;

static PyMethodDef methods[] = {
  {"_power_spectrum_alloc", py_power_spectrum_alloc, METH_VARARGS,
   "_power_spectrum_alloc(filename); allocate a new PowerSpectrum"},
  {"_power_spectrum_call", py_power_spectrum_call, METH_VARARGS,
   "_power_spectrum_call(_ps, k); return P(k)"},
  
  {"_growth_D", py_growth_D, METH_VARARGS, "growth factor D"},
  {"_growth_f", py_growth_f, METH_VARARGS, "growth rate f"},

  {NULL, NULL, 0, NULL}
};


static struct PyModuleDef module = {
  PyModuleDef_HEAD_INIT,
  "_lsslognormal", // name of this module
  "A package for lognormal mock (large-scale structure cosmology)", // Doc String
  -1,
  methods
};

PyMODINIT_FUNC
PyInit__lsslognormal(void) {
  //py_power_spectrum_module_init();    
  
  return PyModule_Create(&module);
}
