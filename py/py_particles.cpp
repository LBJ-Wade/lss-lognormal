#include <iostream>
#include <vector>

#include "particle.h"
#include "error.h"
#include "py_particles.h"
#include "py_assert.h"

static void py_particles_free(PyObject *obj);

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

using namespace std;

PyMODINIT_FUNC
py_particles_module_init()
{
  import_array();

  return NULL;
}

PyObject* py_particles_alloc(PyObject* self, PyObject* args)
{
  // _power_particles_alloc()
  // Create a new particles object
  vector<Particle>* const v= new vector<Particle>();
  
  return PyCapsule_New(v, "_Particles", py_particles_free);
}

void py_particles_free(PyObject *obj)
{
  // Delete the particles object
  // Called automatically by Python
  
  vector<Particle>* const v=
    (vector<Particle>*) PyCapsule_GetPointer(obj, "_Particles");
  py_assert_void(v);

  delete v;
}

PyObject* py_particles_asarray(PyObject* self, PyObject* args)
{
  PyObject *py_particles;

  if(!PyArg_ParseTuple(args, "O", &py_particles)) {
    return NULL;
  }

  vector<Particle>* const v=
    (vector<Particle>*) PyCapsule_GetPointer(py_particles, "_Particles");
  py_assert_ptr(v);

  const int nd= 2;
  npy_intp np= v->size();  
  npy_intp dims[]= {np, 3};

  //if(np > 0)
  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, v->data());

  //npy_intp dims[]= {nc, 3};
  //npy_intp strides[]= {(npy_intp) (sizeof(double)*nc*ncz),
  //(npy_intp) (sizeof(double)*ncz),
  //(npy_intp) (sizeof(double))};
  //return PyArray_New(&PyArray_Type, nd, dims, NPY_FLOAT_TYPE, strides,
  //grid->fx, 0, 0, 0);
  //Py_RETURN_NONE;
}
