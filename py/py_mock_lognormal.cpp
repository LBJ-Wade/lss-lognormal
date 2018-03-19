#include <vector>
#include "power_spectrum.h"
#include "mock_lognormal.h"
#include "py_mock_lognormal.h"
#include "py_assert.h"

using namespace std;

PyObject* py_mock_lognormal_generate(PyObject* self, PyObject* args)
{
  // _mock_lognormal_generate(ps,
  //                          nc, boxsize, np, seed,
  //                          growth, bias, f, fix_amplitude,
  //                          particles)

  // Power spectrum is assumed to be at z=0 and converted to RSD
  // P_output = growth^2 (b + f*mu^2)^2 P(k)
  //
  // ps: PowerSpectrum
  // 
  PyObject *py_ps, *py_particles;
  int nc, np, seed, fix_amplitude;
  double boxsize, growth, bias, f;
  
  if(!PyArg_ParseTuple(args, "OidiidddiO",
		       &py_ps,
		       &nc, &boxsize, &np, &seed,
		       &growth, &bias, &f, &fix_amplitude,
		       &py_particles)) {
    return NULL;
  }

  PowerSpectrum* const ps=
    (PowerSpectrum*) PyCapsule_GetPointer(py_ps, "_PowerSpectrum");
  py_assert_ptr(ps);

  vector<Particle>* const particles=
    (vector<Particle>*) PyCapsule_GetPointer(py_particles, "_Particles");
  py_assert_ptr(particles);

  Grid* const grid= new Grid(nc);

  set_power_spectrum(ps, boxsize, growth, bias, f, grid);

  compute_gaussian_power_spectrum(grid); // P(k) => P_g(k)

  generate_delta_k(seed, fix_amplitude, grid); // P_g(k) => delta_g(k)

  compute_lognormal_density(grid); // delta_g(k) => delta(x) ~ exp(delta_g(x))

  generate_particles(grid, seed + 101, np, particles);

  delete grid;

  Py_RETURN_NONE;
}
  

PyObject* py_mock_lognormal_generate_octant(PyObject* self, PyObject* args)
{
  // _mock_lognormal_generate(ps,
  //                          nc, boxsize, nbar, seed,
  //                          growth, bias, f, fix_amplitude,
  //                          particles)

  // Power spectrum is assumed to be at z=0 and converted to RSD
  // P_output = growth^2 (b + f*mu^2)^2 P(k)
  //
  // ps: PowerSpectrum
  // 
  PyObject *py_ps, *py_particles;
  int nc, seed, fix_amplitude;
  double nbar, boxsize, growth, bias, f;
  double r_min, r_max;
  
  if(!PyArg_ParseTuple(args, "OiddidddddiO",
		       &py_ps,
		       &nc, &boxsize, &nbar, &seed,
		       &r_min, &r_max,
		       &growth, &bias, &f, &fix_amplitude,
		       &py_particles)) {
    return NULL;
  }

  PowerSpectrum* const ps=
    (PowerSpectrum*) PyCapsule_GetPointer(py_ps, "_PowerSpectrum");
  py_assert_ptr(ps);

  vector<Particle>* const particles=
    (vector<Particle>*) PyCapsule_GetPointer(py_particles, "_Particles");
  py_assert_ptr(particles);

  Grid* const grid= new Grid(nc);

  set_power_spectrum(ps, boxsize, growth, bias, f, grid);

  compute_gaussian_power_spectrum(grid); // P(k) => P_g(k)

  generate_delta_k(seed, fix_amplitude, grid); // P_g(k) => delta_g(k)

  compute_lognormal_density(grid); // delta_g(k) => delta(x) ~ exp(delta_g(x))

  generate_particles_octant(grid, seed + 101, nbar, r_min, r_max,particles);

  delete grid;

  Py_RETURN_NONE;
}

PyObject* py_random_generate_octant(PyObject* self, PyObject* args)
{
  // _random_generate_ocatant(nc, boxsize, nbar, seed,
  //                          growth, bias, f, fix_amplitude,
  //                          particles)
  PyObject *py_particles;
  int nc, seed;
  double nbar, boxsize;
  double r_min, r_max;
  
  if(!PyArg_ParseTuple(args, "iddiddO",
		       &nc, &boxsize, &nbar, &seed,
		       &r_min, &r_max,
		       &py_particles)) {
    return NULL;
  }

  vector<Particle>* const particles=
    (vector<Particle>*) PyCapsule_GetPointer(py_particles, "_Particles");
  py_assert_ptr(particles);

  generate_randoms_octant(nc, boxsize, seed + 20001, nbar, r_min, r_max,particles);

  Py_RETURN_NONE;
}
  

  
