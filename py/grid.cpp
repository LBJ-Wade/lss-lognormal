#include <iostream>
#include <cassert>
#include <cstring>
#include <cmath>
#include "grid.h"

using namespace std;

Grid::Grid(const int nc_, const double boxsize_) :
  nc(nc_), ncz(2*(nc_/2+1)), boxsize(boxsize_)
{
  const size_t ngrid= (size_t) nc*nc*ncz;
  fx= (double*) fftw_malloc(sizeof(double)*ngrid);
  assert(fx);
  fk= reinterpret_cast<fftw_complex*>(fx);
  
  plan_forward= fftw_plan_dft_r2c_3d(nc, nc, nc, fx, fk,
				     FFTW_ESTIMATE);

  plan_inverse= fftw_plan_dft_c2r_3d(nc, nc, nc, fk, fx,
				     FFTW_ESTIMATE);

  cerr << "allocated " << ngrid << " grids\n";
}

Grid::~Grid()
{
  fftw_free(fx);
  fftw_destroy_plan(plan_forward);
  fftw_destroy_plan(plan_inverse);
}

void Grid::fft_forward()
{
  // fk= V/nc^3 sum_x f(x) exp(-ik x)
  
  assert(mode == fft_mode_x);
  fftw_execute(plan_forward);
  mode = fft_mode_k;

  const double fac= pow(boxsize/nc, 3.0);
  const size_t ngrid= (size_t) nc*nc*2*(nc/2 + 1);
  for(int i=0; i<ngrid; ++i)
    fx[i] *= fac;
}

void Grid::fft_inverse()
{
  // fx = 1/V sum_k f(k) exp(ik x)
  assert(mode == fft_mode_k);
  fftw_execute(plan_inverse);
  mode= fft_mode_x;

  const double fac= 1.0/(boxsize*boxsize*boxsize);
  const size_t ngrid= (size_t) nc*nc*2*(nc/2 + 1);

  for(int i=0; i<ngrid; ++i)
    fx[i] *= fac;
}

void Grid::clear()
{
  memset(fx, 0, sizeof(double)*(nc*nc*ncz));
}
