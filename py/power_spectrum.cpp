//
// Linear power spectrum
//
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include "power_spectrum.h"
#include "error.h"

using namespace std;

PowerSpectrum::PowerSpectrum(const char filename[]) :
  log_k(0), log_P(0), k_min(0.0), k_max(0.0)
{
  load(filename);
}

PowerSpectrum::~PowerSpectrum()
{
  free(log_k);
  gsl_interp_accel_free(acc);
  gsl_interp_free(interp);
}

void PowerSpectrum::load(const char filename[])
{
  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    throw ErrorFileNotFound();
  }

  vector<double> v_k, v_P;
  double k, P;

  char line[1024];
  
  while(fgets(line, 1024, fp)) {
    if(line[0] == '#')
      continue;

    int ret = sscanf(line, "%lg %lg", &k, &P);

    if(ret != 2) {
      cerr << "Error: unable to parse line, " << line;
      throw ErrorIO();
    }

    v_k.push_back(k);
    v_P.push_back(P);
  }

  int ret= fclose(fp); assert(ret == 0);
  
  // Allocate ps->log_k, ps->log_P and fill the arrays
  const size_t nlines= v_k.size(); assert(v_P.size() == nlines);
  
  log_k= (double*) malloc(2*nlines*sizeof(double));
  assert(log_k);
  log_P= log_k + nlines;

  for(size_t j=0; j<nlines; j++) {
    log_k[j]= log(v_k[j]);
    log_P[j]= log(v_P[j]);
  }
  
  interp= gsl_interp_alloc(gsl_interp_cspline, nlines);
  acc= gsl_interp_accel_alloc();

  const size_t n_required= gsl_interp_min_size(interp);
  if(nlines < n_required) {
    fprintf(stderr, "Error: Not enough power spectrum data points for cubic spline; %lu data points < %lu required\n", nlines, n_required);
    throw RuntimeError();
  }

  gsl_interp_init(interp, log_k, log_P, nlines);

  k_min= exp(log_k[0]);
  k_max= exp(log_k[nlines-1]);
}

double PowerSpectrum::P(const double k) const
{
  if(k <= k_min || k > k_max)
    return 0.0;
  
  double log_Pk=
    gsl_interp_eval(interp, log_k, log_P, log(k), acc);
  
  return exp(log_Pk);
}

