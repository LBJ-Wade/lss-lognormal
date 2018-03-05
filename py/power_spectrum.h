#ifndef POWER_SPECTRUM_H
#define POWER_SPECTRUM_H 1

#include <vector>
#include <gsl/gsl_spline.h>

class PowerSpectrum {
 public:
  PowerSpectrum(const char filename[]);
  ~PowerSpectrum();
  void load(const char filename[]);
  double P(const double k) const;
  
  double *log_k, *log_P;
  double k_min, k_max;

 private:
  gsl_interp* interp;
  gsl_interp_accel* acc;
};

#endif
