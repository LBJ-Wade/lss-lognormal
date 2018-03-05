#ifndef GRID_H
#define GRID_H 1

#include <fftw3.h>

enum FFTMode {fft_mode_unknown, fft_mode_x, fft_mode_k};

class Grid {
 public:
  Grid(const int nc_);
  ~Grid();

  void fft_forward();
  void fft_inverse();
  void clear();

  double* fx;
  fftw_complex* fk;
  const int nc;
  double boxsize;
  FFTMode mode;
  double x0_box[3];


 private:
  const int ncz;
  fftw_plan plan_forward;
  fftw_plan plan_inverse;
};

#endif

