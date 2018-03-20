#include <iostream>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <valarray>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "power_spectrum.h"
#include "particle.h"
#include "grid.h"
#include "mock_lognormal.h"


using namespace std;

inline double nbar_r(const double r)
{
  const double dr = r - 2118.308;
  return 1.64954511e-03*exp(-1.13301267e-03*dr + 2.62878572e-07*dr*dr);
}

//
// Step 1. Set power spectrum grid P(k)
//
void set_power_spectrum(PowerSpectrum const * const ps,
			const double boxsize,
			const double d,
			const double b,
			const double f,
			Grid* const grid)
{
  // Input:
  //     ps (PowerSpectrum): target P(k)
  //                         i.e. the power spectrum of the resulting mock
  //     boxsize (double):   size of the output mock on a side
  //     d (double):         growth factor D
  //     b (double):         bias
  //     f (double):         growth rate f = dln D(a)/dln a
  // Output:
  //     grid (Grid):        grid of P(k_vector)
  grid->clear();
  grid->boxsize= boxsize;
  const int nc= grid->nc;
  const int nckz= nc/2 + 1;

  fftw_complex* const pk= (fftw_complex*) grid->fx; 

  const double fac= 2.0*M_PI/boxsize;
  
  for(int ix=0; ix<nc; ++ix) {
   double kx= ix <= nc/2 ? fac*ix : fac*(ix - nc);
   if(2*ix == nc) continue;
   for(int iy=0; iy<nc; ++iy) {
    double ky= iy <= nc/2 ? fac*iy : fac*(iy - nc);
    if(2*iy == nc) continue;
    
    int iz0= (ix == 0 && iy == 0);

    for(int iz=iz0; iz<nc/2; ++iz) {
      double kz= fac*iz;
      
      double k= sqrt(kx*kx + ky*ky + kz*kz);
      double mu= kz/k;
      double fac= d*(b + f*mu*mu);
      double fac2= fac*fac;
      
      size_t index= (ix*nc + iy)*nckz + iz;
      pk[index][0]= fac2*ps->P(k);
      pk[index][1]= 0.0;
    }
   }
  }

  pk[0][0]= pk[0][1]= 0.0;

  grid->mode= fft_mode_k;
}

//
// Step 2: compute gaussian P(k)
//
void compute_gaussian_power_spectrum(Grid* const grid)
{
  // Input:
  //     grid (Grid): grid of lognormal P(k)
  // Output:
  //     grid (Grid): grid of Gaussian P_g(k)
  
  // Convert the P(k) grid to P_g(k) grid where P_g(k), 
  // `Gaussian power spectrum`, is the power spectrum of delta_g such that
  // the lognormal field have the expected P(k)
  // 1 + delta(x) ~ exp(delta_g(x))
  
  // P(k) => xi(r)
  grid->fft_inverse();

  //
  // Convert non-linear xi(r) to gaussian xi_g(r)
  // xi(r) => xi_g(r) = log(1 + xi(r))
  //
  double* xi= grid->fx;
  size_t nc= grid->nc;
  size_t ncz= 2*(nc/2 + 1);
  
  for(size_t ix=0; ix<nc; ++ix) {
    for(size_t iy=0; iy<nc; ++iy) {
      for(size_t iz=0; iz<nc; ++iz) {
	size_t index= (ix*nc + iy)*ncz + iz;
	assert(1.0 + xi[index] > 0.0);
	xi[index]= log(1.0 + xi[index]);
      }
    }
  }

  // xi_g(r) => P_g(r)
  grid->fft_forward();
}

//
// Step 3: generate random realisation of delta_g(k) based on P_g
//
void generate_delta_k(const unsigned long seed,
		      const bool fix_amplitude,
		      Grid* const grid) 
{
  // Convert P(k) grid to random delta(k) grid such that
  // <|delta(k)|^2> = P(k)
  //
  // input:  grid as P(k)
  // output: grid as delta(k)

  const int nc= grid->nc;
  const size_t nckz= nc/2 + 1;
  const double boxsize= grid->boxsize;
  const double vol= boxsize*boxsize*boxsize;
  fftw_complex* const fk= (fftw_complex*) grid->fx;
  
  // P(k) = 1/V <delta(k) delta^*(k)
  valarray<double> v_corr(1.0, nc);

  gsl_rng* rng= gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(rng, seed);

  size_t negative= 0;
  double P_min= 0.0;

  
  for(int ix=0; ix<nc; ++ix) {
   if(2*ix == nc) continue;
   //double kx= ix <= nc/2 ? ix : ix - nc;
   
   for(int iy=0; iy<nc; ++iy) {
    if(2*iy == nc) continue;
    //double ky= iy <= nc/2 ? iy : iy - nc;
    
    int iz0= (ix == 0 && iy == 0);
    for(int iz=iz0; iz<nc/2; ++iz) {
      //double kz= iz;
      //double k= sqrt(kx*kx + ky*ky + kz*kz);

      size_t index= (ix*nc + iy)*nckz + iz;

      double ampl=1.0;
      if(!fix_amplitude) {
	do
	  ampl = gsl_rng_uniform(rng);
	while(ampl == 0.0);
    
	ampl= -log(ampl);
      }

      double phase= gsl_rng_uniform(rng)*2*M_PI;
      
      double delta2= vol*fk[index][0];

      double delta_k_mag= 0.0;
      if(fk[index][0] < P_min)
	P_min= fk[index][0];
	  
      if( fk[index][0] > 0.0)
	delta_k_mag= sqrt(ampl*delta2);
      else
	negative++;

      fk[index][0]= delta_k_mag*cos(phase);
      fk[index][1]= delta_k_mag*sin(phase);
	  
    }
   }
  }

  gsl_rng_free(rng);

  fprintf(stderr, "P_min= %e\n", P_min);
  fprintf(stderr, "negative P(k): %zu\n", negative);

  //
  // reality condition delta(-k) = delta(k)^*
  //
  for(int ix=0; ix<nc; ++ix) {
    if(2*ix == nc) continue;
    int iix= ix == 0 ? 0 :  nc - ix;
    assert(0 <= iix && iix < nc);
    for(int iy=0; iy<nc; ++iy) {
      if(2*iy == nc) continue;
      int iiy= iy == 0 ? 0 : nc - iy;

      size_t index= (ix*nc + iy)*nckz;    // index of k
      size_t iindex= (iix*nc + iiy)*nckz; // index of -k

      fk[iindex][0]= fk[index][0];
      fk[iindex][1]= -fk[index][1];
    }
  }
}

//
// Convert delta_g(k) to delta(k)
//
void compute_lognormal_density(Grid* const grid)
{
  // input: delta_g(k)
  // output: 1 + delta_target \prop exp(delta_g(x))

  assert(grid->mode == fft_mode_k);
  grid->fft_inverse();
  
  assert(grid->mode == fft_mode_x);
  size_t nc= grid->nc;
  size_t ncz= 2*(nc/2 + 1);
  double* const delta= grid->fx;
  long double nsum= 0.0;
  long double n2sum= 0.0;
  double nmax= 0.0;
  double dmax= 0.0;
  long double d2sum= 0.0;
  
  for(size_t ix=0; ix<nc; ++ix) {
   for(size_t iy=0; iy<nc; ++iy) {
    for(size_t iz=0; iz<nc; ++iz) {
      size_t index= (ix*nc + iy)*ncz + iz;
	  
      double d= delta[index];
      double n= exp(d); // lognormal density
      delta[index] = n;
      if(d > dmax) dmax= d;
      d2sum += d*d;
      
      // lognormal density (unnormalised)
      nsum += n;
      n2sum += n*n;
      if(n > nmax) nmax= n;
    }
   }
  }
  
  cerr << "dmax= " << dmax << endl;
  cerr << "d rms= " << sqrt(d2sum/(nc*nc*nc)) << endl;
  // normalization
  const double nbar= nsum/(nc*nc*nc);
  cerr << "nbar= " << nbar << endl;
  cerr << "n2sum= " << n2sum << endl;

  // Normalise density to n(x)/nbar = 1 + delta(x)
  for(size_t ix=0; ix<nc; ++ix) {
    for(size_t iy=0; iy<nc; ++iy) {
      for(size_t iz=0; iz<nc; ++iz) {
	size_t index= (ix*nc + iy)*ncz + iz;

	delta[index] /= nbar;
      }
    }
  }

  cerr << "Lognormal rms: " << sqrt(n2sum/(nbar*nbar*nc*nc*nc) - 1.0) << endl;
  //cerr << "nmax " << nmax/nbar << endl;
  //return nmax/nbar; // max(1 + delta)
}

//
// step 4: generate particles
//
void generate_particles(Grid const * const grid,
			const unsigned long seed, const size_t np,
			vector<Particle>* const v)
{
  // Samples random number of particles from each grid (Poisson sampling)
  //
  // Input:
  //     seed (unsigned long): random seed of poisson sampling
  //     
  //     grid of 1 + delta(x)
  // Output:
  //     v (vector<Particle>): mock particles are added to this vector
  assert(grid->mode == fft_mode_x);

  gsl_rng* rng= gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(rng, seed);

  const size_t nc= grid->nc;
  const size_t ncz= 2*(nc/2 + 1);

  // mean number density over all box
  // number of particles per cell
  const double num_bar= static_cast<double>(np)/(nc*nc*nc);
  const double dx= grid->boxsize/nc;

  double const * const n= grid->fx;

  // Allocate memory for output particles
  size_t n_alloc= static_cast<size_t>(np + 10.0*sqrt(static_cast<double>(np)));
  v->reserve(v->size() + n_alloc);

  Particle p;
  //int ix[3];

  size_t count= 0;
  size_t count_negative= 0;
  long double num_total_mean= 0;

  for(size_t ix=0; ix<nc; ++ix) {
    for(size_t iy=0; iy<nc; ++iy) {
      for(size_t iz=0; iz<nc; ++iz) {
	size_t index= (ix*nc + iy)*ncz + iz;

	// mean number of particles in the cell
	double num_grid= n[index]*num_bar;
	num_total_mean += num_grid;

	// random realisation of number of particles in cell
	int num= gsl_ran_poisson(rng, num_grid);

	if(n[index] < 0.0)
	  count_negative += num;
	
	for(int i=0; i<num; ++i) {
	  p.x[0]= (ix + gsl_rng_uniform(rng))*dx;
	  p.x[1]= (iy + gsl_rng_uniform(rng))*dx;
	  p.x[2]= (iz + gsl_rng_uniform(rng))*dx;
	  v->push_back(p);
	}
	count += num;
      }
    }
  }
    
  gsl_rng_free(rng);

  fprintf(stderr, "np= %zu, count= %zu\n", np, count);
  fprintf(stderr, "num_total_mean= %Lf\n", num_total_mean);
  if(count_negative > 0)
    fprintf(stderr, "negative density= %zu / %zu\n", count_negative, count);
}

//
// step 4a: generate particles in a octant shell
//
void generate_particles_octant(Grid const * const grid,
			       const unsigned long seed,
			       const double nbar,
			       const double r_min, const double r_max,
			       vector<Particle>* const v)
{
  // Samples random number of particles from each grid (Poisson sampling)
  //
  // Input:
  //     seed (unsigned long): random seed of poisson sampling
  //     
  //     grid of 1 + delta(x)
  // Output:
  //     v (vector<Particle>): mock particles are added to this vector
  assert(grid->mode == fft_mode_x);

  cerr << "start printing mock in an octant shell\n";
  
  gsl_rng* rng= gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(rng, seed);

  const size_t nc= grid->nc;
  const size_t ncz= 2*(nc/2 + 1);

  // mean number density over all box
  // number of particles per cell
  const double dx= grid->boxsize/nc;
  const double num_bar= nbar*dx*dx*dx; //static_cast<double>(np)/(nc*nc*nc);
  
  double const * const n= grid->fx;


  Particle p;

  size_t count= 0;
  size_t count_negative= 0;

  const double ir_max = r_max/dx;
  const double ir2_max = (r_max/dx)*(r_max/dx);
  const double ir2_min = (r_min/dx)*(r_min/dx);
  const double r2_min = r_min*r_min;
  const double r2_max = r_max*r_max;

  const double vol = 0.5*M_PI/3.0*(r_max*r_max*r_max - r_min*r_min*r_min);
  const double np= nbar*vol;
  
  // Allocate memory for output particles
  size_t n_alloc= static_cast<size_t>(np + 10.0*sqrt(np));
  v->reserve(v->size() + n_alloc);


  // Loop over all grids
  for(size_t ix=0; ix<nc; ++ix) {
    if(ix > ir_max) break; // All the following grids are beyond r_max
    for(size_t iy=0; iy<nc; ++iy) {
      if(ix*ix + iy*iy > ir2_max) break;
      for(size_t iz=0; iz<nc; ++iz) {
	if(ix*ix + iy*iy + iz*iz > ir2_max)
	  break;
	else if((ix + 1)*(ix + 1) + (iy + 1)*(iy + 1) + (iz + 1)*(iz + 1)
		< ir2_min)
	  continue;
	
	size_t index= (ix*nc + iy)*ncz + iz;

	// mean number of particles in the cell
	double num_grid= n[index]*num_bar;

	// random realisation of number of particles in cell
	int num= gsl_ran_poisson(rng, num_grid);

	if(n[index] < 0.0)
	  count_negative += num;
	
	for(int i=0; i<num; ++i) {
	  p.x[0]= (ix + gsl_rng_uniform(rng))*dx;
	  p.x[1]= (iy + gsl_rng_uniform(rng))*dx;
	  p.x[2]= (iz + gsl_rng_uniform(rng))*dx;
	  double r2=  p.x[0]*p.x[0] + p.x[1]* p.x[1] + p.x[2]* p.x[2];
	  if(r2_min <= r2 && r2 < r2_max)
	    v->push_back(p);
	}
	count += num;
      }
    }
  }
    
  gsl_rng_free(rng);

  fprintf(stderr, "np= %.1lf, count= %zu\n", np, count);
  // this 'count' is before r2_min <= r2 < r2_max cut
  if(count_negative > 0)
    fprintf(stderr, "negative density= %zu / %zu\n", count_negative, count);
}

void generate_randoms_octant(const size_t nc, const double boxsize,
			     const unsigned long seed,
			     const double rand_ratio,
			     const double r_min, const double r_max,
			     vector<Particle>* const v)
{
  // Samples random number of particles from each grid (Poisson sampling)
  //
  // Args:
  //     seed (unsigned long): random seed of poisson sampling
  //     
  //     grid of 1 + delta(x)
  // Output:
  //     v (vector<Particle>): random particles are added to this vector

  gsl_rng* rng= gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(rng, seed);

  // mean number density over all box
  // number of particles per cell
  const double dx= boxsize/nc;
  const double vol= dx*dx*dx;

  
  Particle p;

  size_t count= 0;

  const double ir_max = r_max/dx;
  const double ir2_min = (r_min/dx)*(r_min/dx);
  const double ir2_max = (r_max/dx)*(r_max/dx);
  const double r2_min = r_min*r_min;
  const double r2_max = r_max*r_max;

  //const double vol = 0.5*M_PI/3.0*(r_max*r_max*r_max - r_min*r_min*r_min);
  //const double np= nbar*vol;
  
  // Allocate memory for output particles
  //size_t n_alloc= static_cast<size_t>(np + 10.0*sqrt(np));
  //v->reserve(v->size() + n_alloc);


  // Loop over all grids
  for(size_t ix=0; ix<nc; ++ix) {
    if(ix > ir_max) break; // All the following grids are beyond r_max
    for(size_t iy=0; iy<nc; ++iy) {
      if(ix*ix + iy*iy > ir2_max) break;
      for(size_t iz=0; iz<nc; ++iz) {
	if(ix*ix + iy*iy + iz*iz > ir2_max)
	  break;
	else if((ix + 1)*(ix + 1) + (iy + 1)*(iy + 1) + (iz + 1)*(iz + 1)
		< ir2_min) // grid does not intersect with r_min < r
	  continue;

	const double r_corner
	  = dx*sqrt(static_cast<double>(ix*ix +iy*iy + iz*iz));
	const double nbar_max = nbar_r(r_corner);
	const double num_bar = rand_ratio*nbar_max*vol;
	
	// random realisation of number of particles in cell
	int num= gsl_ran_poisson(rng, num_bar);

	for(int i=0; i<num; ++i) {
	  p.x[0]= (ix + gsl_rng_uniform(rng))*dx;
	  p.x[1]= (iy + gsl_rng_uniform(rng))*dx;
	  p.x[2]= (iz + gsl_rng_uniform(rng))*dx;
	  double r2=  p.x[0]*p.x[0] + p.x[1]* p.x[1] + p.x[2]* p.x[2];
	  if(r2_min <= r2 && r2 < r2_max) {
	    double r= sqrt(r2);
	    double prob= gsl_rng_uniform(rng);
	    double nbar= nbar_r(r);
	    if(prob < nbar/nbar_max) {	    
	      v->push_back(p);
	    }
	  }
	}
	count += num;
      }
    }
  }
    
  gsl_rng_free(rng);

  //fprintf(stderr, "count= %zu %zu\n", count, v->size());
  // this 'count' is before r2_min <= r2 < r2_max cut
}
