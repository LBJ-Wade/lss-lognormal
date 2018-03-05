#ifndef MOCK_LOGNORMAL_H
#define MOCK_LOGNORMAL_H 1

#include <vector>
#include "power_spectrum.h"
#include "particle.h"
#include "grid.h"

void set_power_spectrum(PowerSpectrum const * const ps,
			const double boxsize,
			const double d,
			const double b,
			const double f,
			Grid* const grid);

void compute_gaussian_power_spectrum(Grid* const grid);

void compute_lognormal_density(Grid* const grid);

void generate_delta_k(const unsigned long seed,
		      const bool fix_amplitude,
		      Grid* const grid);

void generate_particles(Grid const * const grid,
			const unsigned long seed, const size_t np,
			std::vector<Particle>* const v);

#endif
