
void generate_delta_k(const unsigned long seed,
		      const bool fix_amplitude,
		      Grid* const grid) 
{
  // Convert P(k) grid to random delta(k) grid such that
  // <|delta(k)|^2> = P(k)
  //
  // input grid: P(k)
  // output grid: delta_k(k)

  assert(grid->mode == fft_mode_k);
  
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
   double kx= ix <= nc/2 ? ix : ix - nc;
   
   for(int iy=0; iy<nc; ++iy) {
    if(2*iy == nc) continue;
    double ky= iy <= nc/2 ? iy : iy - nc;
    
    int iz0= (ix == 0 && iy == 0);
    for(int iz=iz0; iz<nc/2; ++iz) {
      double kz= iz;
      double k= sqrt(kx*kx + ky*ky + kz*kz);

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
  // reality condition
  //
  for(int ix=0; ix<nc; ++ix) {
    if(2*ix == nc) continue;
    int iix= ix == 0 ? 0 :  nc - ix;
    assert(0 <= iix && iix < nc);
    for(int iy=0; iy<nc; ++iy) {
      if(2*iy == nc) continue;
      int iiy= iy == 0 ? 0 : nc - iy;

      size_t index= (ix*nc + iy)*nckz;
      size_t iindex= (iix*nc + iiy)*nckz;

      fk[iindex][0]= fk[index][0];
      fk[iindex][1]= -fk[index][1];
    }
  }
}
