/*
 * Copyright (c) 2010 Centre National de la Recherche Scientifique.
 * written by Nathanael Schaeffer (CNRS, LGIT, Grenoble, France).
 * 
 * nathanael.schaeffer@ujf-grenoble.fr
 * 
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software. You can use,
 * modify and/or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 * 
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * 
 */

/** \internal \file sht_legendre.c
 \brief Compute legendre polynomials and associated functions.
 
 - Normalization of the associated functions for various spherical harmonics conventions are possible, with or without Condon-Shortley phase.
 - When computing the derivatives (with respect to colatitude theta), there are no singularities.
 - written by Nathanael Schaeffer / LGIT,CNRS, with some ideas and code from the GSL 1.13 and Numerical Recipies.

 The associated legendre functions are computed using the following recurrence formula :
 \f[  Y_l^m(x) = a_l^m \, x \  Y_{l-1}^m(x) + b_l^m \ Y_{l-2}^m(x)  \f]
 with starting values :
 \f[  Y_m^m(x) = a_m^m \ (1-x^2)^{m/2}  \f]
 and
 \f[  Y_{m+1}^m(x) = a_{m+1}^m \ Y_m^m(x)  \f]
 where \f$ a_l^m \f$ and \f$ b_l^m \f$ are coefficients that depend on the normalization, and which are
 precomputed once for all by \ref legendre_precomp.
*/

// range checking for legendre functions. Comment out for slightly faster legendre function generation.
//#define LEG_RANGE_CHECK

#if SHT_VERBOSE > 1
  #define LEG_RANGE_CHECK
#endif

/// computes sin(t)^n from cos(t). ie returns (1-x^2)^(n/2), with x = cos(t)
inline double sint_pow_n(double cost, int n)
{
	double val = 1.0;
	double s2 = (1.-cost)*(1.+cost);		// sin(t)^2 = 1 - cos(t)^2
	if (n&1) val *= sqrt(s2);	// = sin(t)
	n >>= 1;
	while(n>0) {
		if (n&1) val *= s2;
		n >>= 1;
		s2 *= s2;
	}
	return val;		// = sint(t)^n
}

/// Returns the value of a legendre polynomial of degree l and order im*MRES, noramalized for spherical harmonics, using recurrence.
/// Requires a previous call to \ref legendre_precomp().
/// Output compatible with the GSL function gsl_sf_legendre_sphPlm(l, m, x)
double legendre_sphPlm(const int l, const int im, double x)
{
	double *al;
	long int i,m;
	double ymm, ymmp1;

	m = im*MRES;
#ifdef LEG_RANGE_CHECK
	if ( (l>LMAX) || (l<m) || (im>MMAX) ) shtns_runerr("argument out of range in legendre_sphPlm");
#endif

	al = alm[im];
	ymm = al[0] * sint_pow_n(x, m);		// l=m
	if (l==m) return ((double) ymm);

	ymmp1 = ymm * al[1] * x;				// l=m+1
	al+=2;
	if (l == m+1) return ((double) ymmp1);
	
	for (i=m+2; i<l; i+=2) {
		ymm   = al[1]*x*ymmp1 + al[0]*ymm;
		ymmp1 = al[3]*x*ymm + al[2]*ymmp1;
		al+=4;
	}
	if (i==l) {
		ymmp1 = al[1]*x*ymmp1 + al[0]*ymm;
	}
	return ((double) ymmp1);
}

/// Compute values of legendre polynomials noramalized for spherical harmonics,
/// for a range of l=m..lmax, at given m and x, using recurrence.
/// Requires a previous call to \ref legendre_precomp().
/// Output compatible with the GSL function gsl_sf_legendre_sphPlm_array(lmax, m, x, yl)
/// \param lmax maximum degree computed, \param im = m/MRES with m the SH order, \param x argument, x=cos(theta).
/// \param[out] yl is a double array of size (lmax-m+1) filled with the values.
void legendre_sphPlm_array(const int lmax, const int im, double x, double *yl)
{
	double *al;
	long int l,m;
	double ymm, ymmp1;

	m = im*MRES;
#ifdef LEG_RANGE_CHECK
	if ((lmax > LMAX)||(lmax < m)||(im>MMAX)) shtns_runerr("argument out of range in legendre_sphPlm_array");
#endif

	al = alm[im];
	ymm = al[0] * sint_pow_n(x, m);	// l=m
	yl[0] = ymm;
	if (lmax==m) return;		// done.

	ymmp1 = ymm * al[1] * x;		// l=m+1
	yl[1] = ymmp1;
	al+=2;
	if (lmax==m+1) return;		// done.

	yl -= m;			// shift pointer
	for (l=m+2; l<lmax; l+=2) {
		ymm   = al[1]*x*ymmp1 + al[0]*ymm;
		ymmp1 = al[3]*x*ymm   + al[2]*ymmp1;
		yl[l] = ymm;
		yl[l+1] = ymmp1;
		al+=4;
	}
	if (l==lmax) {
		yl[l] = al[1]*x*ymmp1 + al[0]*ymm;
	}
	return;
}


/// Compute values of a legendre polynomial normalized for spherical harmonics derivatives, for a range of l=m..lmax, using recurrence.
/// Requires a previous call to \ref legendre_precomp(). Output is not directly compatible with GSL :
/// - if m=0 : returns ylm(x)  and  d(ylm)/d(theta) = -sin(theta)*d(ylm)/dx
/// - if m>0 : returns ylm(x)/sin(theta)  and  d(ylm)/d(theta).
/// This way, there are no singularities, everything is well defined for x=[-1,1], for any m.
/// \param lmax maximum degree computed, \param im = m/MRES with m the SH order, \param x argument, x=cos(theta).
/// \param sint = sqrt(1-x^2) to avoid recomputation of sqrt.
/// \param[out] yl is a double array of size (lmax-m+1) filled with the values (divided by sin(theta) if m>0)
/// \param[out] dyl is a double array of size (lmax-m+1) filled with the theta-derivatives.
void legendre_sphPlm_deriv_array(const int lmax, const int im, const double x, const double sint, double *yl, double *dyl)
{
	double *al;
	long int l,m;
	double st, y0, y1, dy0, dy1;

	m = im*MRES;
#ifdef LEG_RANGE_CHECK
	if ((lmax > LMAX)||(lmax < m)||(im>MMAX)) shtns_runerr("argument out of range in legendre_sphPlm_deriv_array");
#endif
	al = alm[im];

	st = sint;
	y0 = al[0];
	dy0 = 0.0;
	if (im>0) {		// m > 0
		l = m-1;			// compute  sin(theta)^(m-1)
		do {
			if (l&1) y0 *= st;
			st *= st;
		} while(l >>= 1);
		dy0 = x*m*y0;
		st = sint*sint;		// st = sin(theta)^2 is used in the recurrence for m>0
	}
	yl[0] = y0;		// l=m
	dyl[0] = dy0;
	if (lmax==m) return;		// done.

  #ifndef _GCC_VEC_
	y1 = al[1] * x * y0;
	dy1 = al[1]*( x*dy0 - st*y0 );
	yl[1] = y1;		// l=m+1
	dyl[1] = dy1;
	if (lmax==m+1) return;		// done.
	al+=2;

	yl -= m;	dyl -= m;			// shift pointers
	for (l=m+2; l<lmax; l+=2) {
		y0 = al[1]*x*y1 + al[0]*y0;
		dy0 = al[1]*(x*dy1 - y1*st) + al[0]*dy0;
		yl[l] = y0;		dyl[l] = dy0;
		y1 = al[3]*x*y0 + al[2]*y1;
		dy1 = al[3]*(x*dy0 - y0*st) + al[2]*dy1;
		yl[l+1] = y1;		dyl[l+1] = dy1;
		al+=4;
	}
	if (l==lmax) {
		yl[l] = al[1]*x*y1 + al[0]*y0;
		dyl[l] = al[1]*(x*dy1 - y1*st) + al[0]*dy0;
	}
  #else
  	v2d vy1 = vset(dy0, y0);
  	v2d vx = vdup(x);
  	v2d vst = vset(st, 0.0);

  	v2d vty = vxchg(vy1);		// swap
  	vty = vdup(al[1]) * ( vx*vy1 - vst*vty );

	yl[1] = vlo_to_dbl(vty);		// l=m+1
	dyl[1] = vhi_to_dbl(vty);
	al+=2;

	yl -= m;	dyl -= m;			// shift pointers
	for (l=m+2; l<=lmax; l++) {		// vectorized loop : 4* and 2+ (instead of 7* and 3+)
		v2d vy0 = vy1;
		vy1 = vty;
		vty = vxchg(vty);		// swap
		vty = vdup(al[0])*vy0  +  vdup(al[1]) * (vx*vy1 - vst*vty);
		yl[l] = vlo_to_dbl(vty);		dyl[l] = vhi_to_dbl(vty);
		al+=2;
	}
  #endif
}

/// Precompute constants for the recursive generation of Legendre associated functions, with given normalization.
/// this function is called by \ref shtns_set_size, and assumes up-to-date values in \ref shtns.
/// For the same conventions as GSL, use \c legendre_precomp(sht_orthonormal,1);
/// \param[in] norm : normalization of the associated legendre functions (\ref shtns_norm).
/// \param[in] with_cs_phase : Condon-Shortley phase (-1)^m is included (1) or not (0)
/// \param[in] mpos_renorm : Optional renormalization for m>0.
///  1.0 (no renormalization) is the "complex" convention, while 0.5 leads to the "real" convention (with FFTW).
void legendre_precomp(enum shtns_norm norm, int with_cs_phase, double mpos_renorm)
{
	double *dl0;
	long int im, m, l, lm;
	double t1, t2;

#if SHT_VERBOSE > 1
	printf("        > Condon-Shortley phase = %d, normalization = %d\n", with_cs_phase, norm);
#endif

	if (with_cs_phase != 0) with_cs_phase = 1;		// force to 1 if !=0

	im = ((MMAX+2)>>1)*2;		// alloc memory for arrays.
	alm = (double **) malloc( im * sizeof(double *) + (2*NLM)*sizeof(double) );
	al0 = (double *) (alm + im);
	dlm = (double **) malloc( im * sizeof(double *) + (4*NLM)*sizeof(double) );
	dl0 = (double *) (dlm + im);

/// - Precompute the factors alm and blm of the recurrence relation :
  if (norm == sht_schmidt) {
	for (im=0, lm=0; im<=MMAX; im++) {		/// <b> For Schmidt semi-normalized </b>
		m = im*MRES;
		alm[im] = al0 + lm;
		t2 = sqrt(2*m+1);
		al0[lm] = 1.0/t2;		/// starting value will be divided by \f$ \sqrt{2m+1} \f$ 
		al0[lm+1] = t2;		// l=m+1
		lm+=2;
		for (l=m+2; l<=LMAX; l++) {
			t1 = sqrt((l+m)*(l-m));
			al0[lm+1] = (2*l-1)/t1;		/// \f[  a_l^m = \frac{2l-1}{\sqrt{(l+m)(l-m)}}  \f]
			al0[lm] = - t2/t1;			/// \f[  b_l^m = -\sqrt{\frac{(l-1+m)(l-1-m)}{(l+m)(l-m)}}  \f]
			t2 = t1;	lm+=2;
		}
	}
  } else {
	for (im=0, lm=0; im<=MMAX; im++) {		/// <b> For orthonormal or 4pi-normalized </b>
		m = im*MRES;
		alm[im] = al0 + lm;
		t2 = 2*m+1;
		al0[lm] = 1.0;		// will contain the starting value.
		al0[lm+1] = sqrt(2*m+3);		// l=m+1
		lm+=2;
		for (l=m+2; l<=LMAX; l++) {
			t1 = (l+m)*(l-m);
			al0[lm+1] = sqrt(((2*l+1)*(2*l-1))/t1);			/// \f[  a_l^m = \sqrt{\frac{(2l+1)(2l-1)}{(l+m)(l-m)}}  \f]
			al0[lm] = - sqrt(((2*l+1)*t2)/((2*l-3)*t1));	/// \f[  b_l^m = -\sqrt{\frac{2l+1}{2l-3}\,\frac{(l-1+m)(l-1-m)}{(l+m)(l-m)}}  \f]
			t2 = t1;	lm+=2;
		}
	}
  }

/// - Compute and store the prefactor (independant of x) of the starting value for the recurrence :
/// \f[  Y_m^m(x) = Y_0^0 \ \sqrt{ \prod_{k=1}^{m} \frac{2k+1}{2k} } \ \ (-1)^m \ (1-x^2)^{m/2}  \f]
	if ((norm == sht_fourpi)||(norm == sht_schmidt)) {
		t1 = 1.0;
		al0[0] = t1;		/// \f$ Y_0^0 = 1 \f$  for Schmidt or 4pi-normalized 
	} else {
		t1 = 0.25L/M_PIl;
		al0[0] = sqrt(t1);		/// \f$ Y_0^0 = 1/\sqrt{4\pi} \f$ for orthonormal
	}
	t1 *= mpos_renorm;		// renormalization for m>0
	for (im=1, m=0; im<=MMAX; im++) {
		while(m<im*MRES) {
			m++;
			t1 *= ((double)m + 0.5)/m;	// t1 *= (m+0.5)/m;
		}
		t2 = sqrt(t1);
		if ( m & with_cs_phase ) t2 = -t2;		/// optional \f$ (-1)^m \f$ Condon-Shortley phase.
		alm[im][0] *= t2;
	}

/// - Compute and store coefficients for computation of derivative also.
	for (im=1, lm=0; im<=MMAX; im++) {
		m = im*MRES;
		t2 = 1.0/m;
		dlm[im] = dl0 + lm;
		double *al = alm[im];
			dl0[lm] = al[0];	// ymm initial value.
		for (l=m+1; l<=LMAX; l++) {
			t1 = (l+m)*(l-m);
			if (norm != sht_schmidt) t1 *= (2*l+1.)/(2*l-1.);
			dl0[lm] = al[0];	dl0[lm+1] = al[1];		// ylm recursion coeffs
			dl0[lm+2] = -t2*sqrt(t1);
			dl0[lm+3] = t2*l;
			lm+=4;	al+=2;
		}
	}
}

/// returns the value of the Legendre Polynomial of degree l.
/// l is arbitrary, a direct recurrence relation is used, and a previous call to legendre_precomp() is not required.
double legendre_Pl(const int l, double x)
{
	long int i;
	double p, p0, p1;

	if ((l==0)||(x==1.0)) return ( 1. );
	if (x==-1.0) return ( (l&1) ? -1. : 1. );

	p0 = 1.0;		/// \f$  P_0 = 1  \f$
	p1 = x;			/// \f$  P_1 = x  \f$
	for (i=2; i<=l; i++) {		 /// recurrence : \f[  l P_l = (2l-1) x P_{l-1} - (l-1) P_{l-2}	 \f] (works ok up to l=100000)
		p = (x*(2*i-1)*p1 - (i-1)*p0)/i;
		p0 = p1;
		p1 = p;
	}
	return ((double) p1);
}


/// Generates the abscissa and weights for a Gauss-Legendre quadrature.
/// Newton method from initial Guess to find the zeros of the Legendre Polynome
/// \param x = abscissa, \param w = weights, \param n points.
/// \note Reference:  Numerical Recipes, Cornell press.
void gauss_nodes(long double *x, long double *w, int n)
{
	long int i,l,m;
	long double z, z1, p1, p2, p3, pp;
	long double eps;

	if (sizeof(eps) > 8) {
		eps = 1.1e-19;		// desired precision, minimum = 1.0842e-19 (long double i387)
	} else {
		eps = 2.3e-16;		// desired precision, minimum = 2.2204e-16 (double)
	}

	m = (n+1)/2;
	for (i=1;i<=m;i++) {
		z = cos(M_PI*(i-0.25)/(n+0.5));	// initial guess
		do {
			p1 = z;		// P_1
			p2 = 1.0;	// P_0
			for(l=2;l<=n;l++) {		 // recurrence : l P_l = (2l-1) z P_{l-1} - (l-1) P_{l-2}	(works ok up to l=100000)
				p3 = p2;
				p2 = p1;
				p1 = ((2*l-1)*z*p2 - (l-1)*p3)/l;		// The Legendre polynomial...
			}
			pp = - n*(z*p1-p2)/((1.-z)*(1.+z));	// ... and its derivative.
			z1 = z;
			z = z-p1/pp;
		} while ( fabsl(z-z1) > eps );
		x[i-1] = z;		// Build up the abscissas.
		w[i-1] = 2.0/(((1.-z)*(1.+z))*(pp*pp));		// Build up the weights.
		x[n-i] = -z;
		w[n-i] = w[i-1];
	}

// as we started with initial guesses, we should check if the gauss points are actually unique.
	for (i=m-1; i>0; i--) {
		if (((double) x[i]) == ((double) x[i-1])) shtns_runerr("bad gauss points");
	}

#if SHT_VERBOSE > 1
// test integral to compute :
	z = 0;
	for (i=0;i<m;i++) {
		z += w[i]*x[i]*x[i];
	}
	printf("          Gauss quadrature for 3/2.x^2 = %Lg (should be 1.0) error = %Lg\n",z*3.,z*3.-1.0);
#endif
}
