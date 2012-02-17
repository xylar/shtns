/*
 * Copyright (c) 2010-2012 Centre National de la Recherche Scientifique.
 * written by Nathanael Schaeffer (CNRS, ISTerre, Grenoble, France).
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
 - written by Nathanael Schaeffer / CNRS, with some ideas and code from the GSL 1.13 and Numerical Recipies.

 The associated legendre functions are computed using the following recurrence formula :
 \f[  Y_l^m(x) = a_l^m \, x \  Y_{l-1}^m(x) + b_l^m \ Y_{l-2}^m(x)  \f]
 with starting values :
 \f[  Y_m^m(x) = a_m^m \ (1-x^2)^{m/2}  \f]
 and
 \f[  Y_{m+1}^m(x) = a_{m+1}^m \ Y_m^m(x)  \f]
 where \f$ a_l^m \f$ and \f$ b_l^m \f$ are coefficients that depend on the normalization, and which are
 precomputed once for all by \ref legendre_precomp.
*/

#if SHT_VERBOSE > 1
  #define LEG_RANGE_CHECK
#endif

/// \internal high precision version of \ref a_sint_pow_n
long double a_sint_pow_n_hp(long double val, long double cost, long int n)
{
	long double s2 = (1.-cost)*(1.+cost);		// sin(t)^2 = 1 - cos(t)^2
	long int k = n >> 7;
	if (sizeof(s2) > 8) k = 0;		// enough accuracy, we do not bother.

	if (n&1) val *= sqrtl(s2);	// = sin(t)
	do {
		if (n&2) val *= s2;
		n >>= 1;
		s2 *= s2;
	} while(n > k);
	n >>= 1;
	while(n > 0) {		// take care of very large power n
		n--;
		val *= s2;		
	}
	return val;		// = sint(t)^n
}

/// \internal computes val.sin(t)^n from cos(t). ie returns val.(1-x^2)^(n/2), with x = cos(t)
/// works with very large values of n (up to 2700).
double a_sint_pow_n(double val, double cost, int n)
{
	double s2 = (1.-cost)*(1.+cost);		// sin(t)^2 = 1 - cos(t)^2
	int k = n >> 7;

	if (n&1) val *= sqrt(s2);	// = sin(t)
	do {
		if (n&2) val *= s2;
		n >>= 1;
		s2 *= s2;
	} while(n > k);
	n >>= 1;
	while(n > 0) {		// take care of very large power n
		n--;
		val *= s2;		
	}
	return val;		// = sint(t)^n
}


/// \internal Returns the value of a legendre polynomial of degree l and order im*MRES, noramalized for spherical harmonics, using recurrence.
/// Requires a previous call to \ref legendre_precomp().
/// Output compatible with the GSL function gsl_sf_legendre_sphPlm(l, m, x)
double legendre_sphPlm(shtns_cfg shtns, const int l, const int im, double x)
{
	double *al;
	long int i,m;
	long double ymm, ymmp1;

	m = im*MRES;
#ifdef LEG_RANGE_CHECK
	if ( (l>LMAX) || (l<m) || (im>MMAX) ) shtns_runerr("argument out of range in legendre_sphPlm");
#endif

	al = shtns->alm[im];
	ymm = al[0];		// l=m
	if (m>0) ymm *= SHT_LEG_SCALEF;
	ymmp1 = ymm;
	if (l==m) goto done;

	ymmp1 *= al[1] * x;				// l=m+1
	al+=2;
	if (l == m+1) goto done;

	for (i=m+2; i<l; i+=2) {
		ymm   = al[1]*x*ymmp1 + al[0]*ymm;
		ymmp1 = al[3]*x*ymm + al[2]*ymmp1;
		al+=4;
	}
	if (i==l) {
		ymmp1 = al[1]*x*ymmp1 + al[0]*ymm;
	}
done:
	if (m>0) ymmp1 *= a_sint_pow_n_hp(1.0/SHT_LEG_SCALEF, x, m);
	return ((double) ymmp1);
}



/// \internal high precision version of \ref legendre_sphPlm_array
void legendre_sphPlm_array_hp(shtns_cfg shtns, const int lmax, const int im, const double cost, double *yl)
{
	double *al;
	long int l,m;
	int rescale = 0;		// flag for rescale.
	long double ymm, ymmp1, x;

	m = im*MRES;
#ifdef LEG_RANGE_CHECK
	if ((lmax > LMAX)||(lmax < m)||(im>MMAX)) shtns_runerr("argument out of range in legendre_sphPlm_array");
#endif

	x = cost;
	al = shtns->alm[im];
	yl -= m;			// shift pointer
	ymm = al[0];	// l=m
	if (m>0) {
		if ((lmax <= SHT_L_RESCALE) || (sizeof(ymm) > 8)) {
			ymm = a_sint_pow_n_hp(ymm, x, m);
		} else {
			rescale = 1;
			ymm *= SHT_LEG_SCALEF;
		}
	}
	yl[m] = ymm;
	if (lmax==m) goto done;		// done.

	ymmp1 = ymm * al[1] * x;		// l=m+1
	yl[m+1] = ymmp1;
	al+=2;
	if (lmax==m+1) goto done;		// done.

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
done:
	if (rescale != 0) {
		ymm = a_sint_pow_n_hp(1.0/SHT_LEG_SCALEF, x, m);
		for (l=m; l<=lmax; l++) {		// rescale.
			yl[l] *= ymm;
		}
	}
	return;
}

/// \internal Compute values of legendre polynomials noramalized for spherical harmonics,
/// for a range of l=m..lmax, at given m and x, using recurrence.
/// Requires a previous call to \ref legendre_precomp().
/// Output compatible with the GSL function gsl_sf_legendre_sphPlm_array(lmax, m, x, yl)
/// \param lmax maximum degree computed, \param im = m/MRES with m the SH order, \param x argument, x=cos(theta).
/// \param[out] yl is a double array of size (lmax-m+1) filled with the values.
void legendre_sphPlm_array(shtns_cfg shtns, const int lmax, const int im, const double cost, double *yl)
{
	double *al;
	int l,m;
	double ymm, ymmp1, x;

	if ((im>0) && (lmax > SHT_L_RESCALE)) {		// switch to high precision code.
		legendre_sphPlm_array_hp(shtns, lmax, im, cost, yl);
		return;
	}

	m = im*MRES;
#ifdef LEG_RANGE_CHECK
	if ((lmax > LMAX)||(lmax < m)||(im>MMAX)) shtns_runerr("argument out of range in legendre_sphPlm_array");
#endif

	x = cost;
	al = shtns->alm[im];
	yl -= m;			// shift pointer
	ymm = al[0];	// l=m
	if (m>0) ymm = a_sint_pow_n(ymm, x, m);
	yl[m] = ymm;
	if (lmax==m) return;

	ymmp1 = ymm * al[1] * x;		// l=m+1
	yl[m+1] = ymmp1;
	al+=2;
	if (lmax==m+1) return;

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
}

/// \internal high precision version of \ref legendre_sphPlm_deriv_array
void legendre_sphPlm_deriv_array_hp(shtns_cfg shtns, const int lmax, const int im, const double cost, const double sint, double *yl, double *dyl)
{
	double *al;
	long int l,m;
	int rescale = 0;		// flag for rescale.
	long double x, st, y0, y1, dy0, dy1;

	x = cost;
	m = im*MRES;
#ifdef LEG_RANGE_CHECK
	if ((lmax > LMAX)||(lmax < m)||(im>MMAX)) shtns_runerr("argument out of range in legendre_sphPlm_deriv_array");
#endif
	al = shtns->alm[im];
	yl -= m;	dyl -= m;			// shift pointers

	st = sint;
	y0 = al[0];
	dy0 = 0.0;
	if (m>0) {		// m > 0
		l = m-1;
		if ((lmax <= SHT_L_RESCALE) || (sizeof(y0) > 8)) {
			if (l&1) {
				y0 = a_sint_pow_n_hp(y0, x, l-1) * sint;		// avoid computation of sqrt
			} else  y0 = a_sint_pow_n_hp(y0, x, l);
		} else {
			rescale = 1;
			y0 *= SHT_LEG_SCALEF;
		}
		dy0 = x*m*y0;
		st *= st;			// st = sin(theta)^2 is used in the recurrence for m>0
	}
	yl[m] = y0;		// l=m
	dyl[m] = dy0;
	if (lmax==m) goto done;		// done.

	y1 = al[1] * x * y0;
	dy1 = al[1]*( x*dy0 - st*y0 );
	yl[m+1] = y1;		// l=m+1
	dyl[m+1] = dy1;
	if (lmax==m+1) goto done;		// done.
	al+=2;

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
done:
	if (rescale != 0) {
		l = m-1;			// compute  sin(theta)^(m-1)
		if (l&1) {
			y0 = a_sint_pow_n_hp(1.0/SHT_LEG_SCALEF, x, l-1) * sint;		// avoid computation of sqrt
		} else  y0 = a_sint_pow_n_hp(1.0/SHT_LEG_SCALEF, x, l);
		for (l=m; l<=lmax; l++) {
			yl[l] *= y0;		dyl[l] *= y0;		// rescale
		}
	}
	return;
}

/// \internal Compute values of a legendre polynomial normalized for spherical harmonics derivatives, for a range of l=m..lmax, using recurrence.
/// Requires a previous call to \ref legendre_precomp(). Output is not directly compatible with GSL :
/// - if m=0 : returns ylm(x)  and  d(ylm)/d(theta) = -sin(theta)*d(ylm)/dx
/// - if m>0 : returns ylm(x)/sin(theta)  and  d(ylm)/d(theta).
/// This way, there are no singularities, everything is well defined for x=[-1,1], for any m.
/// \param lmax maximum degree computed, \param im = m/MRES with m the SH order, \param x argument, x=cos(theta).
/// \param sint = sqrt(1-x^2) to avoid recomputation of sqrt.
/// \param[out] yl is a double array of size (lmax-m+1) filled with the values (divided by sin(theta) if m>0)
/// \param[out] dyl is a double array of size (lmax-m+1) filled with the theta-derivatives.
void legendre_sphPlm_deriv_array(shtns_cfg shtns, const int lmax, const int im, const double cost, const double sint, double *yl, double *dyl)
{
	double *al;
	int l,m;
	double x, st, y0, y1, dy0, dy1;

	if ((im>0) && (lmax > SHT_L_RESCALE)) {		// switch to high precision code.
		legendre_sphPlm_deriv_array_hp(shtns, lmax, im, cost, sint, yl, dyl);
		return;
	}

	x = cost;
	m = im*MRES;
#ifdef LEG_RANGE_CHECK
	if ((lmax > LMAX)||(lmax < m)||(im>MMAX)) shtns_runerr("argument out of range in legendre_sphPlm_deriv_array");
#endif
	al = shtns->alm[im];
	yl -= m;	dyl -= m;			// shift pointers

	st = sint;
	y0 = al[0];
	dy0 = 0.0;
	if (m>0) {		// m > 0
		l = m-1;
		if (l&1) {
			y0 = a_sint_pow_n(y0, x, l-1) * sint;		// avoid computation of sqrt
		} else  y0 = a_sint_pow_n(y0, x, l);
		dy0 = x*m*y0;
		st *= st;			// st = sin(theta)^2 is used in the recurrence for m>0
	}
	yl[m] = y0;		// l=m
	dyl[m] = dy0;
	if (lmax==m) return;		// done.

	y1 = al[1] * x * y0;
	dy1 = al[1]*( x*dy0 - st*y0 );
	yl[m+1] = y1;		// l=m+1
	dyl[m+1] = dy1;
	if (lmax==m+1) return;		// done.
	al+=2;

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
}

/// \internal Precompute constants for the recursive generation of Legendre associated functions, with given normalization.
/// this function is called by \ref shtns_set_size, and assumes up-to-date values in \ref shtns.
/// For the same conventions as GSL, use \c legendre_precomp(sht_orthonormal,1);
/// \param[in] norm : normalization of the associated legendre functions (\ref shtns_norm).
/// \param[in] with_cs_phase : Condon-Shortley phase (-1)^m is included (1) or not (0)
/// \param[in] mpos_renorm : Optional renormalization for m>0.
///  1.0 (no renormalization) is the "complex" convention, while 0.5 leads to the "real" convention (with FFTW).
void legendre_precomp(shtns_cfg shtns, enum shtns_norm norm, int with_cs_phase, double mpos_renorm)
{
	double *al0, *bl0;
	double **alm, **blm;
	long int im, m, l, lm;
	double t1, t2;

#if SHT_VERBOSE > 1
	printf("        > Condon-Shortley phase = %d, normalization = %d\n", with_cs_phase, norm);
#endif

	if (with_cs_phase != 0) with_cs_phase = 1;		// force to 1 if !=0

	im = (MMAX+1)*sizeof(double*) + (MIN_ALIGNMENT-1);		// alloc memory for arrays + sse alignement.
	alm = (double **) malloc( im + (2*NLM)*sizeof(double) );
	al0 = (double *) PTR_ALIGN( alm + (MMAX+1) );
	bl0 = al0;		blm = alm;		// by default analysis recurrence is the same
	if ((norm == sht_schmidt) || (mpos_renorm != 1.0)) {
		blm = (double **) malloc( im + (2*NLM)*sizeof(double) );
		bl0 = (double *) PTR_ALIGN( blm + (MMAX+1) );
	}

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

/// - Compute analysis recurrence coefficients if necessary
	if ((norm == sht_schmidt) || (mpos_renorm != 1.0)) {
		im = ((MMAX+2)>>1)*2;
		blm = (double **) malloc( im * sizeof(double *) + (2*NLM)*sizeof(double) );
		bl0 = (double *) (blm + im);
		for (lm=0; lm<2*NLM; lm++) bl0[lm] = al0[lm];		// copy
		for (im=0, lm=0; im<=MMAX; im++) {
			m = im*MRES;
			blm[im] = bl0 + lm;
			double c0 = 1.0;
			if (m>0) c0 = 1.0/mpos_renorm;
			if (norm == sht_schmidt) {
				bl0[lm+1] *= (2*m+3)/(2*m+1.);
				c0 *= 2*m+1;
			}
			bl0[lm] *= c0;
			lm+=2;
			for (l=m+2; l<=LMAX; l++) {
				if (norm == sht_schmidt) {
					bl0[lm] *= (2*l+1)/(2*l-3.);
					bl0[lm+1] *= (2*l+1)/(2*l-1.);
				}
				lm+=2;
			}
		}
	}

	shtns->alm = alm;		shtns->blm = blm;
}

/// \internal returns the value of the Legendre Polynomial of degree l.
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


/// \internal Generates the abscissa and weights for a Gauss-Legendre quadrature.
/// Newton method from initial Guess to find the zeros of the Legendre Polynome
/// \param x = abscissa, \param w = weights, \param n points.
/// \note Reference:  Numerical Recipes, Cornell press.
void gauss_nodes(long double *x, long double *w, int n)
{
	long int i,l,m, k;
	long double z, z1, p1, p2, p3, pp;
	long double eps;

	long double l_1[n-1];		// cache for 1/l on stack.
	for (l=2; l<=n; l++) l_1[l-2] = 1.0L/l;		// precompute 1/l

	eps = 2.3e-16;		// desired precision, minimum = 2.2204e-16 (double)
	if (sizeof(eps) > 8)	eps = 1.1e-19;		// desired precision, minimum = 1.0842e-19 (long double i387)

	m = (n+1)/2;
	for (i=1;i<=m;i++) {
		k=10;		// maximum Newton iteration count to prevent infinite loop.
		z = (1.0 - (n-1.)/(8.*n*n*n)) * cos(M_PI*(i-0.25)/(n+0.5));	// initial guess
		do {
			p1 = z;		// P_1
			p2 = 1.0;	// P_0
			for(l=2;l<=n;l++) {		 // recurrence : l P_l = (2l-1) z P_{l-1} - (l-1) P_{l-2}	(works ok up to l=100000)
				p3 = p2;
				p2 = p1;
				p1 = ((2*l-1)*z*p2 - (l-1)*p3)*l_1[l-2];		// The Legendre polynomial...
			}
			pp = ((1.-z)*(1.+z))/(n*(p2-z*p1));			// ... and its inverse derivative.
			z1 = z;
			z -= p1*pp;		// Newton's method
		} while (( fabsl(z-z1) > (z1+z)*0.5*eps ) && (--k > 0));
		x[i-1] = z;		// Build up the abscissas.
		w[i-1] = 2.0*pp*pp/((1.-z)*(1.+z));		// Build up the weights.
		x[n-i] = -z;
		w[n-i] = w[i-1];
	}
	if (n&1) x[n/2] = 0.0;		// exactly zero.

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
