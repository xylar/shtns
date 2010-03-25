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

// it appears that long double precision does not help the precision of the SHT (if the SHT itself is computed in double precision)
//#define LEG_LONG_DOUBLE

#ifdef LEG_LONG_DOUBLE
  #define LEG_FLOAT_TYPE long double
  #define LEG_SQRT(a) sqrtl(a)
#else
  #define LEG_FLOAT_TYPE double
  #define LEG_SQRT(a) sqrt(a)
#endif

#if SHT_VERBOSE > 1
  #define LEG_RANGE_CHECK
#endif

/// computes sin(t)^n from cos(t). ie returns (1-x^2)^(n/2), with x = cos(t)
inline LEG_FLOAT_TYPE sint_pow_n(LEG_FLOAT_TYPE cost, int n)
{
	LEG_FLOAT_TYPE val = 1.0;
	LEG_FLOAT_TYPE s2 = (1.-cost)*(1.+cost);		// sin(t)^2 = 1 - cos(t)^2
	if (n&1) val *= LEG_SQRT(s2);	// = sin(t)
	n >>= 1;
	while(n>0) {
		if (n&1) val *= s2;
		n >>= 1;
		s2 *= s2;
	}
	return val;		// = sint(t)^n
}

int *mmidx;				///< index array (size MMAX+1)
LEG_FLOAT_TYPE *alm;	///< coefficient list for recurrence (size 2*NLM)

/// Returns the value of a legendre polynomial of degree l and order im*MRES, noramalized for spherical harmonics, using recurrence.
/// Requires a previous call to \ref legendre_precomp().
/// Output compatible with the GSL function gsl_sf_legendre_sphPlm(l, m, x)
double legendre_sphPlm(const int l, const int im, double x)
{
	long int i,m,lm;
	LEG_FLOAT_TYPE ymm, ymmp1;

	m = im*MRES;
#ifdef LEG_RANGE_CHECK
	if ( (l>LMAX) || (l<m) || (im>MMAX) ) shtns_runerr("argument out of range in legendre_sphPlm");
#endif

	lm = mmidx[im];
	ymm = alm[lm] * sint_pow_n(x, m);		// l=m
	if (l==m) return ((double) ymm);

	ymmp1 = ymm * alm[lm+1] * x;				// l=m+1
	lm+=2;
	if (l == m+1) return ((double) ymmp1);
	
	for (i=m+2; i<l; i+=2) {
		ymm   = alm[lm+1]*x*ymmp1 + alm[lm]*ymm;
		ymmp1 = alm[lm+3]*x*ymm + alm[lm+2]*ymmp1;
		lm+=4;
	}
	if (i==l) {
		ymmp1 = alm[lm+1]*x*ymmp1 + alm[lm]*ymm;
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
	long int l,m,lm;
	LEG_FLOAT_TYPE ymm, ymmp1;

	m = im*MRES;
#ifdef LEG_RANGE_CHECK
	if ((lmax > LMAX)||(lmax < m)||(im>MMAX)) shtns_runerr("argument out of range in legendre_sphPlm_array");
#endif

	lm = mmidx[im];
	ymm = alm[lm] * sint_pow_n(x, m);	// l=m
	yl[0] = ymm;
	if (lmax==m) return;		// done.

	ymmp1 = ymm * alm[lm+1] * x;		// l=m+1
	yl[1] = ymmp1;
	lm+=2;
	if (lmax==m+1) return;		// done.

	yl -= m;			// shift pointer
	for (l=m+2; l<lmax; l+=2) {
		ymm   = alm[lm+1]*x*ymmp1 + alm[lm]*ymm;
		ymmp1 = alm[lm+3]*x*ymm + alm[lm+2]*ymmp1;
		yl[l] = ymm;
		yl[l+1] = ymmp1;
		lm+=4;
	}
	if (l==lmax) {
		yl[l] = alm[lm+1]*x*ymmp1 + alm[lm]*ymm;
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
	long int l,m,lm;
	LEG_FLOAT_TYPE st, y0, y1, dy0, dy1;

	m = im*MRES;
#ifdef LEG_RANGE_CHECK
	if ((lmax > LMAX)||(lmax < m)||(im>MMAX)) shtns_runerr("argument out of range in legendre_sphPlm_deriv_array");
#endif
	lm = mmidx[im];

	st = sint;
	y0 = alm[lm];
	dy0 = 0.0;
	if (im>0) {		// m > 0
		l = m-1;			// compute  sin(theta)^(m-1)
		while(l>0) {
			if (l&1) y0 *= st;
			l >>= 1;
			st *= st;
		}
		dy0 = x*m*y0;
		st = sint*sint;		// st = sin(theta)^2 is used in the recurrence for m>0
	}
	yl[0] = y0;		// l=m
	dyl[0] = dy0;
	if (lmax==m) return;		// done.

	y1 = alm[lm+1] * x * y0;
	dy1 = alm[lm+1]*( x*dy0 - st*y0 );
	yl[1] = y1;		// l=m+1
	dyl[1] = dy1;
	lm+=2;
	if (lmax==m+1) return;		// done.

	yl -= m;	dyl -= m;			// shift pointers
	for (l=m+2; l<lmax; l+=2) {
		y0 = alm[lm+1]*x*y1 + alm[lm]*y0;
		dy0 = alm[lm+1]*(x*dy1 - y1*st) + alm[lm]*dy0;
		yl[l] = y0;		dyl[l] = dy0;
		y1 = alm[lm+3]*x*y0 + alm[lm+2]*y1;
		dy1 = alm[lm+3]*(x*dy0 - y0*st) + alm[lm+2]*dy1;
		yl[l+1] = y1;		dyl[l+1] = dy1;
		lm+=4;
	}
	if (l==lmax) {
		yl[l] = alm[lm+1]*x*y1 + alm[lm]*y0;
		dyl[l] = alm[lm+1]*(x*dy1 - y1*st) + alm[lm]*dy0;
	}

/*	// Simple loop, without temporary variables.
	yl -= m;	dyl -= m;			// shift pointers
	for (l=m+2; l<=lmax; l++) {
		yl[l] = alm[lm+1]*x*yl[l-1] + alm[lm]*yl[l-2];
		dyl[l] = alm[lm+1]*(x*dyl[l-1] - yl[l-1]*st) + alm[lm]*dyl[l-2];
		lm+=2;
	}

/*	// Alternative for evaluating the derivative (not better)
	for (l=m+2; l<=lmax; l++) {
//		dyl/dx = - (l * x * y[l] - c1 * (l+m) * y[l-1]) / (1-x^2);
//		=> dyl/dtheta = (l * x * y[l] - c1 * (l+m) * y[l-1]) / sqrt(1-x^2);
		if (m==0) {
			const double c1 = sqrt( (2.*l+1.)/(2.*l-1.) );
			dyl[l] = l*(x*yl[l] - c1*yl[l-1])/st;
			if ( 1.-fabs(x) < 1.e-15 ) dyl[l] = 0.0;	// -l*(l+1)/2 *sin(theta)
		} else {
			const double c1 = sqrt(((2.*l+1.)/(2.*l-1.)) * ((double)(l-m)/(double)(l+m)));
			dyl[l] = l*x*yl[l] - c1*(l+m)*yl[l-1];
		}
	}
*/

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
	long int im, m, l, lm;
	LEG_FLOAT_TYPE t1, t2;

#if SHT_VERBOSE > 1
	printf("        > using custom fast recurrence for legendre polynomials\n");
	printf("        > Condon-Shortley phase = %d, normalization = %d\n", with_cs_phase, norm);
#endif

	if (with_cs_phase != 0) with_cs_phase = 1;		// force to 1 if !=0

	mmidx = (int *) malloc( (MMAX+1) * sizeof(int) );
	alm = (LEG_FLOAT_TYPE *) malloc( 2*NLM * sizeof(LEG_FLOAT_TYPE) );

/// - Precompute the factors alm and blm of the recurrence relation :
  if (norm == sht_schmidt) {
	for (im=0, lm=0; im<=MMAX; im++) {		/// <b> For Schmidt semi-normalized </b>
		m = im*MRES;
		mmidx[im] = lm;
		t2 = LEG_SQRT(2*m+1);
		alm[lm] = 1.0/t2;		/// starting value will be divided by \f$ \sqrt{2m+1} \f$ 
		alm[lm+1] = t2;		// l=m+1
		lm+=2;
		for (l=m+2; l<=LMAX; l++) {
			t1 = LEG_SQRT((l+m)*(l-m));
			alm[lm+1] = (2*l-1)/t1;		/// \f[  a_l^m = \frac{2l-1}{\sqrt{(l+m)(l-m)}}  \f]
			alm[lm] = - t2/t1;			/// \f[  b_l^m = -\sqrt{\frac{(l-1+m)(l-1-m)}{(l+m)(l-m)}}  \f]
			t2 = t1;	lm+=2;
		}
	}
  } else {
	for (im=0, lm=0; im<=MMAX; im++) {		/// <b> For orthonormal or 4pi-normalized </b>
		m = im*MRES;
		mmidx[im] = lm;
		t2 = 2*m+1;
		alm[lm] = 1.0;		// will contain the starting value.
		alm[lm+1] = LEG_SQRT(2*m+3);		// l=m+1
		lm+=2;
		for (l=m+2; l<=LMAX; l++) {
			t1 = (l+m)*(l-m);
			alm[lm+1] = LEG_SQRT(((2*l+1)*(2*l-1))/t1);			/// \f[  a_l^m = \sqrt{\frac{(2l+1)(2l-1)}{(l+m)(l-m)}}  \f]
			alm[lm] = - LEG_SQRT(((2*l+1)*t2)/((2*l-3)*t1));	/// \f[  b_l^m = -\sqrt{\frac{2l+1}{2l-3}\,\frac{(l-1+m)(l-1-m)}{(l+m)(l-m)}}  \f]
			t2 = t1;	lm+=2;
		}
	}
  }

/// - Compute and store the prefactor (independant of x) of the starting value for the recurrence :
/// \f[  Y_m^m(x) = Y_0^0 \ \sqrt{ \prod_{k=1}^{m} \frac{2k+1}{2k} } \ \ (-1)^m \ (1-x^2)^{m/2}  \f]
	if ((norm == sht_fourpi)||(norm == sht_schmidt)) {
		t1 = 1.0;
		alm[0] = t1;		/// \f$ Y_0^0 = 1 \f$  for Schmidt or 4pi-normalized 
	} else {
		t1 = 0.25L/M_PIl;
		alm[0] = LEG_SQRT(t1);		/// \f$ Y_0^0 = 1/\sqrt{4\pi} \f$ for orthonormal
	}
	t1 *= mpos_renorm;		// renormalization for m>0
	for (im=1, m=0; im<=MMAX; im++) {
		while(m<im*MRES) {
			m++;
			t1 *= ((LEG_FLOAT_TYPE)m + 0.5)/m;	// t1 *= (m+0.5)/m;
		}
		t2 = LEG_SQRT(t1);
		if ( m & with_cs_phase ) t2 = -t2;		/// optional \f$ (-1)^m \f$ Condon-Shortley phase.
		alm[mmidx[im]] *= t2;
	}
}

/// returns the value of the Legendre Polynomial of degree l.
/// l is arbitrary, a direct recurrence relation is used, and a previous call to legendre_precomp() is not required.
double legendre_Pl(const int l, double x)
{
	long int i;
	LEG_FLOAT_TYPE p1, p2, p3;

	if ((l==0)||(x==1.0)) return ( 1. );
	if (x==-1.0) return ( (l&1) ? -1. : 1. );

	p2 = 1.0;		/// \f$  P_0 = 1  \f$
	p1 = x;			/// \f$  P_1 = x  \f$
	for (i=2; i<=l; i++) {		 /// recurrence : \f[  l P_l = (2l-1) x P_{l-1} - (l-1) P_{l-2}	 \f] (works ok up to l=100000)
		p3 = p2;
		p2 = p1;
		p1 = (x*(2*i-1)*p2 - (i-1)*p3)/i;
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
