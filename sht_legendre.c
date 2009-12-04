/** \file sht_legendre.c
 \brief Compute legendre polynomials and associated functions.
 The normalization of the associated functions is for spherical harmonics, and the Condon-Shortley phase is included.
 When computing the derivatives (with respect to colatitude theta), there are no singularities.
 written by Nathanael Schaeffer / LGIT,CNRS, with some ideas and code from the GSL 1.13 and Numerical Recipies.
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
LEG_FLOAT_TYPE *alm;	///< coefficient list for recursion (size NLM)

/// Returns the value of a legendre polynomial of degree l and order m, noramalized for spherical harmonics, using recursion.
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
/// for a range of l=m..lmax, at given m and x, using recursion.
/// Requires a previous call to \ref legendre_precomp().
/// Output compatible with the GSL function gsl_sf_legendre_sphPlm_array(lmax, m, x, yl)
/// \param lmax maximum degree computed, \param m order, \param x argument.
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


/// Compute values of a legendre polynomial normalized for spherical harmonics derivatives, for a range of l=m..lmax, using recursion.
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
		st = sint*sint;		// st = sin(theta)^2 is used in the recursion for m>0
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

/// precompute constants for the legendre recursion, of size specified by \ref shtns_set_size
void legendre_precomp()
{
	long int im, m, l, lm;
	LEG_FLOAT_TYPE t1, t2;

#if SHT_VERBOSE > 1
	printf("        > using custom fast recursion for legendre polynomials\n");
#endif

	mmidx = (int *) malloc( (MMAX+1) * sizeof(int) );
	alm = (LEG_FLOAT_TYPE *) malloc( 2*NLM * sizeof(LEG_FLOAT_TYPE) );

/// precompute the factors alm and blm of the recurrence relation :
/// y(l,m) = alm[l,m]*x*y(l-1,m) - blm[l,m]*y(l-2,m)
	for (im=0, lm=0; im<=MMAX; im++) {
		m = im*MRES;
		mmidx[im] = lm;
		if (m < LMAX) {	// l=m+1
			alm[lm] = 0.0;		// will contain the starting value.
			alm[lm+1] = LEG_SQRT(m+m+3);
			t2 = (m+m+1);	lm+=2;
		}
		for (l=m+2; l<=LMAX; l++) {
			t1 = (l+m)*(l-m);
			alm[lm+1] = LEG_SQRT((l+l+1)*(l+l-1)/t1);
			alm[lm] = - LEG_SQRT(((l+l+1)*t2)/((l+l-3)*t1));
			t2 = t1;	lm+=2;
		}
	}

/// Starting value for recursion :
/// Y_m^m(x) = sqrt( (2m+1)/(4pi m) gamma(m+1/2)/gamma(m) ) (-1)^m (1-x^2)^(m/2) / pi^(1/4)
// alternate method : direct recursive computation, using the following properties:
//	gamma(x+1) = x*gamma(x)   et   gamma(1.5)/gamma(1) = sqrt(pi)/2
#ifdef LEG_LONG_DOUBLE
	t1 = 0.25L/M_PIl;
#else
	t1 = 0.25/M_PI;
#endif
	alm[0] = LEG_SQRT(t1);          // Y00 = 1/sqrt(4.pi)
	for (im=1, m=0; im<=MMAX; im++) {
		while(m<im*MRES) {
			m++;
			t1 *= ((LEG_FLOAT_TYPE)m + 0.5)/m;	// t1 *= (m+0.5)/m;
		}
		t2 = (m&1) ? -1.0 : 1.0;	// (-1)^m  Condon-Shortley phase.
		alm[mmidx[im]] = t2 * LEG_SQRT(t1);
	}
}

/// returns the value of the Legendre Polynomial of degree l.
/// Does not require a previous call to legendre_precomp(), l is arbitrary.
double legendre_Pl(const int l, double x)
{
	long int i;
	LEG_FLOAT_TYPE p1, p2, p3;

	if ((l==0)||(x==1.0)) return ( 1. );
	if (x==-1.0) return ( (l&1) ? -1. : 1. );

	p2 = 1.0;		// P_0
	p1 = x;			// P_1
	for (i=2; i<=l; i++) {		 // recurrence : l P_l = (2l-1) z P_{l-1} - (l-1) P_{l-2}	(works ok up to l=100000)
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
