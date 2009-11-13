/**************************************************
 * Spherical Harmonics (pre)computation           *
 *    some ideas and code taken from the GSL 1.13 *
 *    written by Nathanael Schaeffer / LGIT,CNRS  *
 **************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


//#define SHT_LEGENDRE_GSL
//#define LEG_LONG_DOUBLE

#ifdef LEG_LONG_DOUBLE
  #define LEG_FLOAT_TYPE long double
  #define LEG_SQRT(a) sqrtl(a)
#else
  #define LEG_FLOAT_TYPE double
  #define LEG_SQRT(a) sqrt(a)
#endif

#ifdef SHT_LEGENDRE_GSL

#include <gsl/gsl_sf_legendre.h>

double legendre_sphPlm(int l, int m, double x) {
	return gsl_sf_legendre_sphPlm(l, m, x);
}

int legendre_sphPlm_array(int lmax, int m, double x, double *yl) {
	return gsl_sf_legendre_sphPlm_array(lmax, m, x, yl); 
}

int legendre_sphPlm_deriv_array(int lmax, int m, double x, double *yl, double *dyl)
{
	double st,stmin;
	int res, l;

	st = sqrt((1.-x)*(1.+x));
	if (m==0) {
		res = gsl_sf_legendre_sphPlm_deriv_array(lmax, m, x, yl, dyl);
		for (l=m; l<=lmax; l++) {
			dyl[l-m] *= -st;		// multiply by - sin(theta) to compute dyl/dtheta
		}
		return res;
	}
	
	if (st > 0.0) {
		res = gsl_sf_legendre_sphPlm_deriv_array(lmax, m, x, yl, dyl);
		for (l=m; l<=lmax; l++) {
			dyl[l-m] *= -st;		// multiply by - sin(theta) to compute dyl/dtheta
			yl[l-m] *= 1./st;
		}
		return res;
	}
	
	if (m >= 2) {
		for (l=m; l<=lmax; l++) {	// spherical harmonics for m>1 are zero near the poles.
			yl[l-m] = 0.0;
			dyl[l-m] = 0.0;
		}
		return 0.0;		// success.
	} else {	// m=1 here
		stmin = 1.e-6;
	// gsl function diverges for m=1 and sint=0 => use (bad) approximation
		printf("bad approximation used for m=1 spherical harmonics neer the poles !\n");
		res=gsl_sf_legendre_sphPlm_array(lmax, m, sqrt((1.-stmin)*(1.+stmin)) *((x<0.)? -1.:1.), dyl);
		for (l=m; l<=lmax; l++) {		// d(Pl1)/dt |(t=0) = Pl1(epsilon)/sin(epsilon)
			dyl[l-m] *= 1./stmin;
			yl[l-m] = dyl[l-m];
		}
		return res;
	}
}

void legendre_precomp()
{
#if SHT_VERBOSE > 0
	printf("        > using GSL for legendre polynomials\n");
#endif
	// nothing to do when using the gsl.
}

/// returns the value of the Legendre Polynomial of degree l
double legendre_Pl(const int l, double x)
{
	return gsl_sf_legendre_Pl(l, x);
}

#else

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

LEG_FLOAT_TYPE *alm;
LEG_FLOAT_TYPE *blm;
int *mmidx;

/// compute value of a legendre polynomial, using recursion.
/// (requires a previous call to legendre_precompute())
double legendre_sphPlm(const int l, const int m, double x)
{
	LEG_FLOAT_TYPE ymm, ymmp1;
	int i,lm;

	if ( (l<m) || (l>LMAX) || (m>MMAX*MRES) || (m%MRES) ) return 0.0;		// out of range.

	if ( ((m>0)&&((x==1.)||(x==-1.))) )  return 0.0;

	lm = mmidx[m/MRES];
	ymm = alm[lm] * sint_pow_n(x, m);		// l=m
	lm++;
	if (l==m) return ((double) ymm);

	ymmp1 = ymm * alm[lm] * x;				// l=m+1
	lm++;
	if (l == m+1) return ((double) ymmp1);
	
	for (i=m+2; i<l; i+=2) {
		ymm   = alm[lm]*x*ymmp1 + blm[lm]*ymm;
		ymmp1 = alm[lm+1]*x*ymm + blm[lm+1]*ymmp1;
		lm+=2;
	}
	if (i==l) {
		ymmp1 = alm[lm]*x*ymmp1 + blm[lm]*ymm;
	}
	return ((double) ymmp1);
}

/// compute value of a legendre polynomial for a range of l, using recursion.
/// requires a previous call to legendre_precompute()
int legendre_sphPlm_array(const int lmax, const int m, double x, double *yl)
{
	LEG_FLOAT_TYPE ymm, ymmp1;
	int l,lm;

	if ((lmax < m)||(lmax > LMAX)||(m>MMAX*MRES)||(m % MRES)) return -1;

	if ((m>0)&&((x==1.)||(x==-1.))) {
		for(l=m; l<=lmax; l++) {
			yl[l-m] = 0.0;
		}
		return 0;
	}

	lm = mmidx[m/MRES];
	ymm = alm[lm] * sint_pow_n(x, m);	// l=m
	yl[0] = ymm;
	lm++;
	if (lmax==m) return 0;		// done.

	ymmp1 = ymm * alm[lm] * x;		// l=m+1
	yl[1] = ymmp1;
	lm++;
	if (lmax==m+1) return 0;		// done.

	for (l=m+2; l<lmax; l+=2) {
		ymm   = alm[lm]*x*ymmp1 + blm[lm]*ymm;
		ymmp1 = alm[lm+1]*x*ymm + blm[lm+1]*ymmp1;
		yl[l-m] = ymm;
		yl[l-m+1] = ymmp1;
		lm+=2;
	}
	if (l==lmax) {
		yl[l-m] = alm[lm]*x*ymmp1 + blm[lm]*ymm;
	}
	return 0;
}

/// compute value of a legendre polynomial for a range of l, using recursion.
/// requires a previous call to legendre_precompute()
/// for m=0 : returns ylm(x)  and  d(ylm)/d(theta)
/// for m>0 : returns ylm(x)/sin(theta)  and  d(ylm)/d(theta)
int legendre_sphPlm_deriv_array(const int lmax, const int m, double x, double *yl, double *dyl)
{
	double st;
	int l,lm;

	if ((lmax < m)||(lmax > LMAX)||(m>MMAX*MRES)||(m % MRES)) return -1;

	if ((m>1)&&((x==1.)||(x==-1.))) {
		for(l=m; l<=lmax; l++) {
			yl[l-m] = 0.0;	dyl[l-m] = 0.0;
		}
		return 0;
	}

	lm = mmidx[m/MRES];
	st = (1.-x)*(1.+x);		// st = sin(theta)^2 is used in the recursion for m>0
	
	if (m==0) {
		st = sqrt(st);		// we need the square root for the recursion.
		yl[0] = alm[lm];	// l=m
		dyl[0] = 0.0;
		lm++;
	} else {		// m > 0
		double stm1 = alm[lm];
		double st2 = st;
		l = m-1;			// compute  sin(theta)^(m-1)
		if (l&1) stm1 *= sqrt(st2);
		l >>= 1;
		while(l>0) {
			if (l&1) stm1 *= st2;
			l >>= 1;
			st2 *= st2;
		}
		yl[0] = stm1;	// l=m
		dyl[0] = x*m*stm1;
		lm++;
	}

	if (lmax==m) return 0;		// done.

	yl[1] = alm[lm] * x * yl[0];		// l=m+1
	dyl[1] = alm[lm]*( x*dyl[0] - st*yl[0] );
	lm++;

	for (l=m+2; l<=lmax; l++) {
		yl[l-m] = alm[lm]*x*yl[l-m-1] + blm[lm]*yl[l-m-2];
		dyl[l-m] = alm[lm]*(x*dyl[l-m-1] - yl[l-m-1]*st) + blm[lm]*dyl[l-m-2];
		lm++;
	}

	return 0;
}


/// precompute constants for the legendre recursion, of size specified by \ref shtns_set_size
void legendre_precomp()
{
	LEG_FLOAT_TYPE t1, t2;
	int im, m, l, lm;

#if SHT_VERBOSE > 0
	printf("        > using custom fast recursion for legendre polynomials\n");
#endif

	alm = (LEG_FLOAT_TYPE *) malloc( 2*NLM * sizeof(LEG_FLOAT_TYPE) );
	blm = alm + NLM;
	mmidx = (int *) malloc( (MMAX+1) * sizeof(int) );

// precompute the factors of the recurrence relation :
// y(l,m) = alm[l,m]*x*y(l-1,m) - blm[l,m]*y(l-2,m)
	for (im=0, lm=0; im<=MMAX; im++) {
		m = im*MRES;
		mmidx[im] = lm;
		//l=m;
			alm[lm] = 0.0;	blm[lm] = 0.0;	lm++;
		if (m < LMAX) {	// l=m+1
			alm[lm] = LEG_SQRT(m+m+3);
			blm[lm] = 0.0;
			t2 = (m+m+1);	lm++;
		}
		for (l=m+2; l<=LMAX; l++) {
			t1 = (l+m)*(l-m);
			alm[lm] = LEG_SQRT((l+l+1)*(l+l-1)/t1);
			blm[lm] = - LEG_SQRT(((l+l+1)*t2)/((l+l-3)*t1));
			t2 = t1;	lm++;
		}
	}

// Starting value for recursion.
// Y_m^m(x) = sqrt( (2m+1)/(4pi m) gamma(m+1/2)/gamma(m) ) (-1)^m (1-x^2)^(m/2) / pi^(1/4)
// alternate method : direct recursive computation, using the following properties:
//	gamma(x+1) = x*gamma(x)   et   gamma(1.5)/gamma(1) = sqrt(pi)/2
	alm[0] = 0.5/LEG_SQRT(M_PI);          // Y00 = 1/sqrt(4pi)
	m=1;	t1 = 1.0;
	for (im=1; im<=MMAX; im++) {
		while(m<im*MRES) {
			t1 *= (m+0.5)/m;
			m++;
		}
		t2 = (m&1) ? -1.0 : 1.0;	// (-1)^m
		alm[mmidx[im]] = t2 * LEG_SQRT( (2.0+1.0/m)/(8.0*M_PI) * t1 );
	}
}

/// returns the value of the Legendre Polynomial of degree l
double legendre_Pl(const int l, double x)
{
	double p1, p2, p3;
	int i;

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

/* test code

int main()
{
	double x,v1, v2, dmax;
	int n=30;
	int mmax=15;
	int mres=1;
	int lmax=50;
	int i,im,m,l;
	
	legendre_precomp(lmax, mmax, mres);
	dmax = 0.;
	
	for (i=0;i<n;i++) {
		x = cos((i+0.5)*M_PI/n);
		for (im=0; im<=mmax; im++) {
			m=im*mres;
			for (l=m; l<=lmax; l+=mres) {
				v1 = gsl_sf_legendre_sphPlm(l,m,x);
				v2 = legendre_sphPlm(l,m,x);
//				v1 = legendre_sphPlm(l,m,x+1.e-16);
				if ((fabs((v2 - v1)/v1) > dmax)&&(fabs(v1)>1e-16)) {
					printf("m=%d l=%d x=%f   gsl=%e   me=%e   diff=%e  err=%e\n",m,l,x,v1,v2,v2-v1,(v2-v1)/v1);
					dmax = fabs((v2-v1)/v1);
				}
			}
		} 
	}
	printf("max error = %e\n",dmax);
}
*/
#endif



/// Generates the abscissa and weights for a Gauss-Legendre quadrature.
/// Newton method from initial Guess to find the zeros of the Legendre Polynome
/// \param x = abscissa, \param w = weights, \param n points.
/// \note Reference:  Numerical Recipes, Cornell press.
void gauss_nodes(long double *x, long double *w, int n)
{
	long double z, z1, p1, p2, p3, pp, eps;
	long int i,l,m;
	long double pi = M_PI;

	eps = 1.1e-19;	// desired precision, minimum = 1.0842e-19 (long double)

	m = (n+1)/2;
	for (i=1;i<=m;i++) {
		z = cosl(pi*((long double)i-0.25)/((long double)n+0.5));	// initial guess
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

