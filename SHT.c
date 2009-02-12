/**
///////////////////////////////////////////////
// SHT : Spherical Harmonic Transform
//   requires SHT.h for size parameters.
//////////////////////////////////////////////
*/

#include <complex.h>
#include <math.h>
// FFTW la derivee d/dx = ik	(pas de moins !)
#include <fftw3.h>
// GSL for Legendre functions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>

	// cycle counter from FFTW
	#include "cycle.h"


// SHT.h : parameter for SHT (sizes : LMAX, NLAT, MMAX, MRES, NPHI)
#ifndef LMAX
# include "SHT.h"
#endif

/*
#if MMAX == 0
  #undef MRES
  #define MRES 1
#endif
*/

#ifndef MRES
	int _mres_ = -1;	// if no MRES is defined, use run-time MRES instead.
	int _nlm_ = 0;
	#define MRES _mres_
#endif

/** ACCESS TO SPHERICAL HARMONICS COMPONENTS **
  NLM : total number of (l,m) components.
  LM(l,m) : macro returning array index for given l and m.
  LiM(l,im) : macro returning array index for given l and im.
  LM_LOOP( action ) : macro that performs "action" for every (l,m), with l and lm set, but neither m nor im.
  el[NLM], l2[NLM], l_2[NLM] : floating point arrays containing l, l*(l+1) and 1/(l*(l+1))
**/
// NLM : total number of Spherical Harmonic coefficients.
#define NLM ( (MMAX+1)*(LMAX+1) - MRES*(MMAX*(MMAX+1))/2 )
// LM(l,m) : index in the Spherical Harmonic coefficient array [ (l,m) space ]
#define LiM(l,im) ( (im*(2*LMAX+2 -MRES*(im+1)))/2 + l )
#define LM(l,m) ( (m*(2*LMAX+2 -(m+MRES)))/(2*MRES) + l )

// LM_LOOP : loop over all (l,im) and perform "action"  : l and lm are defined. (but NOT m and im)
//  double-loop version : assumes contiguous l-storage with 1 stride. (not compatible with even/odd packed)
#define LM_LOOP( action ) for (im=0, lm=0; im<=MMAX; im++) { for (l=im*MRES; l<=LMAX; l++, lm++) { action } }
//  single-loop + lookup array : no assumption on l-storage is made.
//#define LM_LOOP( action ) for (lm=0; lm<NLM; lm++) { l=el[lm]; { action } }


// half the theta length for even/odd decomposition. (NLAT+1)/2 allows odd NLAT.
#ifndef NLAT
 #ifndef NLAT_2
  #error "NLAT or NLAT_2 must be defined"
 #endif
 #ifndef _SHT_EO_
 #define NLAT (2*NLAT_2)
 #else
  #define NLAT NLAT_2
 #endif
#endif

#ifndef NLAT_2
 #ifndef NLAT
  #error "NLAT or NLAT_2 must be defined"
 #endif
 #ifndef _SHT_EO_
  #define NLAT_2 ((NLAT+1)/2)
 #else
  #define NLAT_2 NLAT
 #endif
#endif

#ifndef M_PI
# define M_PI 3.1415926535897932384626433832795
#endif

// useful values for some basic spherical harmonic representations
// Y00_1 = 1/Y00 = spherical harmonic representation of 1 (l=0,m=0)
#define Y00_1 sqrt(4.*M_PI)
// Y10_ct = spherical harmonic representation of cos(theta) (l=1,m=0)
#define Y10_ct sqrt(4.*M_PI/3.)

// TEST FOR PARAMETERS AT COMPILE TIME
#if NLAT_2 <= LMAX/2
	#error "NLAT_2 must be greater than LMAX/2 !"
#endif
#if LMAX < MMAX*MRES
	#error "MMAX*MRES should not exceed LMAX !"
#endif
#if NPHI <= 2*MMAX
	#error "NPHI and MMAX must conform to the sampling condition NPHI > 2*MMAX !"
#endif


struct DtDp {		// theta and phi derivatives stored together.
	double t, p;
};

const double pi = M_PI;
double ct[NLAT], st[NLAT], st_1[NLAT];	// cos(theta), sin(theta), 1/sin(theta);
#if MRES != _mres_
	double el[NLM], l2[NLM], l_2[NLM];	// l, l(l+1) and 1/(l(l+1))
#else
	double *el, *l2, *l_2;		// l, l(l+1) and 1/(l(l+1))
#endif

long int tm[MMAX+1];	// start theta value for SH (polar optimization : near the poles the legendre polynomials go to zero for high m's)

double* ylm[MMAX+1];		// matrix for inverse transform (synthesis)
struct DtDp* dylm[MMAX+1];	// theta and phi derivative of Ylm matrix
double* zlm[MMAX+1];		// matrix for direct transform (analysis)
struct DtDp* dzlm[MMAX+1];

#ifdef SHT_DCT
double* ylm_dct[MMAX+1];	// matrix for inverse transform (synthesis) using dct.
double* ykm_dct[MMAX+1];
struct DtDp* dylm_dct[MMAX+1];	// theta and phi derivative of Ylm matrix
struct DtDp* dykm_dct[MMAX+1];
double* zlm_dct[MMAX+1];		// matrix for direct transform (analysis)
struct DtDp* dzlm_dct[MMAX+1];
#endif

fftw_plan ifft, fft;	// plans for FFTW.
fftw_plan idct, dct, dctm0;
fftw_plan idctm, dctm;	// time_sht special.
unsigned fftw_plan_mode = FFTW_PATIENT;		// defines the default FFTW planner mode.

/*
// compute non-linear terms in m-spectral space and l-physical space, with FFTW conventions, ie :
// y = x0 + sum(m=1..MMAX) [ xm.exp(i.m.phi) + xm*.exp(-i.m.phi) ]
// this should be faster for small m's (MMAX <= 8), and requires a legendre transform before.
// No aliasing problems.
void NLspec(complex double *x, complex double *y, complex double *nl)
{
	long int m,i,l;

	m=0;	// 2MMAX+1 terms
		for (l=0;l<NLAT;l++)
			nl[m*NLAT+l] = creal(x[0*NLAT+l]) * creal(y[0*NLAT+l]);	// 1 term
		for (i=1; i<=MMAX; i++) {		// 2(MMAX)
			for (l=0;l<NLAT;l++)
//				(double) nl[m*NLAT+l] += 2* creal(x[i*NLAT+l]*conj(y[i*NLAT+l]));	//xi*conj(yi) + conj(xi)*yi
				(double) nl[m*NLAT+l] += 2*( creal(x[i*NLAT+l])*creal(y[i*NLAT+l]) + cimag(x[i*NLAT+l])*cimag(y[i*NLAT+l]) );
		}
	for (m=1; m<=MMAX; m++) {	// total 2(MMAX)+1-m terms
		for (l=0;l<NLAT;l++)
			nl[m*NLAT+l] = creal(x[0*NLAT+l])*y[m*NLAT+l] + x[m*NLAT+l]*creal(y[0*NLAT+l]);	// 2 terms
		for (i=1; i<m; i++) {		// 3(m-1) terms
			for (l=0;l<NLAT;l++)
				nl[m*NLAT+l] += conj(x[i*NLAT+l])*y[(i+m)*NLAT+l] + x[(i+m)*NLAT+l]*conj(y[i*NLAT+l]) + x[i*NLAT+l]*y[(m-i)*NLAT+l];
		}
		for (i=m; i<=MMAX-m; i++) {	// 2(MMAX-2m+1) terms
			for (l=0;l<NLAT;l++)
				nl[m*NLAT+l] += conj(x[i*NLAT+l])*y[(i+m)*NLAT+l] + x[(i+m)*NLAT+l]*conj(y[i*NLAT+l]);
		}
	}
}
*/

/// two pseudo-fft functions for compatibility of data with FFTW and MMAX=0
inline void ifft_m0_c2r(complex double *BrF, double *Br)		// in-place only
{
	long int i;
	for (i=1; i<NLAT; i++)		// zero is already in-place.
//		Br[i] = creal(BrF[i]);
		Br[i] = Br[2*i];
}

inline void fft_m0_r2c(double *Br, complex double *BrF)		//in-place only
{
	long int i;

/* algorithm with temp array */
	double tmp[NLAT];

	for (i=1; i<NLAT; i++)
		tmp[i] = Br[i];		// copy to temp array
	for (i=1; i<NLAT; i++)
//		(double) BrF[i] = tmp[i];	// copy back, interleaved.	*/
		Br[2*i] = tmp[i];	// copy back, interleaved.	*/

/* simple reverse loop algorithm (BAD performance on x86)
	for (i=NLAT-1; i>0; i--)		// zero is already in-place.
		(double) BrF[i] = Br[i];		// backward loop sucks...	*/

/* first upper half can be moved safely in upward direction. (improves things a bit)
	for (i=(NLAT+1)/2; i<NLAT; i++)
		(double) BrF[i] = Br[i];	// upper half can be moved safely.
	for (i=(NLAT-1)/2; i>0; i--)
		(double) BrF[i] = Br[i];		// backward loop sucks...	*/

/* upper halves are moved in upward direction, recursiveley. (may improve, or may not...)
	long int i0 = NLAT;
	while(i0 > 1) {
		//printf("i0=%d ",i0);
		for (i=(i0+1)/2; i<i0; i++)
			(double) BrF[i] = Br[i];	// upper half can be moved safely.
		i0 = (i0+1)/2;
	}	//	*/
}

inline void fft_m0_r2eo(double *Br, double *reo)
{
	long int i;
		for (i=0; i<NLAT/2; i++) {	// compute symmetric and antisymmetric parts.
			reo[2*i]   = Br[i] + Br[NLAT-(i+1)];
			reo[2*i+1] = Br[i] - Br[NLAT-(i+1)];
		}
//		if (i < NLAT_2) {	// NLAT is odd : special equator handling
		if (NLAT & 1) {		// NLAT is odd : special equator handling
			reo[2*i] = Br[i];	reo[2*i+1] = 0.0;
		}
}

/**
	SHT FUNCTIONS
**/

// run-time truncation at LMAX and MMAX
#ifndef LTR
  #define LTR LMAX
#endif
#ifndef MTR
  #define MTR MMAX
#endif
/////////////////////////////////////////////////////
//   Scalar Spherical Harmonics Transform
// input  : ShF = spatial/fourrier data : complex double array of size NLAT*(NPHI/2+1) or double array of size NLAT*(NPHI/2+1)*2
// output : Slm = spherical harmonics coefficients : complex double array of size NLM

void spat_to_SH(complex double *BrF, complex double *Qlm)
{
	#ifndef _SHT_EO_
		#include "SHT/spat_to_SH.c"
	#else
		#include "SHT/spat_to_SHo.c"
	#endif
}

#ifdef SHT_DCT
void spat_to_SH_dct(complex double *BrF, complex double *Qlm)
{
	complex double *Sl;		// virtual pointers for given im
	double *zl;
	long int k,im,m,l;
  #if NPHI > 1
  #else
  #endif

  #if NPHI > 1
	#define BR0 BrF

	fftw_execute_dft_r2c(fft,(double *) BrF, BrF);
	if (MRES & 1) {		// odd m's are present
		for (im=1; im<=MMAX; im+=2) {
			for (k=0; k<NLAT; k++)	BrF[im*NLAT + k] *= st[k];	// corection beforce DCT
		}
	}
  #else
	#define BR0 ((double *)BrF)
//	fft_m0_r2c((double *) BrF, BrF);
  #endif
	fftw_execute_r2r(dct,(double *) BrF, (double *) BrF);		// DCT

	im=0;	m=0;
		Sl = &Qlm[LiM(0,im)];		// virtual pointer for l=0 and im
		zl = zlm_dct[im];
		for (l=m; l<LMAX; l+=2) {		// l has parity of m
			Sl[l] = 0.0;	Sl[l+1] = 0.0;
			for (k=l; k<NLAT; k+=2) {		// for m=0, zl coeff with k<l are zeros.
				(double) Sl[l]   += (double) BR0[k]   * zl[k];
				(double) Sl[l+1] += (double) BR0[k+1] * zl[k+1];
			}
			zl += NLAT;
		}
		if ((LMAX & 1) == 0) {	// if (l == LMAX)  <=>  if ((LMAX & 1) == 0) for m=0
			Sl[l] = 0.0;
			for (k=l; k<NLAT; k+=2) {		// for m=0, DCT coeff with k<l are zeros.
				(double) Sl[l]   += (double) BR0[k]   * zl[k];
			}
		}
		BrF += NLAT;
	for (im=1; im<=MMAX; im++) {
		m=im*MRES;
		Sl = &Qlm[LiM(0,im)];		// virtual pointer for l=0 and im
		zl = zlm_dct[im];
		for (l=m; l<LMAX; l+=2) {		// l has parity of m
			Sl[l] = 0.0;	Sl[l+1] = 0.0;
			for (k=0; k<NLAT; k+=2) {
				Sl[l]   += BrF[k]   * zl[k];
				Sl[l+1] += BrF[k+1] * zl[k+1];
			}
			zl += NLAT;
		}
		if (l == LMAX) {
			Sl[l] = 0.0;
			for (k=0; k<NLAT; k+=2) {
				Sl[l]   += BrF[k]   * zl[k];
			}
		}
		BrF += NLAT;
	}
}
#endif

/////////////////////////////////////////////////////
//   Scalar inverse Spherical Harmonics Transform
// input  : Qlm = spherical harmonics coefficients : complex double array of size NLM [unmodified]
// output : BrF = spatial/fourrier data : complex double array of size NLAT*(NPHI/2+1) or double array of size NLAT*(NPHI/2+1)*2

void SH_to_spat(complex double *Qlm, complex double *BrF)
{
	#ifndef _SHT_EO_
		#include "SHT/SH_to_spat.c"
	#else
		#include "SHT/SHo_to_spat.c"
	#endif
}

#ifdef SHT_DCT
void SH_to_spat_dct(complex double *Qlm, complex double *BrF)
{
  #include "SHT/dct_SH_to_spat.c"
}
#endif

//void SH_to_grad_spat(complex double *Slm, complex double *BtF, complex double *BpF)
// "grad_spat" and "sph_to_spat" are the same : so we alias them.
//#include "SHT/S_to_spat.c"
#define SH_to_grad_spat(S,Gt,Gp) SHsph_to_spat(S, Gt, Gp)

/////////////////////////////////////////////////////
//   Spheroidal/Toroidal to (theta,phi) components inverse Spherical Harmonics Transform
// input  : Slm,Tlm = spherical harmonics coefficients of Spheroidal and Toroidal scalars : 
//          complex double array of size NLM [unmodified]
// output : BtF, BpF = theta, and phi vector components, spatial/fourrier data : 
//          complex double array of size NLAT*(NPHI/2+1) or double array of size NLAT*(NPHI/2+1)*2
void SHsphtor_to_spat(complex double *Slm, complex double *Tlm, complex double *BtF, complex double *BpF)
{
#ifndef _SHT_EO_
  #include "SHT/SHst_to_spat.c"
#else
  #include "SHT/SHost_to_spat.c"
#endif
}

#ifdef SHT_DCT
void SHsphtor_to_spat_dct(complex double *Slm, complex double *Tlm, complex double *BtF, complex double *BpF)
{
  #include "SHT/dct_SHst_to_spat.c"
}
#endif

void SHsph_to_spat(complex double *Slm, complex double *BtF, complex double *BpF)
{
#include "SHT/SHs_to_spat.c"
}

void SHtor_to_spat(complex double *Tlm, complex double *BtF, complex double *BpF)
{
#include "SHT/SHt_to_spat.c"
}

void spat_to_SHsphtor(complex double *BtF, complex double *BpF, complex double *Slm, complex double *Tlm)
{
#ifndef _SHT_EO_
  #include "SHT/spat_to_SHst.c"
#else
  #include "SHT/spat_to_SHost.c"
#endif
}

double SH_to_point(complex double *Qlm, double cost, double phi)
{
	double yl[LMAX+1];
	complex double eimp;
	double vr;
	complex double *Ql;
	long int l,m,im;
	
	vr = 0.0;
	m=0;
		gsl_sf_legendre_sphPlm_array(LTR, m, cost, &yl[m]);
		for (l=m; l<=LTR; l++)
			vr += yl[l] * creal( Qlm[l] );
	for (im=1; im<=MTR; im++) {
		m = im*MRES;
		gsl_sf_legendre_sphPlm_array(LTR, m, cost, &yl[m]);
		eimp = 2.*(cos(m*phi) + I*sin(m*phi));
		Ql = &Qlm[LiM(0,im)];	// virtual pointer for l=0 and im
		for (l=m; l<=LTR; l++)
			vr += yl[l] * creal( Ql[l]*eimp );
	}
	return vr;
}

void SHqst_to_point(complex double *Qlm, complex double *Slm, complex double *Tlm, double cost, double phi,
					   double *vr, double *vt, double *vp)
{
	double yl[LMAX+1];
	double dtyl[LMAX+1];
	complex double eimp;
	double sint, vst, vtt, vsp, vtp;
	complex double *Ql, *Sl, *Tl;
	long int l,m,im;
	
	sint = sqrt(1.0 - cost*cost);
	vst = 0.; vtt = 0.; vsp = 0.; vtp =0.; *vr = 0.;
	m=0;
		gsl_sf_legendre_sphPlm_deriv_array(LTR, m, cost, &yl[m], &dtyl[m]);
		for (l=m; l<=LTR; l++) {
			*vr += yl[l] * creal( Qlm[l] );
			vst += dtyl[l] * creal( Slm[l] );
			vtt += dtyl[l] * creal( Tlm[l] );
		}
	for (im=1; im<=MTR; im++) {
		m = im*MRES;
		gsl_sf_legendre_sphPlm_deriv_array(LTR, m, cost, &yl[m], &dtyl[m]);
		eimp = 2.*(cos(m*phi) + I*sin(m*phi));
		Ql = &Qlm[LiM(0,im)];	Sl = &Slm[LiM(0,im)];	Tl = &Tlm[LiM(0,im)];
		for (l=m; l<=LTR; l++) {
			*vr += yl[l] * creal( Ql[l]*eimp );
			vst += dtyl[l] * creal(Sl[l]*eimp);
			vtt += dtyl[l] * creal(Tl[l]*eimp);
			vsp += (yl[l] *m) * creal(I*Sl[l]*eimp);
			vtp += (yl[l] *m) * creal(I*Tl[l]*eimp);
		}
	}
	*vt = vtp/sint + (-sint*vst);	// Bt = I.m/sint *T  + dS/dt
	*vp = vsp/sint - (-sint*vtt);	// Bp = I.m/sint *S  - dT/dt
}


/**
	INITIALIZATION FUNCTIONS
**/

void runerr(const char * error_text)
{
	printf("*** Run-time error : %s\n",error_text);
	exit(1);
}

// compute number of modes for spherical harmonic description.
inline long int nlm_calc(long int lmax, long int mmax, long int mres)
{
	if (mmax*mres > lmax) mmax = lmax/mres;
	return( (mmax+1)*(lmax+1) - mres*(mmax*(mmax+1))/2 );	// this is wrong if lmax < mmax*mres
/*
	long int im,l,lm;
	lm=0;
	for (im=0; im<=mmax; im++) {	// perform the loop : this gives the exact result !
		for (l=im*mres; l<=lmax; l++)	lm++;
//		if (im*mres <= lmax) lm += lmax+1-im*mres;
	}
	return lm;
*/
}

void EqualSpaceNodes(double *x, double *w, int n)
{
	int j;
	long double f;
	long double pi = M_PI;

// cos theta of latidunal points (equaly spaced in theta)
#ifndef _SHT_EO_
	f = pi/(n-1.0);
#else
	f = pi/(2.*n-1.0);
#endif
	for (j=0; j<n; j++) {
		x[j] = cosl(f*j);
	}

// weights
	for (j=0; j<n; j++) {
		w[j] = 0.0;	// unsupported yet...
	}
	printf("          ! warning : only synthesis (inverse transform) supported so far for this grid !\n");
}

// Generates the abscissa and weights for a Gauss-Legendre quadrature.
// Newton method from initial Guess to find the zeros of the Legendre Polynome
// x = abscissa, w = weights, n points.
// Reference:  Numerical Recipes, Cornell press.
void GaussNodes(long double *x, long double *w, int n)
{
	long double z, z1, p1, p2, p3, pp, eps;
	long int i,j,m;
	long double pi = M_PI;

//	eps = 1.0e-15;	// desired precision, minimum = 2.2204e-16 (double)
	eps = 1.1e-19;	// desired precision, minimum = 1.0842e-19 (long double)

#ifdef _SHT_EO_
	n *=2;
	printf(" Even/odd separation, n=%d\n",n);
#endif
	m = (n+1)/2;
	for (i=1;i<=m;i++) {
		z = cosl(pi*((long double)i-0.25)/((long double)n+0.5));
		z1 = z+1.;
		while ( fabsl(z-z1) > eps )
		{
			p1 = 1.0;
			p2 = 0.0;
			for(j=1;j<=n;j++) {
				p3 = p2;
				p2 = p1;
				p1 = ((2*j-1)*z*p2-(j-1)*p3)/j;	// The Legendre polynomial...
			}
			pp = ((long double)n)*(z*p1-p2)/(z*z-1.0);                       // ... and its derivative.
			z1 = z;
			z = z1-p1/pp;
		}
		x[i-1] = z;		// Build up the abscissas.
		w[i-1] = 2.0/((1-z*z)*(pp*pp));		// Build up the weights.
#ifndef _SHT_EO_
		x[n-i] = -z;
		w[n-i] = w[i-1];
#endif
	}

// as we started with initial guesses, we should check if the gauss points are actually unique.
	for (i=m-1; i>0; i--) {
		if (x[i] == x[i-1]) runerr("bad gauss points\n");
	}

#ifdef _SH_DEBUG_
// test integral to compute :
	z = 0;
	for (i=0;i<m;i++) {
		z += w[i]*x[i]*x[i];
	}
	printf("          Gauss quadrature for 3/2.x^2 = %Lg (should be 1.0) error = %Lg\n",z*3.,z*3.-1.0);
#endif
}

// initialize FFTs using FFTW. stride = NLAT, (contiguous l)
void planFFT()
{
	complex double *ShF;
	double *Sh;
	int nfft = NPHI;
	int ncplx = NPHI/2 +1;
	int nreal;
	int ndct = NLAT;
	fftw_r2r_kind r2r_kind;
	fftw_iodim dims, hdims[2];
	
	nreal = 2*ncplx;
	
// Allocate dummy Spatial Fields.
	ShF = (complex double *) fftw_malloc(ncplx * NLAT * sizeof(complex double));
	Sh = (double *) ShF;

 #if NPHI > 1
	printf("[FFTW] Mmax=%d, Nphi=%d\n",MMAX,NPHI);

	if (NPHI <= 2*MMAX) runerr("[FFTW] the sampling condition Nphi > 2*Mmax is not met.");
	if (NPHI < 3*MMAX) printf("       ! Warning : 2/3 rule for anti-aliasing not met !\n");
	
// IFFT : unnormalized.
	ifft = fftw_plan_many_dft_c2r(1, &nfft, NLAT, ShF, &ncplx, NLAT, 1, Sh, &nreal, NLAT, 1, fftw_plan_mode);
	if (ifft == NULL)
		runerr("[FFTW] ifft planning failed !");

// FFT : must be normalized.
	fft = fftw_plan_many_dft_r2c(1, &nfft, NLAT, Sh, &nreal, NLAT, 1, ShF, &ncplx, NLAT, 1, fftw_plan_mode);
	if (fft == NULL)
		runerr("[FFTW] fft planning failed !");
 #else
	printf("          => no fft required for NPHI=1.\n");
 #endif

#ifdef SHT_DCT
/* LATITUDINAL DCT (THETA) */
/*	r2r_kind = FFTW_REDFT10;
	dct = fftw_plan_many_r2r(1, &ndct, MMAX+1, Sh, &ndct, 2, 2*NLAT, Sh, &ndct, 2, 2*NLAT, &r2r_kind, fftw_plan_mode );
	dct2 = fftw_plan_many_r2r(1, &ndct, MMAX+1, Sh+1, &ndct, 2, 2*NLAT, Sh+1, &ndct, 2, 2*NLAT, &r2r_kind, fftw_plan_mode );
	r2r_kind = FFTW_REDFT01;
	idct = fftw_plan_many_r2r(1, &ndct, MMAX+1, Sh, &ndct, 2, 2*NLAT, Sh, &ndct, 2, 2*NLAT, &r2r_kind, fftw_plan_mode );
	idct2 = fftw_plan_many_r2r(1, &ndct, MMAX+1, Sh+1, &ndct, 2, 2*NLAT, Sh+1, &ndct, 2, 2*NLAT, &r2r_kind, fftw_plan_mode );
*/

	dims.n = NLAT;	dims.is = 2;	dims.os = 2;		// real and imaginary part.
	hdims[0].n = MMAX+1;	hdims[0].is = 2*NLAT; 	hdims[0].os = 2*NLAT;
	hdims[1].n = 2;			hdims[1].is = 1; 	hdims[1].os = 1;
  #if NPHI > 1
	r2r_kind = FFTW_REDFT10;
	dct = fftw_plan_guru_r2r(1, &dims, 2, hdims, Sh, Sh, &r2r_kind, fftw_plan_mode );
	r2r_kind = FFTW_REDFT01;
	idct = fftw_plan_guru_r2r(1, &dims, 2, hdims, Sh, Sh, &r2r_kind, fftw_plan_mode );
	if ((dct == NULL)||(idct == NULL))
		runerr("[FFTW] dct planning failed !");
/// M=0 DCT
	r2r_kind = FFTW_REDFT10;
	dctm0 = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, 1, NLAT, Sh, &ndct, 1, NLAT, &r2r_kind, fftw_plan_mode );
	if (dctm0 == NULL)
		runerr("[FFTW] dctm0 planning failed !");
  #else
	r2r_kind = FFTW_REDFT10;
	dct = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, 1, NLAT, Sh, &ndct, 1, NLAT, &r2r_kind, fftw_plan_mode );
	dctm0 = dct;
	r2r_kind = FFTW_REDFT01;
	idct = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, 1, NLAT, Sh, &ndct, 1, NLAT, &r2r_kind, fftw_plan_mode );
	if ((dct == NULL)||(idct == NULL))
		runerr("[FFTW] dct planning failed !");
  #endif

/// TIME SHT SPECIAL
	hdims[0].n = 1;
	r2r_kind = FFTW_REDFT10;
	dctm = fftw_plan_guru_r2r(1, &dims, 2, hdims, Sh, Sh, &r2r_kind, fftw_plan_mode );
	r2r_kind = FFTW_REDFT01;
	idctm = fftw_plan_guru_r2r(1, &dims, 2, hdims, Sh, Sh, &r2r_kind, fftw_plan_mode );
	if ((dctm == NULL)||(idctm == NULL))
		runerr("[FFTW] dct planning failed !");

/// TIME SHT SPECIAL

//	dct_norm = 1.0/(2*NLAT);
#endif

//	fft_norm = 1.0/nfft;
	fftw_free(ShF);
//	printf("       done.\n");
}

// Perform some optimization on the SHT matrices.
void OptimizeMatrices(double eps)
{
	double *yg;
	int im,m,l,it;

/// POLAR OPTIMIZATION : analyzing coefficients, some can be safely neglected.
	for (im=0;im<=MMAX;im++) {
		m = im*MRES;
		tm[im] = NLAT_2;
		for (l=m;l<=LMAX;l++) {
			it=0;
			while( fabs(ylm[im][it*(LMAX-m+1) + (l-m)]) < eps ) { it++; }
			if (tm[im] > it) tm[im] = it;
		}
	}
	if (eps > 0.0) {
		printf("          polar optimization threshold = %e\n",eps);
#ifdef _SH_DEBUG_
		printf("          tm[im]=");
		for (im=0;im<=MMAX;im++)
			printf(" %d",tm[im]);
		printf("\n");
#endif
	}

/// Compression of dylm and dzlm for m=0, as .p is 0
	im=0;	m=0;
		yg = (double *) dylm[im];
		for (it=0; it<NLAT_2; it++) {
			for (l=m; l<=LMAX; l++)
				yg[it*(LMAX-m+1) + (l-m)] = dylm[im][it*(LMAX-m+1) + (l-m)].t;
		}
		yg = (double *) dzlm[im];
		for (l=m; l<LMAX; l+=2) {
			for (it=0; it<NLAT_2; it++) {
				yg[(l-m)*NLAT_2 + it*2] = dzlm[im][(l-m)*NLAT_2 + it*2].t;
				yg[(l-m)*NLAT_2 + it*2+1] = dzlm[im][(l-m)*NLAT_2 + it*2+1].t;
			}
		}
		if (l==LMAX) {		// last l is stored right away, without interleaving.
			for (it=0; it<NLAT_2; it++) {
				yg[(l-m)*NLAT_2 + it] = dzlm[im][(l-m)*NLAT_2 + it].t;
			}
		}
}

/** initialize SH transform.
input : eps = polar optimization threshold : polar coefficients below that threshold are neglected (for high ms).
        eps is the value under wich the polar values of the Legendre Polynomials Plm are neglected, leading to increased performance (a few percent).
	0 = no polar optimization;  1.e-14 = VERY safe;  1.e-10 = safe;  1.e-6 = aggresive.
*/
#ifndef SHT_DCT
void init_SH(double eps)
{
	double dtylm[LMAX+1];		// temp storage for derivative : d(P_l^m(x))/dx
	double iylm_fft_norm = 2.0*M_PI/NPHI;	// normation FFT pour zlm
	double *dyl;
	double t,tmax;
	long int it,im,m,l;
	long double xg[NLAT], wg[NLAT];	// gauss points and weights.

#ifdef _SHT_EO_
	iylm_fft_norm *= 2.0;	// normation must be multiplied by 2.
#endif
	printf("[init_SH] Lmax=%d, Nlat=%d, Mres=%d, Mmax*Mres=%d, Nlm=%d\n",LMAX,NLAT,MRES,MMAX*MRES,NLM);
	if (MMAX*MRES > LMAX) runerr("[init_SH] MMAX*MRES should not exceed LMAX");
	if ((NLAT_2)*2 <= LMAX) runerr("[init_SH] NLAT_2*2 should be at least LMAX+1");
	if (MRES <= 0) runerr("[init_SH] MRES must be > 0");
	if (2*NLAT <= 3*LMAX) printf("          ! Warning : anti-aliasing condition in theta direction not met.\n");

#ifdef SHT_POLES
	printf("          => using Equaly Spaced Nodes including poles\n");
	EqualSpaceNodes(ct,wg,NLAT);		// equaly-spaced points and weights.
#else
	printf("          => using Gauss Nodes\n");
	GaussNodes(xg,wg,NLAT);	// generate gauss nodes and weights : ct = ]1,-1[ = cos(theta)
#endif
	for (it=0; it<NLAT; it++) {
		ct[it] = xg[it];
		st[it] = sqrtl(1.0 - xg[it]*xg[it]);
		st_1[it] = 1.0/sqrtl(1.0 - xg[it]*xg[it]);
	}

#ifdef _SH_DEBUG_
	printf(" NLAT=%d, NLAT_2=%d\n",NLAT,NLAT_2);
// TEST if gauss points are ok.
	tmax = 0.0;
	for (it = 0; it<NLAT_2; it++) {
		t = gsl_sf_legendre_sphPlm(NLAT, 0, ct[it]);
		if (t>tmax) tmax = t;
//		printf("i=%d, x=%12.12g, p=%12.12g\n",i,ct[i],t);
	}
	printf("          max zero at Gauss node for Plm[l=LMAX+1,m=0] : %g\n",tmax);
	if (NLAT_2 < 100) {
		printf("          Gauss nodes :");
		for (it=0;it<NLAT_2; it++)
			printf(" %g",ct[it]);
		printf("\n");
	}
#endif

// Allocate legendre functions lookup tables.
	ylm[0] = (double *) fftw_malloc(sizeof(double)* NLM*NLAT_2);
	dylm[0] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* NLM*NLAT_2);
	zlm[0] = (double *) fftw_malloc(sizeof(double)* NLM*NLAT_2);
	dzlm[0] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* NLM*NLAT_2);
	for (im=0; im<MMAX; im++) {
		m = im*MRES;
		ylm[im+1] = ylm[im] + NLAT_2*(LMAX+1 -m);
		dylm[im+1] = dylm[im] + NLAT_2*(LMAX+1 -m);
		zlm[im+1] = zlm[im] + NLAT_2*(LMAX+1 -m);
		dzlm[im+1] = dzlm[im] + NLAT_2*(LMAX+1 -m);
	}
#ifdef _SH_DEBUG_
	printf("          Memory used for Ylm and Zlm matrices = %.3f Mb\n",6.0*sizeof(double)*NLM*NLAT_2/(1024.*1024.));
#endif

// Even/Odd symmetry : ylm is even or odd across equator, as l-m is even or odd => only NLAT_2 points required.
// for synthesis (inverse transform)
	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
//		ylm[im] = (double *) fftw_malloc(sizeof(double)* (LMAX+1-m)*NLAT_2);
//		dylm[im] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* (LMAX+1-m)*NLAT_2);
		for (it=0; it<NLAT_2; it++) {
			if ((m==1)&&(st[it]==0.)) {		// gsl function diverges for m=1 and sint=0 => use approximation
				gsl_sf_legendre_sphPlm_array(LMAX, m, ct[it], ylm[im] + it*(LMAX-m+1));
				gsl_sf_legendre_sphPlm_array(LMAX, m, sqrt(sqrt(ct[it+1])), dtylm);
				for (l=m; l<=LMAX; l++) {		// d(Pl1)/dt |(t=0) = Pl1(epsilon)/sin(epsilon)
					dylm[im][it*(LMAX-m+1) + (l-m)].t = dtylm[l-m] / sqrt(1. - sqrt(ct[it+1]));
					dylm[im][it*(LMAX-m+1) + (l-m)].p = dylm[im][it*(LMAX-m+1) + (l-m)].t;
				}
			} else {
				gsl_sf_legendre_sphPlm_deriv_array(LMAX, m, ct[it], ylm[im] + it*(LMAX-m+1), dtylm);	// fixed im legendre functions lookup table.
				for (l=m; l<=LMAX; l++) {
					dylm[im][it*(LMAX-m+1) + (l-m)].t = -st[it] *dtylm[l-m];	// d(Plm(cost))/dt = -sin(t).d(Plm(x))/dx
					dylm[im][it*(LMAX-m+1) + (l-m)].p = ylm[im][it*(LMAX-m+1) + (l-m)] *m/st[it];	// 1/sint(t) dYlm/dphi
					if (st[it]==0.) dylm[im][it*(LMAX-m+1) + (l-m)].p = 0.0;
				}
			}
		}
	}
	
// for analysis (decomposition, direct transform) : transpose and multiply by gauss weight and other normalizations.
// interleave l and l+1 : this stores data in the way it will be read.
	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
//		zlm[im] = (double *) fftw_malloc(sizeof(double)* (LMAX+1-m)*NLAT_2);
//		dzlm[im] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* (LMAX+1-m)*NLAT_2);
		for (it=0;it<NLAT_2;it++) {
			for (l=m;l<LMAX;l+=2) {
				zlm[im][(l-m)*NLAT_2 + it*2]    =  ylm[im][it*(LMAX-m+1) + (l-m)]   * wg[it] *iylm_fft_norm;
				zlm[im][(l-m)*NLAT_2 + it*2 +1] =  ylm[im][it*(LMAX-m+1) + (l+1-m)] * wg[it] *iylm_fft_norm;
				dzlm[im][(l-m)*NLAT_2 + it*2].t = dylm[im][it*(LMAX-m+1) + (l-m)].t * wg[it] *iylm_fft_norm /(l*(l+1));
				dzlm[im][(l-m)*NLAT_2 + it*2].p = dylm[im][it*(LMAX-m+1) + (l-m)].p * wg[it] *iylm_fft_norm /(l*(l+1));
				dzlm[im][(l-m)*NLAT_2 + it*2+1].t = dylm[im][it*(LMAX-m+1) + (l+1-m)].t * wg[it] *iylm_fft_norm /((l+1)*(l+2));
				dzlm[im][(l-m)*NLAT_2 + it*2+1].p = dylm[im][it*(LMAX-m+1) + (l+1-m)].p * wg[it] *iylm_fft_norm /((l+1)*(l+2));
				if (l == 0) {		// les derivees sont nulles pour l=0 (=> m=0)
					dzlm[im][(l-m)*NLAT_2 + it*2].t = 0.0;
					dzlm[im][(l-m)*NLAT_2 + it*2].p = 0.0;
				}
			}
			if (l==LMAX) {		// last l is stored right away, without interleaving.
				zlm[im][(l-m)*NLAT_2 + it]    =  ylm[im][it*(LMAX-m+1) + (l-m)]   * wg[it] *iylm_fft_norm;
				dzlm[im][(l-m)*NLAT_2 + it].t = dylm[im][it*(LMAX-m+1) + (l-m)].t * wg[it] *iylm_fft_norm /(l*(l+1));
				dzlm[im][(l-m)*NLAT_2 + it].p = dylm[im][it*(LMAX-m+1) + (l-m)].p * wg[it] *iylm_fft_norm /(l*(l+1));
			}
		}
	}

	OptimizeMatrices(eps);

 #if NPHI > 1
	planFFT();		// initialize fftw
 #else
	printf("          => no fft required for NPHI=1.\n");
 #endif

// Additional arrays :
#if MRES == _mres_
	_nlm_ = nlm_calc(LMAX, MMAX, MRES);	// store "precomputed" NLM
	#undef NLM
	#define NLM _nlm_
	el = (double *) fftw_malloc(3*NLM * sizeof(double));	// NLM defined at runtime.
	l2 = el + NLM;	l_2 = el + 2*NLM;
#endif
	it = 0;
	for (im=0;im<=MMAX;im++) {
		for (l=im*MRES;l<=LMAX;l++) {
			el[it] = l;	l2[it] = l*(l+1.0);	l_2[it] = 1.0/(l*(l+1.0));
			it++;
		}
	}
	l_2[0] = 0.0;	// undefined for l=0 => replace with 0.
}

#else

void init_SH(double eps)
{
	double Z[2*NLAT_2], dZt[2*NLAT_2], dZp[2*NLAT_2];		// equally spaced theta points.
	double *yk, *yg, *dygt, *dygp;		// temp storage for Plm(xg)
	struct DtDp *dyg, *dyk;
	double dtylm[LMAX+1];		// temp storage for derivative : d(P_l^m(x))/dx
	double iylm_fft_norm = 2.0*M_PI/NPHI;	// normation FFT pour zlm
	long int it,im,m,l, lm, dlm;
	long double t,tsum;
	long double *cktg;		// temp storage for cos(k*tg);
	long double xg[NLAT], wg[NLAT], sg[NLAT_2];//, sg_1[NLAT_2], sg_2[NLAT_2];	// gauss points and weights.

#ifdef _SHT_EO_
	iylm_fft_norm *= 2.0;	// normation must be multiplied by 2.
#endif
	printf("[init_SH] Lmax=%d, Nlat=%d, Mres=%d, Mmax*Mres=%d, Nlm=%d\n",LMAX,NLAT,MRES,MMAX*MRES,NLM);
	if (MMAX*MRES > LMAX) runerr("[init_SH] MMAX*MRES should not exceed LMAX");
	if ((NLAT_2)*2 <= LMAX+1) runerr("[init_SH] NLAT_2*2 should be at least LMAX+2 (DCT)");
	if (MRES <= 0) runerr("[init_SH] MRES must be > 0");
	if (NLAT & 1) runerr("[init_SH] NLAT must be even");
	if (2*NLAT <= 3*LMAX) printf("          ! Warning : anti-aliasing condition in theta direction not met.\n");
	
	printf("          => using Equaly Spaced Nodes with DCT acceleration\n");
	for (it=0; it<NLAT; it++) {	// Chebychev points : equaly spaced but skipping poles.
		long double th = M_PI*(it+0.5)/NLAT;
		ct[it] = cosl(th);
		st[it] = sinl(th);
		st_1[it] = 1.0/sinl(th);
	}

#ifdef _SH_DEBUG_
	printf(" NLAT=%d, NLAT_2=%d\n",NLAT,NLAT_2);
	if (NLAT_2 < 100) {
		printf("          DCT nodes :");
		tsum = 0.0;
		for (it=0;it<NLAT_2; it++) {
			printf(" %g",ct[it]);
			t = fabs(ct[it]*ct[it] + st[it]*st[it] -1.0);
			if (t > tsum) tsum=t;
		}
	}
	printf("\n max st^2 + ct^2 -1 = %lg\n",tsum);
#endif

/// Allocate legendre functions lookup tables.
	ylm[0] = (double *) fftw_malloc(sizeof(double)* NLM*NLAT_2);
	dylm[0] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* NLM*NLAT_2);
	zlm[0] = (double *) fftw_malloc(sizeof(double)* NLM*NLAT_2);
	dzlm[0] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* NLM*NLAT_2);
	for (im=0; im<MMAX; im++) {
		m = im*MRES;
		ylm[im+1] = ylm[im] + NLAT_2*(LMAX+1 -m);
		dylm[im+1] = dylm[im] + NLAT_2*(LMAX+1 -m);
		zlm[im+1] = zlm[im] + NLAT_2*(LMAX+1 -m);
		dzlm[im+1] = dzlm[im] + NLAT_2*(LMAX+1 -m);
	}
	for(im=0, lm=0, dlm=0; im<=MMAX; im++) {	// how much memory to allocate for ylm dct
		m = im*MRES;
		for (l=m;l<=LMAX;l+=2) {
			lm += l+2;	//(l+2 - (m&1));
			dlm += l+2;
		}
	}
#ifdef _SH_DEBUG_
	printf("          Memory used for Ylm and Zlm matrices = %.3f Mb\n",6.0*sizeof(double)*NLM*NLAT_2/(1024.*1024.));
	printf("          Memory used for Ylm_dct matrices = %.3f Mb\n",sizeof(double)*(lm + 2*dlm)/(1024.*1024.));
#endif
	ylm_dct[0] = (double *) fftw_malloc(sizeof(double)* lm);
	dylm_dct[0] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* dlm);
	zlm_dct[0] = (double *) fftw_malloc( sizeof(double)* NLM*(2*NLAT_2 +1) );	// quantite a revoir...
	ykm_dct[0] = (double *) fftw_malloc(sizeof(double)* NLM*(LMAX+1));
	dykm_dct[0] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* NLM*(LMAX+1));
	for (im=0; im<MMAX; im++) {
		m = im*MRES;
		for (l=m, lm=0, dlm=0; l<=LMAX; l+=2) {
			lm += l+2;	//(l+2 - (m&1));
			dlm += l+2;
		}
		ylm_dct[im+1] = ylm_dct[im] + lm;
		dylm_dct[im+1] = dylm_dct[im] + dlm;
		zlm_dct[im+1] = zlm_dct[im] + ((LMAX-m+2)/2)*NLAT;
		ykm_dct[im+1] = ykm_dct[im] + (LMAX+1)*(LMAX+1-m);
		dykm_dct[im+1] = dykm_dct[im] + (LMAX+1)*(LMAX+1-m);
	}

	dct = fftw_plan_r2r_1d( 2*NLAT_2, Z, Z, FFTW_REDFT10, FFTW_ESTIMATE );	// quick and dirty dct.
	idct = fftw_plan_r2r_1d( 2*NLAT_2, Z, Z, FFTW_REDFT01, FFTW_ESTIMATE );	// quick and dirty idct.

// Even/Odd symmetry : ylm is even or odd across equator, as l-m is even or odd => only NLAT_2 points required.
/// for synthesis (inverse transform)
	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
		for (it=0; it<NLAT_2; it++) {
			gsl_sf_legendre_sphPlm_deriv_array(LMAX, m, ct[it], ylm[im] + it*(LMAX-m+1), dtylm);	// fixed im legendre functions lookup table.
			for (l=m; l<=LMAX; l++) {
				dylm[im][it*(LMAX-m+1) + (l-m)].t = - dtylm[l-m];	// d(Plm(cos(t)))/dt = -sin(t) d(Plm(x))/dx
				dylm[im][it*(LMAX-m+1) + (l-m)].p = ylm[im][it*(LMAX-m+1) + (l-m)] *m/st[it];	// 1/sint(t) dYlm/dphi
				if (st[it]==0.) dylm[im][it*(LMAX-m+1) + (l-m)].p = 0.0;
			}
		}
	// go to DCT space
		yg = ylm_dct[im];
		dyg = dylm_dct[im];
		yk = ykm_dct[im];
		dyk = dykm_dct[im];
		for (it=0;it<=LMAX;it+=2) {
			for(l=m; l<=LMAX; l++) {
				yk[(it/2)*(LMAX+1-m) + (l-m)] = 0.0;
				dyk[(it/2)*(LMAX+1-m) + (l-m)].t = 0.0;
				dyk[(it/2)*(LMAX+1-m) + (l-m)].p = 0.0;
			}
		}
		for (l=m; l<=LMAX; l++) {
			if (m & 1) {	// m odd
				for (it=0; it<NLAT_2; it++) {
					Z[it] = ylm[im][it*(LMAX-m+1) + (l-m)] / st[it];	// P[l-1](x)	/st
					dZt[it] = dylm[im][it*(LMAX-m+1) + (l-m)].t * st[it];	// P[l](x)	*1
					dZp[it] = dylm[im][it*(LMAX-m+1) + (l-m)].p;		// P[l-1](x)	*1
				}
			} else {	// m even
				for (it=0; it<NLAT_2; it++) {
					Z[it] = ylm[im][it*(LMAX-m+1) + (l-m)];		// P[l](x)	*1
					dZt[it] = dylm[im][it*(LMAX-m+1) + (l-m)].t;	// P[l-1](x)	/st
//					dZp[it] = dylm[im][it*(LMAX-m+1) + (l-m)].p / st[it];	// P[l-2](x)	/st
					dZp[it] = ylm[im][it*(LMAX-m+1) + (l-m)] * m/(st[it]*st[it]);	// P[l-2](x)	/st
//					dZp[it] = ylm[im][it*(LMAX-m+1) + (l-m)] * m;	// P[l](x)	*st
				}
			}
			if ((l-m)&1) {	// odd
				for (it=NLAT_2; it<2*NLAT_2; it++) {
					Z[it]   = - Z[2*NLAT_2-it-1];	// reconstruct even/odd part
					dZt[it] =   dZt[2*NLAT_2-it-1];
					dZp[it] = - dZp[2*NLAT_2-it-1];
				}
			} else {	// even
				for (it=NLAT_2; it<2*NLAT_2; it++) {
					Z[it] =     Z[2*NLAT_2-it-1];	// reconstruct even/odd part
					dZt[it] = - dZt[2*NLAT_2-it-1];
					dZp[it] =   dZp[2*NLAT_2-it-1];
				}
			}
			fftw_execute(dct);
			fftw_execute_r2r(dct, dZt, dZt);
			fftw_execute_r2r(dct, dZp, dZp);
#ifdef _SH_DEBUG_
			if (LMAX <= 12) {
				printf("\nl=%d, m=%d ::\t", l,m);
				for(it=0;it<2*NLAT_2;it++) printf("%f ",Z[it]);
				printf("\n     dYt ::\t", l,m);
				for(it=0;it<2*NLAT_2;it++) printf("%f ",dZt[it]);
				printf("\n     dYp ::\t", l,m);
				for(it=0;it<2*NLAT_2;it++) printf("%f ",dZp[it]);
			}
#endif
			for (it=(l-m)&1; it<=l; it+=2) {
				yg[it] = Z[it]/(2*NLAT);	// store non-zero coeffs.
				dyg[it].p = dZp[it]/(2*NLAT);
				yk[(it/2)*(LMAX+1-m) + (l-m)] = Z[it]/(2*NLAT);	// and transpose
				dyk[(it/2)*(LMAX+1-m) + (l-m)].p = dZp[it]/(2*NLAT);
			}
			for (it=(l+1-m)&1; it<=l; it+=2) {
				dyg[it].t = dZt[it]/(2*NLAT);	// store non-zero coeffs.
				dyk[(it/2)*(LMAX+1-m) + (l-m)].t = dZt[it]/(2*NLAT);
			}
			if ((l-m)&1) {
				yg += l+1;	//(l+1 - (m&1));	// l-m odd : go to next line of storage.
				dyg += l+1;
			}
		}
		for (it=0; it<NLAT_2; it++) {		// do corrections for dylm.t non-DCT array.
			for (l=m; l<=LMAX; l++) {
				dylm[im][it*(LMAX-m+1) + (l-m)].t *= st[it];	// d(Plm(cos(t)))/dt = -sin(t) d(Plm(x))/dx
			}
		}

		// Compact the coefficients for improved cache efficiency.
		yg = ykm_dct[im];
		dyg = dykm_dct[im];
		for (it=0; it<= LMAX; it+=2) {
			l = (it < m) ? m : it+(m&1);
			while (l<=LMAX) {
				yg[0] = yk[(it/2)*(LMAX+1-m) + (l-m)];
				l++;	yg++;
			}
			l = (it < m) ? m : it-(m&1);
			while (l<=LMAX) {
				dyg[0].t = dyk[(it/2)*(LMAX+1-m) + (l-m)].t;
				dyg[0].p = dyk[(it/2)*(LMAX+1-m) + (l-m)].p;
				l++;	dyg++;

			}
		}
		if (im == 0) {		// compact m=0 dylm because .p = 0 :
			dyg = dykm_dct[im];
			yg = (double *) dykm_dct[im];
			for (it=0; it< LMAX; it+=2) {
				for (l=it; l<=LMAX; l++) {
					yg[0] = dyg[0].t;
					yg++;	dyg++;
				}
			}
			dyg = dylm_dct[im];
			yg = (double *) dylm_dct[im];
			for (l=0; l<=LMAX; l+=2) {
				for (it=0; it<=l; it+=2) {
					yg[it] = dyg[it].t;
					yg[it+1] = dyg[it+1].t;
				}
				yg += l+2;	dyg += l+2;
			}
		}
	}
/// for analysis (decomposition, direct transform) : use gauss-legendre quadrature for dct components
	cktg = (long double *) malloc( sizeof(long double) * (2*NLAT_2)*NLAT_2);
	GaussNodes(xg,wg,NLAT);	// generate gauss nodes and weights : xg = ]1,-1[ = cos(theta)
	for (it=0; it<NLAT_2; it++) {
//		long double thg = acosl(xg[it]);
		sg[it] = sqrtl(1. - xg[it]*xg[it]);
//		sg[it] = sinl(thg);
//		sg_1[it] = 1./sqrtl(1. - xg[it]*xg[it]);
//		sg_2[it] = 1./(1. - xg[it]*xg[it]);
//		for (l=0; l<NLAT; l++) cktg[l*NLAT_2 + it] = cosl(thg*l) *wg[it];	// precompute
		long double T0,T1,T2;		// precompute cos(k*tg) using Chebychev recurrence
		T0 = 1.0;	T1 = xg[it];
		cktg[0*NLAT_2 + it] = T0*wg[it];
		cktg[1*NLAT_2 + it] = T1*wg[it];
		for (l=2; l<(2*NLAT_2); l++) {
			T2 = 2.*xg[it]*T1 - T0;
			cktg[l*NLAT_2 + it] = T2*wg[it];
			T0 = T1;	T1 = T2;
		}
	}
#ifdef _SH_DEBUG_
	for (it=0, tsum=0.0; it<NLAT_2; it++) {
		t = fabsl(xg[it]*xg[it] + sg[it]*sg[it] -1.0);
		if (t > tsum) tsum=t;
	}
	printf("          Gauss nodes : max st^2 + ct^2 -1 = %Lg\n",tsum);
	for (it=0, tsum=0.0; it<NLAT_2; it++) {
		t = fabsl( sg[it] - sinl(acosl(xg[it])) );
		if (t > tsum) tsum=t;
	}
	printf("          Gauss nodes : max sg - sin(tg) = %Lg\n",tsum);
	for (it=0, tsum=0.0; it<NLAT_2; it++) {
		t = fabsl(cktg[(NLAT-1)*NLAT_2 +it] - cosl(acosl(xg[it])*(NLAT-1))*wg[it]);
//		printf("%Lg  ",t);
		if (t > tsum) tsum=t;
	}
	printf("          Gauss nodes : max of Tn(x)-cos(n*x) = %Lg\n",tsum);
#endif

// we need the legendre functions lookup tables for gauss points also !
	yg = (double *) malloc( sizeof(double) * 3*(LMAX+1)*NLAT_2);
	dygt = yg + (LMAX+1)*NLAT_2;
	dygp = yg + 2*(LMAX+1)*NLAT_2;

	inline void calc_Zlm_dct(int l,int m)
	{
		int k0,k1, k,it;
		long double sum, dtsum, dpsum;
		double min,max;
		min=1e32;	max=0.0;

		k0 = (l-m)&1;	k1 = 1-k0;
		for (k=0; k<(2*NLAT_2); k++) {  Z[k] = 0.0;  dZt[k] = 0.0;  dZp[k] = 0.0; }
		for (k=k0; k<(2*NLAT_2); k+=2) {
			sum = 0.0;	dpsum = 0.0;
			for (it=0; it<NLAT_2; it++) {
				  sum += cktg[k*NLAT_2 + it] *   yg[it*(LMAX+1) + (l-m)];
				dpsum += cktg[k*NLAT_2 + it] * dygp[it*(LMAX+1) + (l-m)];
				if (fabs(cktg[k*NLAT_2 + it] * yg[it*(LMAX+1) + (l-m)]) < min) min = fabs(cktg[k*NLAT_2 + it] * yg[it*(LMAX+1) + (l-m)]);
				if (fabs(cktg[k*NLAT_2 + it] * yg[it*(LMAX+1) + (l-m)]) > max) max = fabs(cktg[k*NLAT_2 + it] * yg[it*(LMAX+1) + (l-m)]);
			}
			Z[k] = sum * iylm_fft_norm/NLAT_2;
			dZp[k] = dpsum * iylm_fft_norm/(NLAT_2 *l*(l+1));
			if (l==0) { dZp[k] = 0.0; }
		}
#ifdef _SH_DEBUG_
		if (max/min > 1.e14) {
			printf("\nl=%d, m=%d :: min=%g, max=%g, ratio=%g\t",l,m,min,max,max/min);
		}
#endif
		for (k=k1; k<(2*NLAT_2); k+=2) {
			dtsum = 0.0;
			for (it=0; it<NLAT_2; it++) {
				dtsum += cktg[k*NLAT_2 + it] * dygt[it*(LMAX+1) + (l-m)];
				if (fabsl(cktg[k*NLAT_2 + it] * dygt[it*(LMAX+1) + (l-m)]) < min) min = fabs(cktg[k*NLAT_2 + it] * dygt[it*(LMAX+1) + (l-m)]);
				if (fabsl(cktg[k*NLAT_2 + it] * dygt[it*(LMAX+1) + (l-m)]) > max) max = fabs(cktg[k*NLAT_2 + it] * dygt[it*(LMAX+1) + (l-m)]);
			}
			dZt[k] = dtsum * iylm_fft_norm/(NLAT_2 *l*(l+1));
			if (l==0) { dZt[k] = 0.0; }
		}
#ifdef _SH_DEBUG_
		if (max/min > 1.e14) {
			printf("\nl=%d, m=%d :: (d/dt) min=%g, max=%g, ratio=%g\t",l,m,min,max,max/min);
		}
		if (LMAX <= 12) {
			printf("\nl=%d, m=%d ::\t",l,m);
			for (k=0; k<(2*NLAT_2); k++) printf("%f ",Z[k]);
			printf("\n       dZt ::\t");
			for (k=0; k<(2*NLAT_2); k++) printf("%f ",dZt[k]);
			printf("\n       dZp ::\t");
			for (k=0; k<(2*NLAT_2); k++) printf("%f ",dZp[k]);
		}
		for (k=0, sum=0.0; k<l; k++) if (Z[k]*Z[k] > sum) sum = Z[k]*Z[k];
//		printf("\nmax Z[k] for k<l is : %g",sqrt(sum));
#endif
		if (k0==0) {
			for (k=0;k<(2*NLAT_2);k++) zlm_dct[m/MRES][((l-m)>>1)*(2*NLAT_2) +k] = 0.0;
		}
		for (k=k0; k<(2*NLAT_2); k+=2) {
			if (k == 0) {
				zlm_dct[m/MRES][((l-m)>>1)*(2*NLAT_2) +k] = Z[k]*0.5;		// store zlm_dct
			} else {
				zlm_dct[m/MRES][((l-m)>>1)*(2*NLAT_2) +k] = Z[k];		// store zlm_dct
			}
		}
		// prepare idct :
		fftw_execute_r2r(idct, Z, Z);	fftw_execute_r2r(idct, dZt, dZt);	fftw_execute_r2r(idct, dZp, dZp);
		if (m & 1) {	// m odd
			for (it=0; it<NLAT; it++) Z[it] *= st[it];
		} else {	// m even
			for (it=0; it<NLAT; it++) { dZt[it] *= st[it];		dZp[it] *= st[it]; }
		}
	}

	// zlm in DCT space
	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
		printf("computing weights m=%d\r",m);	fflush(stdout);
		for (it=0; it<NLAT_2; it++) {	// compute Plm's at gauss nodes.
			gsl_sf_legendre_sphPlm_deriv_array( LMAX, m, xg[it], &yg[it*(LMAX+1)], &dygt[it*(LMAX+1)] );
			if (m & 1) {	// m odd
				for (l=m; l<=LMAX; l++) {
					dygp[it*(LMAX+1) + (l-m)] = m * yg[it*(LMAX+1) + (l-m)] / sg[it];
					yg[it*(LMAX+1) + (l-m)] /= sg[it];	// Plm/sin(t) = P[l-1](cost)
					dygt[it*(LMAX+1) + (l-m)] *= -sg[it];	// -(dPlm/dx)*sin(t) = dPlm/dt = P[l](cost)
				}
			} else {	// m even
				for (l=m; l<=LMAX; l++) {
					// Plm = P[l](cost)
					dygt[it*(LMAX+1) + (l-m)] *= -1.;	// -dPlm/dx = P[l-1](cost) = 1/sint.dPlm/dt
//					dygp[it*(LMAX+1) + (l-m)] = m * yg[it*(LMAX+1) + (l-m)] * sg_2[it];	// P[l-2](cost)
					dygp[it*(LMAX+1) + (l-m)] = m * yg[it*(LMAX+1) + (l-m)] /(sg[it]*sg[it]);	// P[l-2](cost)
//					dygp[it*(LMAX+1) + (l-m)] = m * yg[it*(LMAX+1) + (l-m)];	// P[l](cost)
				}
			}
		}

		for (l=m;l<LMAX;l+=2) {
			calc_Zlm_dct(l,m);
			for (it=0; it<NLAT_2; it++) {
				zlm[im][(l-m)*NLAT_2 + it*2] = Z[it];
				dzlm[im][(l-m)*NLAT_2 + it*2].p = dZp[it];
				dzlm[im][(l-m)*NLAT_2 + it*2].t = dZt[it];
			}
			calc_Zlm_dct(l+1,m);
			for (it=0; it<NLAT_2; it++) {
				zlm[im][(l-m)*NLAT_2 + it*2 +1] = Z[it];
				dzlm[im][(l-m)*NLAT_2 + it*2 +1].p = dZp[it];
				dzlm[im][(l-m)*NLAT_2 + it*2 +1].t = dZt[it];
			}
		}
		if (l==LMAX) {
			calc_Zlm_dct(l,m);
			for (it=0; it<NLAT_2; it++) {
				zlm[im][(l-m)*NLAT_2 + it] =  Z[it];
				dzlm[im][(l-m)*NLAT_2 + it].p = dZp[it];
				dzlm[im][(l-m)*NLAT_2 + it].t = dZt[it];
			}
		}
	}
	free(yg);	free(cktg);
	fftw_destroy_plan(idct);	fftw_destroy_plan(dct);

	OptimizeMatrices(eps);

	planFFT();		// initialize fftw

// Additional arrays :
#if MRES == _mres_
	_nlm_ = nlm_calc(LMAX, MMAX, MRES);	// store "precomputed" NLM
	#undef NLM
	#define NLM _nlm_
	el = (double *) fftw_malloc(3*NLM * sizeof(double));	// NLM defined at runtime.
	l2 = el + NLM;	l_2 = el + 2*NLM;
#endif
	it = 0;
	for (im=0;im<=MMAX;im++) {
		for (l=im*MRES;l<=LMAX;l++) {
			el[it] = l;	l2[it] = l*(l+1.0);	l_2[it] = 1.0/(l*(l+1.0));
			it++;
		}
	}
	l_2[0] = 0.0;	// undefined for l=0 => replace with 0.
}
#endif
