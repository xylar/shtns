/********************************************************************
 * SHTns : Spherical Harmonic Transform for numerical simulations.  *
 *    written by Nathanael Schaeffer / LGIT,CNRS                    *
 ********************************************************************/

/// \file SHT.c main source file for SHTns.

/* FLAGS */
//#define SHT_AXISYM
//#define SHT_NO_DCT
#define SHT_NLAT_EVEN
//#define SHT_EO

/// 0:no output, 1:output info to stdout, 2:more output (debug info)
#define SHT_VERBOSE 2

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
// FFTW la derivee d/dx = ik	(pas de moins !)
#include <fftw3.h>
// GSL for Legendre functions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>

// cycle counter from FFTW
#include "cycle.h"

/// Minimum performance improve for DCT in sht_auto mode. If not atained, we switch back to gauss.
#define MIN_PERF_IMPROVE_DCT 0.95
/// Try to enforce at least this accuracy for DCT in sht_auto mode.
#define MIN_ACCURACY_DCT 1.e-8

// supported types of sht's
enum shtns_type {
	sht_gauss,	// use gaussian grid and quadrature. highest accuracy.
	sht_auto,	// use a regular grid if dct is faster with goog accuracy, otherwise defaults to gauss.
	sht_reg_fast,	// use fastest algorithm, on a regular grid, mixing dct and regular quadrature.
	sht_reg_dct,	// use pure dct algorithm, on a regular grid.
	sht_quick_init,	// gauss grid, with minimum init time (useful for pre/post-processing).
	sht_reg_poles	// use a synthesis only algo including poles, not suitable for computations.
};

#define SHT_NATIVE_LAYOUT 0
#define SHT_THETA_CONTIGUOUS 256
#define SHT_PHI_CONTIGUOUS 256*2
//#define SHT_CUSTOM_LAYOUT 256*3

/*
struct SHTdef {
	long int nlm;
	long int lmax,nlat,nlat_2;
	long int mmax,mres,nphi;
	long int mtr_dct;

	double *ct, *st, *st_1;		// cos(theta), sin(theta), 1/sin(theta);
	double *el, *l2, *l_2;		// l, l(l+1) and 1/(l(l+1))
	int *li;

	long int *lmidx;		// (virtual) index in SH array of given im.
	long int *tm;		// start theta value for SH (polar optimization : near the poles the legendre polynomials go to zero for high m's)

	double** ylm;		// matrix for inverse transform (synthesis)
	double** zlm;		// matrix for direct transform (analysis)
	double** ykm_dct;	// matrix for inverse transform (synthesis) using dct.
	double* zlm_dct0;	// matrix for direct transform (analysis), only m=0

	fftw_plan ifft, fft;	// plans for FFTW.
	fftw_plan idct, dctm0;
}

struct VSHTdef {
	long int nlm;
	long int lmax,nlat,nlat_2;
	long int mmax,mres,nphi;
	long int mtr_dct;

	double *ct, *st, *st_1;		// cos(theta), sin(theta), 1/sin(theta);
	double *el, *l2, *l_2;		// l, l(l+1) and 1/(l(l+1))
	int *li;

	long int *lmidx;		// (virtual) index in SH array of given im.
	long int *tm;		// start theta value for SH (polar optimization : near the poles the legendre polynomials go to zero for high m's)

	struct DtDp** dylm;	// theta and phi derivative of Ylm matrix
	struct DtDp** dzlm;
	struct DtDp** dykm_dct;	// theta and phi derivative of Ylm matrix

	fftw_plan ifft, fft;	// plans for FFTW.
	fftw_plan idct;

	spat_to_SH
	SH_to_spat
}
*/


long int SHT_FFT = 0;	///< How to perform fft : 0=no fft, 1=in-place, 2=out-of-place.
// parameter for SHT (sizes : LMAX, NLAT, MMAX, MRES, NPHI)
long int LMAX = -1;	// maximum degree (LMAX) of spherical harmonics.
long int NLAT, NLAT_2;	// number of spatial points in Theta direction (latitude) and half of it (using (NLAT+1)/2 allows odd NLAT.)
#ifndef SHT_AXISYM
  long int MMAX,MRES;	// maximum order (MMAX*MRES) of spherical harmonics. MRES is the periodicity along the phi axis.
  long int NPHI;	// number of spatial points in Phi direction (longitude)
#else
  #define MMAX 0
  #define NPHI 1
  #define MRES 1
#endif
long int NLM = 0;	// total number of (l,m) spherical harmonics components.
long int MTR_DCT = -1;	// m truncation for dct. -1 means no dct at all.

// number of double that have to be allocated for a spatial field (includes reserved space)
#define NSPAT_ALLOC (NLAT*(NPHI/2+1)*2)

#ifndef M_PI
# define M_PI 3.1415926535897932384626433832795
#endif


struct DtDp {		// theta and phi derivatives stored together.
	double t, p;
};

double *ct, *st, *st_1;		// cos(theta), sin(theta), 1/sin(theta);
double *el, *l2, *l_2;		// l, l(l+1) and 1/(l(l+1))
int *li;

long int *lmidx;		// (virtual) index in SH array of given im.
// LM(l,m) : index in the Spherical Harmonic coefficient array [ (l,m) space ]
#define LiM(l,im) ( lmidx[im] + l )

long int *tm;		// start theta value for SH (polar optimization : near the poles the legendre polynomials go to zero for high m's)
double** ylm;		// matrix for inverse transform (synthesis)
struct DtDp** dylm;	// theta and phi derivative of Ylm matrix
double** zlm;		// matrix for direct transform (analysis)
struct DtDp** dzlm;

double** ykm_dct;	// matrix for inverse transform (synthesis) using dct.
struct DtDp** dykm_dct;	// theta and phi derivative of Ylm matrix
double* zlm_dct0;	// matrix for direct transform (analysis), only m=0
double* dzlm_dct0;

fftw_plan ifft, fft;	// plans for FFTW.
fftw_plan idct, dct_m0;			// (I)DCT for NPHI>1
fftw_plan idct_r1, dct_r1;		// (I)DCT for axisymmetric case, NPHI=1
unsigned fftw_plan_mode = FFTW_EXHAUSTIVE;		// defines the default FFTW planner mode.

#define SSE __attribute__((aligned (16)))

/*
	SHT FUNCTIONS
*/

// truncation at LMAX and MMAX
#define LTR LMAX
#define MTR MMAX
#undef SHT_VAR_LTR

/////////////////////////////////////////////////////
///   Scalar Spherical Harmonics Transform : spatial field BrF is converted to SH representation Qlm.
/// \param[in] Vr = spatial data : double array of size NLAT*(NPHI/2+1)*2
/// \param[out] Qlm = spherical harmonics coefficients : complex double array of size NLM
void spat_to_SH(double *Vr, complex double *Qlm)
{
	#include "SHT/spat_to_SH.c"
}

/////////////////////////////////////////////////////
///   Backward Scalar Spherical Harmonics Transform : SH representation Qlm is converted to spatial field BrF.
/// \param[in] Qlm = spherical harmonics coefficients : complex double array of size NLM [unmodified]
/// \param[out] Vr = spatial data : double array of size NLAT*(NPHI/2+1)*2
void SH_to_spat(complex double *Qlm, double *Vr)
{
	#include "SHT/SH_to_spat.c"
}

//void SH_to_grad_spat(complex double *Slm, complex double *BtF, complex double *BpF)
// "grad_spat" and "sph_to_spat" are the same : so we alias them.
//#include "SHT/S_to_spat.c"
#define SH_to_grad_spat(S,Gt,Gp) SHsph_to_spat(S, Gt, Gp)

/////////////////////////////////////////////////////
/// Backward Vector Spherical Harmonics Transform : Spheroidal/Toroidal to (theta,phi) vector components.
/// \param[in] Slm/Tlm : SH array of Spheroidal/Toroidal scalar. complex size NLM
/// \param[out] Vt/Vp : theta/phi components of vector field. double array of size NLAT*(NPHI/2+1)*2
void SHsphtor_to_spat(complex double *Slm, complex double *Tlm, double *Vt, double *Vp)
{
	#include "SHT/SHst_to_spat.c"
}

void SHsph_to_spat(complex double *Slm, double *Vt, double *Vp)
{
	#include "SHT/SHs_to_spat.c"
}

void SHtor_to_spat(complex double *Tlm, double *Vt, double *Vp)
{
	#include "SHT/SHt_to_spat.c"
}

/////////////////////////////////////////////////////
/// Forward Vector Spherical Harmonics Transform : (theta,phi) vector field components to Spheroidal/Toroidal scalars.
/// \param[in] Vt/Vp : theta/phi components of vector field. double array of size NLAT*(NPHI/2+1)*2
/// \param[out] Slm/Tlm : SH array of Spheroidal/Toroidal scalar. complex size NLM
void spat_to_SHsphtor(double *Vt, double *Vp, complex double *Slm, complex double *Tlm)
{
	#include "SHT/spat_to_SHst.c"
}

/// Evaluate scalar SH representation \b Qlm at physical point defined by \b cost = cos(theta) and \b phi
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

/// Evaluate vector SH representation \b Qlm at physical point defined by \b cost = cos(theta) and \b phi
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


/*
	SHT FUNCTIONS with variable LTR truncation.
*/

// truncation at function parameter LTR
#undef LTR
#define SHT_VAR_LTR

/// spatial field Vr is converted to SH representation Qlm with maximum degree LTR
void spat_to_SH_l(double *Vr, complex double *Qlm, int LTR)
{
	#include "SHT/spat_to_SH.c"
}

/////////////////////////////////////////////////////
//   Scalar inverse Spherical Harmonics Transform
// input  : Qlm = spherical harmonics coefficients : complex double array of size NLM [unmodified]
// output : Vr = spatial data : double array of size NLAT*(NPHI/2+1)*2

void SH_to_spat_l(complex double *Qlm, double *Vr, int LTR)
{
	#include "SHT/SH_to_spat.c"
}

//void SH_to_grad_spat(complex double *Slm, complex double *BtF, complex double *BpF)
// "grad_spat" and "sph_to_spat" are the same : so we alias them.
//#include "SHT/S_to_spat.c"
#define SH_to_grad_spat_l(S,Gt,Gp,ltr) SHsph_to_spat(S, Gt, Gp, ltr)

/////////////////////////////////////////////////////
//   Spheroidal/Toroidal to (theta,phi) components inverse Spherical Harmonics Transform
// input  : Slm,Tlm = spherical harmonics coefficients of Spheroidal and Toroidal scalars : 
//          complex double array of size NLM [unmodified]
// output : Vt, Vp = theta, and phi vector components, spatial data : 
//          double array of size NLAT*(NPHI/2+1)*2
void SHsphtor_to_spat_l(complex double *Slm, complex double *Tlm, double *Vt, double *Vp, int LTR)
{
	#include "SHT/SHst_to_spat.c"
}

void SHsph_to_spat_l(complex double *Slm, double *Vt, double *Vp, int LTR)
{
#include "SHT/SHs_to_spat.c"
}

void SHtor_to_spat_l(complex double *Tlm, double *Vt, double *Vp, int LTR)
{
#include "SHT/SHt_to_spat.c"
}

void spat_to_SHsphtor_l(double *Vt, double *Vp, complex double *Slm, complex double *Tlm, int LTR)
{
	#include "SHT/spat_to_SHst.c"
}

#undef MTR
#undef SHT_VAR_LTR

/*
	INITIALIZATION FUNCTIONS
*/

void runerr(const char * error_text)
{
	printf("*** Run-time error : %s\n",error_text);
	exit(1);
}

void alloc_SHTarrays()
{
	long int im,m;

	ct = (double *) fftw_malloc(sizeof(double) * NLAT*3);
	st = ct + NLAT;		st_1 = ct + 2*NLAT;
	tm = (long int *) fftw_malloc(sizeof(long int) * (MMAX+1)*2);
	lmidx =  tm + (MMAX+1);
	ylm = (double **) fftw_malloc( sizeof(double *) * (MMAX+1)*3 );
	zlm = ylm + (MMAX+1);		ykm_dct = ylm + (MMAX+1)*2;
	dylm = (struct DtDp **) fftw_malloc( sizeof(struct DtDp *) * (MMAX+1)*3);
	dzlm = dylm + (MMAX+1);		dykm_dct = dylm + (MMAX+1)*2;

	el = (double *) fftw_malloc( 3*NLM*sizeof(double) + NLM*sizeof(int) );	// NLM defined at runtime.
	l2 = el + NLM;	l_2 = el + 2*NLM;
	li = (int *) (el + 3*NLM);

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
#if SHT_VERBOSE > 1
	printf("          Memory used for Ylm and Zlm matrices = %.3f Mb x2\n",3.0*sizeof(double)*NLM*NLAT_2/(1024.*1024.));
#endif
}

/// compute number of spherical harmonics modes (l,m) for given size parameters. Does not require a previous call to init_SH
/*! \code return (mmax+1)*(lmax+1) - mres*(mmax*(mmax+1))/2; \endcode */
long int nlm_calc(long int lmax, long int mmax, long int mres)
{
	if (mmax*mres > lmax) mmax = lmax/mres;
	return( (mmax+1)*(lmax+1) - mres*(mmax*(mmax+1))/2 );	// this is wrong if lmax < mmax*mres
}

/*
long int nlm_calc_eo(long int lmax, long int mmax, long int mres) {
	long int im,l,lm;
	for (im=0, lm=0; im<=mmax; im++) {
		if (im*mres <= lmax) lm += (lmax+2-im*mres)/2;
	}
	return lm;
}
*/

/// Generates the abscissa and weights for a Gauss-Legendre quadrature.
/// Newton method from initial Guess to find the zeros of the Legendre Polynome
/// \param x = abscissa, \param w = weights, \param n points.
/// \note Reference:  Numerical Recipes, Cornell press.
void GaussNodes(long double *x, long double *w, int n)
{
	long double z, z1, p1, p2, p3, pp, eps;
	long int i,j,m;
	long double pi = M_PI;

	eps = 1.1e-19;	// desired precision, minimum = 1.0842e-19 (long double)

	m = (n+1)/2;
	for (i=1;i<=m;i++) {
		z = cosl(pi*((long double)i-0.25)/((long double)n+0.5));
		do {
			p1 = 1.0;
			p2 = 0.0;
			for(j=1;j<=n;j++) {
				p3 = p2;
				p2 = p1;
				p1 = ((2*j-1)*z*p2-(j-1)*p3)/j;	// The Legendre polynomial...
			}
			pp = ((long double)n)*(z*p1-p2)/(z*z-1.0);                       // ... and its derivative.
			z1 = z;
			z = z-p1/pp;
		} while ( fabsl(z-z1) > eps );
		x[i-1] = z;		// Build up the abscissas.
		w[i-1] = 2.0/((1-z*z)*(pp*pp));		// Build up the weights.
		x[n-i] = -z;
		w[n-i] = w[i-1];
	}

// as we started with initial guesses, we should check if the gauss points are actually unique.
	for (i=m-1; i>0; i--) {
		if (x[i] == x[i-1]) runerr("[SHTns] bad gauss points");
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

void EqualPolarGrid()
{
	int j;
	long double f;
	long double pi = M_PI;

#if SHT_VERBOSE > 0
	printf("        => using Equaly Spaced Nodes including poles\n");
#endif
// cos theta of latidunal points (equaly spaced in theta)
	f = pi/(NLAT-1.0);
	for (j=0; j<NLAT; j++) {
		ct[j] = cosl(f*j);
		st[j] = sinl(f*j);
		st_1[j] = 1.0/sinl(f*j);
	}
#if SHT_VERBOSE > 0
	printf("     !! Warning : only synthesis (inverse transform) supported so far for this grid !\n");
#endif
}


/// initialize FFTs using FFTW. stride = NLAT, (contiguous l)
/// \param[in] theta_inc,phi_inc are the increments to go from one data value to the next in theta and phi direction respectively.
/// \param[in] phi_embed is the size of array in which the nphi elements are embedded (if phi_embed > (NPHI/2+1)*2, in-place fft may be used)
void planFFT(int theta_inc, int phi_inc, int phi_embed)
{
	complex double *ShF;
	double *Sh;
	int nfft, ncplx, nreal;

	if (NPHI <= 2*MMAX) runerr("the sampling condition Nphi > 2*Mmax is not met.");

  if (NPHI>1) {
	SHT_FFT = 1;		// yes, do some fft
	nfft = NPHI;
	ncplx = NPHI/2 +1;
	nreal = phi_embed;
	if ((theta_inc != 1)||(phi_inc != NLAT)||(nreal < 2*ncplx)) {
		SHT_FFT = 2;		// we need to do the fft out-of-place.
	}

#if SHT_VERBOSE > 0
	printf("        using FFTW : Mmax=%d, Nphi=%d  (data layout : phi_inc=%d, theta_inc=%d, phi_embed=%d)\n",MMAX,NPHI,phi_inc,theta_inc,phi_embed);
	if (NPHI < 3*MMAX) printf("     !! Warning : 2/3 rule for anti-aliasing not met !\n");
	if (SHT_FFT > 1) printf("        ** out-of-place fft **\n");
#endif

// Allocate dummy Spatial Fields.
	if (SHT_FFT > 1) {
		ShF = (complex double *) fftw_malloc(ncplx * NLAT * sizeof(complex double));
		Sh = (double *) fftw_malloc(ncplx * NLAT * sizeof(complex double));
	} else {
		ShF = (complex double *) fftw_malloc(ncplx * NLAT * sizeof(complex double));
		Sh = (double *) ShF;
	}

// IFFT : unnormalized.  FFT : must be normalized.
		ifft = fftw_plan_many_dft_c2r(1, &nfft, NLAT, ShF, &ncplx, NLAT, 1, Sh, &nreal, phi_inc, theta_inc, fftw_plan_mode);
		if (ifft == NULL) runerr("[FFTW] ifft planning failed !");
		fft = fftw_plan_many_dft_r2c(1, &nfft, NLAT, Sh, &nreal, phi_inc, theta_inc, ShF, &ncplx, NLAT, 1, fftw_plan_mode);
		if (fft == NULL) runerr("[FFTW] fft planning failed !");

#if SHT_VERBOSE > 1
	printf(" *** fft plan :\n");
	fftw_print_plan(fft);
	printf("\n *** ifft plan :\n");
	fftw_print_plan(ifft);
	printf("\n");
#endif

//	fft_norm = 1.0/nfft;
	if (SHT_FFT > 1) fftw_free(Sh);
	fftw_free(ShF);
  } else {
	if (theta_inc != 1) runerr("only contiguous spatial data is supported for NPHI=1");
#if SHT_VERBOSE > 0
	printf("        => no fft for NPHI=1.\n");
#endif
	SHT_FFT = 0;	// no fft.
  }
	dct_m0 = NULL;	idct = NULL;		// set dct plans to uninitialized.
	dct_r1 = NULL;	idct_r1 = NULL;
}

/// initialize DCTs using FFTW. Must be called if MTR_DCT is changed.
void planDCT()
{
	double *Sh;
	int ndct = NLAT;
	fftw_r2r_kind r2r_kind;
	fftw_iodim dims, hdims[2];
	double Sh0[NLAT] SSE;				// temp storage on the stack, aligned.
	
// real NPHI=1, allocate only once since it does not change.
	if ((dct_r1 == NULL)||(idct_r1 == NULL)) {
		Sh = (double *) fftw_malloc( NLAT * sizeof(double) );
		if (dct_r1 == NULL) {
			r2r_kind = FFTW_REDFT10;
			dct_r1 = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, 1, NLAT, Sh, &ndct, 1, NLAT, &r2r_kind, fftw_plan_mode);
		}
		if (idct_r1 == NULL) {
			r2r_kind = FFTW_REDFT01;
			idct_r1 = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, 1, NLAT, Sh, &ndct, 1, NLAT, &r2r_kind, fftw_plan_mode);
		}
		fftw_free(Sh);
		if ((dct_r1 == NULL)||(idct_r1 == NULL))
			runerr("[FFTW] (i)dct_r1 planning failed !");
#if SHT_VERBOSE > 1
			printf(" *** idct_r1 plan :\n");		fftw_print_plan(idct_r1);
			printf("\n *** dct_r1 plan :\n");	fftw_print_plan(dct_r1);	printf("\n");
#endif
	}

#ifndef SHT_AXISYM
	if (idct != NULL) fftw_destroy_plan(idct);
	// Allocate dummy Spatial Fields.
	Sh = (double *) fftw_malloc((NPHI/2 +1) * NLAT*2 * sizeof(double));

	dims.n = NLAT;	dims.is = 2;	dims.os = 2;		// real and imaginary part.
	hdims[0].n = MTR_DCT+1;	hdims[0].is = 2*NLAT; 	hdims[0].os = 2*NLAT;
	hdims[1].n = 2;			hdims[1].is = 1; 	hdims[1].os = 1;

	if (NPHI>1) {		// complex data for NPHI>1, recompute as it does depend on MTR_DCT
		r2r_kind = FFTW_REDFT01;
		idct = fftw_plan_guru_r2r(1, &dims, 2, hdims, Sh, Sh, &r2r_kind, fftw_plan_mode);
		if (idct == NULL)
			runerr("[FFTW] idct planning failed !");
#if SHT_VERBOSE > 1
			printf(" *** idct plan :\n");		fftw_print_plan(idct);	printf("\n");
#endif
		if (dct_m0 == NULL) {
			r2r_kind = FFTW_REDFT10;
//			dct_m0 = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, 2, 2*NLAT, Sh, &ndct, 2, 2*NLAT, &r2r_kind, fftw_plan_mode);
			dct_m0 = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, 2, 2*NLAT, Sh0, &ndct, 1, NLAT, &r2r_kind, fftw_plan_mode);	// out-of-place.
			if (dct_m0 == NULL)
				runerr("[FFTW] dct_m0 planning failed !");
#if SHT_VERBOSE > 1
				printf(" *** dct_m0 plan :\n");		fftw_print_plan(dct_m0);	printf("\n");
#endif
		}
	} else {	// NPHI==1
		if (dct_m0 == NULL) {
			r2r_kind = FFTW_REDFT10;
			dct_m0 = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, 1, NLAT, Sh0, &ndct, 1, NLAT, &r2r_kind, fftw_plan_mode);	// out-of-place.
			if (dct_m0 == NULL)
				runerr("[FFTW] dct_m0 planning failed !");
#if SHT_VERBOSE > 1
				printf(" *** dct_m0 plan :\n");		fftw_print_plan(dct_m0);	printf("\n");
#endif
		}
	}
	fftw_free(Sh);
#endif
}

/// SET MTR_DCT and updates fftw_plan for DCT's
void Set_MTR_DCT(int m)
{
	if ((zlm_dct0 == NULL)||(m == MTR_DCT)) return;
	if ( m < 0 ) {	// don't use dct
		MTR_DCT = -1;
	} else {
		if (m>MMAX) m=MMAX;
		MTR_DCT = m;
		planDCT();
	}
}

int Get_MTR_DCT() {
	return MTR_DCT;
}

/// TIMINGS
	double get_time(int m, int nloop, double *Sh, double *Th, complex double *Slm, complex double *Tlm)
	{
		int i;
		ticks tik0, tik1;

		Set_MTR_DCT(m);
			SH_to_spat(Tlm,Sh);			// caching...
			SHsphtor_to_spat(Slm,Tlm,Sh,Th);
		tik0 = getticks();
		for (i=0; i<nloop; i++) {
			SH_to_spat(Tlm,Sh);
			SHsphtor_to_spat(Slm,Tlm,Sh,Th);
		}
		tik1 = getticks();
	#if SHT_VERBOSE > 1
		printf("m=%d - ticks : %.3f\n", m*MRES, elapsed(tik1,tik0)/(nloop*NLM*NLAT));
	#endif
		return elapsed(tik1,tik0);
	}

double Find_Optimal_SHT()
{
	complex double *Slm, *Tlm;
	double *Sh, *Th;
	int m, i, minc, nloop;
	double t, tsum, tsum0;

	Sh = (double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Th = (double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Slm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Tlm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);

	t = 1.0 / (RAND_MAX/2);		// some random data to play with.
	for (i=0;i<NLM;i++) {
		Slm[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
		Tlm[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
	}

        minc = MMAX/20 + 1;             // don't test every single m.
        nloop = 10;                     // number of loops to get time.
        if (NLM*NLAT > 1024*1024)
                nloop = 1 + (10*1024*1024)/(NLM*NLAT);          // don't use too many loops for large transforms.
#if SHT_VERBOSE > 1
	printf("\nminc = %d, nloop = %d\n",minc,nloop);
#endif

	m = -1;
	tsum0 = get_time(m, nloop, Sh, Th, Slm, Tlm);
	tsum=tsum0;	i = m;
	t = get_time(m, nloop, Sh, Th, Slm, Tlm);
	if (t<tsum0) tsum0 = t;		// recheck m=-1
	for (m=0; m<=MMAX; m+=minc) {
		t = get_time(m, nloop, Sh, Th, Slm, Tlm);
		if ((m==-1)&&(t<tsum0)) tsum0 = t;
		if (t < tsum) {	tsum = t;	i = m; }
	}
	m=-1;	// recheck m=-1;
		t = get_time(m, nloop, Sh, Th, Slm, Tlm);
		if (t<tsum0) tsum0 = t;
		if (t < tsum) { tsum = t;	i = m; }
		t = get_time(m, nloop, Sh, Th, Slm, Tlm);	// twice
		if (t<tsum0) tsum0 = t;
		if (t < tsum) { tsum = t;	i = m; }

	free(Tlm);	free(Slm);	fftw_free(Th);		fftw_free(Sh);

	Set_MTR_DCT(i);
	return(tsum/tsum0);	// returns ratio of "optimal" time over "no_dct" time
}


/// Perform some optimization on the SHT matrices.
void OptimizeMatrices(double eps)
{
	double *yg;
	int im,m,l,it;

/// POLAR OPTIMIZATION : analyzing coefficients, some can be safely neglected.
	if (eps > 0.0) {
		for (im=0;im<=MMAX;im++) {
			m = im*MRES;
			tm[im] = NLAT_2;
			for (l=m;l<=LMAX;l++) {
				it=0;
				while( fabs(ylm[im][it*(LMAX-m+1) + (l-m)]) < eps ) { it++; }
				if (tm[im] > it) tm[im] = it;
			}
		}
#if SHT_VERBOSE > 0
		printf("        + polar optimization threshold = %.1e\n",eps);
#endif
#if SHT_VERBOSE > 1
		printf("          tm[im]=");
		for (im=0;im<=MMAX;im++)
			printf(" %d",tm[im]);
		printf("\n");
#endif
	} else {
		for (im=0;im<=MMAX;im++)	tm[im] = 0;
	}

/// Compression of dylm and dzlm for m=0, as .p is 0
	im=0;	m=0;
	yg = (double *) dylm[im];
	for (it=0; it<NLAT_2; it++) {
		for (l=m; l<=LMAX; l++)
			yg[it*(LMAX-m+1) + (l-m)] = dylm[im][it*(LMAX-m+1) + (l-m)].t;
	}
	yg = (double *) dzlm[im];
	if (yg != NULL) {	// for sht_reg_poles there is no dzlm defined.
		for (l=1; l<LMAX-1; l+=2) {		// l=0 is zero, so we start at l=1.
			for (it=0; it<NLAT_2; it++) {
				yg[(l-m-1)*NLAT_2 + it*2] = dzlm[im][(l-m+1)*NLAT_2 + it*2].t;		// l+1
				yg[(l-m-1)*NLAT_2 + it*2+1] = dzlm[im][(l-m-1)*NLAT_2 + it*2+1].t;	// l
			}
		}
		if (l==LMAX-1) {
			for (it=0; it<NLAT_2; it++) {
				yg[(l-m-1)*NLAT_2 + it*2] = dzlm[im][(l-m+1)*NLAT_2 + it].t;		// l+1
				yg[(l-m-1)*NLAT_2 + it*2+1] = dzlm[im][(l-m-1)*NLAT_2 + it*2+1].t;	// l
			}
		}
		if (l==LMAX) {		// last l is stored right away, without interleaving.
			for (it=0; it<NLAT_2; it++) {
				yg[(l-m-1)*NLAT_2 + it] = dzlm[im][(l-m-1)*NLAT_2 + it*2+1].t;		// l (odd)
			}
		}
	}
}



void init_SH_synth()
{
	double dtylm[LMAX+1];		// temp storage for derivative : d(P_l^m(x))/dx
	long int it,im,m,l;

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
}

void init_SH_gauss()
{
	double dtylm[LMAX+1];		// temp storage for derivative : d(P_l^m(x))/dx
	double iylm_fft_norm = 2.0*M_PI/NPHI;	// normation FFT pour zlm
	double t,tmax;
	long int it,im,m,l;
	long double xg[NLAT], wg[NLAT];	// gauss points and weights.

#if SHT_VERBOSE > 0
	printf("        => using Gauss Nodes\n");
#endif
	GaussNodes(xg,wg,NLAT);	// generate gauss nodes and weights : ct = ]1,-1[ = cos(theta)
	for (it=0; it<NLAT; it++) {
		ct[it] = xg[it];
		st[it] = sqrtl(1.0 - xg[it]*xg[it]);
		st_1[it] = 1.0/sqrtl(1.0 - xg[it]*xg[it]);
	}

#if SHT_VERBOSE > 1
	printf(" NLAT=%d, NLAT_2=%d\n",NLAT,NLAT_2);
// TEST if gauss points are ok.
	tmax = 0.0;
	for (it = 0; it<NLAT_2; it++) {
		t = gsl_sf_legendre_sphPlm(NLAT, 0, ct[it]);
		if (t>tmax) tmax = t;
//		printf("i=%d, x=%12.12g, p=%12.12g\n",it,ct[it],t);
	}
	printf("          max zero at Gauss node for Plm[l=NLAT,m=0] : %g\n",tmax);
	if (NLAT_2 < 100) {
		printf("          Gauss nodes :");
		for (it=0;it<NLAT_2; it++)
			printf(" %g",ct[it]);
		printf("\n");
	}
#endif

	init_SH_synth();
	
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
}

/// Clenshaw algorithm to evaluate a partial Chebyshev series dct[] up to degree n at any ct = cos(theta)
/// works for n>=2. see http://en.wikipedia.org/wiki/Clenshaw_algorithm
inline eval_dct(double *val, double *dct, long int n, double ct)
{
	double t1, t2, tmp;
	double ct2 = 2.*ct;
	long int i;

	t2 = dct[n];
	t1 = dct[n-1] + ct2*t2;
	for (i=n-2; i>0; i--) {
		tmp = t1;
		t1 = dct[i] + ct2*t1 - t2;
		t2 = tmp;
	}
	*val = dct[0] + ct*t1 - t2;
}

inline eval_dct_cplx(complex double *val, complex double *dct, long int n, double ct)
{
	complex double t1, t2, tmp;
	double ct2 = 2.*ct;
	long int i;

	t2 = dct[n];
	t1 = dct[n-1] + ct2*t2;
	for (i=n-2; i>0; i--) {
		tmp = t1;
		t1 = dct[i] + ct2*t1 - t2;
		t2 = tmp;
	}
	*val = dct[0] + ct*t1 - t2;
}

void init_SH_dct(int analysis)
{
	fftw_plan dct, idct;
	double *yk, *yk0, *dyk0, *yg;		// temp storage
	struct DtDp *dyg, *dyk;
	double dtylm[LMAX+1];		// temp storage for derivative : d(P_l^m(x))/dx
	double iylm_fft_norm = 2.0*M_PI/NPHI;	// FFT normation for zlm
	long int it,im,m,l;
	long int sk, dsk;
	double Z[2*NLAT_2], dZt[2*NLAT_2], dZp[2*NLAT_2];		// equally spaced theta points.
	double is1[NLAT];		// tabulate values for integrals.

#if SHT_VERBOSE > 0
	printf("        => using Equaly Spaced Nodes with DCT acceleration\n");
#endif
	if ((NLAT_2)*2 <= LMAX+1) runerr("[SHTns] NLAT_2*2 should be at least LMAX+2 (DCT)");
	if (NLAT & 1) runerr("[SHTns] NLAT must be even (DCT)");
	for (it=0; it<NLAT; it++) {	// Chebychev points : equaly spaced but skipping poles.
		long double th = M_PI*(it+0.5)/NLAT;
		ct[it] = cosl(th);
		st[it] = sinl(th);
		st_1[it] = 1.0/sinl(th);
	}
#if SHT_VERBOSE > 1
	{
	double tsum, t;
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
	}
#endif

#define KMAX (LMAX+1)

	for(im=0, sk=0, dsk=0; im<=MMAX; im++) {	// how much memory to allocate for ykm_dct ?
		m = im*MRES;
		for (it=0; it<= KMAX; it+=2) {
			l = (it < m) ? m : it-(m&1);
			sk += LMAX+1 - l;
		}
		for (it=0; it<= KMAX; it+=2) {
			l = (it-2 < m) ? m : it-2+(m&1);
			dsk += LMAX+1 - l;
		}
	}
	for (l=0, it=0; l<=LMAX; l+=2)	// how much memory for zlm_dct0 ?
		it += (2*NLAT_2 -l);
	for (l=1, im=0; l<=LMAX; l+=2)	// how much memory for dzlm_dct0 ?
		im += (2*NLAT_2 -l+1);

#if SHT_VERBOSE > 1
	printf("          Memory used for Ykm_dct matrices = %.3f Mb\n",sizeof(double)*(sk + 2.*dsk + it)/(1024.*1024.));
#endif
	ykm_dct[0] = (double *) fftw_malloc(sizeof(double)* sk);
	dykm_dct[0] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* dsk);
	zlm_dct0 = (double *) fftw_malloc( sizeof(double)* it );
	dzlm_dct0 = (double *) fftw_malloc( sizeof(double)* im );
	for (im=0; im<MMAX; im++) {
		m = im*MRES;
		for (it=0, sk=0; it<= KMAX; it+=2) {
			l = (it < m) ? m : it-(m&1);
			sk += LMAX+1 - l;
		}
		for (it=0, dsk=0; it<= KMAX; it+=2) {
			l = (it-2 < m) ? m : it-2+(m&1);
			dsk += LMAX+1 - l;
		}
		ykm_dct[im+1] = ykm_dct[im] + sk;
		dykm_dct[im+1] = dykm_dct[im] + dsk;
	}

	dct = fftw_plan_r2r_1d( 2*NLAT_2, Z, Z, FFTW_REDFT10, FFTW_ESTIMATE );	// quick and dirty dct.
	idct = fftw_plan_r2r_1d( 2*NLAT_2, Z, Z, FFTW_REDFT01, FFTW_ESTIMATE );	// quick and dirty idct.

#if SHT_VERBOSE > 1
	ticks tik0, tik1;
	tik0 = getticks();
#endif

// precomputation for scalar product of Chebychev polynomials.
	for(it=0; it<NLAT; it++)
		is1[it] = 1./(1. - 4.*it*it);

// Even/Odd symmetry : ylm is even or odd across equator, as l-m is even or odd => only NLAT_2 points required.
/// for synthesis (inverse transform)
	// temp memory for ykm_dct.
	yk = (double *) malloc( sizeof(double) * (KMAX+1)*(LMAX+1) );
	dyk = (struct DtDp *) malloc( sizeof(struct DtDp)* (KMAX+1)*(LMAX+1) );
	if (analysis) {
		yk0 = (double *) malloc( sizeof(double) * (LMAX/2+1)*(2*NLAT_2) * 2 );		// temp for zlm_dct0
		dyk0 = yk0 + (LMAX/2+1)*(2*NLAT_2);
	}

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
		for (it=0;it<=KMAX;it+=2) {
			for(l=m; l<=LMAX; l++) {
				yk[(it/2)*(LMAX+1-m) + (l-m)] = 0.0;
				dyk[(it/2)*(LMAX+1-m) + (l-m)].t = 0.0;
				dyk[(it/2)*(LMAX+1-m) + (l-m)].p = 0.0;
			}
		}
		for (l=m; l<=LMAX; l++) {
			if (m & 1) {	// m odd
				for (it=0; it<NLAT_2; it++) {
					Z[it] = ylm[im][it*(LMAX-m+1) + (l-m)] * st[it];	// P[l+1](x)	*st
//					Z[it] = ylm[im][it*(LMAX-m+1) + (l-m)] / st[it];	// P[l-1](x)	/st
					dZt[it] = dylm[im][it*(LMAX-m+1) + (l-m)].t * st[it];	// P[l](x)	*1
					dZp[it] = dylm[im][it*(LMAX-m+1) + (l-m)].p;		// P[l-1](x)	*1
				}
			} else {	// m even
				for (it=0; it<NLAT_2; it++) {
					Z[it] = ylm[im][it*(LMAX-m+1) + (l-m)];		// P[l](x)	*1
//					dZt[it] = dylm[im][it*(LMAX-m+1) + (l-m)].t;	// P[l-1](x)	/st
//					dZp[it] = ylm[im][it*(LMAX-m+1) + (l-m)] * m/(st[it]*st[it]);	// P[l-2](x)	/st
					dZt[it] = dylm[im][it*(LMAX-m+1) + (l-m)].t *st[it]*st[it];	// P[l+1](x)	*st
					dZp[it] = ylm[im][it*(LMAX-m+1) + (l-m)] * m;	// P[l](x)	*st
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
#if SHT_VERBOSE > 1
			if (LMAX <= 12) {
				printf("\nl=%d, m=%d ::\t", l,m);
				for(it=0;it<2*NLAT_2;it++) printf("%e ",Z[it]/(2*NLAT));
				printf("\n     dYt ::\t");
				for(it=0;it<2*NLAT_2;it++) printf("%e ",dZt[it]/(2*NLAT));
				printf("\n     dYp ::\t");
				for(it=0;it<2*NLAT_2;it++) printf("%e ",dZp[it]/(2*NLAT));
			}
#endif
			for (it=(l-m)&1; it<=l+1; it+=2) {
				yk[(it/2)*(LMAX+1-m) + (l-m)] = Z[it]/(2*NLAT);	// and transpose
				dyk[(it/2)*(LMAX+1-m) + (l-m)].p = dZp[it]/(2*NLAT);
			}
			for (it=(l+1-m)&1; it<=l+1; it+=2) {
				dyk[(it/2)*(LMAX+1-m) + (l-m)].t = dZt[it]/(2*NLAT);
			}
		}
		for (it=0; it<NLAT_2; it++) {		// do corrections for dylm.t non-DCT array.
			for (l=m; l<=LMAX; l++) {
				dylm[im][it*(LMAX-m+1) + (l-m)].t *= st[it];	// d(Plm(cos(t)))/dt = -sin(t) d(Plm(x))/dx
			}
		}
		
	/* compute analysis coefficients (fast way)
	 * Wklm = int(Tk*Ylm) = int(Tk.sum(i,a_ilm*Ti)) = sum(i, a_ilm* int(Tk*Ti)) = sum(i, a_ilm*Jik)
	 * with Jik = int(Tk*Ti) = 1/(1-(k-i)^2) + 1/(1-(k+i)^2)
	*/
		if (analysis) {
	#if SHT_VERBOSE > 0
		printf("computing weights m=%d\r",m);	fflush(stdout);
	#endif
		for (l=m; l<=LMAX; l++) {
			long int k0,k1, k,i,d;
			double Jik, yy, dyy;

			k0 = (l-m)&1;	k1 = 1-k0;
			for(k=0; k<NLAT; k++) {	Z[k] = 0.0;		dZt[k] = 0.0;	dZp[k] = 0.0; }
			for (i=k0; i<=l+1; i+=2) {		// i+k even
				yy = yk[(i/2)*(LMAX+1-m) + (l-m)] * iylm_fft_norm/NLAT_2;
				dyy = dyk[(i/2)*(LMAX+1-m) + (l-m)].p * iylm_fft_norm/(NLAT_2 *l*(l+1));
				if (i==0) {	yy*=0.5;	dyy*=0.5; }
				for (k=k0; k<NLAT; k+=2) {
					d = (k<i) ? i-k : k-i;		// d=|i-k|
					Jik = is1[(i+k)/2] + is1[d/2];
					Z[k] += yy * Jik;
					if (m&1) dZp[k] += dyy * Jik;
				}
			}
			if (l != 0) {
				for (i=k1; i<=l+1; i+=2) {		// i+k even
					yy = dyk[(i/2)*(LMAX+1-m) + (l-m)].t * iylm_fft_norm/(NLAT_2 *l*(l+1));
					if (i==0) yy*=0.5;
					for (k=k1; k<NLAT; k+=2) {
						d = (k<i) ? i-k : k-i;		// d=|i-k|
						Jik = is1[(i+k)/2] + is1[d/2];
						dZt[k] += yy * Jik;
					}
				}
			}
#if SHT_VERBOSE > 1
		if (LMAX <= 12) {
			printf("\nl=%d, m=%d ::\t",l,m);
			for (k=0; k<(2*NLAT_2); k++) printf("%f ",Z[k]);
			printf("\n       dZt ::\t");
			for (k=0; k<(2*NLAT_2); k++) printf("%f ",dZt[k]);
			if (m&1) {
				printf("\n       dZp ::\t");
				for (k=0; k<(2*NLAT_2); k++) printf("%f ",dZp[k]);
			}
		}
#endif

			if (m == 0) {		// we store zlm in dct space for m=0
				if (k0==0) 	{
					yk0[((l-m)>>1)*(2*NLAT_2)] = Z[0]*0.5;         // store zlm_dct (k=0)
					for (k=1; k<(2*NLAT_2); k++) yk0[((l-m)>>1)*(2*NLAT_2) +k] = 0.0;		// zero out.
					k0=2;
				}
				for (k=k0; k<(2*NLAT_2); k+=2)
					yk0[((l-m)>>1)*(2*NLAT_2) +k] = Z[k];             // store zlm_dct
					
				if (l>0) {
					if (k1==0) 	{
						dyk0[((l-1-m)>>1)*(2*NLAT_2)] = dZt[0]*0.5;         // store dzlm_dct (k=0)
						for (k=1; k<(2*NLAT_2); k++) dyk0[((l-1-m)>>1)*(2*NLAT_2) +k] = 0.0;		// zero out.
						k1=2;
					}
					for (k=k1; k<(2*NLAT_2); k+=2)
						dyk0[((l-1-m)>>1)*(2*NLAT_2) +k] = dZt[k];             // store dzlm_dct
				}
			}

			fftw_execute_r2r(idct, Z, Z);	fftw_execute_r2r(idct, dZt, dZt);
			if (m == 0) {
				for (it=0; it<NLAT; it++) { dZp[it] = 0.0; 	dZt[it] *= st_1[it]; }
			} else if (m & 1) {	// m odd
				fftw_execute_r2r(idct, dZp, dZp);
				for (it=0; it<NLAT; it++) {	Z[it] *= st_1[it]; }
			} else {	// m even
				for (it=0; it<NLAT; it++) { dZp[it] = Z[it]*m/(l*(l+1)*st[it]); 	dZt[it] *= st_1[it]; }
			}
			
			sk = (l-m)&1;
			if ((sk == 0)&&(l == LMAX)) {
				for (it=0; it<NLAT_2; it++) {
					zlm[im][(l-m)*NLAT_2 + it] =  Z[it];
					dzlm[im][(l-m)*NLAT_2 + it].p = dZp[it];
					dzlm[im][(l-m)*NLAT_2 + it].t = dZt[it];
				}
			} else {
				for (it=0; it<NLAT_2; it++) {
					zlm[im][(l-m-sk)*NLAT_2 + it*2 +sk] = Z[it];
					dzlm[im][(l-m-sk)*NLAT_2 + it*2 +sk].p = dZp[it];
					dzlm[im][(l-m-sk)*NLAT_2 + it*2 +sk].t = dZt[it];
				}
			}
		}
		}

		// Compact the coefficients for improved cache efficiency.
		yg = ykm_dct[im];
		dyg = dykm_dct[im];
		for (it=0; it<= KMAX; it+=2) {
			l = (it < m) ? m : it-(m&1);
			while (l<=LMAX) {
				yg[0] = yk[(it/2)*(LMAX+1-m) + (l-m)];
				l++;	yg++;
			}
		}
		for (it=0; it<= KMAX; it+=2) {
			l = (it-2 < m) ? m : it-2+(m&1);
			while (l<=LMAX) {
				dyg[0].t = dyk[(it/2)*(LMAX+1-m) + (l-m)].t;
				dyg[0].p = dyk[(it/2)*(LMAX+1-m) + (l-m)].p;
				l++;	dyg++;
			}
		}
		if (im == 0) {		// compact m=0 dylm because .p = 0 :
			dyg = dykm_dct[im];
			yg = (double *) dykm_dct[im];
			for (it=0; it<= KMAX; it+=2) {
				for (l=it-2; l<=LMAX; l++) {
					if (l>=0) {
						yg[0] = dyg[0].t;
						yg++;	dyg++;
					}
				}
			}
		}
	}
	
	// compact yk to zlm_dct0
	if (analysis) {
		yg = zlm_dct0;
		for (l=0; l<=LMAX; l+=2) {
			for (it=l; it<2*NLAT_2; it++) {	// for m=0, zl coeff with i<l are zeros.
				*yg = yk0[it];
				yg++;
			}
			yk0 += 2*NLAT_2;
		}
		yg = dzlm_dct0;
		for (l=1; l<=LMAX; l+=2) {
			for (it=l-1; it<2*NLAT_2; it++) {	// for m=0, dzl coeff with i<l-1 are zeros.
				*yg = dyk0[it];
				yg++;
			}
			dyk0 += 2*NLAT_2;
		}
		free(yk0 - (2*NLAT_2)*(LMAX/2+1));
	}

#if SHT_VERBOSE > 1
	tik1 = getticks();
	printf("\n    ticks : %.3f\n", elapsed(tik1,tik0)/(NLM*NLAT*(MMAX+1)));
#endif
	free(dyk);	free(yk);
	fftw_destroy_plan(idct);	fftw_destroy_plan(dct);
}


/// return the max error for a back-and-forth SHT transform.
/// \param[in] disp : 1 = print more informations about accuracy, 0 = silent.
double SHT_error()
{
	complex double *Tlm0, *Slm0, *Tlm, *Slm;
	double *Sh, *Th;
	double t, tmax, n2,  err;
	long int i, jj, nlm_cplx;
	
	srand( time(NULL) );	// init random numbers.
	
	Tlm0 = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Slm0 = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Slm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Tlm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Sh = (double *) fftw_malloc( NSPAT_ALLOC * sizeof(double) );
	Th = (double *) fftw_malloc( NSPAT_ALLOC * sizeof(double) );

// m = nphi/2 is also real if nphi is even.
	nlm_cplx = ( MMAX*2 == NPHI ) ? LiM(MRES*MMAX,MMAX) : NLM;
	t = 1.0 / (RAND_MAX/2);
	for (i=0; i<NLM; i++) {
		if ((i<=LMAX)||(i>=nlm_cplx)) {		// m=0 or m*2=nphi : real random data
			Slm0[i] = t*((double) (rand() - RAND_MAX/2));
			Tlm0[i] = t*((double) (rand() - RAND_MAX/2));
		} else {							// m>0 : complex random data
			Slm0[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
			Tlm0[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
		}
	}

	SH_to_spat(Slm0,Sh);		// scalar SHT
	spat_to_SH(Sh, Slm);
	for (i=0, tmax=0., n2=0.; i<NLM; i++) {		// compute error
		t = cabs(Slm[i] - Slm0[i]);
		n2 += t*t;
		if (t>tmax) { tmax = t; jj = i; }
	}
	err = tmax;
#if SHT_VERBOSE > 1
	printf("        scalar SH - poloidal   rms error = %.3g  max error = %.3g for l=%d,lm=%d\n",sqrt(n2/NLM),tmax,li[jj],jj);
#endif

	Slm0[0] = 0.0; 	Tlm0[0] = 0.0;		// l=0, m=0 n'a pas de signification sph/tor
	SHsphtor_to_spat(Slm0, Tlm0, Sh, Th);		// vector SHT
	spat_to_SHsphtor(Sh, Th, Slm, Tlm);
	for (i=0, tmax=0., n2=0.; i<NLM; i++) {		// compute error
		t = cabs(Slm[i] - Slm0[i]);
		n2 += t*t;
		if (t>tmax) { tmax = t; jj = i; }
	}
	if (tmax > err) err = tmax;
#if SHT_VERBOSE > 1
	printf("        vector SH - spheroidal rms error = %.3g  max error = %.3g for l=%d,lm=%d\n",sqrt(n2/NLM),tmax,li[jj],jj);
#endif
	for (i=0, tmax=0., n2=0.; i<NLM; i++) {		// compute error
		t = cabs(Tlm[i] - Tlm0[i]);
		n2 += t*t;
		if (t>tmax) { tmax = t; jj = i; }
	}
	if (tmax > err) err = tmax;
#if SHT_VERBOSE > 1
	printf("                  - toroidal   rms error = %.3g  max error = %.3g for l=%d,lm=%d\n",sqrt(n2/NLM),tmax,li[jj],jj);
#endif
	return(err);		// return max error.
}

#ifndef _HGID_
  #define _HGID_ "unknown"
#endif


/*! Initialization of Spherical Harmonic transforms (backward and forward, vector and scalar, ...) of given size.
 * <b>This function must be called after shtns_geometry and before any SH transform.</b> and sets all global variables.
 * \param lmax : maximum SH degree that we want to describe.
 * \param mmax : number of azimutal wave numbers.
 * \param mres : \c 2.pi/mres is the azimutal periodicity. \c mmax*mres is the maximum SH order.
 * \param nlat,nphi : respectively the number of latitudinal and longitudinal grid points.
 * \param flags allows to choose the type of transform (see \ref shtns_type) and the spatial data layout (see \ref spat)
 * \param eps : polar optimization threshold : polar values of Legendre Polynomials below that threshold are neglected (for high m), leading to increased performance (a few percents)
 *  0 = no polar optimization;  1.e-14 = VERY safe;  1.e-10 = safe;  1.e-6 = aggresive, but still good accuracy.
*/
void shtns_init(enum shtns_type flags, double eps, int lmax, int mmax, int mres, int nlat, int nphi)
{
	double t;
	int im,m,l,lm;
	int theta_inc, phi_inc, phi_embed;

	if (lmax < 1) runerr("[SHTns] lmax must be larger than 1");

	switch (flags & 0xFFFF00) {
		case SHT_NATIVE_LAYOUT : 	theta_inc=1;  phi_inc=nlat;  phi_embed=2*(nphi/2+1);  break;
		case SHT_THETA_CONTIGUOUS :	theta_inc=1;  phi_inc=nlat;  phi_embed=nphi;  break;
		default :
		case SHT_PHI_CONTIGUOUS :	phi_inc=1;  theta_inc=nphi;  phi_embed=nphi;  break;
	}
	flags = flags & 255;	// clear higher bits.

	// copy to global variables.
#ifdef SHT_AXISYM
	if ((mmax != MMAX)||(nphi != NPHI)) runerr("[SHTns] axisymmetric version : only Mmax=0 and Nphi=1 allowed");
#else
	MMAX = mmax;	MRES = mres;	NPHI = nphi;
#endif
	LMAX = lmax;	NLAT_2 = (nlat+1)/2;
	NLAT = nlat;
	NLM = nlm_calc(LMAX,MMAX,MRES);
  #if SHT_VERBOSE > 0
	printf("[SHTns] build " __DATE__ ", " __TIME__ ", id: " _HGID_ "\n");
	printf("        Lmax=%d, Nlat=%d, Mres=%d, Mmax*Mres=%d, Nlm=%d\n",LMAX,NLAT,MRES,MMAX*MRES,NLM);
  #endif
	if (MMAX*MRES > LMAX) runerr("[SHTns] MMAX*MRES should not exceed LMAX");
	if ((NLAT_2)*2 <= LMAX) runerr("[SHTns] NLAT_2*2 should be at least LMAX+1");
	if (MRES <= 0) runerr("[SHTns] MRES must be > 0");
#ifdef SHT_NLAT_EVEN
	if ((NLAT & 1)&&(flags != sht_reg_poles)) runerr("[SHTns] NLAT must be even.");
#endif

	alloc_SHTarrays();	// allocate dynamic arrays
	planFFT(theta_inc, phi_inc, phi_embed);		// initialize fftw
	zlm_dct0 = NULL;	// used as a flag.
  #if SHT_VERBOSE > 0
	if (2*NLAT <= 3*LMAX) printf("     !! Warning : anti-aliasing condition in theta direction not met.\n");
  #endif
// Additional arrays init :
	for (im=0, lm=0; im<=MMAX; im++) {
		m = im*MRES;
		lmidx[im] = lm -m;		// virtual pointer for l=0
		for (l=im*MRES;l<=LMAX;l++) {
			el[lm] = l;	l2[lm] = l*(l+1.0);	l_2[lm] = 1.0/(l*(l+1.0));
			li[lm] = l;
			lm++;
		}
	}
	l_2[0] = 0.0;	// undefined for l=0 => replace with 0.

#ifndef SHT_NO_DCT
	if (flags == sht_reg_dct) {	// pure dct.
		init_SH_dct(1);
		OptimizeMatrices(0.0);
		Set_MTR_DCT(MMAX);
	}
	if ((flags == sht_auto)||(flags == sht_reg_fast))
	{
		init_SH_dct(1);
		OptimizeMatrices(eps);
  #if SHT_VERBOSE > 0
		printf("finding optimal MTR_DCT ...\r");	fflush(stdout);
  #endif
		t = Find_Optimal_SHT();
  #if SHT_VERBOSE > 0
		printf("        + optimal MTR_DCT = %d  (%.1f%% performance gain)\n", MTR_DCT*MRES, 100.*(1/t-1));
  #endif
		if (t > MIN_PERF_IMPROVE_DCT) {
			Set_MTR_DCT(-1);		// turn off DCT.
		} else {
			t = SHT_error();
			if (t > MIN_ACCURACY_DCT) {
  #if SHT_VERBOSE > 0
				printf("     !! Not enough accuracy (%.3g) => turning off DCT.\n",t);
  #endif
				Set_MTR_DCT(-1);		// turn off DCT.
			}
		}
  #if SHT_VERBOSE < 2
		if (MTR_DCT == -1) {			// free memory used by DCT and disables DCT.
			fftw_free(zlm_dct0);	fftw_free(dykm_dct[0]);	fftw_free(ykm_dct[0]);		// free now useless arrays.
			zlm_dct0 = NULL;			// this completely disables DCT.
			if (idct != NULL) fftw_destroy_plan(idct);	// free unused dct plans
			if (dctm0 != NULL) fftw_destroy_plan(dctm0);
			if (flags == sht_auto) {
				flags = sht_gauss;		// switch to gauss grid, even better accuracy.
	#if SHT_VERBOSE > 0
				printf("        => switching back to Gauss Grid for higher accuracy.\n");
	#endif
			}
		}
  #endif
	}
	if (flags == sht_gauss)
#else
	if ((flags == sht_gauss)||(flags == sht_auto)||(flags == sht_reg_fast)||(flags == sht_reg_dct))
#endif
	{
		MTR_DCT = -1;		// we do not use DCT !!!
		init_SH_gauss();
		OptimizeMatrices(eps);
	}
	if (flags == sht_reg_poles)
	{
		fftw_plan_mode = FFTW_ESTIMATE;		// quick fftw init
		MTR_DCT = -1;		// we do not use DCT !!!
		fftw_free(dzlm[0]);	fftw_free(zlm[0]);	// no inverse transform.
		dzlm[0] = NULL;		zlm[0] = NULL;		// mark as unused.
		EqualPolarGrid();
		init_SH_synth();
		OptimizeMatrices(eps);
	}
	if (flags == sht_quick_init)
	{
		fftw_plan_mode = FFTW_ESTIMATE;		// quick fftw init
		MTR_DCT = -1;		// we do not use DCT !!!
		init_SH_gauss();
		OptimizeMatrices(eps);
	}

	if ((flags != sht_reg_poles)&&(flags != sht_quick_init)) {
		t = SHT_error();		// compute SHT accuracy.
  #if SHT_VERBOSE > 0
		printf("        + SHT accuracy = %.3g\n",t);
  #endif
  #if SHT_VERBOSE < 2
		if (t > 1.e-3) runerr("[SHTns] bad SHT accuracy");		// stop if something went wrong (but not in debug mode)
  #endif
	}
}


/*
	SHT FUNCTIONS for m=0 only (axisymmetric)
*/

#define SHT_AXISYM
// truncation at LMAX and MMAX
#define LTR LMAX
#define MTR MMAX
#undef SHT_VAR_LTR

void spat_to_SH_m0(double *Vr, complex double *Qlm)
{
	#include "SHT/spat_to_SH.c"
}

void SH_to_spat_m0(complex double *Qlm, double *Vr)
{
	#include "SHT/SH_to_spat.c"
}

void SHsphtor_to_spat_m0(complex double *Slm, complex double *Tlm, double *Vt, double *Vp)
{
	#include "SHT/SHst_to_spat.c"
}

void SHsph_to_spat_m0(complex double *Slm, double *Vt)
{
	#include "SHT/SHs_to_spat.c"
}

void SHtor_to_spat_m0(complex double *Tlm, double *Vp)
{
	#include "SHT/SHt_to_spat.c"
}

void spat_to_SHsphtor_m0(double *Vt, double *Vp, complex double *Slm, complex double *Tlm)
{
	#include "SHT/spat_to_SHst.c"
}

#undef SHT_AXISYM

/*
	SHT FUNCTIONS with assumed equatorial symmetry
*/

#define SHT_NO_DCT
#define LTR LMAX
#define MTR MMAX
#undef SHT_VAR_LTR

void SHeo_to_spat(complex double *Qlm, double *Vr, int parity)
{
	#include "SHT/SHe_to_spat.c"
}

void spat_to_SHeo(double *Vr, complex double *Qlm, int parity)
{
	#include "SHT/spat_to_SHe.c"
}

void SHeo_sphtor_to_spat(complex double *Slm, complex double *Tlm, double *Vt, double *Vp, int parity)
{
	#include "SHT/SHest_to_spat.c"
}

void spat_to_SHeo_sphtor(double *Vt, double *Vp, complex double *Slm, complex double *Tlm, int parity)
{
	#include "SHT/spat_to_SHest.c"
}
