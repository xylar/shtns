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


// SHT.h : parameter for SHT (sizes : LMAX, NLAT, MMAX, MRES, NPHI)
#ifndef LMAX
# include "SHT.h"
#endif

#ifndef MRES
# define MRES 1
#endif

/** ACCESS TO SPHERICAL HARMONICS COMPONENTS **
  NLM : total number of (l,m) components.
  LM(l,m) : macro returning array index for given l and m.
  LiM(l,im) : macro returning array index for given l and im.
  LM_LOOP( action ) : macro that performs "action" for every (l,m), with l and lm set, but neither m nor im.
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
#if MRES == 0
	#error "MRES cannot be 0 !"
#endif



struct DtDp {		// theta and phi derivatives stored together.
	double t, p;
};

const double pi = M_PI;
double ct[NLAT], st[NLAT], st_1[NLAT];	// cos(theta), sin(theta), 1/sin(theta);
double el[NLM], l2[NLM], l_2[NLM];		// l, l(l+1) and 1/(l(l+1))

long int tm[MMAX+1];	// start theta value for SH (polar optimization : near the poles the legendre polynomials go to zero for high m's)

double* ylm[MMAX+1];		// matrix for inverse transform (synthesis)
struct DtDp* dylm[MMAX+1];	// theta and phi derivative of Ylm matrix
double* zlm[MMAX+1];		// matrix for direct transform (analysis)
struct DtDp* dzlm[MMAX+1];

fftw_plan ifft, fft;	// plans for FFTW.
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
		gsl_sf_legendre_sphPlm_array(LMAX, m, cost, &yl[m]);
		for (l=m; l<=LMAX; l++)
			vr += yl[l] * creal( Qlm[l] );
	for (im=1; im<=MMAX; im++) {
		m = im*MRES;
		gsl_sf_legendre_sphPlm_array(LMAX, m, cost, &yl[m]);
		eimp = cos(m*phi) + I*sin(m*phi);
		Ql = &Qlm[LiM(0,im)];	// virtual pointer for l=0 and im
		for (l=m; l<=LMAX; l++)
			vr += yl[l] * creal( Ql[l]*eimp );
	}
	return vr;
}

double SHqst_to_point(complex double *Qlm, complex double *Slm, complex double *Tlm, double cost, double phi,
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
		gsl_sf_legendre_sphPlm_deriv_array(LMAX, m, cost, &yl[m], &dtyl[m]);
		for (l=m; l<=LMAX; l++) {
			*vr += yl[l] * creal( Qlm[l] );
			vst += (-sint *dtyl[l]) * creal( Slm[l] );
			vtt += (-sint *dtyl[l]) * creal( Tlm[l] );
		}
	for (im=1; im<=MMAX; im++) {
		m = im*MRES;
		gsl_sf_legendre_sphPlm_deriv_array(LMAX, m, cost, &yl[m], &dtyl[m]);
		eimp = cos(m*phi) + I*sin(m*phi);
		Ql = &Qlm[LiM(0,im)];	Sl = &Slm[LiM(0,im)];	Tl = &Tlm[LiM(0,im)];
		for (l=m; l<=LMAX; l++) {
			*vr += yl[l] * creal( Ql[l]*eimp );
			vst += (-sint *dtyl[l]) * creal(Sl[l]*eimp);
			vtt += (-sint *dtyl[l]) * creal(Tl[l]*eimp);
			vsp += (yl[l] *m/sint) * creal(I*Sl[l]*eimp);
			vtp += (yl[l] *m/sint) * creal(I*Tl[l]*eimp);
		}
	}
	*vt = vtp + vst;	// Bt = I.m/sint *T  + dS/dt
	*vp = vsp - vtt;	// Bp = I.m/sint *S  - dT/dt
}


/*
	INITIALIZATION FUNCTIONS
*/

void runerr(const char * error_text)
{
	printf("*** Run-time error : %s\n",error_text);
	exit(1);
}

void EqualSpaceNodes(double *x, double *w, int n)
{
	double f;
	int j;

// cos theta of latidunal points (equaly spaced in theta)
#ifndef _SHT_EO_
	f = pi/(n-1.0);
#else
	f = pi/(2.*n-1.0);
#endif
	for (j=0; j<n; j++) {
		x[j] = cos(f*j);
	}

// weights
	for (j=0; j<n; j++) {
		w[j] = 0.0;	// unsupported yet...
	}
	printf("          ! warning : only synthesis (inverse transform) supported so far for equaly spaced grid !\n");
}


// Generates the abscissa and weights for a Gauss-Legendre quadrature.
// Newton method from initial Guess to find the zeros of the Legendre Polynome
// x = abscissa, w = weights, n points.
// Reference:  Numerical Recipes, Cornell press.
void GaussNodes(double *x, double *w, int n)
{
	double z, z1, p1, p2, p3, pp, eps;
	long int i,j,m;

	eps = 1.0e-15;	// desired precision, minimum = 2.2204e-16 (double)

#ifdef _SHT_EO_
	n *=2;
#endif
	m = (n+1)/2;
	for (i=1;i<=m;i++) {
		z = cos(pi*((double)i-0.25)/((double)n+0.5));
		z1 = z+1.;
		while ( fabs(z-z1) > eps )
		{
			p1 = 1.0;
			p2 = 0.0;
			for(j=1;j<=n;j++) {
				p3 = p2;
				p2 = p1;
				p1 = ((2*j-1)*z*p2-(j-1)*p3)/j;	// The Legendre polynomial...
			}
			pp = ((double)n)*(z*p1-p2)/(z*z-1.0);                       // ... and its derivative.
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
	for (i=m; i>0; i--) {
		if (x[i] == x[i-1]) runerr("bad gauss points\n");
	}

#ifdef _SH_DEBUG_
// test integral to compute :
	z = 0;
	for (i=0;i<m;i++) {
		z += w[i]*x[i]*x[i];
	}
	printf("          quadrature for 3/2.x^2 = %g (should be 1.0) error = %g\n",z*3.0,z*3.0-1.0);
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
	
	nreal = 2*ncplx;
	
// Allocate dummy Spatial Fields.
	ShF = (complex double *) fftw_malloc(ncplx * NLAT * sizeof(complex double));
	Sh = (double *) ShF;

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

//	fft_norm = 1.0/nfft;
	fftw_free(ShF);
//	printf("       done.\n");
}

/** initialize SH transform.
input : eps = polar optimization threshold : polar coefficients below that threshold are neglected (for high ms).
        eps is the value under wich the polar values of the Legendre Polynomials Plm are neglected, leading to increased performance (a few percent).
	0 = no polar optimization;  1.e-14 = VERY safe;  1.e-10 = safe;  1.e-6 = aggresive.
*/
void init_SH(double eps)
{
	double wg[NLAT];	// gauss points and weights.
	double dtylm[LMAX+1];		// temp storage for derivative : d(P_l^m(x))/dx
	double iylm_fft_norm = 2.0*pi/NPHI;	// normation FFT pour zlm
	double t,tmax;
	long int it,im,m,l;

#ifdef _SHT_EO_
	iylm_fft_norm *= 2.0;	// normation must be multiplied by 2.
#endif
	printf("[init_SH] Lmax=%d, Nlat=%d, Mres=%d, Mmax*Mres=%d, Nlm=%d\n",LMAX,NLAT,MRES,MMAX*MRES,NLM);
	if (MMAX*MRES > LMAX) runerr("[init_SH] MMAX*MRES should not exceed LMAX");
	if (NLAT <= LMAX) runerr("[init_SH] NLAT should be at least LMAX+1");

#ifdef SHT_EQUAL
	printf("          => using Equaly Spaced Nodes\n");
	EqualSpaceNodes(ct,wg,NLAT);		// equaly-spaced points and weights.
#else
	printf("          => using Gauss Nodes\n");
	GaussNodes(ct,wg,NLAT);	// generate gauss nodes and weights : ct = ]1,-1[ = cos(theta)
#endif
	for (it=0; it<NLAT; it++) {
		st[it] = sqrt(1.0 - ct[it]*ct[it]);
		st_1[it] = 1.0/sqrt(1.0 - ct[it]*ct[it]);
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
	printf("          Gauss nodes :");
	for (it=0;it<NLAT_2; it++)
		printf(" %g",ct[it]);
	printf("\n");
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
					dylm[im][it*(LMAX-m+1) + (l-m)].t = -st[it] * dtylm[l-m];	// d(Plm(cos(t)))/dt = -sin(t) d(Plm(x))/dx
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

// POLAR OPTIMIZATION : analyzing coefficients, some can be safely neglected.
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

 #if NPHI > 1
	planFFT();		// initialize fftw
 #else
	printf("          => no fft required for NPHI=1.\n");
 #endif

// Additional arrays :
	it = 0;
	for (im=0;im<=MMAX;im++) {
		for (l=im*MRES;l<=LMAX;l++) {
			el[it] = l;	l2[it] = l*(l+1.0);	l_2[it] = 1.0/(l*(l+1.0));
			it++;
		}
	}
	l_2[0] = 0.0;	// undefined for l=0 => replace by 0.
}
