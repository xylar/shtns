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

// NLM : total number of Spherical Harmonic coefficients.
#define NLM ( (MMAX+1)*(LMAX+1) - MRES*(MMAX*(MMAX+1))/2 )
// LM(l,m) : index in the Spherical Harmonic coefficient array [ (l,m) space ]
#define LiM(l,im) ( (im*(2*LMAX+2 -MRES*(im+1)))/2 + l )
#define LM(l,m) ( (m*(2*LMAX+2 -(m+MRES)))/(2*MRES) + l )

#ifndef M_PI
# define M_PI 3.1415926535897932384626433832795
#endif

// useful values for some basic spherical harmonic representations
// Y00_1 = 1/Y00 = spherical harmonic representation of 1 (l=0,m=0)
#define Y00_1 sqrt(4.*M_PI)
// Y10_ct = spherical harmonic representation of cos(theta) (l=1,m=0)
#define Y10_ct sqrt(4.*M_PI/3.)

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

/////////////////////////////////////////////////////
//   Scalar Spherical Harmonics Transform
// input  : ShF = spatial/fourrier data : complex double array of size NLAT*(NPHI/2+1) or double array of size NLAT*(NPHI/2+1)*2
// output : Slm = spherical harmonics coefficients : complex double array of size NLM
void spat_to_SH(complex double *ShF, complex double *Slm)
{
	complex double fpm[NLAT];	// symmetric (even) and anti-symmetric (odd) parts, interleaved.
	complex double *Sl;		// virtual pointers for given im
	double *zl;
	long int i,im,m,l;

	fftw_execute_dft_r2c(fft,(double *) ShF, ShF);

	im = 0;
		m=im*MRES;
		for (i=0;i<NLAT/2;i++) {	// compute symmetric and antisymmetric parts. m=0 : everything is REAL
			(double) fpm[2*i] = (double) ShF[i] + (double) ShF[NLAT-(i+1)];
			(double) fpm[2*i+1] = (double) ShF[i] - (double) ShF[NLAT-(i+1)];
		}
		l=m;
		Sl = &Slm[LiM(0,im)];		// virtual pointer for l=0 and im
		zl = zlm[im];
		ShF += NLAT;
		while (l<LMAX) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
			Sl[l] = 0.0;	Sl[l+1] = 0.0;		// Slm[LiM(l,im)] = 0.0;	Slm[LiM(l+1,im)] = 0.0;
			for (i=0;i<NLAT;i+=2) {
				(double) Sl[l] += (double) fpm[i] * zl[i];		// Slm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
				(double) Sl[l+1] += (double) fpm[i+1] * zl[i+1];	// Slm[LiM(l+1,im)] += zlm[im][(l+1-m)*NLAT/2 + i] * fm[i];
			}
			l+=2;
			zl += NLAT;
		}
		if (l==LMAX) {
			Sl[l] = 0.0;	// Slm[LiM(l,im)] = 0.0;
			for (i=tm[im];i<NLAT/2;i++) {	// polar optimization
				(double) Sl[l] += zl[i] * (double) fpm[2*i];	// Slm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
			}
		}
	for (im=1;im<=MMAX;im++) {
		m=im*MRES;
		for (i=tm[im];i<NLAT/2;i++) {	// compute symmetric and antisymmetric parts.
			fpm[2*i] = ShF[i] + ShF[NLAT-(i+1)];
			fpm[2*i+1] = ShF[i] - ShF[NLAT-(i+1)];
		}
		l=m;
		Sl = &Slm[LiM(0,im)];		// virtual pointer for l=0 and im
		zl = zlm[im];
		ShF += NLAT;
		while (l<LMAX) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
			Sl[l] = 0.0;	Sl[l+1] = 0.0;		// Slm[LiM(l,im)] = 0.0;	Slm[LiM(l+1,im)] = 0.0;
			for (i=tm[im]*2;i<NLAT;i+=2) {	// tm[im] : polar optimization
				Sl[l] += fpm[i] * zl[i];		// Slm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
				Sl[l+1] += fpm[i+1] * zl[i+1];	// Slm[LiM(l+1,im)] += zlm[im][(l+1-m)*NLAT/2 + i] * fm[i];
			}
			l+=2;
			zl += NLAT;
		}
		if (l==LMAX) {
			Sl[l] = 0.0;	// Slm[LiM(l,im)] = 0.0;
			for (i=tm[im];i<NLAT/2;i++) {	// polar optimization
				Sl[l] += zl[i] * fpm[2*i];	// Slm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
			}
		}
	}
}

/////////////////////////////////////////////////////
//   Scalar inverse Spherical Harmonics Transform
// input  : Slm = spherical harmonics coefficients : complex double array of size NLM [unmodified]
// output : ShF = spatial/fourrier data : complex double array of size NLAT*(NPHI/2+1) or double array of size NLAT*(NPHI/2+1)*2
void SH_to_spat(complex double *Slm, complex double *ShF)
{
	complex double fe, fo;		// even and odd parts
	complex double *Sl;
	double *yl;
	long int i,im,m,l;

	im = 0;
		m = im*MRES;
		Sl = &Slm[LiM(0,im)];	// virtual pointer for l=0 and im
		i=0;
		yl = ylm[im] + i*(LMAX-m+1) -m;
		while (i<NLAT/2) {	// ops : NLAT/2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
			fe = 0.0;	fo = 0.0;
			while (l<LMAX) {	// compute even and odd parts
				(double) fe += yl[l] * (double) Sl[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Slm[LiM(l,im)];
				(double) fo += yl[l+1] * (double) Sl[l+1];	// fo += ylm[im][i*(LMAX-m+1) + (l+1-m)] * Slm[LiM(l+1,im)];
				l+=2;
			}
			if (l==LMAX) {
				(double) fe += yl[l] * (double) Sl[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Slm[LiM(l,im)];
			}
			ShF[i] = fe + fo;
			ShF[NLAT-(i+1)] = fe - fo;
			i++;
			yl += (LMAX-m+1);
		}
		if (i<(NLAT+1)/2) {	// NLAT impair, fe only
			l=m;
			fe = 0.0;
			while (l<LMAX) {	// compute even and odd parts
				(double) fe += yl[l] * (double) Sl[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Slm[LiM(l,im)];
				l+=2;
			}
			if (l==LMAX) {
				(double) fe += yl[l] * (double) Sl[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Slm[LiM(l,im)];
			}
			ShF[i] = fe;
			i++;
			yl += (LMAX-m+1);
		}
		ShF += NLAT;
	for (im=1; im<=MMAX; im++) {
		m = im*MRES;
		Sl = &Slm[LiM(0,im)];	// virtual pointer for l=0 and im
		i=0;
		while (i<tm[im]) {	// polar optimization
			ShF[i] = 0.0;
			ShF[NLAT-tm[im] + i] = 0.0;	// south pole zeroes <=> ShF[im*NLAT + NLAT-(i+1)] = 0.0;
			i++;
		}
		yl = ylm[im] + i*(LMAX-m+1) -m;
		while (i<NLAT/2) {	// ops : NLAT/2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
			fe = 0.0;	fo = 0.0;
			while (l<LMAX) {	// compute even and odd parts
				fe += yl[l] * Sl[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Slm[LiM(l,im)];
				fo += yl[l+1] * Sl[l+1];	// fo += ylm[im][i*(LMAX-m+1) + (l+1-m)] * Slm[LiM(l+1,im)];
				l+=2;
			}
			if (l==LMAX) {
				fe += yl[l] * Sl[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Slm[LiM(l,im)];
			}
			ShF[i] = fe + fo;
			ShF[NLAT-(i+1)] = fe - fo;
			i++;
			yl += (LMAX-m+1);
		}
		ShF += NLAT;
	}
	for(im=MMAX+1; im<=NPHI/2; im++) {	// padding for high m's
		for (i=0;i<NLAT;i++)
			ShF[i] = 0.0;
		ShF += NLAT;
	}

	ShF -= NLAT*(NPHI/2+1);		// restore original pointer
	fftw_execute_dft_c2r(ifft, ShF, (double *) ShF);
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
#include "SHT/ST_to_spat.c"
}

void SHsph_to_spat(complex double *Slm, complex double *BtF, complex double *BpF)
{
#include "SHT/S_to_spat.c"
}

void SHtor_to_spat(complex double *Tlm, complex double *BtF, complex double *BpF)
{
#include "SHT/T_to_spat.c"
}

void spat_to_SHsphtor(complex double *BtF, complex double *BpF, complex double *Slm, complex double *Tlm)
{
	complex double teo[NLAT], peo[NLAT];	// theta and phi even and odd parts
	complex double *Sl, *Tl;		// virtual pointers for given im
	struct DtDp *dzl;
	long int i,im,m,l;

	fftw_execute_dft_r2c(fft,(double *) BtF, BtF);
	fftw_execute_dft_r2c(fft,(double *) BpF, BpF);

	im = 0;		// dzl.p = 0.0 : and evrything is REAL
		m=im*MRES;
		for (i=0;i<NLAT/2;i++) {	// compute symmetric and antisymmetric parts.
			(double) teo[2*i] = (double) BtF[i] + (double) BtF[NLAT-(i+1)];
			(double) teo[2*i+1] = (double) BtF[i] - (double) BtF[NLAT-(i+1)];
			(double) peo[2*i] = (double) BpF[i] + (double) BpF[NLAT-(i+1)];
			(double) peo[2*i+1] = (double) BpF[i] - (double) BpF[NLAT-(i+1)];
		}
		l=m;
		Sl = &Slm[LiM(0,im)];		// virtual pointer for l=0 and im
		Tl = &Tlm[LiM(0,im)];
		dzl = dzlm[im];
		BtF += NLAT;	BpF += NLAT;
		while (l<LMAX) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
			Sl[l] = 0.0;	Sl[l+1] = 0.0;		// Slm[LiM(l,im)] = 0.0;	Slm[LiM(l+1,im)] = 0.0;
			Tl[l] = 0.0;	Tl[l+1] = 0.0;
			for (i=0;i<NLAT;i+=2) {
				(double) Sl[l] += dzl[i].t * (double) teo[i+1];
				(double) Tl[l] -= dzl[i].t * (double) peo[i+1];
				
				(double) Sl[l+1] += dzl[i+1].t * (double) teo[i];
				(double) Tl[l+1] -= dzl[i+1].t * (double) peo[i];
			}
			l+=2;
			dzl += NLAT;
		}
		if (l==LMAX) {
			Sl[l] = 0.0;	Tl[l] = 0.0;
			for (i=0;i<NLAT/2;i++) {
				(double) Sl[l] += dzl[i].t * (double) teo[2*i+1];
				(double) Tl[l] -= dzl[i].t * (double) peo[2*i+1];
			}
		}
	for (im=1;im<=MMAX;im++) {
		m=im*MRES;
		for (i=tm[im];i<NLAT/2;i++) {	// compute symmetric and antisymmetric parts.
			teo[2*i] = BtF[i] + BtF[NLAT-(i+1)];
			teo[2*i+1] = BtF[i] - BtF[NLAT-(i+1)];
			peo[2*i] = BpF[i] + BpF[NLAT-(i+1)];
			peo[2*i+1] = BpF[i] - BpF[NLAT-(i+1)];
		}
		l=m;
		Sl = &Slm[LiM(0,im)];		// virtual pointer for l=0 and im
		Tl = &Tlm[LiM(0,im)];
		dzl = dzlm[im];
		BtF += NLAT;	BpF += NLAT;
		while (l<LMAX) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
			Sl[l] = 0.0;	Sl[l+1] = 0.0;		// Slm[LiM(l,im)] = 0.0;	Slm[LiM(l+1,im)] = 0.0;
			Tl[l] = 0.0;	Tl[l+1] = 0.0;
			for (i=tm[im]*2;i<NLAT;i+=2) {	// tm[im] : polar optimization
				Sl[l] += dzl[i].t *teo[i+1] - dzl[i].p *peo[i]*I;		// ref: these E. Dormy p 72.
				Tl[l] -= dzl[i].t *peo[i+1] + dzl[i].p *teo[i]*I;
				
				Sl[l+1] += dzl[i+1].t *teo[i] - dzl[i+1].p *peo[i+1]*I;
				Tl[l+1] -= dzl[i+1].t *peo[i] + dzl[i+1].p *teo[i+1]*I;
			}
			l+=2;
			dzl += NLAT;
		}
		if (l==LMAX) {
			Sl[l] = 0.0;	Tl[l] = 0.0;
			for (i=tm[im];i<NLAT/2;i++) {	// polar optimization
				Sl[l] += dzl[i].t *teo[2*i+1] - dzl[i].p *peo[2*i]*I;
				Tl[l] -= dzl[i].t *peo[2*i+1] + dzl[i].p *teo[2*i]*I;
			}
		}
	}
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
	f = pi/(n-1.0);
	for (j=0; j<n; j++) {
		x[j] = cos(f*j);
	}

// weights
	for (j=0; j<n; j++) {
		w[j] = 0.0;	// unsupported yet...
	}
	printf("          warning! only synthesis (inverse transform) supported so far for equaly spaced grid)\n");
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

	m = (n+1)/2;
	for (i=1;i<=m;i++) {
		z = cos(pi*((double)i-0.25)/((double)n+0.5));
		z1 = z+1;
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
		x[i-1] = -z;		// Build up the abscissas.
		x[n-i] = z;
		w[i-1] = 2.0/((1-z*z)*(pp*pp));		// Build up the weights.
		w[n-i] = w[i-1];
	}

// as we started with initial guesses, we should check if the gauss points are actually unique.
	for (i=n; i>0; i--) {
		if (x[i] == x[i-1]) runerr("bad gauss points\n");
	}

#ifdef _SH_DEBUG_
// test integral to compute :
	z = 0;
	for (i=0;i<n;i++) {
		z += w[i]*x[i]*x[i];
	}
	printf("          quadrature for 3/2.x^2 = %g (should be 1.0) error = %g\n",z*1.5,z*1.5-1.0);
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
	printf("       done.\n");
}

/** initialize SH transform.
input : eps = polar optimization threshold : polar coefficients below that threshold are neglected (for high ms).
        eps is the value under wich the polar values of the Legendre Polynomials Plm are neglected, leading to increased performance (a few percent).
	0 = no polar optimization;  1.e-14 = VERY safe;  1.e-10 = safe;  1.e-6 = aggresive.
*/
void init_SH(double eps)
{
	double xg[NLAT], wg[NLAT];	// gauss points and weights.
	double dtylm[LMAX+1];		// temp storage for derivative : d(P_l^m(x))/dx
	double iylm_fft_norm = 2.0*pi/NPHI;	// normation FFT pour zlm
	double t,tmax;
	long int it,im,m,l;

	printf("[init_SH] Lmax=%d, Nlat=%d, Mres=%d, Mmax*Mres=%d, Nlm=%d\n",LMAX,NLAT,MRES,MMAX*MRES,NLM);
	if (MMAX*MRES > LMAX) runerr("[init_SH] MMAX*MRES should not exceed LMAX");
	if (NLAT <= LMAX) runerr("[init_SH] NLAT should be at least LMAX+1");
	
#ifdef SHT_EQUAL
	printf("          => using Equaly Spaced Nodes\n");
	EqualSpaceNodes(ct,wg,NLAT);		// equaly-spaced points and weights.
	for (it=0;it<NLAT; it++) {
		st[it] = sqrt(1.0 - ct[it]*ct[it]);
		st_1[it] = 1.0/sqrt(1.0 - ct[it]*ct[it]);
	}
#else
	printf("          => using Gauss Nodes\n");
	GaussNodes(xg,wg,NLAT);	// generate gauss nodes and weights : xg = ]-1,1[ = cos(theta)
	for (it=0; it<NLAT; it++) {
		ct[it] = xg[NLAT-1-it];			// on prend theta = ]0,pi/2[ => cos(theta) = ]1,0[
		st[it] = sqrt(1.0 - ct[it]*ct[it]);
		st_1[it] = 1.0/sqrt(1.0 - ct[it]*ct[it]);
	}
#endif

#ifdef _SH_DEBUG_
// TEST if gauss points are ok.
	tmax = 0.0;
	for (it = 0; it<NLAT/2; it++) {
		t = gsl_sf_legendre_sphPlm(NLAT, 0, ct[it]);
		if (t>tmax) tmax = t;
//		printf("i=%d, x=%12.12g, p=%12.12g\n",i,ct[i],t);
	}
	printf("          max zero at Gauss node for Plm[l=LMAX+1,m=0] : %g\n",tmax);
	printf("          Gauss nodes :");
	for (it=0;it<NLAT/2; it++)
		printf(" %g",ct[it]);
	printf("\n");
#endif

// Allocate legendre functions lookup tables.
	ylm[0] = (double *) fftw_malloc(sizeof(double)* NLM*NLAT/2);
	dylm[0] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* NLM*NLAT/2);
	zlm[0] = (double *) fftw_malloc(sizeof(double)* NLM*NLAT/2);
	dzlm[0] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* NLM*NLAT/2);
	for (im=0; im<MMAX; im++) {
		m = im*MRES;
		ylm[im+1] = ylm[im] + NLAT/2*(LMAX+1 -m);
		dylm[im+1] = dylm[im] + NLAT/2*(LMAX+1 -m);
		zlm[im+1] = zlm[im] + NLAT/2*(LMAX+1 -m);
		dzlm[im+1] = dzlm[im] + NLAT/2*(LMAX+1 -m);
	}

// Even/Odd symmetry : ylm is even or odd across equator, as l-m is even or odd => only NLAT/2 points required.
// for synthesis (inverse transform)
	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
//		ylm[im] = (double *) fftw_malloc(sizeof(double)* (LMAX+1-m)*NLAT/2);
//		dylm[im] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* (LMAX+1-m)*NLAT/2);
		for (it=0;it<NLAT/2;it++) {
#ifdef SHT_EQUAL
			if ((m==1)&&(st[it]==0.)) {
				printf("special case : ylm undef (set to zero) for m=1, at the pole.\n");
				for (l=m; l<=LMAX; l++) {
					ylm[im][it*(LMAX-m+1) + (l-m)] = 0.0;
					dylm[im][it*(LMAX-m+1) + (l-m)].t = 0.0;
					dylm[im][it*(LMAX-m+1) + (l-m)].p = 0.0;
				}
			} else {
#endif
			gsl_sf_legendre_sphPlm_deriv_array(LMAX, m, ct[it], ylm[im] + it*(LMAX-m+1), dtylm);	// fixed im legendre functions lookup table.
			for (l=m; l<=LMAX; l++) {
				dylm[im][it*(LMAX-m+1) + (l-m)].t = -st[it] * dtylm[l-m];	// d(Plm(cos(t)))/dt = -sin(t) d(Plm(x))/dx
				dylm[im][it*(LMAX-m+1) + (l-m)].p = ylm[im][it*(LMAX-m+1) + (l-m)] *m/st[it];	// 1/sint(t) dYlm/dphi
#ifdef SHT_EQUAL
				if (st[it]==0.) dylm[im][it*(LMAX-m+1) + (l-m)].p = 0.0;
			}
#endif
			}
		}
	}
	
// for analysis (decomposition, direct transform) : transpose and multiply by gauss weight and other normalizations.
// interleave l and l+1 : this stores data in the way it will be read.
	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
//		zlm[im] = (double *) fftw_malloc(sizeof(double)* (LMAX+1-m)*NLAT/2);
//		dzlm[im] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* (LMAX+1-m)*NLAT/2);
		for (it=0;it<NLAT/2;it++) {
			for (l=m;l<LMAX;l+=2) {
				zlm[im][(l-m)*NLAT/2 + it*2] = ylm[im][it*(LMAX-m+1) + (l-m)] * wg[it] *iylm_fft_norm;
				zlm[im][(l-m)*NLAT/2 + it*2 +1] = ylm[im][it*(LMAX-m+1) + (l+1-m)] * wg[it] *iylm_fft_norm;
				dzlm[im][(l-m)*NLAT/2 + it*2].t = dylm[im][it*(LMAX-m+1) + (l-m)].t * wg[it] *iylm_fft_norm /(l*(l+1));
				dzlm[im][(l-m)*NLAT/2 + it*2].p = dylm[im][it*(LMAX-m+1) + (l-m)].p * wg[it] *iylm_fft_norm /(l*(l+1));
				dzlm[im][(l-m)*NLAT/2 + it*2+1].t = dylm[im][it*(LMAX-m+1) + (l+1-m)].t * wg[it] *iylm_fft_norm /((l+1)*(l+2));
				dzlm[im][(l-m)*NLAT/2 + it*2+1].p = dylm[im][it*(LMAX-m+1) + (l+1-m)].p * wg[it] *iylm_fft_norm /((l+1)*(l+2));
				if (l == 0) {		// les derivees sont nulles pour l=0 (=> m=0)
					dzlm[im][(l-m)*NLAT/2 + it*2].t = 0.0;
					dzlm[im][(l-m)*NLAT/2 + it*2].p = 0.0;
				}
			}
			if (l==LMAX) {		// last l is stored right away, without interleaving.
				zlm[im][(l-m)*NLAT/2 + it] = ylm[im][it*(LMAX-m+1) + (l-m)] * wg[it] *iylm_fft_norm;
				dzlm[im][(l-m)*NLAT/2 + it].t = dylm[im][it*(LMAX-m+1) + (l-m)].t * wg[it] *iylm_fft_norm /(l*(l+1));
				dzlm[im][(l-m)*NLAT/2 + it].p = dylm[im][it*(LMAX-m+1) + (l-m)].p * wg[it] *iylm_fft_norm /(l*(l+1));
			}
		}
	}

// POLAR OPTIMIZATION : analyzing coefficients, some can be safely neglected.
	for (im=0;im<=MMAX;im++) {
		m = im*MRES;
		tm[im] = 1.0*NLAT;
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

	planFFT();		// initialize fftw

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
