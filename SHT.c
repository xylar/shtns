///////////////////////////////////////////////
// SHT : Spherical Harmonic Transform
//   requires SHT.h for size parameters.
//////////////////////////////////////////////

#include <complex.h>
#include <math.h>
// FFTW la derivee d/dx = ik	(pas de moins !)
#include <fftw3.h>
// GSL for Legendre functions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>


// SHT.h : parameter for SHT (sizes : LMAX, NLAT, MMAX, MRES, NPHI)
#include "SHT.h"
// LMMAX : total number of Spherical Harmonic coefficients.
//#define LMMAX ( ((MMAX+1)*(2*LMAX+2-MMAX))/2 )
#define LMMAX ( (MMAX+1)*(LMAX+1) - MRES* MMAX*(MMAX+1)/2 )
// LM(l,m) : index in the Spherical Harmonic coefficient array [ (l,m) space ]
//#define LM(l,m) ((m*(2*LMAX +1-m))/2 + l)
#define LM(l,m) ( (m/MRES)*(2*LMAX+3 -m)/2 + l-m )
#define LiM(l,im) ( im*(2*LMAX+3 -(im+2)*MRES)/2 + l )

double pi = atan(1)*4;

double* ylm[MMAX+1];	// matrix for direct transform
double* iylm[MMAX+1];	// matrix for inverse transform
int tm[MMAX+1];		// start theta value for SH (polar optimization : near the poles the legendre polynomials go to zero for high m's)

fftw_plan ifft, fft;	// plans for FFTW.
	
unsigned fftw_plan_mode = FFTW_PATIENT;		// defines the default FFTW planner mode.



// ShF : spatial/fourrier data : complex double array of size NLAT*(NPHI/2+1) or double array of size NLAT*(NPHI/2+1)*2
// Slm : spherical harmonics coefficients : complex double array of size LMMAX
void spat_to_SH(complex double *ShF, complex double *Slm)
{
	complex double fp[NLAT/2];	// symmetric (even) part
	complex double fm[NLAT/2];	// anti-symmetric (odd) part
	int i,im,m,l;

	fftw_execute(fft);
/*
// UNOPTIMIZED
	for (im=0;im<=MMAX;im++) {
		m=im*MRES;
		for (l=m;l<=LMAX;l++) {		// ops : 2*NLAT*(LMAX-m+1)
			Slm[LM(l,m)] = 0.0;
			for (i=tm[im];i<NLAT-tm[im];i++) {
				Slm[LM(l,m)] += iylm[im][(l-m)*NLAT + i] * ShF[im*NLAT + i];
			}
		}
	}
*/
// OPTIMIZED.
	for (im=0;im<=MMAX;im++) {
		m=im*MRES;
		for (i=tm[im];i<NLAT/2;i++) {	// compute symmetric and antisymmetric parts.
			fp[i] = ShF[im*NLAT + i] + ShF[im*NLAT + NLAT-(i+1)];
			fm[i] = ShF[im*NLAT + i] - ShF[im*NLAT + NLAT-(i+1)];
		}
		l=m;
		while (l<LMAX) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
			Slm[LiM(l,im)] = 0.0;	Slm[LiM(l+1,im)] = 0.0;
			for (i=tm[im];i<NLAT/2;i++) {	// tm[im] : polar optimization
				Slm[LiM(l,im)] += iylm[im][(l-m)*NLAT + i] * fp[i];
				Slm[LiM(l+1,im)] += iylm[im][(l+1-m)*NLAT + i] * fm[i];
			}
			l+=2;
		}
		if (l==LMAX) {
			Slm[LiM(l,im)] = 0.0;
			for (i=tm[im];i<NLAT/2;i++) {	// polar optimization
				Slm[LiM(l,im)] += iylm[im][(l-m)*NLAT + i] * fp[i];
			}
		}
	}
}

void SH_to_spat(complex double *Slm, complex double *ShF)
{
	complex double fp[NLAT/2];	// symmetric (even) part
	complex double fm[NLAT/2];	// anti-symmetric (odd) part
	int i,im,m,l;

/*
// UNOPTIMIZED
	for (im=0;im<=MMAX;im++) {
		m = im*MRES;
		for (i=0;i<NLAT;i++) {
			ShF[im*NLAT +i] = 0.0;
		}
		for (i=tm[im];i<NLAT-tm[im];i++) {	// ops : 2*(lmax-m+1)*NLAT
			for(l=m;l<=LMAX;l++) {
				ShF[im*NLAT +i] += ylm[im][i*(LMAX-m+1) + (l-m)] * Slm[LM(l,m)];
			}
		}
	}
*/
// OPTIMIZED
	for (im=0;im<=MMAX;im++) {
		m = im*MRES;
		i=0;
		while (i<tm[im]) {	// polar optimization
			ShF[im*NLAT + i] = 0.0;
//			ShF[im*NLAT + NLAT-(i+1)] = 0.0;	// south pole zeroes
			ShF[im*NLAT + NLAT-tm[im] + i] = 0.0;	// south pole zeroes
			i++;
		}
		while (i<NLAT/2) {	// ops : NLAT/2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			fp[i] = 0.0;	fm[i] = 0.0;
			l=m;
			while (l<LMAX) {	// compute even and odd parts
				fp[i] += ylm[im][i*(LMAX-m+1) + (l-m)] * Slm[LiM(l,im)];
				fm[i] += ylm[im][i*(LMAX-m+1) + (l+1-m)] * Slm[LiM(l+1,im)];
				l+=2;
			}
			if (l==LMAX) {
				fp[i] += ylm[im][i*(LMAX-m+1) + (l-m)] * Slm[LiM(l,im)];
			}
			ShF[im*NLAT +i] = fp[i] + fm[i];
			ShF[im*NLAT + NLAT-(i+1)] = fp[i] - fm[i];
			i++;
		}
	}

	for (im=MMAX+1;im<=NPHI/2;im++) {		// padding for high m's
		for (i=0;i<NLAT;i++) {
			ShF[im*NLAT +i] = 0.0;
		}
	}

	fftw_execute(ifft);
}




void runerr(const char * error_text)
{
	printf("*** Run-time error : %s\n",error_text);
	exit(1);
}

// Generates the abscissa and weights for a Gauss-Legendre quadrature.
// Newton method from initial Guess to find the zeros of the Legendre Polynome
// x = abscissa, w = weights, n points.
// Reference:  Numerical Recipes, Cornell press.
void Gauss(double *x, double *w, int n)
{
	double z, z1, p1, p2, p3, pp, eps;
	int i,j,m;

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
		if (x[i] == x[i-1]) runerr("bas gauss points\n");
	}
}

// initialize FFTs using FFTW. stride = NLAT, (contiguous l)
void planFFT()
{
	int nfft = NPHI;
	int ncplx = NPHI/2 +1;
	int nreal;
	
	nreal = 2*ncplx;
	
// Allocate Spatial Fields.
	ShF = (complex double *) fftw_malloc(ncplx * NLAT * sizeof(complex double));
	Sh = (double *) ShF;	// alias for inplace.

	printf("[FFTW] Mmax=%d, Nphi=%d\n",MMAX,NPHI);

	if (NPHI < 2*MMAX) runerr("[FFTW] the condition Nphi >= 2*Mmax is not met.");
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
	printf("       done.\n");
}

// initialize SH transform.
void init_SH()
{
	double xg[NLAT], wg[NLAT];	// gauss points and weights.
	double eps = POLAR_OPT_THRESHOLD;	// eps : polar coefficients below that threshold are neglected (for high ms)
	double iylm_fft_norm = 2.0*pi/NPHI;	// normation FFT pour iylm
	double t,tmax;
	int it,im,m,l;

	printf("[init_SH] Lmax=%d, Nlat=%d, Mres=%d, Mmax*Mres=%d, LMmax=%d\n",LMAX,NLAT,MRES,MMAX*MRES,LMMAX);
	if (MMAX*MRES > LMAX) runerr("[init_SH] MMAX*MRES should not exceed LMAX");
	if (NLAT <= LMAX) runerr("[init_SH] NLAT should be at least LMAX+1");
	
	Gauss(xg,wg,NLAT);	// generate gauss nodes and weights

#ifdef _SH_DEBUG_
// TEST if gauss points are ok.
	tmax = 0.0;
	for (it = 0; it<NLAT; it++) {
		t = gsl_sf_legendre_sphPlm(NLAT, 0, xg[it]);
		if (t>tmax) tmax = t;
//		printf("i=%d, x=%12.12g, p=%12.12g\n",i,xg[i],t);
	}
	printf("          max zero at Gauss node for Plm[l=LMAX+1,m=0] : %g\n",tmax);
#endif

// for synthesis (inverse transform)
	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
		ylm[im] = (double *) malloc(sizeof(double)* (LMAX-m+1)*NLAT);
		for (it=0;it<NLAT;it++) {
//			gsl_sf_legendre_sphPlm_deriv_array(LMAX, m, xg[it], ylm[im], &dylm[it*LMMAX + LM(im,im)]);	// fixed im legendre functions lookup table.
			gsl_sf_legendre_sphPlm_array(LMAX, m, xg[it], ylm[im] + it*(LMAX-m+1));	// fixed im legendre functions lookup table.
		}
	}
	
// for analysis (decomposition, direct transform) : transpose and multiply by gauss weight.
	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
		iylm[im] = (double *) malloc(sizeof(double)* (LMAX-m+1)*NLAT);
		for (it=0;it<NLAT;it++) {
			for (l=m;l<=LMAX;l++) {
				iylm[im][(l-m)*NLAT + it] = ylm[im][it*(LMAX-m+1) + (l-m)] * wg[it] *iylm_fft_norm;
			}
		}
	}

// POLAR OPTIMIZATION : analysing coefficients, some can be safely neglected.
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
}
