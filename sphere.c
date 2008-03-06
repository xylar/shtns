// base flow. Cartesien spectral.

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
// FFTW la derivee d/dx = ik	(pas de moins !)
#include <fftw3.h>
// for fitting...
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>



#define LMAX 79
#define MMAX 0
#define LMMAX ( ((MMAX+1)*(2*LMAX+2-MMAX))/2 )
#define NLAT (LMAX+1)
#define NPHI (MMAX*2)

#define LM(l,m) ((m*(2*LMAX +1-m))/2 + l)

double pi = atan(1)*4;
double *ylm, *dylm, *iylm;

// ((MMAX+1)*(2*LMAX+2-MMAX))/2 -1

fftw_plan ifft, fft;		// plans de FFTW.
double fft_norm;			// normalisation required for the Direct Fourier Transform.

complex double *Slm;
complex double *ShF;
double *Sh;
	
double xg[NLAT];
double wg[NLAT];

double spat[NLAT];
double spec[LMMAX];

unsigned fftw_plan_mode = FFTW_PATIENT;		// defines the default FFTW planner mode.


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
		w[i-1] = 2.0*pi*w[i-1];			// normalize
		w[n-i] = w[i-1];
	}

// as we started with initial guesses, we should check if the gauss points are unique.
	for (i=n; i>0; i--) {
		if (x[i] == x[i-1]) {
			printf("ERROR : bad gauss points !\n");
			exit(1);
		}
	}
}

void write_vect(char *fn, double *vec, int N)
{
	FILE *fp;
	int i;
	
	fp = fopen(fn,"w");
	for (i=0;i<N;i++) {
		fprintf(fp,"%.6g ",vec[i]);
	}
	fclose(fp);
}

void runerr(const char * error_text)
{
	printf("*** Run-time error : %s\n",error_text);
	exit(1);
}

void planFFT()	// les FFTs. STRIDE = 1 (m's ranges en contigus.)
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
	ifft = fftw_plan_many_dft_c2r(1, &nfft, NLAT, ShF, &ncplx, 1, ncplx, Sh, &nreal, 1, nreal, fftw_plan_mode);
	if (ifft == NULL)
		runerr("[FFTW] ifft planning failed !");

// FFT : must be normalized.
	fft = fftw_plan_many_dft_r2c(1, &nfft, NLAT, Sh, &nreal, 1, nreal, ShF, &ncplx, 1, ncplx, fftw_plan_mode);
	if (fft == NULL)
		runerr("[FFTW] fft planning failed !");

	fft_norm = 1.0/nfft;
}

void init_HS()
{
	int i,m;
	
	Gauss(xg,wg,NLAT);	// generate gauss nodes and weights

// for synthesis (inverse transform)
	ylm = (double *) malloc(sizeof(double)* LMMAX*NLAT);
	dylm = (double *) malloc(sizeof(double)* LMMAX*NLAT);
	
	for (i=0;i<NLAT;i++) {
		for (m=0; m<=MMAX; m++) {
			gsl_sf_legendre_sphPlm_deriv_array(LMAX, m, xg[i], &ylm[i*LMMAX + LM(m,m)], &dylm[i*LMMAX + LM(m,m)]);	// fixed m legendre functions lookup table.
		}
	}
	
// for analysis (decomposition, direct transform) : transpose and multiply by gauss weight.
	iylm = (double *) malloc(sizeof(double)* LMMAX*NLAT);
	
	for (i=0;i<NLAT;i++) {
		for (m=0; m<LMMAX; m++) {
			iylm[m*NLAT + i] = ylm[i*LMMAX + m] * wg[i];
		}
	}
}


int main()
{

	double t,tmax;
	int i,m;

	init_HS();

	tmax = 0.0;
	for (i = 0; i<NLAT; i++) {
		t = gsl_sf_legendre_sphPlm(NLAT, 0, xg[i]);
		if (t>tmax) tmax = t;
//		printf("i=%d, x=%12.12g, p=%12.12g\n",i,xg[i],t);
	}
	printf("max zero : %g\n",tmax);

//	planFFT();

// test case :
	for (m=0;m<LMMAX;m++) spec[m] = 0.0;
	spec[0] = 1.0;
	spec[4] = 2.0;
	spec[3] = 0.5;
	spec[10] = 0.7;

	write_vect("spec0",spec,LMMAX);

// SYNTHESE (INVERSE)
	for (i=0;i<NLAT;i++) {
		spat[i] = 0.0;
		for (m=0;m<LMMAX;m++) {
			spat[i] += ylm[i*LMMAX + m] * spec[m];
		}
	}

	write_vect("spat",spat,LMMAX);
	write_vect("xg",xg,NLAT);
	write_vect("wg",wg,NLAT);

// DECOMPOSITION (DIRECT)
	for (m=0;m<LMMAX;m++) {
		spec[m] = 0.0;
		for (i=0;i<NLAT;i++) {
			spec[m] += wg[i]*ylm[i*LMMAX + m] * spat[i];
		}
	}

	write_vect("spec1",spec,LMMAX);

// DECOMPOSITION (OPTIMIZED)
	for (m=0;m<LMMAX;m++) {
		spec[m] = 0.0;
		for (i=0;i<NLAT;i++) {
			spec[m] += iylm[m*NLAT + i] * spat[i];
		}
	}

	write_vect("spec2",spec,LMMAX);

}

