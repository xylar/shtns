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
#define MMAX 32
#define LMMAX ( ((MMAX+1)*(2*LMAX+2-MMAX))/2 )
#define NLAT (LMAX+1)
#define NPHI (MMAX*2)

#define LM(l,m) ((m*(2*LMAX +1-m))/2 + l)

double pi = atan(1)*4;

double* ylm[MMAX+1];
double* iylm[MMAX+1];

// ((MMAX+1)*(2*LMAX+2-MMAX))/2 -1

fftw_plan ifft, fft;		// plans de FFTW.
double fft_norm;			// normalisation required for the Direct Fourier Transform.

complex double *Slm;	// spherical harmonics l,m space
complex double *ShF;	// Fourier space : theta,m
double *Sh;				// real space : theta,phi
	
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

void write_mx(char *fn, double *mx, int N1, int N2)
{
	FILE *fp;
	int i,j;
	
	fp = fopen(fn,"w");
	for (i=0;i<N1;i++) {
		for(j=0;j<N2;j++) {
			fprintf(fp,"%.6g ",mx[i*N2+j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}


void runerr(const char * error_text)
{
	printf("*** Run-time error : %s\n",error_text);
	exit(1);
}

void planFFT()	// les FFTs. STRIDE = NLAT (m's ranges en contigus.)
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

	fft_norm = 1.0/nfft;
}

void init_HS()
{
	int i,m,l;
	
	planFFT();
	Gauss(xg,wg,NLAT);	// generate gauss nodes and weights

// for synthesis (inverse transform)
	for (m=0; m<=MMAX; m++) {
		ylm[m] = (double *) malloc(sizeof(double)* (LMAX-m+1)*NLAT);
		for (i=0;i<NLAT;i++) {
//			gsl_sf_legendre_sphPlm_deriv_array(LMAX, m, xg[i], ylm[m], &dylm[i*LMMAX + LM(m,m)]);	// fixed m legendre functions lookup table.
			gsl_sf_legendre_sphPlm_array(LMAX, m, xg[i], ylm[m] + i*(LMAX-m+1));	// fixed m legendre functions lookup table.
		}
	}
	
// for analysis (decomposition, direct transform) : transpose and multiply by gauss weight.
	for (m=0; m<=MMAX; m++) {
		iylm[m] = (double *) malloc(sizeof(double)* (LMAX-m+1)*NLAT);
		for (i=0;i<NLAT;i++) {
			for (l=m;l<=LMAX;l++) {
				iylm[m][(l-m)*NLAT + i] = ylm[m][i*(LMAX-m+1) + (l-m)] * wg[i] *2.0*pi* fft_norm;
			}
		}
	}
}

	int tm[MMAX+1];


analyse_HS()
{
	int i,m,l;
	double eps= 1.e-6;
	
	for (m=0;m<=MMAX;m++) {
		tm[m] = 1.0*NLAT;
		for (l=m;l<=LMAX;l++) {
			i=0;
			while( ylm[m][i*(LMAX-m+1) + (l-m)] < eps ) { i++; }
			if (tm[m] > i) tm[m] = i;
		}
//		tm[m] = 0;
	}

//	write_vect("tm",tm,MMAX+1);
}

int main()
{

	double t,tmax;
	int i,m,l,jj;

	init_HS();
	analyse_HS();

	tmax = 0.0;
	for (i = 0; i<NLAT; i++) {
		t = gsl_sf_legendre_sphPlm(NLAT, 0, xg[i]);
		if (t>tmax) tmax = t;
//		printf("i=%d, x=%12.12g, p=%12.12g\n",i,xg[i],t);
	}
	printf("max zero : %g\n",tmax);

// test FFT :
	for(i=0;i<NLAT*(NPHI/2+1);i++) {
		ShF[i] = 0;
	}
	ShF[0] = 1.0+I;
	ShF[NLAT] = 2.0+I;
	ShF[NLAT*2] = 3.0+I;
	
	fftw_execute(ifft);
	write_mx("sph",Sh,NPHI,NLAT);
	fftw_execute(fft);
	write_mx("sphF",Sh,NPHI/2+1,2*NLAT);

// test case...
	Slm = (complex double *) malloc(sizeof(complex double)* LMMAX);
	for (i=0;i<LMMAX;i++) {
		Slm[i] = 0.0;
	}
	
	Slm[LM(0,0)] = 1.0;
	Slm[LM(3,3)] = 2.0;
	Slm[LM(3,1)] = 3.0;
	Slm[LM(10,5)] = 4.0;
	Slm[LM(55,30)] = 5.0;
	
for (jj=0;jj<1000;jj++) {

// synthese (inverse legendre)
	for (m=0;m<=MMAX;m++) {
		for (i=0;i<NLAT;i++) {
			ShF[m*NLAT +i] = 0.0;
		}
		for (i=tm[m];i<NLAT-tm[m];i++) {	// optimize
			for(l=m;l<=LMAX;l++) {
				ShF[m*NLAT +i] += ylm[m][i*(LMAX-m+1) + (l-m)] * Slm[LM(l,m)];
			}
		}
	}
	fftw_execute(ifft);
//	write_mx("sph",Sh,NPHI,NLAT);

// analyse (direct legendre)
	fftw_execute(fft);
	for (m=0;m<=MMAX;m++) {
		for (l=m;l<=LMAX;l++) {
			Slm[LM(l,m)] = 0.0;
//			for (i=0;i<NLAT;i++) {
			for (i=tm[m];i<NLAT-tm[m];i++) {
				Slm[LM(l,m)] += iylm[m][(l-m)*NLAT + i] * ShF[m*NLAT + i];
			}
		}
	}
//	write_vect("ylm",Slm,LMMAX*2);

}

}

