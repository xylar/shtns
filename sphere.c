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


// LMAX : maximum degree of Spherical Harmonic
#define LMAX 79
// MMAX : max number of fourrier decomposition (degree = MMAX * MRES)
#define MMAX 32
// MRES : azimutal symmetry
#define MRES 2
// LMMAX : total number of Spherical Harmonic coefficients.
//#define LMMAX ( ((MMAX+1)*(2*LMAX+2-MMAX))/2 )
#define LMMAX ( (MMAX+1)*(LMAX+1) - MRES* MMAX*(MMAX+1)/2 )
// NLAT : number of latitudinal (theta) gauss points, at least LMAX+1, should be EVEN
#define NLAT (LMAX+1)
// NPHI : number of azimutal grid points, at least MMAX*2
#define NPHI (MMAX*2)

// LM(l,m) : index in the Spherical Harmonic coefficient array [ (l,m) space ]
//#define LM(l,m) ((m*(2*LMAX +1-m))/2 + l)
#define LM(l,m) ( (m/MRES)*(2*LMAX+3 -m)/2 + l-m )
#define LiM(l,im) ( im*(2*LMAX+3 -(im+2)*MRES)/2 + l )

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

// as we started with initial guesses, we should check if the gauss points are actually unique.
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

void planFFT()	// les FFTs. STRIDE = NLAT (l ranges en contigus.)
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

int tm[MMAX+1];		// start theta value for spherical harmonics. (because at the pole the legendre polynoms go to zero for high m's

void init_HS()
{
	int it,im,m,l;
	double eps= 1.e-14;	// precision : polar coefficients below that threshold are neglected (for high ms)

	if (MMAX*MRES > LMAX) runerr("[init_HS] resolution mismatch : MMAX*MRES should not exceed LMAX");
	if (NLAT <= LMAX) runerr("[init_HS] resolution mismatch : NLAT should be at least LMAX+1");
	
	planFFT();		// initialize fftw
	Gauss(xg,wg,NLAT);	// generate gauss nodes and weights

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
				iylm[im][(l-m)*NLAT + it] = ylm[im][it*(LMAX-m+1) + (l-m)] * wg[it] *2.0*pi* fft_norm;
			}
		}
	}

// POLAR OPTIMIZATION : analysing coefficients, some can be neglected.
	for (im=0;im<=MMAX;im++) {
		m = im*MRES;
		tm[im] = 1.0*NLAT;
		for (l=m;l<=LMAX;l++) {
			it=0;
			while( ylm[im][it*(LMAX-m+1) + (l-m)] < eps ) { it++; }
			if (tm[im] > it) tm[im] = it;
		}
//		tm[im] = 0;	// this would cancel the optimisation.
	}

//	write_vect("tm",tm,MMAX+1);
	for (im=0;im<=MMAX;im++)
		printf("%d ",tm[im]);

}



int main()
{
	complex double fp[NLAT/2];	// partie symetrique
	complex double fm[NLAT/2];	// partie anti-sym

	double t,tmax;
	int i,im,m,l,jj;

	init_HS();

// TEST if gauss points are ok.
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
	Slm[LM(3,1)] = 3.0;
	Slm[LM(3,3)] = 2.0;
	Slm[LM(10,5)] = 4.0;
	Slm[LM(55,12)] = 5.0;
	
for (jj=0;jj<3000;jj++) {

// synthese (inverse legendre)
	for (im=0;im<=MMAX;im++) {
		m = im*MRES;
/*
// Initial version without symmetric/antisymmetric optimization
		for (i=0;i<NLAT;i++) {
			ShF[im*NLAT +i] = 0.0;
		}
		for (i=tm[im];i<NLAT-tm[im];i++) {	// ops : 2*(lmax-m+1)*NLAT
			for(l=m;l<=LMAX;l++) {
				ShF[im*NLAT +i] += ylm[im][i*(LMAX-m+1) + (l-m)] * Slm[LM(l,m)];
			}
		}
*/

/* EVEN - ODD optimized
		for (i=0;i<NLAT;i++) {
			ShF[im*NLAT +i] = 0.0;
		}
		for (i=tm[im];i<NLAT/2;i++) {		// ops : 3*(lmax-m+1)*NLAT/2	= 25% faster
			for(l=m;l<=LMAX;l+=2) {		// EVEN
				t = ylm[im][i*(LMAX-m+1) + (l-m)] * Slm[LM(l,m)];
				ShF[im*NLAT +i] += t;
				ShF[im*NLAT +NLAT-(i+1)] += t;
			}
			for(l=m+1;l<=LMAX;l+=2) {	// ODD
				t = ylm[im][i*(LMAX-m+1) + (l-m)] * Slm[LM(l,m)];
				ShF[im*NLAT +i] += t;
				ShF[im*NLAT +NLAT-(i+1)] -= t;
			}
		}
*/
		i=0;
		while (i<tm[im]) {	// polar optimization
			ShF[im*NLAT + i] = 0.0;
//			ShF[im*NLAT + NLAT-(i+1)] = 0.0;	// south pole zeroes
			ShF[im*NLAT + NLAT-tm[im] + i] = 0.0;	// south pole zeroes
			i++;
		}
		while (i<NLAT/2) {	// NLAT/2 * [ (lmax-m+1)*2 + 6]	: almost twice as fast.
			fp[i] = 0.0;	fm[i] = 0.0;
			l=m;
			while (l<LMAX) {	// compute even and odd parts
				fp[i] += ylm[im][i*(LMAX-m+1) + (l-m)] * Slm[LM(l,m)];
				fm[i] += ylm[im][i*(LMAX-m+1) + (l+1-m)] * Slm[LM(l+1,m)];
				l+=2;
			}
			if (l==LMAX) {
				fp[i] += ylm[im][i*(LMAX-m+1) + (l-m)] * Slm[LM(l,m)];
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
//	write_mx("sph",Sh,NPHI,NLAT);


// analyse (direct legendre)
	fftw_execute(fft);
	for (im=0;im<=MMAX;im++) {
		m=im*MRES;

/*		UNOPTIMIZED
		for (l=m;l<=LMAX;l++) {
			Slm[LM(l,m)] = 0.0;
//			for (i=0;i<NLAT;i++) {
			for (i=tm[im];i<NLAT-tm[im];i++) {
				Slm[LM(l,m)] += iylm[im][(l-m)*NLAT + i] * ShF[im*NLAT + i];
			}
		}
*/

		for (i=tm[im];i<NLAT/2;i++) {	// compute symmetric and antisymmetric parts.
			fp[i] = ShF[im*NLAT + i] + ShF[im*NLAT + NLAT-(i+1)];
			fm[i] = ShF[im*NLAT + i] - ShF[im*NLAT + NLAT-(i+1)];
		}

		l=m;
		while (l<LMAX) {
			Slm[LM(l,m)] = 0.0;	Slm[LM(l+1,m)] = 0.0;
			for (i=tm[im];i<NLAT/2;i++) {	// tm[im] : polar optimization
				Slm[LM(l,m)] += iylm[im][(l-m)*NLAT + i] * fp[i];
			}
			for (i=tm[im];i<NLAT/2;i++) {	// tm[im] : polar optimization
				Slm[LM(l+1,m)] += iylm[im][(l+1-m)*NLAT + i] * fm[i];
			}
			l+=2;
		}
		if (l==LMAX) {
			Slm[LM(l,m)] = 0.0;
			for (i=tm[im];i<NLAT/2;i++) {	// polar optimization
				Slm[LM(l,m)] += iylm[im][(l-m)*NLAT + i] * fp[i];
			}
		}
/*	alternative : on calcule d'abord les pairs puis les impairs, mais un peu moins bon au niveau du cache...
		for (l=m;l<=LMAX;l+=2) {	// l-m pair (0, 2 ...)
			Slm[LM(l,m)] = 0.0;
//			for (i=0;i<NLAT;i++) {
			for (i=tm[im];i<NLAT/2;i++) {	// polar optimization
				Slm[LM(l,m)] += iylm[im][(l-m)*NLAT + i] * fp[i];
			}
		}
		for (l=m+1;l<=LMAX;l+=2) {	// l-m impair (1, 3, ...)
			Slm[LM(l,m)] = 0.0;
//			for (i=0;i<NLAT;i++) {
			for (i=tm[im];i<NLAT/2;i++) {	// polar optimization
				Slm[LM(l,m)] += iylm[im][(l-m)*NLAT + i] * fm[i];
			}
		}
*/
	}

}

	write_vect("ylm",Slm,LMMAX*2);

}

