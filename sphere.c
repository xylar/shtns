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

#define LM(l,m) ((m*(2*LMAX +1-m))/2 + l)

double pi = atan(1)*4;
double *ylm;
double *dylm;

// ((MMAX+1)*(2*LMAX+2-MMAX))/2 -1

double spat[NLAT];
double spec[LMMAX];


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


int main()
{
	double xg[NLAT];
	double wg[NLAT];

	double t,tmax;
	int i,m;

	Gauss(xg,wg,NLAT);

	tmax = 0.0;
	for (i = 0; i<NLAT; i++) {
		t = gsl_sf_legendre_sphPlm(NLAT, 0, xg[i]);
		if (t>tmax) tmax = t;
//		printf("i=%d, x=%12.12g, p=%12.12g\n",i,xg[i],t);
	}
	printf("max zero : %g\n",tmax);

	ylm = (double *) malloc(sizeof(double)* LMMAX*NLAT);
	dylm = (double *) malloc(sizeof(double)* LMMAX*NLAT);

	m = 0;
	for (i=0;i<NLAT;i++) {
	 	gsl_sf_legendre_sphPlm_deriv_array(LMAX, m, xg[i], &ylm[i*LMMAX], &dylm[i*LMMAX]);	// table des fonctions de legendre à m fixé.
	}

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

}

