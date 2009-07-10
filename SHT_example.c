/** \example SHT_example.c
 \brief An example program that performs backward and forward Spherical Harmonic Transforms using SHTns
**/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
// FFTW la derivee d/dx = ik	(pas de moins !)
#include <fftw3.h>
// GSL for Legendre functions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>



#include "SHT.h"

/// a simple function that writes a vector to a file
void write_vect(char *fn, double *vec, int N);

/// a simple function that writes a matrix to a file.
void write_mx(char *fn, double *mx, int N1, int N2);

// polar optimization threshold
#define POLAR_OPT_THR 0

int main()
{
	complex double *Slm, *Tlm;	// spherical harmonics l,m space
	double *Sh, *Th;		// real space : theta,phi
	long int i,im,lm;

//	init_SH(shtns_type, eps,         lmax, mmax, mres, nlat, nphi);
	init_SH( sht_gauss, POLAR_OPT_THR, 340, 3,    5,   512,  12 );

	Sh = (double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Th = (double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));

	Slm = (complex double *) malloc(sizeof(complex double)* NLM);
	Tlm = (complex double *) malloc(sizeof(complex double)* NLM);


// SH_to_spat
	LM_LOOP( Slm[lm]=0.0;  Tlm[lm] = 0.0; )		/* this is the same as :
													for (lm=0;lm<NLM;lm++) {
														Slm[lm] = 0.0;	Tlm[lm] = 0.0;
													} */

	Slm[LM(1,1)] = 1;				// access to SH coefficient
	Slm[LiM(1,0)] = Y10_ct;
	SH_to_spat(Slm,Sh);

	write_mx("spat",Sh,NPHI,NLAT);
	SHsphtor_to_spat(Slm,Tlm,Sh,Th);		// vector transform
	write_mx("spatt",Sh,NPHI,NLAT);
	write_mx("spatp",Th,NPHI,NLAT);

// spat_to_SH
	for (im=0;im<NPHI;im++) {
		for (i=0;i<NLAT;i++) {
			Sh[im*NLAT+i] = ct[i];			// use cos(theta) array
		}
	}
	spat_to_SH(Sh,Slm);
	write_vect("ylm",Slm,NLM*2);

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

