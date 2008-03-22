// base flow. Cartesien spectral.

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
// FFTW la derivee d/dx = ik	(pas de moins !)
#include <fftw3.h>
// GSL for Legendre functions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>

#include <time.h>


complex double *Slm;	// spherical harmonics l,m space
complex double *ShF, *ThF, *NLF;	// Fourier space : theta,m
double *Sh, *Th, *NL;		// real space : theta,phi (alias of ShF)


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

#include "SHTfast.c"

void write_ylm(char *fn, complex double *y)
{
	FILE *fp;
	int im,m,l;
	
	fp = fopen(fn,"w");
	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
		fprintf(fp,"m=%d ::",m);
		for (l=m;l<=LMAX;l++) {
			fprintf(fp," %.6g,%.6g",creal(y[LM(l,m)]), cimag(y[LM(l,m)]));
		}
		fprintf(fp,"\n");
	}
}

int main()
{
	complex double t1, t2;
	double t,tmax;
	int i,it,im,m,l,jj;
	FILE *fp;

	init_SH();

#ifdef _SH_DEBUG_
	fp = fopen("cheby","w");
	for (it=0;it<NLAT; it++)
		fprintf(fp,"%d\t%f\t%f\t%f\n", it, ct[it], wc[it], wc[it+NLAT]);
#endif
	
	ShF = (complex double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Sh = (double *) ShF;
	ThF = (complex double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Th = (double *) ThF;
	NLF = (complex double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	NL = (double *) NLF;

	for (i=0;i<NLAT*NPHI;i++) {
		Sh[i] =0.0;
	}
	Sh[0] = 1.0;
	Sh[2] = 2.0;
	Sh[3] = 3.0;
	write_mx("DCTa",Sh,NPHI,NLAT);
	fftw_execute_r2r(idct,Sh, Sh);		// DCT of weighted data
	write_mx("DCTb",Sh,NPHI,NLAT);
	
	fftw_execute_r2r(dct,Sh, Sh);		// DCT of weighted data
	write_mx("DCTc",Sh,NPHI,NLAT);
	
	

// test case...
	Slm = (complex double *) malloc(sizeof(complex double)* NLM);
	for (i=0;i<NLM;i++)
		Slm[i] = 0.0;
	
	Slm[LM(0,0)] = 1.0;
//	Slm[LM(1,1)] = 3.0;
	Slm[LM(4,2)] = 2.0+I;
// 	Slm[LM(10,4)] = -4.0-I;
// 	Slm[LM(55,12)] = 5.0-2.0*I;

	write_ylm("ylm0",Slm);
	m = clock();
	for (jj=0;jj<1000;jj++) {

// synthese (inverse legendre)
		SH_to_spat(Slm,ShF);
		SH_to_spat(Slm,ThF);
		for (i=0; i< NLAT*NPHI; i++) {
			ThF[i] *= ShF[i];
		}
// analyse (direct legendre)
		spat_to_SH(ShF,Slm);
	}
	m = clock() - m;
	printf("SHT + iSHT time : %d\n",m);

	write_ylm("ylm",Slm);


	
}

