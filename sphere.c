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

#include "SHT.c"
	
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

int main()
{
	complex double t1, t2;
	double t,tmax;
	int i,im,m,l,jj;

	init_SH();

	t1 = 1.0+2.0*I;
	t2 = 1.0-I;
	printf("test : %f, %f, %f, %f\n",creal(t1),cimag(t1), creal(t2),cimag(t2));

	(double) t1 = 8.0 +I;
	(double) t2 = 8.1;
	printf("test : %f, %f, %f, %f\n",creal(t1),cimag(t1), creal(t2),cimag(t2));

	write_vect("cost",ct,NLAT/2);
	write_vect("sint",st,NLAT/2);

	ShF = (complex double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Sh = (double *) ShF;
	ThF = (complex double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Th = (double *) ThF;
	NLF = (complex double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	NL = (double *) NLF;
/*
// test FFT :
	for(i=0;i<NLAT*(NPHI/2+1);i++) {
		ShF[i] = 0;	ThF[i] = 0;
	}
	ShF[0] = 1.0-I;
	ShF[NLAT] = 1.0-I;
	ShF[2*NLAT] = 2.0+I;
	ShF[5*NLAT] = 3.0+I;
	ShF[8*NLAT] = 5.0-I;
	ThF[0] = 1.0-I;
	ThF[3*NLAT] = 2.0-I;
	ThF[NLAT] = 3.0+I;
	ThF[6*NLAT] = 8.0-I;
	ThF[7*NLAT] = 1.0-3.*I;

	if (MMAX>0) {
		fftw_execute_dft_c2r(ifft,ShF,Sh);
		write_mx("sph",Sh,NPHI,NLAT);
		fftw_execute_dft_r2c(fft,Sh,ShF);
		write_mx("sphF",Sh,NPHI/2+1,2*NLAT);
	}

// compare FFT NL terms to pure spectral ...
	l = clock();
	for(jj=0;jj<10000;jj++) {
	for(i=0;i<NLAT*(NPHI/2+1);i++) {
		ShF[i] = 0;	ThF[i] = 0;
	}
	ShF[0] = 1.0-I;
	ShF[NLAT] = 1.0-I;
	ShF[2*NLAT] = 2.0+I;
	ShF[5*NLAT] = 3.0+I;
	ShF[8*NLAT] = 5.0-I;
	ThF[0] = 1.0-I;
	ThF[3*NLAT] = 2.0-I;
	ThF[NLAT] = 3.0+I;
	ThF[6*NLAT] = 8.0-I;
	ThF[7*NLAT] = 1.0-3.*I;

		fftw_execute_dft_c2r(ifft,ShF,Sh);
		fftw_execute_dft_c2r(ifft,ThF,Th);
		for(i=0;i<NLAT*NPHI;i++) {
			NL[i] = Sh[i] * Th[i];
		}
		fftw_execute_dft_r2c(fft,NL,NLF);
	}
	l = clock() - l;
	printf("fft nl time : %d\n",l);
		write_mx("NLfft",NL,NPHI/2+1,2*NLAT);

	if (MMAX>0) {
		fftw_execute_dft_c2r(ifft,ShF,Sh);
		write_mx("sph",Sh,NPHI,NLAT);
		fftw_execute_dft_r2c(fft,Sh,ShF);
		write_mx("sphF",Sh,NPHI/2+1,2*NLAT);
	}
	m = clock();
	for(jj=0;jj<10000;jj++) {
	for(i=0;i<NLAT*(NPHI/2+1);i++) {
		ShF[i] = 0;	ThF[i] = 0;
	}
	ShF[0] = 1.0-I;
	ShF[NLAT] = 1.0-I;
	ShF[2*NLAT] = 2.0+I;
	ShF[5*NLAT] = 3.0+I;
	ShF[8*NLAT] = 5.0-I;
	ThF[0] = 1.0-I;
	ThF[3*NLAT] = 2.0-I;
	ThF[NLAT] = 3.0+I;
	ThF[6*NLAT] = 8.0-I;
	ThF[7*NLAT] = 1.0-3.*I;
		NLspec(ShF,ThF, NLF);
	}
	m = clock() - m;
	printf("spectral nl time : %d\n",m);
	printf("fft/spectral time : %f\n",(double) l / (double) m);
		write_mx("NLspec",NL,NPHI/2+1,2*NLAT);
/*
// test Ylm :
	im = 0; l=0; m=im*MRES;
	write_vect("y00",&iylm[im][(l-m)*NLAT/2],NLAT/2);
	write_vect("dy00",&idylm[im][(l-m)*NLAT/2].t,NLAT);
	im = 0; l=1; m=im*MRES;
	write_vect("y10",&iylm[im][(l-m)*NLAT/2],NLAT/2);
	write_vect("dty10",&idylm[im][(l-m)*NLAT/2].t,NLAT);
	im = 0; l=LMAX; m=im*MRES;
	write_vect("yLmax0",&iylm[im][(l-m)*NLAT/2],NLAT/2);
	write_vect("dtyLmax0",&idylm[im][(l-m)*NLAT/2].t,NLAT);
	im = MMAX; m=im*MRES; l=m;
	write_vect("ymmax-mmax",&iylm[im][(l-m)*NLAT/2],NLAT/2);
	write_vect("dtymmax-mmax",&idylm[im][(l-m)*NLAT/2].t,NLAT);
	im = 3; l=8; m=im*MRES;
	write_vect("y83",&iylm[im][(l-m)*NLAT/2],NLAT/2);
	write_vect("dty83",&idylm[im][(l-m)*NLAT/2].t,NLAT);
	im = 10; l=35; m=im*MRES;
	write_vect("y3510",&iylm[im][(l-m)*NLAT/2],NLAT/2);
	write_vect("dty3510",&idylm[im][(l-m)*NLAT/2].t,NLAT);

*/
// test case...
	Slm = (complex double *) malloc(sizeof(complex double)* NLM);
	for (i=0;i<NLM;i++) {
		Slm[i] = 0.0;
	}
	
	Slm[LM(0,0)] = 1.0;
	Slm[LM(1,1)] = 3.0;
	Slm[LM(4,2)] = 2.0+I;
// 	Slm[LM(10,4)] = -4.0-I;
// 	Slm[LM(55,12)] = 5.0-2.0*I;

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

	write_vect("ylm",Slm,NLM*2);

}

