/*
 * Copyright (c) 2010-2011 Centre National de la Recherche Scientifique.
 * written by Nathanael Schaeffer (CNRS, ISTerre, Grenoble, France).
 * 
 * nathanael.schaeffer@ujf-grenoble.fr
 * 
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software. You can use,
 * modify and/or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 * 
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * 
 */

/// \file time_SHT.c This program performs some spherical harmonic transforms, and displays timings and accuracy.
/// \c make \c time_SHT to compile, and then run.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>

// cycle counter from FFTW
#include "cycle.h"

#include "shtns.h"

complex double *Slm, *Slm0, *Tlm, *Tlm0, *Qlm;	// spherical harmonics l,m space
complex double *ShF, *ThF, *NLF;	// Fourier space : theta,m
double *Sh, *Th, *NL;		// real space : theta,phi (alias of ShF)

int LMAX,MMAX,MRES,NLM;
int NLAT = 0;
int NPHI = 0;

// number of SH iterations
long int SHT_ITER = 50;		// do 50 iterations by default

	
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

double scal_error(complex double *Slm, complex double *Slm0, int ltr)
{
	long int jj,i, nlm_cplx;
	double tmax,t,n2;

	nlm_cplx = (MMAX*2 == NPHI) ? LiM(MRES*MMAX,MMAX) : NLM;
// compute error :
	tmax = 0;	n2 = 0;		jj=0;
	for (i=0;i<NLM;i++) {
	  if (li[i] <= ltr) {
		if ((i <= LMAX)||(i >= nlm_cplx)) {		// m=0, and 2*m=nphi is real
			Slm[i] = creal(Slm[i]-Slm0[i]);
			t = fabs(creal(Slm[i]));
		} else {
			Slm[i] -= Slm0[i];
			t = cabs(Slm[i]);
		}
		n2 += t*t;
		if (t>tmax) { tmax = t; jj = i; }
	  }
	}
	printf("   => max error = %g (l=%.0f,lm=%d)   rms error = %g",tmax,el[jj],jj,sqrt(n2/NLM));
	if (tmax > 1e-3) {
		if (NLM < 15) {
			printf("\n orig:");
			for (i=0; i<NLM;i++)
				if ((i <= LMAX)||(i >= nlm_cplx)) {		// m=0, and 2*m=nphi is real
					printf("  %g",creal(Slm0[i]));
				} else {
					printf("  %g,%g",creal(Slm0[i]),cimag(Slm0[i]));
				}
			printf("\n diff:");
			for (i=0; i<NLM;i++)
				if ((i <= LMAX)||(i >= nlm_cplx)) {		// m=0, and 2*m=nphi is real
					printf("  %g",creal(Slm[i]));
				} else {
					printf("  %g,%g",creal(Slm[i]),cimag(Slm[i]));
				}
		}
		printf("    **** ERROR ****\n");
	}
	else printf("\n");
	return(tmax);
}

double vect_error(complex double *Slm, complex double *Tlm, complex double *Slm0, complex double *Tlm0, int ltr)
{
	long int jj,i, nlm_cplx;
	double tmax,t,n2;

// compute error :
	tmax = 0;	n2 = 0;		jj=0;
	for (i=0;i<NLM;i++) {
	  if (li[i] <= ltr) {
		if ((i <= LMAX)||(i >= LiM(MRES*(NPHI+1)/2,(NPHI+1)/2))) {
			Slm[i] = creal(Slm[i]-Slm0[i]);
			t = fabs(creal(Slm[i]));
		} else {
			Slm[i] -= Slm0[i];
			t = cabs(Slm[i]);
		}
		n2 += t*t;
		if (t>tmax) { tmax = t; jj = i; }
	  }
	}
	printf("   Spheroidal => max error = %g (l=%.0f,lm=%d)    rms error = %g",tmax,el[jj],jj,sqrt(n2/NLM));
	if (tmax > 1e-3) { printf("    **** ERROR ****\n"); }
		else printf("\n");
//	write_vect("Slm",Slm,NLM*2);

// compute error :
	tmax = 0;	n2 = 0;		jj=0;
	for (i=0;i<NLM;i++) {
	  if (li[i] <= ltr) {
		if ((i <= LMAX)||(i >= LiM(MRES*(NPHI+1)/2,(NPHI+1)/2))) {
			Tlm[i] = creal(Tlm[i]- Tlm0[i]);
			t = fabs(creal(Tlm[i]));
		} else {
			Tlm[i] -= Tlm0[i];
			t = cabs(Tlm[i]);
		}
		n2 += t*t;
		if (t>tmax) { tmax = t; jj = i; }
	  }
	}
	printf("   Toroidal => max error = %g (l=%.0f,lm=%d)    rms error = %g",tmax,el[jj],jj,sqrt(n2/NLM));
	if (tmax > 1e-3) { printf("    **** ERROR ****\n"); }
		else printf("\n");
//	write_vect("Tlm",Tlm,NLM*2);
}

int test_SHT_fly()
{
	long int jj,i, nlm_cplx;
	clock_t tcpu;
	double ts,ta;

	for (i=0;i<NLM;i++) Slm[i] = Slm0[i];	// restore test case...

	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
		SH_to_spat_fly(Slm,Sh);
	}
	tcpu = clock() - tcpu;
	ts = tcpu / (1000.*SHT_ITER);

	tcpu = clock();
	spat_to_SH_fly(Sh,Slm);
	for (jj=1; jj< SHT_ITER; jj++) {
		spat_to_SH_fly(Sh,Tlm);
	}
	tcpu = clock() - tcpu;
	ta = tcpu / (1000.*SHT_ITER);
	printf("   SHT time on-the-fly : \t synthesis = %f ms \t analysis = %f ms\n", ts, ta);

	scal_error(Slm, Slm0, LMAX);
	return (int) tcpu;
}


int test_SHT()
{
	long int jj,i, nlm_cplx;
	clock_t tcpu;
	double ts, ta;

	for (i=0;i<NLM;i++) Slm[i] = Slm0[i];	// restore test case...

	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
		SH_to_spat(Slm,Sh);
	}
	tcpu = clock() - tcpu;
	ts = tcpu / (1000.*SHT_ITER);

	tcpu = clock();
	spat_to_SH(Sh,Slm);
	for (jj=1; jj< SHT_ITER; jj++) {
		spat_to_SH(Sh,Tlm);
	}
	tcpu = clock() - tcpu;
	ta = tcpu / (1000.*SHT_ITER);
	printf("   SHT time : \t synthesis = %f ms \t analysis = %f ms\n", ts, ta);

	scal_error(Slm, Slm0, LMAX);
	return (int) tcpu;
}

int test_SHT_parity(int eo)
{
	long int jj,i, nlm_cplx;
	clock_t tcpu;

	for (i=0;i<NLM;i++) Slm[i] = Slm0[i];	// restore test case...
	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
// synthese (inverse legendre)
		SHeo_to_spat(Slm,Sh,eo);
		spat_to_SHeo(Sh,Slm,eo);
	}
	tcpu = clock() - tcpu;
	printf("   2iSHT + NL + SHT x%d time : %d ms\n", SHT_ITER, (int )tcpu / 1000);

	scal_error(Slm, Slm0, LMAX);
	return (int) tcpu;
}

int test_SHT_l(int ltr)
{
	int jj,i;
	clock_t tcpu;
	double ts, ta;

	for (i=0;i<NLM;i++) Slm[i] = Slm0[i];	// restore test case...

	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
		SH_to_spat_l(Slm,Sh,ltr);
	}
	tcpu = clock() - tcpu;
	ts = tcpu / (1000.*SHT_ITER);

	tcpu = clock();
		spat_to_SH_l(Sh,Slm,ltr);
	for (jj=1; jj< SHT_ITER; jj++) {
		spat_to_SH_l(Sh,Tlm,ltr);
	}
	tcpu = clock() - tcpu;
	ta = tcpu / (1000.*SHT_ITER);
	printf("   SHT time truncated at l=%d : synthesis = %f ms, analysis = %f ms\n", ltr, ts, ta);

	scal_error(Slm, Slm0, ltr);
	return (int) tcpu;
}

int test_SHT_vect_l(int ltr)
{
	int jj,i;
	clock_t tcpu;
	double ts, ta;

	complex double *S2 = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	complex double *T2 = (complex double *) fftw_malloc(sizeof(complex double)* NLM);

	for (i=0;i<NLM;i++) {
		Slm[i] = Slm0[i];	Tlm[i] = Tlm0[i];
	}
	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
		SHsphtor_to_spat_l(Slm,Tlm,Sh,Th,ltr);
	}
	tcpu = clock() - tcpu;
	ts = tcpu / (1000.*SHT_ITER);
	tcpu = clock();
		spat_to_SHsphtor_l(Sh,Th,Slm,Tlm, ltr);
	for (jj=1; jj< SHT_ITER; jj++) {
		spat_to_SHsphtor_l(Sh,Th,S2,T2, ltr);
	}
	tcpu = clock() - tcpu;
	ta = tcpu / (1000.*SHT_ITER);
	printf("   vector SHT time trucated at l=%d : \t synthesis %f ms \t analysis %f ms\n", ltr, ts, ta);

	fftw_free(T2);	fftw_free(S2);
	vect_error(Slm, Tlm, Slm0, Tlm0, ltr);
	return (int) tcpu;
}

int test_SHT_vect()
{
	int jj,i;
	clock_t tcpu;
	double ts, ta;

	complex double *S2 = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	complex double *T2 = (complex double *) fftw_malloc(sizeof(complex double)* NLM);

	for (i=0;i<NLM;i++) {
		Slm[i] = Slm0[i];	Tlm[i] = Tlm0[i];
	}
	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
		SHsphtor_to_spat(Slm,Tlm,Sh,Th);
	}
	tcpu = clock() - tcpu;
	ts = tcpu / (1000.*SHT_ITER);
	tcpu = clock();
		spat_to_SHsphtor(Sh,Th,Slm,Tlm);
	for (jj=1; jj< SHT_ITER; jj++) {
		spat_to_SHsphtor(Sh,Th,S2,T2);
	}
	tcpu = clock() - tcpu;
	ta = tcpu / (1000.*SHT_ITER);
	printf("   vector SHT time : \t synthesis %f ms \t analysis %f ms\n", ts, ta);

	fftw_free(T2);	fftw_free(S2);
	vect_error(Slm, Tlm, Slm0, Tlm0, LMAX);
	return (int) tcpu;
}

int test_SHT_vect_fly()
{
	int jj,i;
	clock_t tcpu;
	double ts;

	for (i=0;i<NLM;i++) {
		Slm[i] = Slm0[i];	Tlm[i] = Tlm0[i];
	}
	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
		SHsphtor_to_spat_fly(Slm,Tlm,Sh,Th);
	}
	tcpu = clock() - tcpu;
	ts = tcpu / (1000.*SHT_ITER);
	printf("   vector SHT time on-the-fly : \t synthesis %f ms\n", ts);

	spat_to_SHsphtor(Sh,Th,Slm,Tlm);
	vect_error(Slm, Tlm, Slm0, Tlm0, LMAX);
	return (int) tcpu;
}

int test_SHT_vect3d()
{
	int jj,i;
	clock_t tcpu;
	double ts, ta;
	
	complex double *Q2 = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	complex double *S2 = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	complex double *T2 = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	
	for (i=0;i<NLM;i++) {
		Slm[i] = Slm0[i];	Tlm[i] = Tlm0[i];	Qlm[i] = Tlm0[i];
	}

	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
		SHqst_to_spat(Qlm,Slm,Tlm,NL,Sh,Th);
	}
	tcpu = clock() - tcpu;
	ts = tcpu / (1000.*SHT_ITER);
	tcpu = clock();
		spat_to_SHqst(NL,Sh,Th,Qlm,Slm,Tlm);
	for (jj=1; jj< SHT_ITER; jj++) {
		spat_to_SHqst(NL,Sh,Th,Q2,S2,T2);
	}
	tcpu = clock() - tcpu;
	ta = tcpu / (1000.*SHT_ITER);
	printf("   3D vector SHT time : \t synthesis %f ms \t analysis %f ms\n", ts, ta);

	vect_error(Slm, Tlm, Slm0, Tlm0, LMAX);
	scal_error(Qlm, Tlm0, LMAX);
	return (int) tcpu;
}

int test_SHT_vect3d_fly()
{
	int jj,i;
	clock_t tcpu;
	double ts;
	
	for (i=0;i<NLM;i++) {
		Slm[i] = Slm0[i];	Tlm[i] = Tlm0[i];	Qlm[i] = Tlm0[i];
	}

	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
		SHqst_to_spat_fly(Qlm,Slm,Tlm,NL,Sh,Th);
	}
	tcpu = clock() - tcpu;
	ts = tcpu / (1000.*SHT_ITER);
	printf("   3D vector SHT time on-the-fly : \t synthesis %f ms\n", ts);	
	spat_to_SHqst(NL,Sh,Th,Qlm,Slm,Tlm);

	vect_error(Slm, Tlm, Slm0, Tlm0, LMAX);
	scal_error(Qlm, Tlm0, LMAX);
	return (int) tcpu;
}


/*
fftw_plan ifft_in, ifft_out;
fftw_plan fft_in, fft_out;
fftw_plan fft_tr, ifft_tr;


// we want to test if in-place is faster than out-of place or not.
init_fft_tests()
{
	complex double *ShF, *Shout;
	double *Sh;
	int nfft, ncplx, nreal;
	unsigned fftw_plan_mode = FFTW_EXHAUSTIVE;		// defines the default FFTW planner mode.

	nfft = NPHI;
	ncplx = NPHI/2 +1;
	nreal = 2*ncplx;

// Allocate dummy Spatial Fields.
	ShF = (complex double *) fftw_malloc(ncplx * NLAT * sizeof(complex double));
	Sh = (double *) ShF;

// IFFT : unnormalized
	ifft_in = fftw_plan_many_dft_c2r(1, &nfft, NLAT, ShF, &ncplx, NLAT, 1, Sh, &nreal, NLAT, 1, fftw_plan_mode);
	if (ifft_in == NULL) printf("ifft_in failed\n");
// FFT : must be normalized.
	fft_in = fftw_plan_many_dft_r2c(1, &nfft, NLAT, Sh, &nreal, NLAT, 1, ShF, &ncplx, NLAT, 1, fftw_plan_mode);
	if (fft_in == NULL) printf("fft_in failed\n");
printf("in-place done\n");
	printf("** ifft in-place :\n");	fftw_print_plan(ifft_in);
	printf("\n** fft in-place :\n");	fftw_print_plan(fft_in);

	Shout = (complex double *) fftw_malloc(ncplx * NLAT * sizeof(complex double));
	ifft_out = fftw_plan_many_dft_c2r(1, &nfft, NLAT, Shout, &ncplx, NLAT, 1, Sh, &nfft, NLAT, 1, fftw_plan_mode);
	if (ifft_out == NULL) printf("ifft_out failed\n");
	fft_out = fftw_plan_many_dft_r2c(1, &nfft, NLAT, Sh, &nfft, NLAT, 1, Shout, &ncplx, NLAT, 1, fftw_plan_mode);
	if (fft_out == NULL) printf("fft_out failed\n");
printf("\nout-of-place done\n");
	printf("** ifft out-of-place :\n");	fftw_print_plan(ifft_out);
	printf("\n** fft out-of-place :\n");	fftw_print_plan(fft_out);

	ifft_tr = fftw_plan_many_dft_c2r(1, &nfft, NLAT, Shout, &ncplx, NLAT, 1, Sh, &nfft, 1, NPHI, fftw_plan_mode);
	if (ifft_out == NULL) printf("ifft_out failed\n");
	fft_tr = fftw_plan_many_dft_r2c(1, &nfft, NLAT, Sh, &nfft, 1, NPHI, Shout, &ncplx, NLAT, 1, fftw_plan_mode);
	if (fft_out == NULL) printf("fft_out failed\n");
printf("\ntranspose done\n");
	printf("** ifft + transpose :\n");	fftw_print_plan(ifft_tr);
	printf("\n** fft + transpose :\n"); fftw_print_plan(fft_tr);

	fftw_free(Shout);	fftw_free(ShF);
}

do_fft_tests()
{
	complex double *Sho;
	int jj;
	clock_t tcpu;

	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
		fftw_execute_dft_c2r(ifft_in, ShF, (double *) ShF);
	}
	tcpu = clock() - tcpu;
	printf("  ifft in-place : %d\n", (int) tcpu);

	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
		fftw_execute_dft_r2c(fft_in, (double *) ShF, ShF);
	}
	tcpu = clock() - tcpu;
	printf("  fft in-place : %d\n", (int) tcpu);

	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
		Sho = (complex double *) fftw_malloc( (NPHI/2+1) * NLAT * sizeof(complex double));
		fftw_execute_dft_c2r(ifft_out, Sho, (double *) ShF);
		fftw_free(Sho);
	}
	tcpu = clock() - tcpu;
	printf("  ifft out-of-place (+malloc) : %d\n", (int) tcpu);

	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
		Sho = (complex double *) fftw_malloc( (NPHI/2+1) * NLAT * sizeof(complex double));
		fftw_execute_dft_r2c(fft_out, (double *) ShF, Sho);
		fftw_free(Sho);
	}
	tcpu = clock() - tcpu;
	printf("  fft out-of-place (+malloc) : %d\n", (int) tcpu);

	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
		Sho = (complex double *) fftw_malloc( (NPHI/2+1) * NLAT * sizeof(complex double));
		fftw_execute_dft_c2r(ifft_tr, Sho, (double *) ShF);
		fftw_free(Sho);
	}
	tcpu = clock() - tcpu;
	printf("  ifft transpose (+malloc) : %d\n", (int) tcpu);

	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
		Sho = (complex double *) fftw_malloc( (NPHI/2+1) * NLAT * sizeof(complex double));
		fftw_execute_dft_r2c(fft_tr, (double *) ShF, Sho);
		fftw_free(Sho);
	}
	tcpu = clock() - tcpu;
	printf("  fft transpose (+malloc) : %d\n", (int) tcpu);

}
*/

void usage()
{
	printf("\nUsage: time_SHT lmax [options] \n");
	printf("        where lmax is the maxiumum spherical harmonic degree.\n");
	printf("** available options :\n");
	printf(" -mmax=<mmax> : defines the maximum spherical harmonic order <mmax>\n");
	printf(" -nphi=<nphi> : defines the number of azimutal (longitude) point\n");
	printf(" -nlat=<nlat> : defines the number of grid points in theta (latitude)\n");
	printf(" -mres=<mres> : the azimutal periodicity (1 for no symmetry; 2 for two-fold symmetry, ...)\n");
	printf(" -polaropt=<thr> : set the threshold for polar optimization. 0 for no polar optimization, 1.e-6 for agressive.\n");
	printf(" -iter=<n> : set the number of back-and-forth transforms to compute timings and errors.\n");
	printf(" -gauss : force gauss grid\n");
	printf(" -reg : force regular grid\n");
	printf(" -oop : force out-of-place transform\n");
	printf(" -transpose : force transpose data (ie phi varies fastest)\n");
	printf(" -nlorder : define non-linear order to be resolved.\n");
}

int main(int argc, char *argv[])
{
	complex double t1, t2;
	double t,tmax,n2;
	int i,im,m,l,jj, m_opt;
	clock_t tcpu;
	ticks tik0, tik1;
	double e0,e1;
	double polaropt = 1.e-8;		// default for polar optimization.
	enum shtns_type shtmode = sht_auto;		// default to "auto" (fastest) mode.
	enum shtns_norm shtnorm = sht_orthonormal;		// default to "orthonormal" SH.
	int layout = SHT_NATIVE_LAYOUT;
	int nlorder = 0;
	char name[20];
	FILE* fw;

	srand( time(NULL) );	// initialise les nombres.

	printf("time_SHT performs some spherical harmonic transforms, and displays timings and accuracy.\n");
	if (argc < 2) {
		usage();	exit(1);
	}

// if some wisdom has been saved, reload it.
	fw = fopen("fftw_wisdom","r");
	if (fw != NULL) {
		printf("> FFTW wisdom found.\n");
		fftw_import_wisdom_from_file(fw);
		fclose(fw);
	}

//	first argument is lmax, and is mandatory.
	sscanf(argv[1],"%lf",&t);	LMAX=t;
	MMAX=-1;	MRES=1;

	for (i=2; i<argc; i++) {		// parse command line
		sscanf(argv[i],"-%[^=]=%lf",name,&t);
		if (strcmp(name,"mmax") == 0) MMAX = t;
		if (strcmp(name,"mres") == 0) MRES = t;
		if (strcmp(name,"nlat") == 0) NLAT = t;
		if (strcmp(name,"nphi") == 0) NPHI = t;
		if (strcmp(name,"polaropt") == 0) polaropt = t;
		if (strcmp(name,"iter") == 0) SHT_ITER = t;
		if (strcmp(name,"gauss") == 0) shtmode = sht_gauss;		// force gauss grid.
		if (strcmp(name,"reg") == 0) shtmode = sht_reg_fast;	// force regular grid.
		if (strcmp(name,"quickinit") == 0) shtmode = sht_quick_init;	// Gauss grid and fast initialization time, but suboptimal fourier transforms.
		if (strcmp(name,"schmidt") == 0) shtnorm = sht_schmidt | SHT_NO_CS_PHASE;
		if (strcmp(name,"4pi") == 0) shtnorm = sht_fourpi | SHT_REAL_NORM;
		if (strcmp(name,"oop") == 0) layout = SHT_THETA_CONTIGUOUS;
		if (strcmp(name,"transpose") == 0) layout = SHT_PHI_CONTIGUOUS;
		if (strcmp(name,"nlorder") == 0) nlorder = t;
	}

	if (MMAX == -1) MMAX=LMAX/MRES;
	NLM = shtns_set_size(LMAX, MMAX, MRES, shtnorm);
	shtns_precompute_auto(shtmode | layout, polaropt, nlorder, &NLAT, &NPHI);
	
// now would be a good time to save fftw's wisdom.
	fw = fopen("fftw_wisdom","w");
	if (fw != NULL) {
		fftw_export_wisdom_to_file(fw);
		fclose(fw);
	}

	m_opt = Get_MTR_DCT();
/*
	t1 = 1.0+2.0*I;
	t2 = 1.0-I;
	printf("test : %f, %f, %f, %f\n",creal(t1),cimag(t1), creal(t2),cimag(t2));

	(double) t1 = 8.0 +I;
	(double) t2 = 8.1;
	printf("test : %f, %f, %f, %f\n",creal(t1),cimag(t1), creal(t2),cimag(t2));
*/
	write_vect("cost",ct,NLAT);
	write_vect("sint",st,NLAT);

	ShF = (complex double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Sh = (double *) ShF;
	ThF = (complex double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Th = (double *) ThF;
	NLF = (complex double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	NL = (double *) NLF;

	Tlm0 = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Slm0 = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Slm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Tlm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Qlm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);

// perform fft tests.
//	init_fft_tests();
//	do_fft_tests();
//	exit(0);

	Set_MTR_DCT(MMAX);

// SH_to_spat
	for (i=0;i<NLM;i++) {
		Slm[i] = 0.0;	Tlm[i] = 0.0;
	}
	for (im=0;im<NPHI;im++) {
		for (i=0;i<NLAT;i++) {
			Sh[im*NLAT+i] = 0.0;
		}
	}
	Slm[LiM(1,0)] = sh10_ct();
	if ((MMAX > 0)&&(MRES==1))
		Slm[LiM(1,1)] = sh11_st();
//	write_vect("ylm0",Slm, NLM*2);
	SH_to_spat(Slm,Sh);
	write_mx("spat",Sh,NPHI,NLAT);
	SHsphtor_to_spat(Slm,Tlm,Sh,Th);
	write_mx("spatt",Sh,NPHI,NLAT);
	write_mx("spatp",Th,NPHI,NLAT);

//	SHqst_to_lat(Slm,Slm,Tlm,ct[0],Sh,Th,Th,NPHI/2,LMAX,MMAX);
//	write_vect("spat_lat", Sh, NPHI/2);

/*	for (i=0;i<(NLAT/2)*NPHI;i++) {
		Sh[i] = 0.0;
	}
	SHeo_to_spat(Slm, Sh, 0);
	write_mx("spate",Sh,NPHI,NLAT/2);
	for (i=0;i<(NLAT/2)*NPHI;i++) {
		Th[i] = 0.0;
	}
	SHeo_to_spat(Slm, Th, 1);
	write_mx("spato",Th,NPHI,NLAT/2);

	for (i=0;i<NLM;i++)
		Slm[i] = 0.0;
	spat_to_SHeo(Sh, Slm, 0);
	write_vect("ylme",Slm, NLM*2);
	for (i=0;i<NLM;i++)
		Tlm[i] = 0.0;
	spat_to_SHeo(Th, Slm, 1);
	write_vect("ylmeo",Slm, NLM*2);
*/
// spat_to_SH
	for (im=0;im<NPHI;im++) {
		for (i=0;i<NLAT;i++) {
			Sh[im*NLAT+i] = ct[i];
		}
	}
	spat_to_SH(Sh,Slm);
	write_vect("ylm",(double *)Slm,NLM*2);

// test case...
	t = 1.0 / (RAND_MAX/2);
	for (i=0;i<NLM;i++) {
		Slm0[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
		Tlm0[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
	}

//	printf("** performing %d scalar SHT with NL evaluation\n", SHT_ITER);
	printf("** performing %d scalar SHT\n", SHT_ITER);
  if (m_opt >= 0) {
	printf(":: OPTIMAL\n");
	Set_MTR_DCT(m_opt);
	test_SHT();
	printf(":: FULL DCT\n");
	Set_MTR_DCT(MMAX);
	test_SHT();
  }
	printf(":: NO DCT\n");
	Set_MTR_DCT(-1);
	test_SHT();
	printf(":: ON THE FLY\n");
	test_SHT_fly();

  if (m_opt >= 0) {
	printf(":: OPTIMAL with LTR\n");
	Set_MTR_DCT(m_opt);
	test_SHT_l(LMAX/2);
	printf(":: FULL DCT with LTR\n");
	Set_MTR_DCT(MMAX);
	test_SHT_l(LMAX/2);
  }
	printf(":: NO DCT with LTR\n");
	Set_MTR_DCT(-1);
	test_SHT_l(LMAX/2);

#define TEST_VECT_SHT
#ifdef TEST_VECT_SHT
	Slm0[LM(0,0)] = 0.0;	// l=0, m=0 n'a pas de signification sph/tor
	Tlm0[LM(0,0)] = 0.0;	// l=0, m=0 n'a pas de signification sph/tor
//	for (i=0;i<NLM;i++) Slm0[i] = 0.0;	// zero out Slm.

	printf("** performing %d vector SHT\n", SHT_ITER);
  if (m_opt >= 0) {
	printf(":: OPTIMAL\n");
	Set_MTR_DCT(m_opt);
	test_SHT_vect();
	printf(":: FULL DCT\n");
	Set_MTR_DCT(MMAX);
	test_SHT_vect();
  }
	printf(":: NO DCT\n");
	Set_MTR_DCT(-1);
	test_SHT_vect();
	printf(":: ON THE FLY\n");
	test_SHT_vect_fly();

  if (m_opt >= 0) {
	printf(":: OPTIMAL with LTR\n");
	Set_MTR_DCT(m_opt);
	test_SHT_vect_l(LMAX/2);
	printf(":: FULL DCT with LTR\n");
	Set_MTR_DCT(MMAX);
	test_SHT_vect_l(LMAX/2);
  }
	printf(":: NO DCT with LTR\n");
	Set_MTR_DCT(-1);
	test_SHT_vect_l(LMAX/2);
	
	printf("** performing %d 3D vector SHT (no DCT) \n", SHT_ITER);
  if (m_opt >= 0) {
	printf(":: OPTIMAL with LTR\n");
	Set_MTR_DCT(m_opt);
	test_SHT_vect3d();
	printf(":: FULL DCT with LTR\n");
	Set_MTR_DCT(MMAX);
	test_SHT_vect3d();
  }
	printf(":: NO DCT\n");
	Set_MTR_DCT(-1);
	test_SHT_vect3d();
	printf(":: ON THE FLY\n");
	test_SHT_vect3d_fly();

#endif

}

