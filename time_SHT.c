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
// cycle counter from FFTW
#include "cycle.h"


complex double *Slm, *Slm0, *Tlm, *Tlm0;	// spherical harmonics l,m space
complex double *ShF, *ThF, *NLF;	// Fourier space : theta,m
double *Sh, *Th, *NL;		// real space : theta,phi (alias of ShF)

// parameters for SHT.c
#define NLAT_2 304
#define LMAX 400
#define NPHI 32
#define MMAX 10
#define MRES 1

#define _SH_DEBUG_
//#define _SHT_EO_
//#define SHT_EQUAL	/* SHT on equal spaced grid + polar points. */
#define SHT_DCT
//#define SHT_DCT_IT_OUT
#include "SHT.c"

// polar optimization threshold
#define POLAR_OPT_THR 1e-6
//#define POLAR_OPT_THR 0
// number of SH iterations
#define SHT_ITER 1000

	
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
	double t,tmax,n2;
	int i,im,m,l,jj, m_opt;
	clock_t tcpu;
	ticks tik0, tik1;
	double e0,e1;

	srand( time(NULL) );	// initialise les nombres.
	init_SH( POLAR_OPT_THR );
	m_opt = MTR_DCT;
	
	if (MMAX > 0) {
		write_mx("yl1",ylm[1],NLAT_2,LMAX);
//		write_mx("dyl1",&dylm[1][0].t,NLAT_2,2*LMAX);
	}

	t1 = 1.0+2.0*I;
	t2 = 1.0-I;
	printf("test : %f, %f, %f, %f\n",creal(t1),cimag(t1), creal(t2),cimag(t2));

	(double) t1 = 8.0 +I;
	(double) t2 = 8.1;
	printf("test : %f, %f, %f, %f\n",creal(t1),cimag(t1), creal(t2),cimag(t2));

	write_vect("cost",ct,NLAT);
	write_vect("sint",st,NLAT);

	ShF = (complex double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Sh = (double *) ShF;
	ThF = (complex double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Th = (double *) ThF;
	NLF = (complex double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	NL = (double *) NLF;

	Tlm0 = (complex double *) malloc(sizeof(complex double)* NLM);
	Slm0 = (complex double *) malloc(sizeof(complex double)* NLM);
	Slm = (complex double *) malloc(sizeof(complex double)* NLM);
	Tlm = (complex double *) malloc(sizeof(complex double)* NLM);

	printf("test\n");
	Set_MTR_DCT(MMAX);
	printf("test\n");
// SH_to_spat
	for (i=0;i<NLM;i++) {
		Slm[i] = 0.0;	Tlm[i] = 0.0;
	}
//	Slm[LiM(1,1)] = 1;
	Slm[LiM(1,0)] = Y10_ct;
	SH_to_spat(Slm,ShF);
	write_mx("spat",Sh,NPHI,NLAT);
#ifdef SHT_DCT
	SH_to_spat_dct(Slm,ShF);
	write_mx("spat_dct",Sh,NPHI,NLAT);
#endif
	SHsphtor_to_spat(Slm,Tlm,ShF,ThF);
	write_mx("spatt",Sh,NPHI,NLAT);
	write_mx("spatp",Th,NPHI,NLAT);
#ifdef SHT_DCT
	SHsphtor_to_spat_dct(Slm,Tlm,ShF,ThF);
	write_mx("spatt_dct",Sh,NPHI,NLAT);
	write_mx("spatp_dct",Th,NPHI,NLAT);
#endif

// spat_to_SH
	for (im=0;im<NPHI;im++) {
		for (i=0;i<NLAT;i++) {
			Sh[im*NLAT+i] = ct[i];
		}
	}
#ifdef SHT_DCT
//	spat_to_SH_dct(ShF,Slm);
//	write_vect("ylm_dct",Slm,NLM*2);		// should be sqrt(4pi/3) if l=1, zero if l!=1.
#else
	spat_to_SH(ShF,Slm);
	write_vect("ylm",Slm,NLM*2);
#endif
//	return(1);

// test case...
	t = 1.0 / (RAND_MAX/2);
	for (i=0;i<NLM;i++) {
		Slm0[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
		Tlm0[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
		Slm[i] = Slm0[i];
	}

/*
/// TIME WITH DIFFERENT MTR_DCT : dichotomy
	double get_time(int m) {
		clock_t tcpu;
		ticks tik0, tik1;
		Set_MTR_DCT(m);
		tcpu = clock();
		tik0 = getticks();
			for (i=0; i<1000; i++) {
				SH_to_spat(Slm,ShF);
//				spat_to_SH(ShF,Slm);
				SHsphtor_to_spat(Slm,Tlm,ShF,ThF);
//				spat_to_SHsphtor(ShF,ThF,Slm,Tlm);
			}
		tik1 = getticks();
		tcpu = clock() - tcpu;
		printf("m=%d - ticks : %.0f, tcpu : %d\n", m, elapsed(tik1,tik0), (int) tcpu);
		return elapsed(tik1,tik0);
	}

	m = -1;
	tmax = get_time(m);	jj = m;
	for (m=-1; m<=MMAX; m++) {
		t = get_time(m);
		if (t < tmax) {
			tmax = t;
			jj = m;
		}
	}
	get_time(-1);
	printf("optimal MTR_DCT = %d\n", jj);
	Set_MTR_DCT(jj);
*/
	
//	Slm[LM(0,0)] = 1.0;
//	Slm[LM(1,1)] = 3.0;
//	Slm[LM(4,2)] = 2.0+I;
// 	Slm[LM(10,4)] = -4.0-I;
// 	Slm[LM(55,12)] = 5.0-2.0*I;


// test case...
	t = 1.0 / (RAND_MAX/2);
	for (i=0;i<NLM;i++) {
		Slm0[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
		Tlm0[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
		Slm[i] = Slm0[i];
	}

/// REGULAR
	printf("[REG] :: performing %d scalar SHT with NL evaluation\n", SHT_ITER);
	Set_MTR_DCT(-1);
	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
// synthese (inverse legendre)
		SH_to_spat(Slm,ShF);
		SH_to_spat(Tlm,ThF);
		for (i=0; i< NLAT*NPHI; i++) {
			ThF[i] *= ShF[i];
		}
// analyse (direct legendre)
		spat_to_SH(ShF,Slm);
	}
	tcpu = clock() - tcpu;
	printf("  2iSHT + NL + SHT x%d time : %d\n", SHT_ITER, (int )tcpu);

// compute error :
	tmax = 0;	n2 = 0;		jj=0;
	for (i=0;i<NLM;i++) {
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
	printf("  => max error = %g (l=%.0f,lm=%d)   rms error = %g\n",tmax,el[jj],jj,sqrt(n2/NLM));
	write_vect("Qlm.reg",Slm,NLM*2);

// restore test case...
	for (i=0;i<NLM;i++) Slm[i] = Slm0[i];
/// DCT
#ifdef SHT_DCT
	printf("[DCT] :: performing %d scalar SHT with NL evaluation\n", SHT_ITER);
#ifdef SHT_DCT_IT_OUT
	printf("         (IT_OUT variant : not working yet...)\n");
#endif
	Set_MTR_DCT(MMAX);
	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
// synthese (inverse legendre)
		SH_to_spat_dct(Slm,ShF);
		SH_to_spat_dct(Tlm,ThF);
		for (i=0; i< NLAT*NPHI; i++) {
			ThF[i] *= ShF[i];
		}
// analyse (direct legendre)
		spat_to_SH(ShF,Slm);
	}
	tcpu = clock() - tcpu;
	printf("  2iSHT + NL + SHT x%d time : %d\n", SHT_ITER, (int )tcpu);

// compute error :
	tmax = 0;	n2 = 0;		jj=0;
	for (i=0;i<NLM;i++) {
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
	printf("  => max error = %g (l=%.0f,lm=%d)   rms error = %g\n",tmax,el[jj],jj,sqrt(n2/NLM));
	write_vect("Qlm.dct",Slm,NLM*2);

// restore test case...
	for (i=0;i<NLM;i++) Slm[i] = Slm0[i];
/// OPT
	printf("[OPT] :: performing %d scalar SHT with NL evaluation\n", SHT_ITER);
	Set_MTR_DCT(m_opt);
	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
// synthese (inverse legendre)
		SH_to_spat(Slm,ShF);
		SH_to_spat(Tlm,ThF);
		for (i=0; i< NLAT*NPHI; i++) {
			ThF[i] *= ShF[i];
		}
// analyse (direct legendre)
		spat_to_SH(ShF,Slm);
	}
	tcpu = clock() - tcpu;
	printf("  2iSHT + NL + SHT x%d time : %d\n", SHT_ITER, (int )tcpu);

// compute error :
	tmax = 0;	n2 = 0;		jj=0;
	for (i=0;i<NLM;i++) {
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
	printf("  => max error = %g (l=%.0f,lm=%d)   rms error = %g\n",tmax,el[jj],jj,sqrt(n2/NLM));
	write_vect("Qlm.dct",Slm,NLM*2);
#endif

#define TEST_VECT_SHT
#ifdef TEST_VECT_SHT
	printf("\n[REG] :: performing %d vector SHT\n", SHT_ITER);
	Set_MTR_DCT(-1);
	Slm0[LM(0,0)] = 0.0;	// l=0, m=0 n'a pas de signification sph/tor
	Tlm0[LM(0,0)] = 0.0;	// l=0, m=0 n'a pas de signification sph/tor
	for (i=0;i<NLM;i++) {
		Tlm0[i] = 0.0;	// zero out Slm.
	}
	for (i=0;i<NLM;i++) {
		Slm[i] = Slm0[i];	Tlm[i] = Tlm0[i];
	}
	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
#ifndef _SHT_EO_
// synthese (inverse legendre)
		SHsphtor_to_spat(Slm,Tlm,ShF,ThF);
// analyse (direct legendre)
		spat_to_SHsphtor(ShF,ThF,Slm,Tlm);
#else
		SHsphtor_to_spat(Slm,Slm,ShF,ThF);
		spat_to_SHsphtor(ShF,ThF,Slm,Slm);
#endif
	}
	tcpu = clock() - tcpu;
	printf("iSHT + SHT x%d time : %d\n", SHT_ITER, (int) tcpu);

//	write_vect("dylm0", dylm_dct[0]

// compute error :
	tmax = 0;	n2 = 0;		jj=0;
	for (i=0;i<NLM;i++) {
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
	printf("  Spheroidal => max error = %g (l=%.0f,lm=%d)    rms error = %g\n",tmax,el[jj],jj,sqrt(n2/NLM));
	write_vect("Slm",Slm,NLM*2);

// compute error :
	tmax = 0;	n2 = 0;		jj=0;
	for (i=0;i<NLM;i++) {
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
	printf("  Toroidal => max error = %g (l=%.0f,lm=%d)    rms error = %g\n",tmax,el[jj],jj,sqrt(n2/NLM));
	write_vect("Tlm",Tlm,NLM*2);
#ifdef SHT_DCT
	printf("[DCT] :: performing %d vector SHT\n", SHT_ITER);
	Set_MTR_DCT(MMAX);
	Slm0[LM(0,0)] = 0.0;	// l=0, m=0 n'a pas de signification sph/tor
	Tlm0[LM(0,0)] = 0.0;	// l=0, m=0 n'a pas de signification sph/tor
	for (i=0;i<NLM;i++) {
		Slm[i] = Slm0[i];	Tlm[i] = Tlm0[i];
	}
	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
#ifndef _SHT_EO_
// synthese (inverse legendre)
		SHsphtor_to_spat_dct(Slm,Tlm,ShF,ThF);
// analyse (direct legendre)
		spat_to_SHsphtor(ShF,ThF,Slm,Tlm);
#else
		SHsphtor_to_spat_dct(Slm,Slm,ShF,ThF);
		spat_to_SHsphtor(ShF,ThF,Slm,Slm);
#endif
	}
	tcpu = clock() - tcpu;
	printf("iSHT + SHT x%d time : %d\n", SHT_ITER, (int) tcpu);

//	write_vect("dylm0", dylm_dct[0]

// compute error :
	tmax = 0;	n2 = 0;		jj=0;
	for (i=0;i<NLM;i++) {
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
	printf("  Spheroidal => max error = %g (l=%.0f,lm=%d)    rms error = %g\n",tmax,el[jj],jj,sqrt(n2/NLM));
	write_vect("Slm",Slm,NLM*2);

// compute error :
	tmax = 0;	n2 = 0;		jj=0;
	for (i=0;i<NLM;i++) {
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
	printf("  Toroidal => max error = %g (l=%.0f,lm=%d)    rms error = %g\n",tmax,el[jj],jj,sqrt(n2/NLM));
	write_vect("Tlm",Tlm,NLM*2);
#endif

	printf("[OPT] :: performing %d vector SHT\n", SHT_ITER);
	Set_MTR_DCT(m_opt);
	Slm0[LM(0,0)] = 0.0;	// l=0, m=0 n'a pas de signification sph/tor
	Tlm0[LM(0,0)] = 0.0;	// l=0, m=0 n'a pas de signification sph/tor
	for (i=0;i<NLM;i++) {
		Slm[i] = Slm0[i];	Tlm[i] = Tlm0[i];
	}
	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
#ifndef _SHT_EO_
// synthese (inverse legendre)
		SHsphtor_to_spat(Slm,Tlm,ShF,ThF);
// analyse (direct legendre)
		spat_to_SHsphtor(ShF,ThF,Slm,Tlm);
#else
		SHsphtor_to_spat(Slm,Slm,ShF,ThF);
		spat_to_SHsphtor(ShF,ThF,Slm,Slm);
#endif
	}
	tcpu = clock() - tcpu;
	printf("iSHT + SHT x%d time : %d\n", SHT_ITER, (int) tcpu);

//	write_vect("dylm0", dylm_dct[0]

// compute error :
	tmax = 0;	n2 = 0;		jj=0;
	for (i=0;i<NLM;i++) {
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
	printf("  Spheroidal => max error = %g (l=%.0f,lm=%d)    rms error = %g\n",tmax,el[jj],jj,sqrt(n2/NLM));
	write_vect("Slm",Slm,NLM*2);

// compute error :
	tmax = 0;	n2 = 0;		jj=0;
	for (i=0;i<NLM;i++) {
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
	printf("  Toroidal => max error = %g (l=%.0f,lm=%d)    rms error = %g\n",tmax,el[jj],jj,sqrt(n2/NLM));
	write_vect("Tlm",Tlm,NLM*2);

#endif
}

