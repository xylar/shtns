// base flow. Cartesien spectral.

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
// cycle counter from FFTW
#include "cycle.h"

#include "SHT.h"

complex double *Slm, *Slm0, *Tlm, *Tlm0;	// spherical harmonics l,m space
complex double *ShF, *ThF, *NLF;	// Fourier space : theta,m
double *Sh, *Th, *NL;		// real space : theta,phi (alias of ShF)

// polar optimization threshold
#define POLAR_OPT_THR 1e-6
//#define POLAR_OPT_THR 0
// number of SH iterations
#define SHT_ITER 50

	
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

int test_SHT()
{
	int jj,i;
	double tmax,t,n2;
	clock_t tcpu;

	for (i=0;i<NLM;i++) Slm[i] = Slm0[i];	// restore test case...
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
	printf("   2iSHT + NL + SHT x%d time : %d\n", SHT_ITER, (int )tcpu);

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
	printf("   => max error = %g (l=%.0f,lm=%d)   rms error = %g\n",tmax,el[jj],jj,sqrt(n2/NLM));
//		write_vect("Qlm.reg",Slm,NLM*2);
	return (int) tcpu;
}

int test_SHT_l(int ltr)
{
	int jj,i;
	double tmax,t,n2;
	clock_t tcpu;

	for (i=0;i<NLM;i++) Slm[i] = Slm0[i];	// restore test case...
	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
// synthese (inverse legendre)
		SH_to_spat_l(Slm,ShF,ltr);
		SH_to_spat_l(Tlm,ThF,ltr);
		for (i=0; i< NLAT*NPHI; i++) {
			ThF[i] *= ShF[i];
		}
// analyse (direct legendre)
		spat_to_SH_l(ShF,Slm,ltr);
	}
	tcpu = clock() - tcpu;
	printf("   2iSHT + NL + SHT x%d with L-truncation at %d. time : %d\n", SHT_ITER, ltr, (int )tcpu);

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
	printf("   => max error = %g (l=%.0f,lm=%d)   rms error = %g\n",tmax,el[jj],jj,sqrt(n2/NLM));
//		write_vect("Qlm.reg",Slm,NLM*2);
	return (int) tcpu;
}



int test_SHT_vect()
{
	int jj,i;
	double tmax,t,n2;
	clock_t tcpu;

	for (i=0;i<NLM;i++) {
		Slm[i] = Slm0[i];	Tlm[i] = Tlm0[i];
	}
	tcpu = clock();
	for (jj=0; jj< SHT_ITER; jj++) {
	#ifndef _SHT_EO_
		SHsphtor_to_spat(Slm,Tlm,ShF,ThF);
		spat_to_SHsphtor(ShF,ThF,Slm,Tlm);
	#else
		SHsphtor_to_spat(Slm,Slm,ShF,ThF);
		spat_to_SHsphtor(ShF,ThF,Slm,Slm);
	#endif
	}
	tcpu = clock() - tcpu;
	printf("   iSHT + SHT x%d time : %d\n", SHT_ITER, (int) tcpu);

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
	printf("   Spheroidal => max error = %g (l=%.0f,lm=%d)    rms error = %g\n",tmax,el[jj],jj,sqrt(n2/NLM));
//	write_vect("Slm",Slm,NLM*2);

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
	printf("   Toroidal => max error = %g (l=%.0f,lm=%d)    rms error = %g\n",tmax,el[jj],jj,sqrt(n2/NLM));
//	write_vect("Tlm",Tlm,NLM*2);
	return (int) tcpu;
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
#ifndef SHT_AXISYM
	//                          ... lmax,mmax,mres, nlat, nphi );
	init_SH( sht_auto, POLAR_OPT_THR, 340, 3, 5,  512, 12 );
#else
	init_SH( sht_auto, POLAR_OPT_THR, 681, 0, 1,  1024, 1 );
#endif
	m_opt = Get_MTR_DCT();

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

	Set_MTR_DCT(MMAX);

// SH_to_spat
	for (i=0;i<NLM;i++) {
		Slm[i] = 0.0;	Tlm[i] = 0.0;
	}
//	Slm[LiM(1,1)] = 1;
	Slm[LiM(1,0)] = Y10_ct;
	SH_to_spat(Slm,ShF);
	write_mx("spat",Sh,NPHI,NLAT);
	SHsphtor_to_spat(Slm,Tlm,ShF,ThF);
	write_mx("spatt",Sh,NPHI,NLAT);
	write_mx("spatp",Th,NPHI,NLAT);

// spat_to_SH
	for (im=0;im<NPHI;im++) {
		for (i=0;i<NLAT;i++) {
			Sh[im*NLAT+i] = ct[i];
		}
	}
	spat_to_SH(ShF,Slm);
	write_vect("ylm",Slm,NLM*2);


// test case...
	t = 1.0 / (RAND_MAX/2);
	for (i=0;i<NLM;i++) {
		Slm0[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
		Tlm0[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
	}

	printf("** performing %d scalar SHT with NL evaluation\n", SHT_ITER);
	printf(":: OPTIMAL\n");
	Set_MTR_DCT(m_opt);
	test_SHT();
	printf(":: FULL DCT\n");
	Set_MTR_DCT(MMAX);
	test_SHT();
	printf(":: NO DCT\n");
	Set_MTR_DCT(-1);
	test_SHT();

	printf(":: OPTIMAL with LTR\n");
	Set_MTR_DCT(m_opt);
	test_SHT_l(LMAX/2);
	printf(":: FULL DCT with LTR\n");
	Set_MTR_DCT(MMAX);
	test_SHT_l(LMAX/2);
	printf(":: NO DCT with LTR\n");
	Set_MTR_DCT(-1);
	test_SHT_l(LMAX/2);

#define TEST_VECT_SHT
#ifdef TEST_VECT_SHT
	Slm0[LM(0,0)] = 0.0;	// l=0, m=0 n'a pas de signification sph/tor
	Tlm0[LM(0,0)] = 0.0;	// l=0, m=0 n'a pas de signification sph/tor
//	for (i=0;i<NLM;i++) Slm0[i] = 0.0;	// zero out Slm.

	printf("** performing %d vector SHT\n", SHT_ITER);
	printf(":: OPTIMAL\n");
	Set_MTR_DCT(m_opt);
	test_SHT_vect();
	printf(":: FULL DCT\n");
	Set_MTR_DCT(MMAX);
	test_SHT_vect();
	printf(":: NO DCT\n");
	Set_MTR_DCT(-1);
	test_SHT_vect();
#endif
}

