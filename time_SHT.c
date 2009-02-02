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
#define NLAT_2 256
#define LMAX 340
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

void spat_to_SHm(long int im, complex double *BrF, complex double *Qlm)
{
	complex double reo[2*NLAT_2];	// symmetric (even) and anti-symmetric (odd) parts, interleaved.
	complex double *Ql;		// virtual pointers for given im
	double *zl;
	long int i,m,l;

// defines how to access even and odd parts of data
	#define re	reo[2*i]
	#define ro	reo[2*i+1]

	BrF += im*NLAT;
	if (im == 0) {		// dzl.p = 0.0 : and evrything is REAL
		m=im*MRES;
 		for (i=0;i<NLAT/2;i++) {	// compute symmetric and antisymmetric parts.
			(double) reo[2*i]   = (double) BrF[i] + (double) BrF[NLAT-(i+1)];
			(double) reo[2*i+1] = (double) BrF[i] - (double) BrF[NLAT-(i+1)];
 		}
 		if (i < NLAT_2) {		// NLAT is odd : special equator handling
			(double) reo[2*i] = (double) BrF[i];	(double) reo[2*i+1] = 0.0;
 		}
		l=m;
		Ql = &Qlm[LiM(0,im)];		// virtual pointer for l=0 and im
		zl = zlm[im];
		while (l<LTR) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
			Ql[l] = 0.0;
			Ql[l+1] = 0.0;
			for (i=0; i < NLAT_2; i++) {
				(double) Ql[l] += (double) re * zl[2*i];		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
				(double) Ql[l+1] += (double) ro * zl[2*i+1];	// Qlm[LiM(l+1,im)] += zlm[im][(l+1-m)*NLAT/2 + i] * fm[i];
			}
			l+=2;
			zl += 2*NLAT_2;
		}
		if (l==LMAX) {
			Ql[l] = 0.0;
			for (i=0;i<NLAT_2;i++) {
				(double) Ql[l] += zl[i] * (double) re;		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
			}
		} else if (l==LTR) {
			Ql[l] = 0.0;
			for (i=0; i < NLAT_2; i++) {
				(double) Ql[l] += (double) re * zl[2*i];		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
			}
		}
	} else {
		m=im*MRES;
 		for (i=tm[im];i<NLAT/2;i++) {	// compute symmetric and antisymmetric parts.
			reo[2*i]   = BrF[i] + BrF[NLAT-(i+1)];
			reo[2*i+1] = BrF[i] - BrF[NLAT-(i+1)];
 		}
 		if (i<NLAT_2) {		// NLAT is odd : special equator handling
			reo[2*i] = BrF[i];		reo[2*i+1] = 0.0;
 		}
		l=m;
		Ql = &Qlm[LiM(0,im)];		// virtual pointer for l=0 and im
		zl = zlm[im];
		while (l<LTR) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
			Ql[l] = 0.0;
			Ql[l+1] = 0.0;
			for (i=tm[im]; i < NLAT_2; i++) {	// tm[im] : polar optimization
				Ql[l]   += re * zl[2*i];		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
				Ql[l+1] += ro * zl[2*i+1];	// Qlm[LiM(l+1,im)] += zlm[im][(l+1-m)*NLAT/2 + i] * fm[i];
			}
			l+=2;
			zl += 2*NLAT_2;
		}
		if (l==LMAX) {
			Ql[l] = 0.0;	// Qlm[LiM(l,im)] = 0.0;
			for (i=tm[im];i<NLAT_2;i++) {	// polar optimization
				Ql[l] += zl[i] * re;	// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
			}
		} else if (l==LTR) {
			Ql[l] = 0.0;
			for (i=tm[im]; i < NLAT_2; i++) {	// tm[im] : polar optimization
				Ql[l]   += re * zl[2*i];		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
			}
		}
	}
	#undef re
	#undef ro
}

void spat_to_SHm_dct(long int im, complex double *ShF, complex double *Slm)
{
	complex double *Sl;		// virtual pointers for given im
	double *zl;
	long int k,m,l;

  #if NPHI > 1
	if (MRES & 1) {		// odd m's are present
		if (im & 1) {
			for (k=0; k<NLAT; k++)	ShF[im*NLAT + k] *= st[k];	// corection beforce DCT
		}
	}
  #endif
	ShF += im*NLAT;
	fftw_execute_r2r(dctm,(double *) ShF, (double *) ShF);		// DCT

	if (im == 0) {
		m=0;
		Sl = &Slm[LiM(0,im)];		// virtual pointer for l=0 and im
		zl = zlm_dct[im];
		for (l=m; l<LMAX; l+=2) {		// l has parity of m
			Sl[l] = 0.0;	Sl[l+1] = 0.0;
			for (k=l; k<=NLAT; k+=2) {		// for m=0, zl coeff with k<l are zeros.
				(double) Sl[l]   += (double) ShF[k]   * zl[k];
				(double) Sl[l+1] += (double) ShF[k+1] * zl[k+1];
			}
			zl += NLAT;
		}
		if ((LMAX & 1) == 0) {	// if (l == LMAX)  <=>  if ((LMAX & 1) == 0) for m=0
			Sl[l] = 0.0;
			for (k=l; k<=NLAT; k+=2) {		// for m=0, DCT coeff with k<l are zeros.
				(double) Sl[l]   += (double) ShF[k]   * zl[k];
			}
		}
	} else {
		m=im*MRES;
		Sl = &Slm[LiM(0,im)];		// virtual pointer for l=0 and im
		zl = zlm_dct[im];
		for (l=m; l<LMAX; l+=2) {		// l has parity of m
			Sl[l] = 0.0;	Sl[l+1] = 0.0;
			for (k=0; k<=NLAT; k+=2) {
				Sl[l]   += ShF[k]   * zl[k];
				Sl[l+1] += ShF[k+1] * zl[k+1];
			}
			zl += NLAT;
		}
		if (l == LMAX) {
			Sl[l] = 0.0;
			for (k=0; k<=NLAT; k+=2) {
				Sl[l]   += ShF[k]   * zl[k];
			}
		}
	}
}


void SHm_to_spat(long int im, complex double *Qlm, complex double *BrF)
{
	complex double fe, fo;		// even and odd parts
	complex double *Ql;
	double *yl;
	long int i,m,l;

	m = im*MRES;
	BrF += im*NLAT;
	if (im == 0) {
		Ql = &Qlm[LiM(0,im)];	// virtual pointer for l=0 and im
		i=0;
		yl  = ylm[im] + i*(LMAX-m+1) -m;
		while (i < NLAT_2) {	// ops : NLAT_2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
			fe = 0.0;	fo = 0.0;
			while (l<LTR) {	// compute even and odd parts
				(double) fe += yl[l] * (double) Ql[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
				(double) fo += yl[l+1] * (double) Ql[l+1];	// fo += ylm[im][i*(LMAX-m+1) + (l+1-m)] * Qlm[LiM(l+1,im)];
				l+=2;
			}
			if (l==LTR) {
				(double) fe += yl[l] * (double) Ql[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
			}
			BrF[i] = fe + fo;
			i++;
			BrF[NLAT-i] = fe - fo;
			yl  += (LMAX-m+1);
		}
	} else {
		Ql = &Qlm[LiM(0,im)];	// virtual pointer for l=0 and im
		i=0;
		while (i<tm[im]) {	// polar optimization
			BrF[i] = 0.0;
			BrF[NLAT-tm[im] + i] = 0.0;	// south pole zeroes <=> BrF[im*NLAT + NLAT-(i+1)] = 0.0;
			i++;
		}
		yl  = ylm[im] + i*(LMAX-m+1) -m;
		while (i < NLAT_2) {	// ops : NLAT_2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
			fe = 0.0;	fo = 0.0;
			while (l<LTR) {	// compute even and odd parts
				fe  += yl[l] * Ql[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
				fo  += yl[l+1] * Ql[l+1];	// fo += ylm[im][i*(LMAX-m+1) + (l+1-m)] * Qlm[LiM(l+1,im)];
				l+=2;
			}
			if (l==LTR) {
				fe  += yl[l] * Ql[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
			}
			BrF[i] = fe + fo;
			i++;
			BrF[NLAT-i] = fe - fo;
			yl  += (LMAX-m+1);
		}
	}
}

void SHm_to_spat_dct(long int im, complex double *Qlm, complex double *BrF)
{
	complex double *Sl;
	double *yl;
	long int it,m,l;

	BrF += NLAT*im;
	if (im == 0) {
		m=0;
		Sl = &Qlm[LiM(0,im)];		// virtual pointer for l=0 and im
		yl = ylm_dct[im];
		for (it=0; it<NLAT; it++)
			BrF[it] = 0.0;		// zero out array (includes DCT padding)
		for (l=m; l<LMAX; l+=2) {
			for (it=0; it<=l; it+=2) {
				(double) BrF[it]   += yl[it] *   (double) Sl[l];
				(double) BrF[it+1] += yl[it+1] * (double) Sl[l+1];
			}
			yl += (l+2 - (m&1));
		}
		if (l==LMAX) {
			for (it=0; it<=l; it+=2)
				(double) BrF[it] += Sl[l] * (double) yl[it];
		}
	} else {
		m=im*MRES;
		Sl = &Qlm[LiM(0,im)];		// virtual pointer for l=0 and im

#ifdef SHT_DCT_IT_OUT
		yl = ykm_dct[im] -m;
		it = 0;
		while (it < m) {
			BrF[it] = 0.0;	BrF[it+1] = 0.0;
			for (l=m; l<LMAX; l+=2) {
				BrF[it]   += Sl[l]   * yl[l];
				BrF[it+1] += Sl[l+1] * yl[l+1];
			}
			it+=2;
			yl += (LMAX-m+1);
		}
		while (it < LMAX) {
			BrF[it] = 0.0;	BrF[it+1] = 0.0;
			for (l=it; l<LMAX; l+=2) {
				BrF[it]   += Sl[l]   * yl[l];
				BrF[it+1] += Sl[l+1] * yl[l+1];
			}
			it+=2;
			yl += (LMAX-m+1);
		}
		while (it < NLAT) {
			BrF[it] = 0.0;	BrF[it+1] = 0.0;
			it+=2;
		}
#else
		yl = ylm_dct[im];
		for (it=0; it<NLAT; it++)
			BrF[it] = 0.0;		// zero out array (includes DCT padding)
		for (l=m; l<LMAX; l+=2) {
			for (it=0; it<=l; it+=2) {
				BrF[it]   += Sl[l]   * yl[it];
				BrF[it+1] += Sl[l+1] * yl[it+1];
			}
			yl += (l+2 - (m&1));
		}
		if (l==LMAX) {
			for (it=0; it<=l; it+=2)
				BrF[it] += Sl[l] * yl[it];
		}
#endif
	}
	fftw_execute_r2r(idctm,(double *) BrF, (double *) BrF);		// iDCT
  #if NPHI>1
	if (MRES & 1) {		// odd m's must be multiplied by sin(theta) which was removed from ylm's
		if (im & 1) {
			for (it=0; it<NLAT; it++)
				BrF[it] *= st[it];
		}
	}
  #endif
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

int main()
{
	complex double t1, t2;
	double t,tmax,n2;
	int i,im,m,l,jj;
	clock_t tcpu;
	ticks tik0, tik1;

	srand( time(NULL) );	// initialise les nombres.
	init_SH( POLAR_OPT_THR );
	
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

// SH_to_spat
	for (i=0;i<NLM;i++) {
		Slm[i] = 0.0;	Tlm[i] = 0.0;
	}
//	Slm[LiM(1,1)] = 1;
	Tlm[LiM(4,1)] = 1.0;
	SH_to_spat(Slm,ShF);
	write_mx("spat",Sh,NPHI,NLAT);
	SH_to_spat_dct(Slm,ShF);
	write_mx("spat_dct",Sh,NPHI,NLAT);

	SHsphtor_to_spat(Slm,Tlm,ShF,ThF);
	write_mx("spatt",Sh,NPHI,NLAT);
	write_mx("spatp",Th,NPHI,NLAT);
	SHsphtor_to_spat_dct(Slm,Tlm,ShF,ThF);
	write_mx("spatt_dct",Sh,NPHI,NLAT);
	write_mx("spatp_dct",Th,NPHI,NLAT);

// spat_to_SH
/*	for (im=0;im<NPHI;im++) {
		for (i=0;i<NLAT;i++) {
			Sh[im*NLAT+i] = ct[i];
		}
	}*/
	spat_to_SH_dct(ShF,Slm);
	write_vect("ylm",Slm,NLM*2);		// should be sqrt(4pi/3) if l=1, zero if l!=1.
//	return(1);

// test case...
	t = 1.0 / (RAND_MAX/2);
	for (i=0;i<NLM;i++) {
		Slm0[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
		Tlm0[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
		Slm[i] = Slm0[i];
	}
	
//	Slm[LM(0,0)] = 1.0;
//	Slm[LM(1,1)] = 3.0;
//	Slm[LM(4,2)] = 2.0+I;
// 	Slm[LM(10,4)] = -4.0-I;
// 	Slm[LM(55,12)] = 5.0-2.0*I;

#ifdef SHT_DCT

/// TEST m by m and compare DCT to normal.
	for (im=0; im <= MMAX; im++) {
		m = im*MRES;
		tcpu = clock();
		tik0 = getticks();
		for (jj=0; jj<SHT_ITER; jj++)	SHm_to_spat(im, Slm, ShF);
		tik1 = getticks();
		tcpu = clock() - tcpu;
		printf("m=%d - iSHT x%d time : %d, ticks : %.0f\n", m, SHT_ITER, ((int )tcpu)/SHT_ITER, elapsed(tik1,tik0)/SHT_ITER);

		tcpu = clock();
		tik0 = getticks();
		for (jj=0; jj<SHT_ITER; jj++)	SHm_to_spat_dct(im, Slm, ShF);
		tik1 = getticks();
		tcpu = clock() - tcpu;
		printf("m=%d - iSHT_dct x%d time : %d, ticks : %.0f\n", m, SHT_ITER, ((int )tcpu)/SHT_ITER, elapsed(tik1,tik0)/SHT_ITER);

		tcpu = clock();
		tik0 = getticks();
		for (jj=0; jj<SHT_ITER; jj++)	spat_to_SHm(im, ShF, Slm);
		tik1 = getticks();
		tcpu = clock() - tcpu;
		printf("m=%d - SHT x%d time : %d, ticks : %.0f\n", m, SHT_ITER, ((int )tcpu)/SHT_ITER, elapsed(tik1,tik0)/SHT_ITER);

/*		tcpu = clock();
		tik0 = getticks();
		for (jj=0; jj<SHT_ITER; jj++)	spat_to_SHm_dct(im, ShF, Slm);
		tik1 = getticks();
		tcpu = clock() - tcpu;
		printf("m=%d - SHT_dct x%d time : %d, ticks : %.0f\n\n", m, SHT_ITER, ((int )tcpu)/SHT_ITER, elapsed(tik1,tik0)/SHT_ITER);*/
	}

#endif

// test case...
	t = 1.0 / (RAND_MAX/2);
	for (i=0;i<NLM;i++) {
		Slm0[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
		Tlm0[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
		Slm[i] = Slm0[i];
	}

/// REGULAR
	printf("[REG] :: performing %d scalar SHT with NL evaluation\n", SHT_ITER);
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
	printf("[DCT] :: performing %d scalar SHT with NL evaluation\n", SHT_ITER);
#ifdef SHT_DCT_IT_OUT
	printf("         (IT_OUT variant : not working yet...)\n");
#endif
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

#define TEST_VECT_SHT
#ifdef TEST_VECT_SHT
	printf("[REG] :: performing %d vector SHT\n", SHT_ITER);
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

	printf("[DCT] :: performing %d vector SHT\n", SHT_ITER);
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
}

