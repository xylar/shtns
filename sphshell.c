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


#include "SHT.c"

complex double* BF[NR];
double* B[NR];
complex double* Blm[NR];

void alloc_fields()
{
	int ir;
// Allocate Spatial Fields.
	for (ir = 0; ir < NR; ir++) {
		BF[ir] = (complex double *) fftw_malloc( 3*NR*(NPHI/2+1)*NLAT * sizeof(complex double));	// 3 components
		B[ir] = (double *) BF;	// alias for inplace.
		Blm[ir] = (complex double *) malloc( 2*NLM * sizeof(complex double));	// Pol/Tor
	}
}

void poltor_to_spat(complex double **Plm, complex double **Tlm, double** Br, double** Bt, double** Bp)
{
	complex double Qlm[NLM];	// l(l+1) * P/r
	complex double Slm[NLM];	// dP/dr + P/r
	int ir,lm;

	for (ir=0; ir<NR; ir++) {
		for (lm=0; lm<NLM; lm++) {		// Solenoidal deduced from radial derivative of Poloidal
			Slm[lm] = Gr[ir].l*Plm[ir-1][lm] + Gr[ir].d*Plm[ir][lm] + Gr[ir].u*Plm[ir+1][lm];
			Qlm[lm] = r_1[ir]*l2[lm] * Plm[ir][lm];
		}
		SH_to_spat(Qlm,(complex double *) Br[ir]);
		SHsphertor_to_spat(Slm, Tlm[ir], (complex double *) Bt[ir], (complex double *) Bp[ir]);
	}
}

void spat_to_PolSphTor(double** Br, double** Bt, double** Bp, complex double **Plm, complex double **Slm, complex double **Tlm)
{
	complex double Q[NLM];	// l(l+1) * P/r
	complex double S[NLM];	// dP/dr + P/r
	complex double T[NLM];
	int ir,lm;

	for (ir=0; ir<NR; ir++) {
		spat_to_SHsphertor((complex double *) Bt[ir], (complex double *) Bp[ir], S[ir], T[ir]);
		spat_to_SH((complex double *) Br[ir], Q);
		for (lm=0; lm<NLM; lm++) {
			Plm[ir][lm] = r[ir]*l_2[lm] * Q[lm];		// poloidal
		}
	}
}


void poltor_to_rot_spat(complex double **Plm, complex double **Tlm, double** Br, double** Bt, double** Bp)
{
	complex double Q[NLM];	// l(l+1) * T/r
	complex double S[NLM];	// dT/dr + T/r
	complex double T[NLM];	//  [ l(l+1)/r^2 - 1/r d2/dr2(r .) ] P
	int ir,lm;

	for (ir=0; ir<NR; ir++) {
		for (lm=0; lm<NLM; lm++) {		// Solenoidal deduced from radial derivative of Poloidal
			Q[lm] = r_1[ir]*l2[lm] * Tlm[ir][lm];
			S[lm] = Gr[ir].l*Tlm[ir-1][lm] + Gr[ir].d*Tlm[ir][lm] + Gr[ir].u*Tlm[ir+1][lm];
			T[lm] = (r_2[ir]*l2[lm] - r_1[ir]*D2r[ir].d)*Plm[ir][lm] - r_1[ir]*D2r[ir].l*Plm[ir-1][lm] - r_1[ir]*D2r[ir].u*Plm[ir+1][lm];
		}
		SH_to_spat(Q,(complex double *) Br[ir]);
		SHsphertor_to_spat(S, T, (complex double *) Bt[ir], (complex double *) Bp[ir]);
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

int main()
{
	double t,tmax;
	int i,im,m,l,jj;

	init_SH();

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

// test Ylm :
	im = 0; l=0; m=im*MRES;
	write_vect("y00",&iylm[im][(l-m)*NLAT/2],NLAT/2);
	im = 0; l=1; m=im*MRES;
	write_vect("y10",&iylm[im][(l-m)*NLAT/2],NLAT/2);
	im = 0; l=LMAX; m=im*MRES;
	write_vect("yLmax0",&iylm[im][(l-m)*NLAT/2],NLAT/2);
	im = MMAX; m=im*MRES; l=m;
	write_vect("ymmax-mmax",&iylm[im][(l-m)*NLAT/2],NLAT/2);
	im = 3; l=8; m=im*MRES;
	write_vect("y83",&iylm[im][(l-m)*NLAT/2],NLAT/2);
	im = 30; l=65; m=im*MRES;
	write_vect("y6530",&iylm[im][(l-m)*NLAT/2],NLAT/2);

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
		SH_to_spat(Slm,ShF);
//		write_mx("sph",Sh,NPHI,NLAT);

// analyse (direct legendre)
		spat_to_SH(ShF,Slm);
	}

	write_vect("ylm",Slm,LMMAX*2);
}
