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
	long int ir;
// Allocate Spatial Fields.
	for (ir = 0; ir < NR; ir++) {
		BF[ir] = (complex double *) fftw_malloc( 3*NR*(NPHI/2+1)*NLAT * sizeof(complex double));	// 3 components
		B[ir] = (double *) BF;	// alias for inplace.
		Blm[ir] = (complex double *) malloc( 2*NLM * sizeof(complex double));	// Pol/Tor
	}
}

void PolTor_to_spat(complex double **Plm, complex double **Tlm, double** Br, double** Bt, double** Bp)
{
	complex double Q[NLM];	// l(l+1) * P/r
	complex double S[NLM];	// dP/dr + P/r = 1/r.d(rP)/dr = Wr(P)
	long int ir,lm;

	for (ir=0; ir<NR; ir++) {
		for (lm=0; lm<NLM; lm++) {		// Solenoidal deduced from radial derivative of Poloidal
			S[lm] = Wr[ir].l*Plm[ir-1][lm] + Wr[ir].d*Plm[ir][lm] + Wr[ir].u*Plm[ir+1][lm];
			Q[lm] = r_1[ir]*l2[lm] * Plm[ir][lm];
		}
		SH_to_spat(Q,(complex double *) Br[ir]);
		SHsphertor_to_spat(S, Tlm[ir], (complex double *) Bt[ir], (complex double *) Bp[ir]);
	}
}

void spat_to_PolSphTor(double** Br, double** Bt, double** Bp, complex double **Plm, complex double **Slm, complex double **Tlm)
{
	complex double Q[NLM];	// l(l+1) * P/r
	long int ir,lm;

	for (ir=0; ir<NR; ir++) {
		spat_to_SHsphertor((complex double *) Bt[ir], (complex double *) Bp[ir], Slm[ir], Tlm[ir]);
		spat_to_SH((complex double *) Br[ir], Q);
		for (lm=0; lm<NLM; lm++) {
			Plm[ir][lm] = r[ir]*l_2[lm] * Q[lm];		// poloidal
		}
	}
}


void PolTor_to_rot_spat(complex double **Plm, complex double **Tlm, double** Br, double** Bt, double** Bp)
{
	complex double Q[NLM];	// l(l+1) * T/r
	complex double S[NLM];	// dT/dr + T/r
	complex double T[NLM];	//  [ l(l+1)/r^2 - 1/r d2/dr2(r .) ] P
	long int ir,lm;

	for (ir=0; ir<NR; ir++) {
		for (lm=0; lm<NLM; lm++) {		// Solenoidal deduced from radial derivative of Poloidal
			Q[lm] = r_1[ir]*l2[lm] * Tlm[ir][lm];
			S[lm] = Wr[ir].l*Tlm[ir-1][lm] + Wr[ir].d*Tlm[ir][lm] + Wr[ir].u*Tlm[ir+1][lm];
			T[lm] = (r_2[ir]*l2[lm] - r_1[ir]*D2r[ir].d)*Plm[ir][lm] - r_1[ir]*D2r[ir].l*Plm[ir-1][lm] - r_1[ir]*D2r[ir].u*Plm[ir+1][lm];
		}
/* Pour Coriolis : ez ^ u
	ez = cos(theta).er - sin(theta).etheta
	cos theta = Y(m=0,l=1) * 2*sqrt(pi/3)		>> peut etre rajouté à Qlm. (=> Br)
	-sin theta = dY(m=0,l=1)/dt * 2*sqrt(pi/3)	>> peut etre rajouté à Slm  (=> Bt)
*/
		// add Background Vorticity for Coriolis Force (l=1, m=0)
		//	(double) Q[1] += Omega0 * 2.0*sqrt(pi/3.0);
		//	(double) S[1] += Omega0 * 2.0*sqrt(pi/3.0);
		SH_to_spat(Q,(complex double *) Br[ir]);
		SHsphertor_to_spat(S, T, (complex double *) Bt[ir], (complex double *) Bp[ir]);
	}
}

// spatial to curl : only for ir = 1 .. NR-1
//	Pol <- Tor
//	Tor <- Q/r - 1/r.d(rS)/dr
void spat_to_rot_PolTor(double** Br, double** Bt, double** Bp, complex double **Plm, complex double **Tlm)
{
	complex double Q[NLM];
	complex double Sl[NLM], Sd[NLM], Su[NLM];
	complex double* St;	// temp pointer.
	long int ir,lm;

	ir = 0;
		spat_to_SHsphertor((complex double *) Bt[ir], (complex double *) Bp[ir], Sd, Su);	// discard Plm[0]
		spat_to_SHsphertor((complex double *) Bt[ir+1], (complex double *) Bp[ir+1], Su, Plm[ir+1]);
	for (ir=1; ir<NR-1; ir++) {
		St = Sl;   Sl = Sd;   Sd = Su;   Su = St;		// rotating buffers.
		spat_to_SHsphertor((complex double *) Bt[ir+1], (complex double *) Bp[ir+1], Su, Plm[ir+1]);
		spat_to_SH((complex double *) Br[ir], Q);
		for (lm=0; lm<NLM; lm++) {
			Tlm[ir][lm] = r_1[ir]*Q[lm] - (Wr[ir].l * Sl[lm] + Wr[ir].d * Sd[lm] + Wr[ir].u * Su[lm]);
		}
	}
}

void NLU(complex double **Plm, complex double **Tlm, double** NLr, double** NLt, double** NLp)
{
	complex double Su[NLM];
	complex double Sw[NLM];
	complex double Tw[NLM];
	double NLtmp[2*(NPHI/2+1)*NLAT];	// should be FFTW allocated
	long int ir,lm;

/* Pour Coriolis : ez ^ u
	ez = cos(theta).er - sin(theta).etheta
	cos theta = Y(m=0,l=1) * 2*sqrt(pi/3)		>> peut etre rajouté à Qlm. (=> Br)
	-sin theta = dY(m=0,l=1)/dt * 2*sqrt(pi/3)	>> peut etre rajouté à Slm  (=> Bt)
*/
	for (ir=0; ir<NR; ir++) {
		for (lm=0; lm<NLM; lm++) {
			Su[lm] = Wr[ir].l*Plm[ir-1][lm] + Wr[ir].d*Plm[ir][lm] + Wr[ir].u*Plm[ir+1][lm];
			Tw[lm] = (r_2[ir]*l2[lm] - r_1[ir]*D2r[ir].d)*Plm[ir][lm] - r_1[ir]*D2r[ir].l*Plm[ir-1][lm] - r_1[ir]*D2r[ir].u*Plm[ir+1][lm];
			Sw[lm] = Wr[ir].l*Tlm[ir-1][lm] + Wr[ir].d*Tlm[ir][lm] + Wr[ir].u*Tlm[ir+1][lm];
		}
//		(double) Sw[1] += Omega0 * 2.0*sqrt(pi/3.0);	// add Background Vorticity for Coriolis Force (l=1, m=0)
		SHsphertor_to_spat(Su, Tlm[ir], (complex double *) Nlt[ir], (complex double *) Nlp[ir]);
		SHsphertor_to_spat(Sw, Tw, (complex double *) NLr[ir], (complex double *) Nltmp);
		for(lm = 0; lm < NPHI*NLAT; lm++) {
			NLt[ir][lm] *= NLr[ir][lm];	NLp[ir][i] *= NLtmp[lm];
		}
		for (lm=0; lm<NLM; lm++) {
			Su[lm] = r_1[ir]*l2[lm] * Plm[ir][lm];
			Sw[lm] = r_1[ir]*l2[lm] * Tlm[ir][lm];
		}
//		(double) Sw[1] += Omega0 * 2.0*sqrt(pi/3.0);	// add Background Vorticity for Coriolis Force (l=1, m=0)
		SH_to_spat(Su,(complex double *) NLr[ir]);
		SH_to_spat(Sw,(complex double *) NLtmp);
		for(lm = 0; lm < NPHI*NLAT; lm++)
			NLr[ir][lm] *= NLtmp[lm];
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

/*
      do ir=1,nr-1
         Aa(ir) = -Ek_Pm*Lra(ir)/2.
         do l=1,lmax
           Ab(l,ir) = -Ek_Pm*(Lrb(ir)-l*(l+1)/r(ir)**2 )/2. +1.0/dtB
         enddo
         Ac(ir) = -Ek_Pm*Lrc(ir)/2.
      enddo
*/

void init_Bmatrix()
{
	double dx_1,dx_2;
	int i,l;

// Boundary conditions
//	r=0 : T=0, P=0  => not computed at i=0.
//	r=1 : T=0, dP/dr= -(l+1)/r P  (insulator) => only P is computed.

	for(i=1; i<NR-1; i++) {
		for (l=0; l<=LMAX; l++) {
			MB[i*(LMAX+1) +l].l =           0.5*eta*Lr[i].l;
			MB[i*(LMAX+1) +l].d = 1.0/dtB + 0.5*eta*(Lr[i].d - r_2[i]*l2[l]);
			MB[i*(LMAX+1) +l].u =           0.5*eta*Lr[i].u;
		}
	}
	i = NR-1;		// CL poloidale : dP/dr = -(l+1)/r P => permet d'approximer d2P/dr2 de maniere discrete avec 2 points.
		dx_1 = 1.0/(r[i]-r[i-1]);	dx_2 = dx_1*dx_1;
		for (l=0; l<=LMAX; l++) {
			MB[i*(LMAX+1) +l].l =           eta * dx_2;
			MB[i*(LMAX+1) +l].d = 1.0/dtB + eta*( -dx_2 - (l+1.0)*r_1[i]*dx_1 -(l+1.0 + 0.5*l2[l])*r_2[i] );
			MB[i*(LMAX+1) +l].u = 0.0;
		}

	for(i=1; i<NR; i++) {
		for (l=0; l<=LMAX; l++) {
			MBp_1[i*(LMAX+1) +l].l =         - MB[i*(LMAX+1) +l].l;
			MBp_1[i*(LMAX+1) +l].d = 2.0/dtB - MB[i*(LMAX+1) +l].d;
			MBp_1[i*(LMAX+1) +l].u =         - MB[i*(LMAX+1) +l].u;
			MBt_1[i*(LMAX+1) +l].l = MBp_1[i*(LMAX+1) +l].l;
			MBt_1[i*(LMAX+1) +l].d = MBp_1[i*(LMAX+1) +l].d;
			MBt_1[i*(LMAX+1) +l].u = MBp_1[i*(LMAX+1) +l].u;
		}
	}
	TriDec(MBt_1, NR-2);
}


int main()
{
	double t,tmax;
	int i,im,m,l,jj;

	init_SH();
	init_rad_sph(rmin, rmax);
	init_Bmatrix();

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

