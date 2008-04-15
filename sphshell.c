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

#define NR 100

#include "grid.c"

struct TriDiagL *MB, *MB_1;

complex double* BF[NR];
double* B[NR];
complex double* Blm[NR];
complex double* RHSlm[NR];

struct BaseField {
	long int mmax,lmax;
	complex double *ylm;
};

double eta = 1.0;	// magnetic diffusivity.
double dtB = 0.001;	// time step for magnetic field.

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
		SHsphtor_to_spat(S, Tlm[ir], (complex double *) Bt[ir], (complex double *) Bp[ir]);
	}
}

void spat_to_PolSphTor(double** Br, double** Bt, double** Bp, complex double **Plm, complex double **Slm, complex double **Tlm)
{
	complex double Q[NLM];	// l(l+1) * P/r
	long int ir,lm;

	for (ir=0; ir<NR; ir++) {
		spat_to_SHsphtor((complex double *) Bt[ir], (complex double *) Bp[ir], Slm[ir], Tlm[ir]);
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

	for (ir=1; ir <= NR-2; ir++) {
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
		SHsphtor_to_spat(S, T, (complex double *) Bt[ir], (complex double *) Bp[ir]);
	}
}

// spatial to curl : only for ir = 1 .. NR-2
//	Pol <- Tor
//	Tor <- Q/r - 1/r.d(rS)/dr
void spat_to_rot_PolTor(double** Br, double** Bt, double** Bp, complex double **Plm, complex double **Tlm)
{
	complex double Q[NLM];		// Q
	complex double S[NLM*3];	// buffers for S.
	complex double *St, *Sl, *Sd, *Su;	// pointers to S (temp, lower, diag, upper)
	long int ir,lm;

	Sl = S;   Sd = S + NLM;   Su = S + 2*NLM;	// init pointers to buffer.

	ir = 0;
		spat_to_SHsphtor((complex double *) Bt[ir], (complex double *) Bp[ir], Sd, Su);	// discard Plm[0]
		spat_to_SHsphtor((complex double *) Bt[ir+1], (complex double *) Bp[ir+1], Su, Plm[ir+1]);
	for (ir=1; ir <= NR-2; ir++) {
		St = Sl;	Sl = Sd;	Sd = Su;	Su = St;		// rotate buffers.
		spat_to_SHsphtor((complex double *) Bt[ir+1], (complex double *) Bp[ir+1], Su, Plm[ir+1]);
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
	for (ir=1; ir <= NR-2; ir++) {
		for (lm=0; lm<NLM; lm++) {
			Su[lm] = Wr[ir].l*Plm[ir-1][lm] + Wr[ir].d*Plm[ir][lm] + Wr[ir].u*Plm[ir+1][lm];
			Tw[lm] = (r_2[ir]*l2[lm] - r_1[ir]*D2r[ir].d)*Plm[ir][lm] - r_1[ir]*D2r[ir].l*Plm[ir-1][lm] - r_1[ir]*D2r[ir].u*Plm[ir+1][lm];
			Sw[lm] = Wr[ir].l*Tlm[ir-1][lm] + Wr[ir].d*Tlm[ir][lm] + Wr[ir].u*Tlm[ir+1][lm];
		}
//		(double) Sw[1] += Omega0 * 2.0*sqrt(pi/3.0);	// add Background Vorticity for Coriolis Force (l=1, m=0)
		SHsphtor_to_spat(Su, Tlm[ir], (complex double *) NLt[ir], (complex double *) NLp[ir]);
		SHsphtor_to_spat(Sw, Tw, (complex double *) NLr[ir], (complex double *) NLtmp);
		for(lm = 0; lm < NPHI*NLAT; lm++) {
			NLt[ir][lm] *= NLr[ir][lm];	NLp[ir][lm] *= NLtmp[lm];
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

	MB = (struct TriDiagL *) malloc( NR* sizeof(struct TriDiagL));
	MB_1 = (struct TriDiagL *) malloc( NR* sizeof(struct TriDiagL));

// Boundary conditions
//	r=0 : T=0, P=0  => not computed at i=0.
//	r=1 : T=0, dP/dr= -(l+1)/r P  (insulator) => only P is computed.

	for(i=1; i<NR-1; i++) {
					MB[i].l =              0.5*eta*Lr[i].l;
		for (l=0; l<=LMAX; l++)	MB[i].d[l] = 1.0/dtB + 0.5*eta*(Lr[i].d - r_2[i]*l*(l+1));
					MB[i].u =              0.5*eta*Lr[i].u;
	}
	i = NR-1;		// CL poloidale : dP/dr = -(l+1)/r P => permet d'approximer d2P/dr2 de maniere discrete avec 2 points.
		dx_1 = 1.0/(r[i]-r[i-1]);	dx_2 = dx_1*dx_1;
					MB[i].l =              eta * dx_2;
		for (l=0; l<=LMAX; l++) MB[i].d[l] = 1.0/dtB + eta*( -dx_2 - (l+1.0)*r_1[i]*dx_1 -(l+1.0 + 0.5*l*(l+1))*r_2[i] );
					MB[i].u = 0.0;

	for(i=1; i<NR; i++) {
					MB_1[i].l =            - MB[i].l;
		for (l=0; l<=LMAX; l++) MB_1[i].d[l] = 2.0/dtB - MB[i].d[l];
					MB_1[i].u =            - MB[i].u;
	}
	TriDec(MB_1, 1, NR-1);		// for Btor : 1 to NR-2, for Bpol : 1 to NR-1
// for poloidal, use : cTriSolve(MB_1, RHS, Plm, 1, NR-1)
// for toroidal, use : cTriSolve(MB_1, RHS, Tlm, 1, NR-2)

/*
	for(i=1; i<NR-1; i++) {
					MB[i].l =              eta*Lr[i].l;
		for (l=0; l<=LMAX; l++)	MB[i].d[l] =           eta*(Lr[i].d - r_2[i]*l2[l]);
					MB[i].u =              eta*Lr[i].u;
	}
	i = NR-1;		// CL poloidale : dP/dr = -(l+1)/r P => permet d'approximer d2P/dr2 de maniere discrete avec 2 points.
		dx_1 = 1.0/(r[i]-r[i-1]);	dx_2 = dx_1*dx_1;
					MB[i].l =    eta * dx_2;
		for (l=0; l<=LMAX; l++) MB[i].d[l] = eta*( -dx_2 - (l+1.0)*r_1[i]*dx_1 -(l+1.0 + 0.5*l2[l])*r_2[i] );
					MB[i].u = 0.0;
*/
}


void alloc_fields()
{
	long int ir;
// Allocate Spatial Fields.
	for (ir = 0; ir < NR; ir++) {
		BF[ir] = (complex double *) fftw_malloc( 3*NR*(NPHI/2+1)*NLAT * sizeof(complex double));	// 3 components
		B[ir] = (double *) BF;	// alias for inplace.
		Blm[ir] = (complex double *) malloc( 2*NLM * sizeof(complex double));	// Pol/Tor
		RHSlm[ir] = (complex double *) malloc( 2*NLM * sizeof(complex double));	// Pol/Tor
	}
}



int main()
{
	double t0,t1;
	int i,im,m,l,jj, it;
	FILE *fp; 

	init_SH();
	init_rad_sph(0.0, 1.0);
	init_Bmatrix();

	alloc_fields();

	for (i=0; i<NR; i++) {
		for (l=0;l<NLM;l++) {
			Blm[i][l] = 0.0;
			RHSlm[i][l] = 0.0;
		}
		t0 = r[i]*pi;	t1 = r[i]*4.493409457909;
//		Blm[i][1] = sin(t0)/(t0*t0) - cos(t0)/t0;
		Blm[i][1] = r[i];
	}



	fp = fopen("prof0","w");
	for (i=0; i<NR; i++) {
		fprintf(fp,"%.6g ",Blm[i][1]);
	}
	fclose(fp);

	for (it=0; it< 10000; it++) {
		t0 = t1;
		t1 = Blm[NR/2][1];
		printf("%g :: tx=%g\n",t1,log(t1/t0)/dtB);

		cTriMul(MB, Blm, RHSlm, 1, NR-1);
//		cTriMul(MB, RHSlm, Blm, 1, NR-1);
		cTriSolve(MB_1, RHSlm, Blm, 1, NR-1);
	}

	fp = fopen("prof1","w");
	for (i=0; i<NR; i++) {
		fprintf(fp,"%.6g ",Blm[i][1]);
	}
	fclose(fp);

}

