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

// number of radial grid points.
#define NR 100
// radial points for inner core (NG = 0 : no inner core)
#define NG 0
#define NU (NR-NG)

#include "grid.c"

// boundary conditions constants.
#define BC_NO_SLIP 0
#define BC_FREE_SLIP 1
#define BC_MAGNETIC 2

struct TriDiagL *MB, *MB_1, *MUt, *MUt_1;
struct PentaDiag *MUp[NU], *MUp_1[NU];

double eta = 1.0;	// magnetic diffusivity.
double nu = 1.0;	// kinematic viscosity.
double dtB = 0.0002;	// time step for magnetic field.
double dtU = 0.0002;	// time step for velocity field.

double Omega0 = 0.0;		// global rotation rate (of outer boundary) => Coriolis force .
double DeltaOmega = 0.0;	// differential rotation (of inner core)

/*
struct VectField {
	double* r,t,p;
};

struct PolTor {
	complex double* p,t;
};

struct VectField U,B,J;
struct PolTor Ulm, Blm, NLb1, NLb2, NLu1, NLu2;
*/

// Magnetic field representation.
double *Br[NR], *Bt[NR], *Bp[NR];	// Spatial representation
complex double *BPlm[NR], *BTlm[NR];	// Spherical Harmonics representation
complex double *NLbp1[NR], *NLbp2[NR], *NLbt1[NR], *NLbt2[NR];	// Adams-Bashforth : 2 non-linear terms stored.
// Velocity field representation.
double *Ur[NR], *Ut[NR], *Up[NR];
complex double *UPlm[NR], *UTlm[NR];
complex double *NLup1[NR], *NLup2[NR], *NLut1[NR], *NLut2[NR];
// spare field J (for current density, or vorticity)
double *Jr[NR], *Jt[NR], *Jp[NR];


// Allocate memory for Magnetic field
void alloc_Bfields()
{
	long int ir;
	for (ir = 0; ir < NR; ir++) {		// shell by shell allocation.
		Br[ir] = (double *) fftw_malloc( 3*(NPHI/2+1)*NLAT * sizeof(complex double));	// 3 components
		Bt[ir] = Br[ir] + 2*(NPHI/2+1)*NLAT;
		Bp[ir] = Bt[ir] + 2*(NPHI/2+1)*NLAT;

		BPlm[ir] = (complex double *) malloc( 6*NLM * sizeof(complex double));	// Pol/Tor
		BTlm[ir] = BPlm[ir] + NLM;
		NLbp1[ir] = BTlm[ir] + NLM;
		NLbt1[ir] = NLbp1[ir] + NLM;
		NLbp2[ir] = NLbt1[ir] + NLM;
		NLbt2[ir] = NLbp2[ir] + NLM;
	}
}

// Allocate memory for Velocity field
void alloc_Ufields()
{
	long int ir;

	// not defined inside the inner core.
	for (ir = 0; ir < NG; ir++) {
		Ur[ir] = NULL;	Ut[ir] = NULL;	Up[ir] = NULL;
		UPlm[ir] = NULL;	UTlm[ir] = NULL;
		NLup1[ir] = NULL;	NLup2[ir] = NULL;	NLut1[ir] = NULL;	NLut2[ir] = NULL;
	}

	for (ir = NG; ir < NR; ir++) {		// shell by shell allocation.
		Ur[ir] = (double *) fftw_malloc( 3*(NPHI/2+1)*NLAT * sizeof(complex double));	// 3 components
		Ut[ir] = Ur[ir] + 2*(NPHI/2+1)*NLAT;
		Up[ir] = Ut[ir] + 2*(NPHI/2+1)*NLAT;

		UPlm[ir] = (complex double *) malloc( 6*NLM * sizeof(complex double));	// Pol/Tor
		UTlm[ir] = UPlm[ir] + NLM;
		NLup1[ir] = UPlm[ir] + 2*NLM;
		NLut1[ir] = UPlm[ir] + 3*NLM;
		NLup2[ir] = UPlm[ir] + 4*NLM;
		NLut2[ir] = UPlm[ir] + 5*NLM;
	}
}

// compute the cartesian coordinates of a vector field at r=0, from its Poloidal value component
inline void Pol_to_cart_spat0(complex double **Plm, double *Vx, double *Vy, double *Vz)
{
	double dx;
	long int lm;

	dx = 2.0*sqrt(3.0) / r[1];
	lm = LM(1,1);	// l=1, m=1
	*Vx = dx * creal( Plm[1][lm] );
	*Vy = dx * cimag( Plm[1][lm] );
	lm = LM(1,0);	// l=1, m=0
	*Vz = dx * (double) Plm[1][lm];
}

// compute the value of the spheroidal scalar at r=0 from the cartesian representation of a vector field at the center.
inline void cart_spat0_to_Sph(double Vx, double Vy, double Vz, complex double *S10, complex double *S11)
{
	double dx;

	dx = 1.0/sqrt(3.0);
	*S10 = Vz * dx;		// l=1, m=0
	*S11 = (Vx + I*Vy) *dx;	// l=1, m=1
}

#define CALC_U(bc)   PolTor_to_spat(UPlm, UTlm, Ur, Ut, Up, NG, NR-1, (bc))
#define CALC_B       PolTor_to_spat(BPlm, BTlm, Br, Bt, Bp, 0, NR-1, BC_MAGNETIC)


// compute spatial field from Poloidal/Toroidal representation.
void PolTor_to_spat(complex double **Plm, complex double **Tlm, double** Vr, double** Vt, double** Vp, long int istart, long int iend, long int BC)
{
	complex double Q[NLM];	// l(l+1) * P/r
	complex double S[NLM];	// dP/dr + P/r = 1/r.d(rP)/dr = Wr(P)
	double dr;
	long int ir,lm,l;

	ir = istart;
	if (r[ir] == 0.0) {	// field at r=0 : S = 2.dP/dr (l=1 seulement) => store Bx, By, Bz (cartesian coordinates)
		Pol_to_cart_spat0(Plm, Vr[ir],Vt[ir],Vp[ir]);
		for (lm=1; lm<NLAT*NPHI; lm++) {	// zero for the rest.
			Vr[ir][lm] = 0.0;	Vt[ir][lm] = 0.0;	Vp[ir][lm] = 0.0;
		}
	} else if (BC == BC_NO_SLIP) {
		// Velocity field (no slip => Ur=0, Ut=0, Up=r.sint.DeltaOmega)
		for (lm=0; lm<NLAT*NPHI; lm++) {
			Vr[ir][lm] = 0.0;	Vt[ir][lm] = 0.0;
		}
		for (lm=0; lm<NPHI; lm++) {
			for (l=0; l<NLAT; l++) Vp[ir][lm*NLAT+l] = r[ir]*st[l]*DeltaOmega;
		}
	} else if (BC == BC_FREE_SLIP) {
		// Velocity field (free slip BC) : P = 0, d(T/r)/dr = 0   => Ur=0
		dr = r[ir+1]/(r[ir]*(r[ir+1]-r[ir]));
		for (lm=0; lm<NLM; lm++) {
			S[lm] = dr * Plm[ir+1][lm];
		}
		SHsphtor_to_spat(S, Tlm[ir], (complex double *) Vt[ir], (complex double *) Vp[ir]);
		for (lm=0; lm<NLAT*NPHI; lm++) {
			Vr[ir][lm] = 0.0;
		}
	}

	for (ir=istart+1; ir<iend; ir++) {
		for (lm=0; lm<NLM; lm++) {		// Solenoidal deduced from radial derivative of Poloidal
			S[lm] = Wr[ir].l*Plm[ir-1][lm] + Wr[ir].d*Plm[ir][lm] + Wr[ir].u*Plm[ir+1][lm];
			Q[lm] = r_1[ir]*l2[lm] * Plm[ir][lm];
		}
		SH_to_spat(Q,(complex double *) Vr[ir]);
		SHsphtor_to_spat(S, Tlm[ir], (complex double *) Vt[ir], (complex double *) Vp[ir]);
	}

	ir = iend;
	if (BC == BC_NO_SLIP) {
		// Velocity field (no slip => U = 0)
		for (lm=0; lm<NLAT*NPHI; lm++) {
			Vr[ir][lm] = 0.0;	Vt[ir][lm] = 0.0;	Vp[ir][lm] = 0.0;
		}
		return;
	} else if (BC == BC_FREE_SLIP) {
		// Velocity field (free slip BC) : P = 0, d(T/r)/dr = 0   => Ur=0
		dr = r[ir-1]/(r[ir]*(r[ir-1]-r[ir]));
		for (lm=0; lm<NLM; lm++) {
			S[lm] = dr * Plm[ir-1][lm];
		}
		SHsphtor_to_spat(S, Tlm[ir], (complex double *) Vt[ir], (complex double *) Vp[ir]);
		for (lm=0; lm<NLAT*NPHI; lm++) {
			Vr[ir][lm] = 0.0;
		}
		return;
	} else if (BC == BC_MAGNETIC) {
		// Magnetic field (insulator BC) : T=0, dP/dr = -(l+1)/r.P => S = -lP/r
		for (lm=0, l=0; lm<NLM; lm++) {
			S[lm] = -l*Plm[ir][lm]*r_1[ir];
			if (l < LMAX) { l++; } else { l=0; }
			Q[lm] = r_1[ir]*l2[lm] * Plm[ir][lm];
		}
		SH_to_spat(Q,(complex double *) Vr[ir]);
		SHsph_to_spat(S, (complex double *) Vt[ir], (complex double *) Vp[ir]);
	}
}

// compute vorticity field from poloidal/toroidal components of velocity [ie poltor_to_rot_spat applied to U]
// IN: Plm, Tlm : pol/tor components of velocity
//     Omega0 : background global rotation (solid body rotation rate) to include in vorticity
//     BC : boundary condition (BC_NO_SLIP or BC_FREE_SLIP)
// OUT: Vr,Vt,Vp : r,theta,phi components of vorticity field
void calc_Vort(complex double **Plm, complex double **Tlm, double** Vr, double** Vt, double** Vp, double Omega0, long int BC)
{
	complex double Q[NLM];	// l(l+1) * T/r
	complex double S[NLM];	// dT/dr + T/r
	complex double T[NLM];	//  [ l(l+1)/r^2 - 1/r d2/dr2(r .) ] P
	double dr;
	long int ir,lm;

	ir = NG;
	if (r[ir] = 0.0) {	// field at r=0 : S = 2.dT/dr (l=1 seulement) => store Vx, Vy, Vz (cartesian coordinates)
		Pol_to_cart_spat0(Tlm, Vr[ir],Vt[ir],Vp[ir]);
		for (lm=1; lm<NLAT*NPHI; lm++) {	// zero for the rest.
			Vr[ir][lm] = 0.0;	Vt[ir][lm] = 0.0;	Vp[ir][lm] = 0.0;
		}
	} else if (BC == BC_NO_SLIP) {
		for (lm=0; lm<NLM; lm++) {		// Solenoidal deduced from radial derivative of Toroidal
			Q[lm] = r_1[ir]*l2[lm] * Tlm[ir][lm];
			S[lm] = Wr[ir].l*Tlm[ir-1][lm] + Wr[ir].d*Tlm[ir][lm] + Wr[ir].u*Tlm[ir+1][lm];
			T[lm] = (r_2[ir]*l2[lm] - r_1[ir]*D2r[ir].d)*Plm[ir][lm] - r_1[ir]*D2r[ir].l*Plm[ir-1][lm] - r_1[ir]*D2r[ir].u*Plm[ir+1][lm];
		}
	}

	for (ir=istart+1; ir < iend; ir++) {
		for (lm=0; lm<NLM; lm++) {		// Solenoidal deduced from radial derivative of Poloidal
			Q[lm] = r_1[ir]*l2[lm] * Tlm[ir][lm];
			S[lm] = Wr[ir].l*Tlm[ir-1][lm] + Wr[ir].d*Tlm[ir][lm] + Wr[ir].u*Tlm[ir+1][lm];
			T[lm] = (r_2[ir]*l2[lm] - r_1[ir]*D2r[ir].d)*Plm[ir][lm] - r_1[ir]*D2r[ir].l*Plm[ir-1][lm] - r_1[ir]*D2r[ir].u*Plm[ir+1][lm];
		}
/* Pour Coriolis : ez ^ u
	ez = cos(theta).er - sin(theta).etheta
	cos theta = Y(m=0,l=1) * 2*sqrt(pi/3)		>> peut etre rajouté à Qlm. (=> Vr)
	-sin theta = dY(m=0,l=1)/dt * 2*sqrt(pi/3)	>> peut etre rajouté à Slm  (=> Vt)
*/
		// add Background Vorticity for Coriolis Force (l=1, m=0)
		(double) Q[1] += Omega0 * 2.0*sqrt(pi/3.0);
		(double) S[1] += Omega0 * 2.0*sqrt(pi/3.0);
		SH_to_spat(Q,(complex double *) Vr[ir]);
		SHsphtor_to_spat(S, T, (complex double *) Vt[ir], (complex double *) Vp[ir]);
	}

	ir = iend;
}

void PolTor_to_rot_spat(complex double **Plm, complex double **Tlm, double** Vr, double** Vt, double** Vp, long int istart, long int iend, long int BC)
// Pol = Tor;  Tor = Lap.Pol
{
	complex double Q[NLM];	// l(l+1) * T/r
	complex double S[NLM];	// dT/dr + T/r
	complex double T[NLM];	//  [ l(l+1)/r^2 - 1/r d2/dr2(r .) ] P
	double dr;
	long int ir,lm;

	ir = istart;
	if (r[ir] = 0.0) {	// field at r=0 : S = 2.dT/dr (l=1 seulement) => store Bx, By, Bz (cartesian coordinates)
		Pol_to_cart_spat0(Tlm, Vr[ir],Vt[ir],Vp[ir]);
		for (lm=1; lm<NLAT*NPHI; lm++) {	// zero for the rest.
			Vr[ir][lm] = 0.0;	Vt[ir][lm] = 0.0;	Vp[ir][lm] = 0.0;
		}
	} else if (BC == BC_NO_SLIP) {
		for (lm=0; lm<NLM; lm++) {		// Solenoidal deduced from radial derivative of Toroidal
			Q[lm] = r_1[ir]*l2[lm] * Tlm[ir][lm];
			S[lm] = Wr[ir].l*Tlm[ir-1][lm] + Wr[ir].d*Tlm[ir][lm] + Wr[ir].u*Tlm[ir+1][lm];
			T[lm] = (r_2[ir]*l2[lm] - r_1[ir]*D2r[ir].d)*Plm[ir][lm] - r_1[ir]*D2r[ir].l*Plm[ir-1][lm] - r_1[ir]*D2r[ir].u*Plm[ir+1][lm];
		}
	}

	for (ir=istart+1; ir < iend; ir++) {
		for (lm=0; lm<NLM; lm++) {		// Solenoidal deduced from radial derivative of Poloidal
			Q[lm] = r_1[ir]*l2[lm] * Tlm[ir][lm];
			S[lm] = Wr[ir].l*Tlm[ir-1][lm] + Wr[ir].d*Tlm[ir][lm] + Wr[ir].u*Tlm[ir+1][lm];
			T[lm] = (r_2[ir]*l2[lm] - r_1[ir]*D2r[ir].d)*Plm[ir][lm] - r_1[ir]*D2r[ir].l*Plm[ir-1][lm] - r_1[ir]*D2r[ir].u*Plm[ir+1][lm];
		}
/* Pour Coriolis : ez ^ u
	ez = cos(theta).er - sin(theta).etheta
	cos theta = Y(m=0,l=1) * 2*sqrt(pi/3)		>> peut etre rajouté à Qlm. (=> Vr)
	-sin theta = dY(m=0,l=1)/dt * 2*sqrt(pi/3)	>> peut etre rajouté à Slm  (=> Vt)
*/
		// add Background Vorticity for Coriolis Force (l=1, m=0)
		//	(double) Q[1] += Omega0 * 2.0*sqrt(pi/3.0);
		//	(double) S[1] += Omega0 * 2.0*sqrt(pi/3.0);
		SH_to_spat(Q,(complex double *) Vr[ir]);
		SHsphtor_to_spat(S, T, (complex double *) Vt[ir], (complex double *) Vp[ir]);
	}

	ir = iend;
}

// spatial to curl : only for ir = 1 .. NR-2
//	Pol <- Tor
//	Tor <- Q/r - 1/r.d(rS)/dr
void spat_to_rot_PolTor(double** Br, double** Bt, double** Bp, complex double **Plm, complex double **Tlm, long int BC)
{
	complex double Q[NLM];		// Q
	complex double S[NLM*3];	// buffers for S.
	complex double *St, *Sl, *Sd, *Su;	// pointers to S (temp, lower, diag, upper)
	long int ir,lm;

	Sl = S;   Sd = S + NLM;   Su = S + 2*NLM;	// init pointers to buffer.

	ir = 0;
		for (lm=0; lm<NLM; lm++) {
			Plm[ir][lm] = 0.0;	Tlm[ir][lm] = 0.0;	Sd[lm] = 0.0;
		}
		lm = LM(1,0);	Sd[lm] = Bp[ir][0] /sqrt(3.0);		// S(l=1,m=0) = Vz / sqrt(3)
		lm = LM(1,1);	Sd[lm] = (Br[ir][0] + I*Bt[ir][0]) /sqrt(3.0);	// S(l=1,m=1) = (Vx+iVy)/sqrt(3)

		spat_to_SHsphtor((complex double *) Bt[ir+1], (complex double *) Bp[ir+1], Su, Plm[ir+1]);

	for (ir=1; ir < NR-2; ir++) {
		St = Sl;	Sl = Sd;	Sd = Su;	Su = St;		// rotate buffers.
		spat_to_SHsphtor((complex double *) Bt[ir+1], (complex double *) Bp[ir+1], Su, Plm[ir+1]);
		spat_to_SH((complex double *) Br[ir], Q);
		for (lm=0; lm<NLM; lm++) {
			Tlm[ir][lm] = r_1[ir]*Q[lm] - (Wr[ir].l * Sl[lm] + Wr[ir].d * Sd[lm] + Wr[ir].u * Su[lm]);
		}
	}

	ir = NR-2;
	if (BC == BC_NO_SLIP) {		// NR-1 : all zero. => P = 0, T=0.
		St = Sl;	Sl = Sd;	Sd = Su;	Su = St;		// rotate buffers.
		spat_to_SH((complex double *) Br[ir], Q);
		for (lm=0; lm<NLM; lm++) {
			Tlm[ir][lm] = r_1[ir]*Q[lm] - (Wr[ir].l * Sl[lm] + Wr[ir].d * Sd[lm]);		// Su=0
		}
		for (lm=0; lm<NLM; lm++) {
			Plm[ir+1][lm] = 0.0;	Tlm[ir+1][lm] = 0.0;	// for paranoia
		}
		return;
	} else if (BC == BC_FREE_SLIP) {
		St = Sl;	Sl = Sd;	Sd = Su;	Su = St;		// rotate buffers.
		spat_to_SHsphtor((complex double *) Bt[ir+1], (complex double *) Bp[ir+1], Su, Plm[ir+1]);
		spat_to_SH((complex double *) Br[ir], Q);
		for (lm=0; lm<NLM; lm++) {
			Tlm[ir][lm] = r_1[ir]*Q[lm] - (Wr[ir].l * Sl[lm] + Wr[ir].d * Sd[lm] + Wr[ir].u * Su[lm]);
		}
		ir = NR-1;	// P = 0
		St = Sl;	Sl = Sd;	Sd = Su;	Su = St;		// rotate buffers.
		spat_to_SH((complex double *) Br[ir], Q);
		for (lm=0; lm<NLM; lm++) {
			Tlm[ir][lm] = r_1[ir]*Q[lm] - (Wr[ir].l * Sl[lm] + Wr[ir].d * Sd[lm]);
		}
		for (lm=0; lm<NLM; lm++) {
			Plm[ir][lm] = 0.0;	// for paranoia
		}
		return;
	} else if (BC == BC_MAGNETIC) {
		St = Sl;	Sl = Sd;	Sd = Su;	Su = St;		// rotate buffers.
		spat_to_SHsphtor((complex double *) Bt[ir+1], (complex double *) Bp[ir+1], Su, Plm[ir+1]);
		spat_to_SH((complex double *) Br[ir], Q);
		for (lm=0; lm<NLM; lm++) {
			Tlm[ir][lm] = r_1[ir]*Q[lm] - (Wr[ir].l * Sl[lm] + Wr[ir].d * Sd[lm] + Wr[ir].u * Su[lm]);
		}
		ir = NR-1;	// T=0, P already set.
		for (lm=0; lm<NLM; lm++) {
			Tlm[ir][lm] = 0.0;	// for paranoia
		}
		return;
	}
}

/*
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
*/

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
	for (ir=NG+1; ir <= NR-2; ir++) {
		for (lm=0; lm<NLM; lm++) {
			Su[lm] = Wr[ir].l*Plm[ir-1][lm] + Wr[ir].d*Plm[ir][lm] + Wr[ir].u*Plm[ir+1][lm];
			Tw[lm] = (r_2[ir]*l2[lm] - r_1[ir]*D2r[ir].d)*Plm[ir][lm] - r_1[ir]*D2r[ir].l*Plm[ir-1][lm] - r_1[ir]*D2r[ir].u*Plm[ir+1][lm];
			Sw[lm] = Wr[ir].l*Tlm[ir-1][lm] + Wr[ir].d*Tlm[ir][lm] + Wr[ir].u*Tlm[ir+1][lm];
		}
		(double) Sw[1] += Omega0 * 2.0*sqrt(pi/3.0);	// add Background Vorticity for Coriolis Force (l=1, m=0)
		SHsphtor_to_spat(Su, Tlm[ir], (complex double *) NLt[ir], (complex double *) NLp[ir]);
		SHsphtor_to_spat(Sw, Tw, (complex double *) NLr[ir], (complex double *) NLtmp);
		for(lm = 0; lm < NPHI*NLAT; lm++) {
			NLt[ir][lm] *= NLr[ir][lm];	NLp[ir][lm] *= NLtmp[lm];
		}
		for (lm=0; lm<NLM; lm++) {
			Su[lm] = r_1[ir]*l2[lm] * Plm[ir][lm];
			Sw[lm] = r_1[ir]*l2[lm] * Tlm[ir][lm];
		}
		(double) Sw[1] += Omega0 * 2.0*sqrt(pi/3.0);	// add Background Vorticity for Coriolis Force (l=1, m=0)
		SH_to_spat(Su,(complex double *) NLr[ir]);
		SH_to_spat(Sw,(complex double *) NLtmp);
		for(lm = 0; lm < NPHI*NLAT; lm++)
			NLr[ir][lm] *= NLtmp[lm];
	}
}

// compute rot(u^B) : spatial B and U must be set.
// B is overwritten with u^B
void induction(double **Br, double **Bt, double **Bp, complex double **NLP, complex double **NLT)
{
	double vr,vt,vp;
	long int ir,lm;

	for (ir=0; ir<NG; ir++) {		// for the inner core : solid body rotation  Up = r.DeltaOmega
		vp = r[ir]*DeltaOmega;
		for (lm=0; lm<NPHI*NLAT; lm++) {
			vr = - vp*Bt[ir][lm];
			vt =   vp*Br[ir][lm];
			Br[ir][lm] = vr;	Bt[ir][lm] = vt;	Bp[ir][lm] = 0.0;
		}
	}
	for (ir=NG; ir<NR; ir++) {
		for (lm=0; lm<NPHI*NLAT; lm++) {
			vr = Ut[ir][lm]*Bp[ir][lm] - Up[ir][lm]*Bt[ir][lm];
			vt = Up[ir][lm]*Br[ir][lm] - Ur[ir][lm]*Bp[ir][lm];
			vp = Ur[ir][lm]*Bt[ir][lm] - Ut[ir][lm]*Br[ir][lm];
			Br[ir][lm] = vr;	Bt[ir][lm] = vt;	Bp[ir][lm] = vp;
		}
	}
/*	ir = NR-1;	// U = 0
		for (lm=0; lm<NPHI*NLAT; lm++) {
			Br[ir][lm] = 0.0;	Bt[ir][lm] = 0.0;	Bp[ir][lm] = 0.0;
		}*/
	spat_to_rot_PolTor(Br, Bt, Bp, NLP, NLT, BC_MAGNETIC);
}

// compute total Energy
double Energy(complex double **Plm, complex double **Tlm)
{
	double E,er;
	long int ir,lm;

	E = 0.0;
	for (ir=0; ir<NR; ir++) {
		er = 0.0;
		for(lm=0; lm<NLM; lm++)
			er += creal(Plm[ir][lm])*creal(Plm[ir][lm]) + cimag(Plm[ir][lm])*cimag(Plm[ir][lm]) +
				creal(Tlm[ir][lm])*creal(Tlm[ir][lm]) + cimag(Tlm[ir][lm])*cimag(Tlm[ir][lm]);
		E += er*dr[ir]*r[ir];
	}
	return E;
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

void write_HS(char *fn, complex double **HS)
{
	FILE *fp;
	int ir,lm;
	
	fp = fopen(fn,"w");
	for (ir=0;ir<NR;ir++) {
		for(lm=0;lm<NLM;lm++) {
			fprintf(fp,"%.6g %.6g ",creal(HS[ir][lm]),cimag(HS[ir][lm]));
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

void write_slice(char *fn, double **v, int im)
{
	FILE *fp;
	int i,j;

	fp = fopen(fn,"w");
		fprintf(fp,"0 ");			// first row = radius
		for(j=0;j<NLAT/2;j++) {
			fprintf(fp,"%.6g ",ct[j]);	// first line = cos(theta)
		}
		for(j=1;j<=NLAT/2;j++) {
			fprintf(fp,"-%.6g ",ct[NLAT/2-j]);	// first line = cos(theta)
		}
	for (i=0;i<NR;i++) {
		fprintf(fp,"\n%.6g ",r[i]);		// first row = radius
		for(j=0;j<NLAT;j++) {
			fprintf(fp,"%.6g ",v[i][im*NLAT + j]);		// data
		}
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
	long int i,l;

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
	i = NR-1;		// CL poloidale : dP/dr = -(l+1)/r P => allows to aproximate d2P/dr2 withe only 2 points.
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
}

// BC : boundary condition, can be BC_NO_SLIP or BC_FREE_SLIP
void init_Umatrix(long int BC)
{
	double dx_1,dx_2;
	double Lp0d[LMAX+1], Lp0u[LMAX+1];	// for r=0 condition, useful for Poloidal.
	long int i,l;

// allocate toroidal matrix
	MUt = (struct TriDiagL *) malloc( NU* sizeof(struct TriDiagL));
	MUt_1 = (struct TriDiagL *) malloc( NU* sizeof(struct TriDiagL));
// allocate poloidal matrix
	MUp[0] = (struct PentaDiag *) malloc( NU*(LMAX+1)* sizeof(struct PentaDiag));
	MUp_1[0] = (struct PentaDiag *) malloc( NU*(LMAX+1)* sizeof(struct PentaDiag));
	for (i=1; i<NU; i++) {
		MUp[i] = MUp[i-1] + (LMAX+1);
		MUp_1[i] = MUp_1[i-1] + (LMAX+1);
	}

// r=0 constraint :
//	T=0, P=0  => not computed at i=0.
//	for l=1 : d2/dr2 P = 0
//	for l>1 : dP/dr = 0
// Boundary conditions (NO SLIP)
//	P=0, T=0, dP/dr = 0  => not computed at boundary
//	with DeltaOmega (couette), at inner core : T(l=1,m=0) = r.DeltaOmega/sqrt(3)
// Boundary conditions (FREE SLIP)
//	P=0, d2/dr2 P = 0, d(T/r)/dr = 0  => only T computed at boundary

/// POLOIDAL
// use toroidal matrix as temporary matrix...
	for(i=NG+1; i<NR-1; i++) {
					MUt[i-NG].l =    Lr[i].l;
		for (l=0; l<=LMAX; l++)	MUt[i-NG].d[l] = (Lr[i].d - r_2[i]*l*(l+1));
					MUt[i-NG].u =    Lr[i].u;
	}
	i = NR-1;		// BC poloidal
		if (BC == BC_FREE_SLIP) {	// free-slip
			dx_2 = -2.0/(r[i]*(r[i]-r[i-1]));	// -2/(r.dr)
		} else {	// default to no-slip
			dx_2 = 1.0/(r[i]-r[i-1]);	dx_2 = 2.0*dx_2*dx_2;	// 2/(dr.dr)
		}
					MUt[i-NG].l =    dx_2;
		for (l=0; l<=LMAX; l++) MUt[i-NG].d[l] = ( -dx_2 - r_2[i]*l*(l+1) );
					MUt[i-NG].u =    0.0;
	i = NG;		// dx_1 for l=1,  dx_2 for l=2
		if (r[i] == 0.0) {	// no inner core, use central condition
			dx_2 = 1.0/(r[i+1]-r[i]);	dx_2 = 2.0*dx_2*dx_2;	// 2/(dr.dr)
			dx_1 = 2.0/(r[i]*(r[i+1]-r[i]));	// 2/(r.dr)
		} else if (BC == BC_FREE_SLIP) {	// free-slip
			dx_2 = 2.0/(r[i]*(r[i+1]-r[i]));	// 2/(r.dr)
			dx_1 = dx_2;	// same for l=1
		} else {	// default to no-slip
			dx_2 = 1.0/(r[i+1]-r[i]);	dx_2 = 2.0*dx_2*dx_2;	// 2/(dr.dr)
			dx_1 = dx_2;	// same for l=1
		}
		for (l=0; l<=LMAX; l++) {	// Lp0d and Lp0u are the coeficient for i = NG
			Lp0d[l] = ( -dx_2 - r_2[i]*l*(l+1) );
			Lp0u[l] = dx_2;
		}
		l=1;
			Lp0d[l] = ( -dx_1 - r_2[i]*l*(l+1) );
			Lp0u[l] = dx_1;

	i=1;
		for (l=0; l<=LMAX; l++) {
			MUp[i][l].l2 = 0.0;
			MUp[i][l].l1 =           0.5*nu*( MUt[i].l*Lp0d[l] + MUt[i].d[l]*MUt[i].l );
			MUp[i][l].d  = 1.0/dtU + 0.5*nu*( MUt[i].l*Lp0u[l] + MUt[i].d[l]*MUt[i].d[l] + MUt[i].u*MUt[i+1].l );
			MUp[i][l].u1 =           0.5*nu*(                    MUt[i].d[l]*MUt[i].u    + MUt[i].u*MUt[i+1].d[l] );
			MUp[i][l].u2 =           0.5*nu*(                                              MUt[i].u*MUt[i+1].u );
		}
	for (i=2; i<NU-1; i++) {
		for (l=0; l<=LMAX; l++) {
			MUp[i][l].l2 =           0.5*nu*( MUt[i].l*MUt[i-1].l );
			MUp[i][l].l1 =           0.5*nu*( MUt[i].l*MUt[i-1].d[l] + MUt[i].d[l]*MUt[i].l );
			MUp[i][l].d  = 1.0/dtU + 0.5*nu*( MUt[i].l*MUt[i-1].u    + MUt[i].d[l]*MUt[i].d[l] + MUt[i].u*MUt[i+1].l );
			MUp[i][l].u1 =           0.5*nu*(                          MUt[i].d[l]*MUt[i].u    + MUt[i].u*MUt[i+1].d[l] );
			MUp[i][l].u2 =           0.5*nu*(                                                    MUt[i].u*MUt[i+1].u );
		}
	}
	for (i=1; i<NU-1; i++) {
		for (l=0; l<=LMAX; l++) {
			MUp_1[i][l].l2 =         - MUp[i][l].l2;
			MUp_1[i][l].l1 =         - MUp[i][l].l1;
			MUp_1[i][l].d  = 2.0/dtU - MUp[i][l].d;
			MUp_1[i][l].u1 =         - MUp[i][l].u1;
			MUp_1[i][l].u2 =         - MUp[i][l].u2;
		}
	}
	PentaDec(MUp_1, 1, NU-2);

/// TOROIDAL
	for(i=NG+1; i<NR-1; i++) {
					MUt[i-NG].l =              0.5*nu*Lr[i].l;
		for (l=0; l<=LMAX; l++)	MUt[i-NG].d[l] = 1.0/dtU + 0.5*nu*(Lr[i].d - r_2[i]*l*(l+1));
					MUt[i-NG].u =              0.5*nu*Lr[i].u;
	}
	i = NR-1;		// BC toroidal (free-slip) : d(T/r)/dr = 0 =>  dT/dr = T/r
		dx_1 = 1.0/(r[i]-r[i-1]);	dx_2 = dx_1*dx_1;
					MUt[i-NG].l =              nu * dx_2;
		for (l=0; l<=LMAX; l++) MUt[i-NG].d[l] = 1.0/dtU + nu*( -dx_2 + r_1[i]*dx_1 + (1.0 - 0.5*l*(l+1))*r_2[i] );
					MUt[i-NG].u = 0.0;

	for(i=1; i<NU; i++) {
					MUt_1[i].l =            - MUt[i].l;
		for (l=0; l<=LMAX; l++) MUt_1[i].d[l] = 2.0/dtU - MUt[i].d[l];
					MUt_1[i].u =            - MUt[i].u;
	}
	TriDec(MUt_1, 1, NU-1);		// for no-slip : 1 to NU-2, for free-slip : 1 to NU-1
// for NO-SLIP, use : cTriSolve(MUt_1, RHS, Tlm, 1, NU-2)
// for FREE-SLIP, use : cTriSolve(MUt_1, RHS, Tlm, 1, NU-1)
}

void step_dyncin(long int nstep)
{
	long int i,l;

	while(nstep > 0) {
		nstep--;
		PolTor_to_spat(BPlm, BTlm, Br, Bt, Bp, BC_MAGNETIC);
		induction(NLbp1, NLbt1);
		cTriMul(MB, BPlm, Btmp, 1, NR-1);
		for (i=0; i<NR; i++) {
			for (l=0;l<NLM;l++) {
				Btmp[i][l] += 1.5*NLbp1[i][l] - 0.5*NLbp2[i][l];
			}
		}
		cTriSolve(MB_1, Btmp, BPlm, 1, NR-1);

		cTriMul(MB, BTlm, Btmp, 1, NR-2);
		for (i=0; i<NR-1; i++) {
			for (l=0;l<NLM;l++) {
				Btmp[i][l] += 1.5*NLbt1[i][l] - 0.5*NLbt2[i][l];
			}
		}
		cTriSolve(MB_1, Btmp, BTlm, 1, NR-2);

		PolTor_to_spat(BPlm, BTlm, Br, Bt, Bp, BC_MAGNETIC);
		induction(NLbp2, NLbt2);
		cTriMul(MB, BPlm, Btmp, 1, NR-1);
		for (i=0; i<NR; i++) {
			for (l=0;l<NLM;l++) {
				Btmp[i][l] += 1.5*NLbp2[i][l] - 0.5*NLbp1[i][l];
			}
		}
		cTriSolve(MB_1, Btmp, BPlm, 1, NR-1);
		cTriMul(MB, BTlm, Btmp, 1, NR-2);
		for (i=0; i<NR-1; i++) {
			for (l=0;l<NLM;l++) {
				Btmp[i][l] += 1.5*NLbt2[i][l] - 0.5*NLbt1[i][l];
			}
		}
		cTriSolve(MB_1, Btmp, BTlm, 1, NR-2);
	}
}

void step_NS(long int nstep)
{
	long int i,l;

	while(nstep > 0) {
		
	}
}

int main (int argc, char *argv[])
{
	double t0,t1,Rm,Rm2,z;
	int i,im,m,l,jj, it, lmtest;
	FILE *fp;
	complex double *Btmp[NR];

	init_SH();
	init_rad_sph(0.0, 1.0);
	init_Bmatrix();

	alloc_Bfields();	alloc_Ufields();
	for (i=0;i<NR;i++)
		Btmp[i] = (complex double *) malloc( NLM * sizeof(complex double));

	sscanf(argv[1],"%lf",&Rm);
	printf("=> Rm = %f\n",Rm);
	// Dudley James : Rmcrit around 54 ... for m=1 // I find Rmcrit = 40 for m=2...
	
	lmtest = LM(1,1);
	for (i=0; i<NR; i++) {
		for (l=0;l<NLM;l++) {
			BPlm[i][l] = 0.0;
			BTlm[i][l] = 0.0;
			UTlm[i][l] = 0.0;
			UPlm[i][l] = 0.0;
		}
//		BPlm[i][1] = r[i];	// dipole
		BTlm[i][lmtest] = r[i]*(1-r[i]);
//		BPlm[i][LM(2,2)] = r[i];
		BPlm[i][lmtest] = r[i];

// WARNING : compared to the article of Dudley & James, the Pol and Tor functions must be divided by r.
// simple roll flows (eq 24)
//		UTlm[i][LM(2,0)] = Rm*r[i]*sin(pi*r[i]);
// 		UPlm[i][LM(2,0)] = 0.14*UTlm[i][LM(2,0)];

// gubins (from Dudley-James)
		UTlm[i][LM(2,0)] = -r[i]*sin(2*pi*r[i]) * tanh(2*pi*(1-r[i]));
		UPlm[i][LM(2,0)] = 0.1* UTlm[i][LM(2,0)];
/*// pekeris (") en j2  (eq 20-21)
		if (r[i] != 0.0) {
			z =  5.7634591968447*r[i];
			UPlm[i][LM(2,2)] = Rm*5.7634591968447 * ((3.0/(z*z*z) -1.0/z)*sin(z) - 3./(z*z) * cos(z));
			UTlm[i][LM(2,2)] = 5.7634591968447 * UPlm[i][LM(2,2)];
		}	*/
//		printf("%f ", UPlm[i][LM(2,2)]);
	}

// representation spatiale de U
// maximum de vitesse :
	PolTor_to_spat(UPlm, UTlm, Ur, Ut, Up, BC_FREE_SLIP);
	t0 = 0.0;
	for(i=0;i<NR; i++) {
		for(im=0;im<NPHI;im++) {
			for(it=0;it<NLAT;it++) {
				t1 = Ur[i][im*NLAT+it]*Ur[i][im*NLAT+it] + Ut[i][im*NLAT+it]*Ut[i][im*NLAT+it] + Up[i][im*NLAT+it]*Up[i][im*NLAT+it];
				if (t1 > t0) t0 = t1;
			}
		}
	}
	Rm2 = sqrt(t0);
	printf(":: max speed = %f, rescaled to Rm = %f\n",Rm2, Rm);
	for (i=0; i<NR; i++) {
		for (l=0;l<NLM;l++) {
			UTlm[i][l] *= Rm/Rm2;
			UPlm[i][l] *= Rm/Rm2;
		}
	}
	PolTor_to_spat(UPlm, UTlm, Ur, Ut, Up, BC_FREE_SLIP);
	write_slice("Ur",Ur,0);
	write_slice("Ut",Ut,0);
	write_slice("Up",Up,0);


	fp = fopen("prof0","w");
	for (i=0; i<NR; i++) {
		fprintf(fp,"%.6g ",BPlm[i][lmtest]);
	}
	fclose(fp);

	write_HS("Bpol0",BPlm);	write_HS("Btor0",BTlm);

	PolTor_to_spat(BPlm, BTlm, Br, Bt, Bp, BC_MAGNETIC);
	induction(NLbp2, NLbt2);

	for (it=0; it< 1000; it++) {
		t0 = t1;
//		t1 = Energy(BPlm, BTlm);
		t1 = BPlm[NR/2][lmtest];
		printf("%g :: tx=%g\n",t1,log(t1/t0)/(2*dtB));

		step_dyncin(1);
	}

	write_HS("Bpol",BPlm);	write_HS("Btor",BTlm);
	write_HS("Upol",UPlm);	write_HS("Utor",UTlm);

	PolTor_to_spat(BPlm, BTlm, Br, Bt, Bp, BC_MAGNETIC);
	write_slice("Br",Br,0);
	write_slice("Bt",Bt,0);
	write_slice("Bp",Bp,0);

	fp = fopen("Pprof","w");
	for (i=0; i<NR; i++) {
		fprintf(fp,"%.6g ",BPlm[i][LM(1,1)]);
	}
	fclose(fp);
	fp = fopen("Tprof","w");
	for (i=0; i<NR; i++) {
		fprintf(fp,"%.6g ",BTlm[i][LM(1,1)]);
	}
	fclose(fp);

	printf("=> Rm = %f (based on max speed)\n",Rm);

}

