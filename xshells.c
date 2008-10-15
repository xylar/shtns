  /////////////////////////////////////////////////////////////////////////
 // XSHELLS : eXtendable Spherical Harmonic Earth-Like Liquid Simulator //
/////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
// FFTW : spatial derivative is d/dx = ik	(no minus sign !)
#include <fftw3.h>
// GSL for Legendre functions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>


#include "SHT.c"

// number of radial grid points.
long int NR,NU;		//  NR: total radial grid points. NU:for velocity field.
long int NG=0;		//  NG: grid points for inner core.

#include "grid.c"

// boundary condition constants.
#define BC_NONE 0
#define BC_NO_SLIP 1
#define BC_FREE_SLIP 2
#define BC_MAGNETIC 3

struct TriDiagL *MB, *MB_1, *MUt, *MUt_1;
struct PentaDiag **MUp, **MUp_1;

double nu, eta;		// viscosity and magnetic diffusivity.
double dtU, dtB;	// time step for navier-stokes and induction equation.
double Omega0;		// global rotation rate (of outer boundary) => Coriolis force .
double DeltaOmega;	// differential rotation (of inner core)

//#define DEB printf("%s:%u pass\n", __FILE__, __LINE__)
#define DEB (0)


struct VectField {
	double **r,**t,**p;
};

struct PolTor {
	complex double **P,**T;
};

struct VectField B, U, W, J, B0;
struct PolTor Blm, Ulm, NLb1, NLb2, NLu1, NLu2;

#include "xshells_io.c"

/// Allocate memory for a dynamic vector field (like U or B) and its pol/tor representation + non-linear term storage.
void alloc_DynamicField(struct PolTor *PT, struct VectField *V, struct PolTor *NL1, struct PolTor *NL2, long int istart, long int iend)
{
	long int ir;

	// alloc radial structure.
	V->r = (double **) malloc( 3*NR * sizeof(double *) );
	V->t = V->r + NR;	V->p = V->r + 2*NR;
	PT->P = (complex double **) malloc( 6*NR * sizeof(complex double *) );
	PT->T = PT->P + NR;
	NL1->P = PT->P + 2*NR;	NL1->T = PT->P + 3*NR;
	NL2->P = PT->P + 4*NR;	NL2->T = PT->P + 5*NR;

	// not defined here, set as NULL pointer.
	for (ir = 0; ir < istart; ir++) {
		PT->P[ir] = NULL;	PT->T[ir] = NULL;
		NL1->P[ir] = NULL;	NL1->T[ir] = NULL;
		NL2->P[ir] = NULL;	NL2->T[ir] = NULL;
		V->r[ir] = NULL;	V->t[ir] = NULL;	V->p[ir] = NULL;
	}
	for (ir = istart; ir <= iend; ir++) {		// shell by shell allocation.
		V->r[ir] = (double *) fftw_malloc( 3*(NPHI/2+1)*NLAT * sizeof(complex double));	// 3 components
		V->t[ir] = V->r[ir] + 2*(NPHI/2+1)*NLAT;
		V->p[ir] = V->t[ir] + 2*(NPHI/2+1)*NLAT;
		PT->P[ir] = (complex double *) malloc( 6*NLM * sizeof(complex double));	// Pol/Tor
		PT->T[ir] = PT->P[ir] + NLM;
		NL1->P[ir] = PT->P[ir] + 2*NLM;
		NL1->T[ir] = PT->P[ir] + 3*NLM;
		NL2->P[ir] = PT->P[ir] + 4*NLM;
		NL2->T[ir] = PT->P[ir] + 5*NLM;
	}
	for (ir = iend+1; ir<NR; ir++) {
		PT->P[ir] = NULL;	PT->T[ir] = NULL;
		NL1->P[ir] = NULL;	NL1->T[ir] = NULL;
		NL2->P[ir] = NULL;	NL2->T[ir] = NULL;
		V->r[ir] = NULL;	V->t[ir] = NULL;	V->p[ir] = NULL;
	}
}

/// Allocate a vector field (not evolved) like J, W, B0, ...
void alloc_VectField(struct VectField *V, long int istart, long int iend)
{
	long int ir;

	// alloc radial structure.
	V->r = (double **) malloc( 3*NR * sizeof(double *) );
	V->t = V->r + NR;	V->p = V->r + 2*NR;
	// not defined here, set as NULL pointer.
	for (ir = 0; ir < istart; ir++) {
		V->r[ir] = NULL;	V->t[ir] = NULL;	V->p[ir] = NULL;
	}
	for (ir = istart; ir <= iend; ir++) {		// shell by shell allocation.
		V->r[ir] = (double *) fftw_malloc( 3*(NPHI/2+1)*NLAT * sizeof(complex double));	// 3 components
		V->t[ir] = V->r[ir] + 2*(NPHI/2+1)*NLAT;
		V->p[ir] = V->t[ir] + 2*(NPHI/2+1)*NLAT;
	}
	for (ir = iend+1; ir<NR; ir++) {
		V->r[ir] = NULL;	V->t[ir] = NULL;	V->p[ir] = NULL;
	}
}



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
	MUt = (struct TriDiagL *) malloc( 2*NU* sizeof(struct TriDiagL));
	MUt_1 = MUt + NU;
// allocate poloidal matrix
	MUp = (struct PentaDiag **) malloc( 2*NU * sizeof(struct PentaDiag *) );
	MUp_1 = MUp + NU;
	MUp[0] = (struct PentaDiag *) malloc( 2*NU*(LMAX+1)* sizeof(struct PentaDiag));
	MUp_1[0] = MUp[0] + NU*(LMAX+1);
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
//	with DeltaOmega (couette), at inner core : T(l=1,m=0) = r.DeltaOmega*Y10_ct
// Boundary conditions (FREE SLIP)
//	P=0, d2/dr2 P = 0, d(T/r)/dr = 0  => only T computed at boundary

/// POLOIDAL
/// (d/dt - nu.Lap)(-Lap.Up) = Toroidal[rot_NL]
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
			MUp[i][l].l1 = -( 1.0/dtU*MUt[i].l    + 0.5*nu*( MUt[i].l*Lp0d[l] + MUt[i].d[l]*MUt[i].l ) );
			MUp[i][l].d  = -( 1.0/dtU*MUt[i].d[l] + 0.5*nu*( MUt[i].l*Lp0u[l] + MUt[i].d[l]*MUt[i].d[l] + MUt[i].u*MUt[i+1].l ) );
			MUp[i][l].u1 = -( 1.0/dtU*MUt[i].u    + 0.5*nu*(                    MUt[i].d[l]*MUt[i].u    + MUt[i].u*MUt[i+1].d[l] ) );
			MUp[i][l].u2 = -(                       0.5*nu*(                                              MUt[i].u*MUt[i+1].u ) );
		}
	for (i=2; i<NU-1; i++) {
		for (l=0; l<=LMAX; l++) {
			MUp[i][l].l2 = -(                       0.5*nu*( MUt[i].l*MUt[i-1].l ) );
			MUp[i][l].l1 = -( 1.0/dtU*MUt[i].l    + 0.5*nu*( MUt[i].l*MUt[i-1].d[l] + MUt[i].d[l]*MUt[i].l ) );
			MUp[i][l].d  = -( 1.0/dtU*MUt[i].d[l] + 0.5*nu*( MUt[i].l*MUt[i-1].u    + MUt[i].d[l]*MUt[i].d[l] + MUt[i].u*MUt[i+1].l ) );
			MUp[i][l].u1 = -( 1.0/dtU*MUt[i].u    + 0.5*nu*(                          MUt[i].d[l]*MUt[i].u    + MUt[i].u*MUt[i+1].d[l] ) );
			MUp[i][l].u2 = -(                       0.5*nu*(                                                    MUt[i].u*MUt[i+1].u ) );
		}
	}

	for (i=1; i<NU-1; i++) {
		for (l=0; l<=LMAX; l++) {
			MUp_1[i][l].l2 =                     - MUp[i][l].l2;
			MUp_1[i][l].l1 = -2.0/dtU*MUt[i].l    - MUp[i][l].l1;
			MUp_1[i][l].d  = -2.0/dtU*MUt[i].d[l] - MUp[i][l].d;
			MUp_1[i][l].u1 = -2.0/dtU*MUt[i].u    - MUp[i][l].u1;
			MUp_1[i][l].u2 =                     - MUp[i][l].u2;
		}
	}
	PentaDec(MUp_1, 1, NU-2);

/// TOROIDAL
/// (d/dt - nu.Lap) Ut = Poloidal(rot_NL)
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



/// compute the cartesian coordinates of a vector field at r=0, from its Poloidal value component
inline void Pol_to_cart_spat0(complex double **Plm, double *Vx, double *Vy, double *Vz)
{
	double dx;
	long int lm;

	dx = 2.0 / (r[1] * Y10_ct);
	lm = LM(1,1);	// l=1, m=1
	*Vx = dx * creal( Plm[1][lm] );
	*Vy = dx * cimag( Plm[1][lm] );
	lm = LM(1,0);	// l=1, m=0
	*Vz = dx * (double) Plm[1][lm];
}

/// compute the value of the spheroidal scalar at r=0 from the cartesian representation of a vector field at the center.
inline void cart_spat0_to_Sph(double Vx, double Vy, double Vz, complex double *S10, complex double *S11)
{
	double dx = Y10_ct;

	*S10 = Vz * dx;		// l=1, m=0
	*S11 = (Vx + I*Vy) *dx;	// l=1, m=1
}

#define CALC_U(bc)   PolTor_to_spat(&Ulm, &U, NG, NR-1, (bc))
#define CALC_B       PolTor_to_spat(&Blm, &B, 0, NR-1, BC_MAGNETIC)

/// compute spatial field from Poloidal/Toroidal representation.
void PolTor_to_spat(struct PolTor *PT, struct VectField *V, long int istart, long int iend, long int BC)
{
	complex double Q[NLM];	// l(l+1) * P/r
	complex double S[NLM];	// dP/dr + P/r = 1/r.d(rP)/dr = Wr(P)
	double dr;
	long int ir,lm,l;

	ir = istart;
	if ((istart == 0)&&(r[ir] == 0.0)) {	// field at r=0 : S = 2.dP/dr (l=1 seulement) => store Bx, By, Bz (cartesian coordinates)
		Pol_to_cart_spat0(PT->P, V->r[ir],V->t[ir],V->p[ir]);
		for (lm=1; lm<NLAT*NPHI; lm++) {	// zero for the rest.
			V->r[ir][lm] = 0.0;	V->t[ir][lm] = 0.0;	V->p[ir][lm] = 0.0;
		}
	} else if (BC == BC_NO_SLIP) {
		// Velocity field (no slip => Ur=0, Ut=0, Up=r.sint.DeltaOmega)
		for (lm=0; lm<NLAT*NPHI; lm++) {
			V->r[ir][lm] = 0.0;	V->t[ir][lm] = 0.0;
		}
		for (lm=0; lm<NPHI; lm++) {
			for (l=0; l<NLAT; l++) V->p[ir][lm*NLAT+l] = r[ir]*st[l]*DeltaOmega;
		}
	} else if (BC == BC_FREE_SLIP) {
		// Velocity field (free slip BC) : P = 0, d(T/r)/dr = 0   => Ur=0
		dr = r[ir+1]/(r[ir]*(r[ir+1]-r[ir]));
		for (lm=0; lm<NLM; lm++) {
			S[lm] = dr * PT->P[ir+1][lm];
		}
		SHsphtor_to_spat(S, PT->T[ir], (complex double *) V->t[ir], (complex double *) V->p[ir]);
		for (lm=0; lm<NLAT*NPHI; lm++) {
			V->r[ir][lm] = 0.0;
		}
	}

	for (ir=istart+1; ir<iend; ir++) {
		for (lm=0; lm<NLM; lm++) {		// Solenoidal deduced from radial derivative of Poloidal
			S[lm] = Wr[ir].l*PT->P[ir-1][lm] + Wr[ir].d*PT->P[ir][lm] + Wr[ir].u*PT->P[ir+1][lm];
			Q[lm] = r_1[ir]*l2[lm] * PT->P[ir][lm];
		}
		SH_to_spat(Q,(complex double *) V->r[ir]);
		SHsphtor_to_spat(S, PT->T[ir], (complex double *) V->t[ir], (complex double *) V->p[ir]);
	}

	ir = iend;
	if (BC == BC_NO_SLIP) {
		// Velocity field (no slip => U = 0)
		for (lm=0; lm<NLAT*NPHI; lm++) {
			V->r[ir][lm] = 0.0;	V->t[ir][lm] = 0.0;	V->p[ir][lm] = 0.0;
		}
		return;
	} else if (BC == BC_FREE_SLIP) {
		// Velocity field (free slip BC) : P = 0, d(T/r)/dr = 0   => Ur=0
		dr = r[ir-1]/(r[ir]*(r[ir-1]-r[ir]));
		for (lm=0; lm<NLM; lm++) {
			S[lm] = dr * PT->P[ir-1][lm];
		}
		SHsphtor_to_spat(S, PT->T[ir], (complex double *) V->t[ir], (complex double *) V->p[ir]);
		for (lm=0; lm<NLAT*NPHI; lm++) {
			V->r[ir][lm] = 0.0;
		}
		return;
	} else if (BC == BC_MAGNETIC) {
		// Magnetic field (insulator BC) : T=0, dP/dr = -(l+1)/r.P => S = -lP/r
		for (lm=0, l=0; lm<NLM; lm++) {
			S[lm] = -l*PT->P[ir][lm]*r_1[ir];
			if (l < LMAX) { l++; } else { l=0; }
			Q[lm] = r_1[ir]*l2[lm] * PT->P[ir][lm];
		}
		SH_to_spat(Q,(complex double *) V->r[ir]);
		SHsph_to_spat(S, (complex double *) V->t[ir], (complex double *) V->p[ir]);
	}
}


/// Compute W = rot(V) + Wz0.ez  from PolTor components.
/// IN: PT  : pol/tor components of vector V
///     Wz0 : background value along ez
/// OUT: W : r,theta,phi components of rot(V) + Wz0.ez
void PolTor_to_rot_spat(struct PolTor *PT, double Wz0, struct VectField *W, long int istart, long int iend, long int BC)
// Pol = Tor;  Tor = Lap.Pol
{
	complex double Q[NLM];	// l(l+1) * T/r
	complex double S[NLM];	// dT/dr + T/r
	complex double T[NLM];	//  [ l(l+1)/r^2 - 1/r d2/dr2(r .) ] P
	double dr;
	long int ir,lm;

	Wz0 *= Y10_ct;		// multiply by representation of cos(theta) in spherical harmonics (l=1,m=0)

	if (istart == 0) {
		ir = istart;
		if (r[ir] = 0.0) {	// field at r=0 : S = 2.dT/dr (l=1 seulement) => store Bx, By, Bz (cartesian coordinates)
			Pol_to_cart_spat0(PT->T, W->r[ir],W->t[ir],W->p[ir]);
			for (lm=1; lm<NLAT*NPHI; lm++) {	// zero for the rest.
				W->r[ir][lm] = 0.0;	W->t[ir][lm] = 0.0;	W->p[ir][lm] = 0.0;
			}
			istart++;
		} else if (BC == BC_NO_SLIP) {
			for (lm=0; lm<NLM; lm++) {		// Solenoidal deduced from radial derivative of Toroidal
				Q[lm] = r_1[ir]*l2[lm] * PT->T[ir][lm];
				S[lm] = Wr[ir].l*PT->T[ir-1][lm] + Wr[ir].d*PT->T[ir][lm] + Wr[ir].u*PT->T[ir+1][lm];
				T[lm] = (r_2[ir]*l2[lm] - r_1[ir]*D2r[ir].d)*PT->P[ir][lm] - r_1[ir]*D2r[ir].l*PT->P[ir-1][lm] - r_1[ir]*D2r[ir].u*PT->P[ir+1][lm];
			}
			istart++;
		}
	}
	for (ir=istart; ir < iend; ir++) {
		for (lm=0; lm<NLM; lm++) {		// Solenoidal deduced from radial derivative of Poloidal
			Q[lm] = r_1[ir]*l2[lm] * PT->T[ir][lm];
			S[lm] = Wr[ir].l*PT->T[ir-1][lm] + Wr[ir].d*PT->T[ir][lm] + Wr[ir].u*PT->T[ir+1][lm];
			T[lm] = (r_2[ir]*l2[lm] - r_1[ir]*D2r[ir].d)*PT->P[ir][lm] - r_1[ir]*D2r[ir].l*PT->P[ir-1][lm] - r_1[ir]*D2r[ir].u*PT->P[ir+1][lm];
		}
/* Pour Coriolis : ez ^ u
	ez = cos(theta).er - sin(theta).etheta
	cos theta = Y(m=0,l=1) * 2*sqrt(pi/3)		>> peut etre rajouté à Qlm. (=> Vr)
	-sin theta = dY(m=0,l=1)/dt * 2*sqrt(pi/3)	>> peut etre rajouté à Slm  (=> Vt)
*/
		// add Background mode (l=1, m=0) [ie Coriolis force (2*Omega0) or Magnetic field (B0z)]
		(double) Q[1] += Wz0;
		(double) S[1] += Wz0;
		SH_to_spat(Q,(complex double *) W->r[ir]);
		SHsphtor_to_spat(S, T, (complex double *) W->t[ir], (complex double *) W->p[ir]);
	}

	ir = iend;
}

/// spatial to curl : only for ir = 1 .. NR-2
///	Pol <- Tor
///	Tor <- Q/r - 1/r.d(rS)/dr
void spat_to_rot_PolTor(struct VectField *V, complex double **Plm, complex double **Tlm, long int istart, long int iend)
{
	complex double Q[NLM];		// Q
	complex double S[NLM*3];	// buffers for S.
	complex double *St, *Sl, *Sd, *Su;	// pointers to S (temp, lower, diag, upper)
	long int ir,lm;

	Sl = S;   Sd = S + NLM;   Su = S + 2*NLM;	// init pointers to buffer.

	ir = istart;
		spat_to_SHsphtor((complex double *) V->t[ir], (complex double *) V->p[ir], Sd, Plm[ir]);
		spat_to_SHsphtor((complex double *) V->t[ir+1], (complex double *) V->p[ir+1], Su, Plm[ir+1]);
		spat_to_SH((complex double *) V->r[ir], Q);
		for (lm=0; lm<NLM; lm++) {
			Tlm[ir][lm] = r_1[ir]*Q[lm] - (Wr[ir].d * Sd[lm] + Wr[ir].u * Su[lm]);		// Sl = 0
		}
	for (ir=istart+1; ir < iend; ir++) {
		St = Sl;	Sl = Sd;	Sd = Su;	Su = St;		// rotate buffers.
		spat_to_SHsphtor((complex double *) V->t[ir+1], (complex double *) V->p[ir+1], Su, Plm[ir+1]);
		spat_to_SH((complex double *) V->r[ir], Q);
		for (lm=0; lm<NLM; lm++) {
			Tlm[ir][lm] = r_1[ir]*Q[lm] - (Wr[ir].l * Sl[lm] + Wr[ir].d * Sd[lm] + Wr[ir].u * Su[lm]);
		}
	}
	ir = iend;
		St = Sl;	Sl = Sd;	Sd = Su;	Su = St;		// rotate buffers.
		spat_to_SH((complex double *) V->r[ir], Q);
		for (lm=0; lm<NLM; lm++) {
			Tlm[ir][lm] = r_1[ir]*Q[lm] - (Wr[ir].l * Sl[lm] + Wr[ir].d * Sd[lm]);		// Su=0
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


/// compute vorticity field from poloidal/toroidal components of velocity [ie poltor_to_rot_spat applied to U]
/// IN: PT : pol/tor components of velocity
///     Om0 : background global rotation (solid body rotation rate) to include in vorticity
/// OUT: W : r,theta,phi components of vorticity field
void calc_Vort(struct PolTor *PT, double Om0, struct VectField *W)
{
	complex double Q[NLM];
	complex double S[NLM];
	complex double T[NLM];
	long int ir,lm;

/* Pour Coriolis : ez ^ u
	ez = cos(theta).er - sin(theta).etheta
	cos theta = Y(m=0,l=1) * 2*sqrt(pi/3)		>> peut etre rajouté à Qlm. (=> Vr)
	-sin theta = dY(m=0,l=1)/dt * 2*sqrt(pi/3)	>> peut etre rajouté à Slm  (=> Vt)
*/
	Om0 = 2.0*Om0 * Y10_ct;	// multiply by representation of cos(theta) in spherical harmonics (l=1,m=0)
//	Plm -= NG;	Tlm -= NG;	Vr -= NG;	Vt -= NG;	Vp -= NG;	// adjust pointers.
	for (ir=NG+1; ir <= NR-2; ir++) {
		for (lm=0; lm<NLM; lm++) {
			T[lm] = (r_2[ir]*l2[lm] - r_1[ir]*D2r[ir].d)*PT->P[ir][lm] - r_1[ir]*D2r[ir].l*PT->P[ir-1][lm] - r_1[ir]*D2r[ir].u*PT->P[ir+1][lm];
			Q[lm] = r_1[ir]*l2[lm] * PT->T[ir][lm];
			S[lm] = Wr[ir].l*PT->T[ir-1][lm] + Wr[ir].d*PT->T[ir][lm] + Wr[ir].u*PT->T[ir+1][lm];
		}
		(double) Q[1] += Om0;		// add Background Vorticity for Coriolis Force (l=1, m=0)
		(double) S[1] += Om0;		// add Background Vorticity for Coriolis Force (l=1, m=0)
		SH_to_spat(Q,(complex double *) W->r[ir]);
		SHsphtor_to_spat(S, T, (complex double *) W->t[ir], (complex double *) W->p[ir]);
	}
}


// *****************************
// ***** NON-LINEAR TERMS ******
// *****************************

/// compute rot(VxW + JxB) and its pol-tor components.
/// output : V, B, J unchanged, W = VxW + JxB
///          NL = PolTor[rot(VxW + JxB)]
void calc_NL_MHD(struct VectField *V, struct VectField *W, struct VectField *B, struct VectField *J, struct PolTor *NL)
{
	double vr,vt,vp;
	long int ir,lm;

	for (ir=NG+1; ir<=NR-2; ir++) {
		for (lm=0; lm<NPHI*NLAT; lm++) {	// catastrophic memory accesses... 4*3 = 12 disjoint arrays.
			vr =  J->t[ir][lm]*B->p[ir][lm] - J->p[ir][lm]*B->t[ir][lm];	// JxB
			vt =  J->p[ir][lm]*B->r[ir][lm] - J->r[ir][lm]*B->p[ir][lm];
			vp =  J->r[ir][lm]*B->t[ir][lm] - J->t[ir][lm]*B->r[ir][lm];

			vr += V->t[ir][lm]*W->p[ir][lm] - V->p[ir][lm]*W->t[ir][lm];	// VxW
			vt += V->p[ir][lm]*W->r[ir][lm] - V->r[ir][lm]*W->p[ir][lm];
			vp += V->r[ir][lm]*W->t[ir][lm] - V->t[ir][lm]*W->r[ir][lm];
			W->r[ir][lm] = vr;	W->t[ir][lm] = vt;	W->p[ir][lm] = vp;
		}
	}
	spat_to_rot_PolTor(W, NL->T, NL->P, NG+1, NR-2);
}

/// compute rot(VxW) and its pol-tor components.
/// output : V unchanged, W = VxW, NL = PolTor[rot(VxW)]
/// U is kept unchanged
void calc_NL_NS(struct VectField *V, struct VectField *W, struct PolTor *NL)
{
	double vr,vt,vp;
	long int ir,lm;

	for (ir=NG+1; ir<=NR-2; ir++) {
		for (lm=0; lm<NPHI*NLAT; lm++) {	// catastrophic memory accesses... 2*3 = 6 disjoint arrays.
			vr = V->t[ir][lm]*W->p[ir][lm] - V->p[ir][lm]*W->t[ir][lm];
			vt = V->p[ir][lm]*W->r[ir][lm] - V->r[ir][lm]*W->p[ir][lm];
			vp = V->r[ir][lm]*W->t[ir][lm] - V->t[ir][lm]*W->r[ir][lm];
			W->r[ir][lm] = vr;	W->t[ir][lm] = vt;	W->p[ir][lm] = vp;
		}
	}
	spat_to_rot_PolTor(W, NL->T, NL->P, NG+1, NR-2);
}

/// compute rot(VxB)
/// IN:  V, B are spatial vector fields (unchanged)
/// OUT: VxB  is the VxB resulting vector field. (Note: VxB can be either of V or B or a third vector)
///      NL is the PolTor component of rot(VxB)
void induction(struct VectField *V, struct VectField *B, struct VectField *VxB, struct PolTor *NL)
{
	double vr,vt,vp,rO;
	long int ir,lm,it;

	for (ir=0; ir<NG; ir++) {		// for the inner core : solid body rotation  Up = r.sint.DeltaOmega
		rO = r[ir]*DeltaOmega;
		it = 0;
		for (lm=0; lm<NPHI*NLAT; lm++) {
			vp = rO*st[it]; it++; if (it>=NLAT) it=0;	// TO BE CHECKED !
			vr = - vp*B->t[ir][lm];
			vt =   vp*B->r[ir][lm];
			VxB->r[ir][lm] = vr;	VxB->t[ir][lm] = vt;	VxB->p[ir][lm] = 0.0;
		}
	}
	for (ir=NG; ir<NR; ir++) {
		for (lm=0; lm<NPHI*NLAT; lm++) {
			vr = V->t[ir][lm]*B->p[ir][lm] - V->p[ir][lm]*B->t[ir][lm];
			vt = V->p[ir][lm]*B->r[ir][lm] - V->r[ir][lm]*B->p[ir][lm];
			vp = V->r[ir][lm]*B->t[ir][lm] - V->t[ir][lm]*B->r[ir][lm];
			VxB->r[ir][lm] = vr;	VxB->t[ir][lm] = vt;	VxB->p[ir][lm] = vp;
		}
	}
/*	ir = NR-1;	// U = 0
		for (lm=0; lm<NPHI*NLAT; lm++) {
			Br[ir][lm] = 0.0;	Bt[ir][lm] = 0.0;	Bp[ir][lm] = 0.0;
		}*/
	spat_to_rot_PolTor(VxB, NL->P, NL->T, 1,NR-2);
}

// **************************
// ***** TIME STEPPING ******
// **************************

/// perform nstep steps of Navier-Stokes equation, with temporary array Alm.
double step_NS(long int nstep, complex double **Alm)
{
	inline step1(struct PolTor *NL, struct PolTor *NLo, complex double **Alm)
	{
		long int i,l;

		CALC_U(BC_NO_SLIP);
		calc_Vort(&Ulm, Omega0, &W);
		calc_NL_NS(&U, &W, NL);
//		calc_NL_MHD(&U, &W, &B, &J, NL);

		cTriMulBC(MUt, Ulm.T +NG, Alm, 1, NU-2);
		for (i=1; i<NU-1; i++) {
			for (l=1;l<NLM;l++) {
				Alm[i][l] += 1.5*NL->T[i+NG][l] - 0.5*NLo->T[i+NG][l];
			}
		}
		cTriSolveBC(MUt_1, Alm, Ulm.T +NG, 1, NU-2);

		cPentaMul(MUp, Ulm.P +NG, Alm, 1, NU-2);
		for (i=1; i<NU-1; i++) {
			for (l=1;l<NLM;l++) {
				Alm[i][l] += 1.5*NL->P[i+NG][l] - 0.5*NLo->P[i+NG][l];
			}
		}
		cPentaSolve(MUp_1, Alm, Ulm.P +NG, 1, NU-2);
	}

	while(nstep > 0) {
		nstep--;

		step1(&NLu1, &NLu2, Alm);
		step1(&NLu2, &NLu1, Alm);
	}

	return 2*nstep*dtU;
}

/// perform nstep steps of inducation equation, with temporary array Alm.
void step_Induction(long int nstep, complex double **Alm)
{
	inline step1(struct PolTor *NL, struct PolTor * NLo, complex double **Btmp)
	{
		long int i,l;

		induction(&U, &B0, &J, NL);
		cTriMul(MB, Blm.P, Btmp, 1, NR-1);
		for (i=0; i<NR; i++) {
			for (l=0;l<NLM;l++) {
				Btmp[i][l] += 1.5*NL->P[i][l] - 0.5*NLo->P[i][l];
			}
		}
		cTriSolve(MB_1, Btmp, Blm.P, 1, NR-1);

		cTriMul(MB, Blm.T, Btmp, 1, NR-2);
		for (i=0; i<NR-1; i++) {
			for (l=0;l<NLM;l++) {
				Btmp[i][l] += 1.5*NL->T[i][l] - 0.5*NLo->T[i][l];
			}
		}
		cTriSolve(MB_1, Btmp, Blm.T, 1, NR-2);
	}

	while(nstep > 0) {
		nstep--;

		step1(&NLb1, &NLb2, Alm);
		step1(&NLb2, &NLb1, Alm);
	}
}


// compute total Energy
double Energy(struct PolTor *PT)
{
	complex double pt;
	double E,er;
	long int ir,lm;

	E = 0.0;
	for (ir=0; ir<NR; ir++) {
		er = 0.0;
		for(lm=0; lm<NLM; lm++) {
			pt = PT->P[ir][lm];
			er += creal(pt)*creal(pt) + cimag(pt)*cimag(pt);
			pt = PT->T[ir][lm];
			er += creal(pt)*creal(pt) + cimag(pt)*cimag(pt);
		}
		E += er*dr[ir]*r[ir];
	}
	return E;
}


void read_Par(char *fname, char *job, long int *iter_max, long int *modulo, double *polar_opt_max)
{
	double tmp;
	int id;
	FILE *fp, *fpw;
	char str[256];
	char name[32];

	printf("[read_Par] reading parameters from '%s'...\n",fname);

	fp = fopen(fname,"r");
	if (fp == NULL) runerr("[read_Par] file not found !");

	fpw = fopen("tmp.par","w");	// keep the parameters for that job.

	// DEFAULT VALUES
	Omega0 = 1.0;		// global rotation is on by default.

	while (!feof(fp))
	{
		fgets(str, 256, fp);
		fputs(str, fpw);		// save file...
//		printf(str);
		if ((str[0] != '#' )&&(strlen(str) > 3)) {	// ignore les lignes commentees (#) et les lignes vides.
			sscanf(str,"%s = %lf",name,&tmp);
//			printf("[%s] = %f\n",name,tmp);
			// PHYSICAL PARAMS
			if (strcmp(name,"Omega0") == 0)		Omega0 = tmp;
			if (strcmp(name,"DeltaOmega") == 0)	DeltaOmega = tmp;
			if (strcmp(name,"nu") == 0)		nu = tmp;
			if (strcmp(name,"eta") == 0)		eta = tmp;
			// NUMERICAL SCHEME
			if (strcmp(name,"NR") == 0)		NR = tmp;
			if (strcmp(name,"NG") == 0)		NG = tmp;
			if (strcmp(name,"dtU") == 0)		dtU = tmp;
			if (strcmp(name,"dtB") == 0)		dtB = tmp;
			// TIME DOMAIN
			if (strcmp(name,"iter_max") == 0)	*iter_max = tmp;
			if (strcmp(name,"modulo") == 0)		*modulo = tmp;
			if (strcmp(name,"job") == 0) 		sscanf(str,"%s = %s",name,job);
			// ALGORITHM TUNING
			if (strcmp(name,"sht_polar_opt_max") == 0)	*polar_opt_max = tmp;
			if (strcmp(name,"fftw_plan_mode") == 0) {
				switch( (int) tmp) {
					case 0: fftw_plan_mode = FFTW_ESTIMATE; break;
					case 1: fftw_plan_mode = FFTW_MEASURE; break;
					case 2: fftw_plan_mode = FFTW_PATIENT; break;	// default
					case 3: fftw_plan_mode = FFTW_EXHAUSTIVE;
				}
			}
		}
	}
	// SOME BASIC COMPUTATIONS
	NU = NR-NG;

	fclose(fp);	fclose(fpw);
	sprintf(str, "%s.%s",fname,job);	rename("tmp.par", str);		// rename file.	
	fflush(stdout);		// when writing to a file.
}

int main (int argc, char *argv[])
{
	double t0,t1,Rm,Rm2,z;
	long int i,im,m,l,jj, it, lmtest;
	long int iter_max, modulo;
	complex double **Alm;		// temp scalar spectral (pol or tor) representation
	FILE *fp;
	char command[100] = "xshells.par";
	char job[40];
	double polar_opt_max = 0.0;	// default SHT optimization.

	printf("[XSHELLS] eXtendable Spherical Harmonic Earth-Like Liquid Simulator\n          by Nathanael Schaeffer / LGIT, build %s, %s\n",__DATE__,__TIME__);

	read_Par(command, job, &iter_max, &modulo, &polar_opt_max);

	init_SH(polar_opt_max);
	init_rad_sph(0.0, 1.0);

	init_Bmatrix();		init_Umatrix(BC_NO_SLIP);

	alloc_DynamicField(&Blm, &B, &NLb1, &NLb2, 0, NR-1);
	alloc_DynamicField(&Ulm, &U, &NLu1, &NLu2, NG, NR-1);
	alloc_VectField(&B0,0,NR-1);	alloc_VectField(&J,0,NR-1);	alloc_VectField(&W,NG,NR-1);

	printf("[Params] job name : %s -- NG=%d, rg=%f\n",job,NG,r[NG]);
	printf("         Ek=%.2e, Ro=%.2e, Pm=%.2e, Re=%.2e\n",nu/Omega0, r[NG]*DeltaOmega/Omega0, nu/eta, r[NG]*DeltaOmega/nu);
	printf("         dtU.Omega=%.2e, dtU.nu.R^2=%.2e, dtB.eta.R^2=%.2e, dtB/dtU=%.2e\n",dtU*Omega0, dtU*nu, dtB*eta, dtB/dtU);

	DEB;

	Alm = (complex double **) malloc( NR * sizeof(complex double *));
	for (i=0;i<NR;i++)
		Alm[i] = (complex double *) malloc( NLM * sizeof(complex double));

	DEB;
	// init B fields.
	for (i=0; i<NR; i++) {
		for (l=0;l<NLM;l++) {
			Blm.P[i][l] = 0.0;
			Blm.T[i][l] = 0.0;
		}
		Blm.P[i][1] = r[i];		// magnetic dipole
	}
	DEB;
	PolTor_to_spat(&Blm, &B0, 0, NR-1, BC_MAGNETIC);	// background magnetic field.
	for (i=0; i<NR; i++) {
		Blm.P[i][1] = 0.0;	// remove dipole.
	}

	DEB;
	// init U fields
	for (i=NG; i<NR; i++) {
		for (l=0;l<NLM;l++) {
			Ulm.T[i][l] = 0.0;
			Ulm.P[i][l] = 0.0;
		}
//		Ulm.T[i][LM(1,0)] = r[i]*DeltaOmega*(1-(r[i]-r[NG])/(r[NR-1]-r[NG])) * Y10_ct;
	}
	i = NG;
		Ulm.T[i][LM(1,0)]  = r[i]*DeltaOmega * Y10_ct;	// rotation differentielle de la graine.

	DEB;

	PolTor_to_spat(&Blm, &B, 0, NR-1, BC_MAGNETIC);
	induction(&U, &B, &B, &NLb2);

	DEB;

	CALC_U(BC_NO_SLIP);
	DEB;
	write_slice("Ur0",U.r,0);	write_slice("Ut0",U.t,0);	write_slice("Up0",U.p,0);
	DEB;
	calc_Vort(&Ulm, Omega0, &W);
	DEB;

	write_slice("Wr0",W.r,0);	write_slice("Wt0",W.t,0);	write_slice("Wp0",W.p,0);

	DEB;

	for (it=0; it< iter_max; it++) {
//		t0 = t1;
//		printf("%g :: tx=%g\n",t1,log(t1/t0)/(2*dtB));
		t0 = creal(Ulm.T[NG+1][1]);
		t1 = creal(Ulm.P[NG+1][2]);
		printf("[it %d] t=%g, P0=%g, T0=%g\n",it,2*it*modulo*dtU,t1, t0);
		if (isnan(t0)) runerr("NaN encountered");

		step_NS(modulo, Alm);
	}

	write_HS("Bpol",Blm.P);	write_HS("Btor",Blm.T);
	write_HS("Upol",Ulm.P);	write_HS("Utor",Ulm.T);

	CALC_U(BC_NO_SLIP);
	write_slice("Ur",U.r,0);	write_slice("Ut",U.t,0);	write_slice("Up",U.p,0);

	calc_Vort(&Ulm, 0.0, &W);		// without background vorticity
	write_slice("Wr",W.r,0);	write_slice("Wt",W.t,0);	write_slice("Wp",W.p,0);

	fp = fopen("Pprof","w");
	for (i=NG; i<NR; i++) {
		fprintf(fp,"%.6g ",Ulm.P[i][LM(2,0)]);
	}
	fclose(fp);
	fp = fopen("Tprof","w");
	for (i=NG; i<NR; i++) {
		fprintf(fp,"%.6g ",Ulm.T[i][LM(1,0)]);
	}
	fclose(fp);

}

