  /////////////////////////////////////////////////////////////////////////
 // XSHELLS : eXtendable Spherical Harmonic Earth-Like Liquid Simulator //
/////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <complex.h>
#include <math.h>
// FFTW : spatial derivative is d/dx = ik	(no minus sign !)
#include <fftw3.h>
// GSL for Legendre functions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>

#include "SHT.c"

double nu, eta;		// viscosity and magnetic diffusivity.
double dtU, dtB;	// time step for navier-stokes and induction equation.
double Omega0;		// global rotation rate (of outer boundary) => Coriolis force .
double DeltaOmega;	// differential rotation (of inner core)

double ftime = 0.0;		// current fluid time.
double t_forcing = 0.0;		// forcing time-scale.

//#define MASK
//#define COIL
//#define _MHD_CURRENT_FREE_SMALL_RM_

#include "grid.c"

struct TriDiagL *MB, *MB_1, *MUt, *MUt_1;
struct PentaDiag **MUp, **MUp_1;

#define DEB printf("%s:%u pass\n", __FILE__, __LINE__)
//#define DEB printf("thread %d :: %s:%u pass\n",omp_get_thread_num(), __FILE__, __LINE__)
//#define DEB (0)

#ifdef MASK
float **mask;

/// Allocate a constant scalar field like mask or constant ... (low precision)
void alloc_float(float ***msk, long int istart, long int iend)
{
	void *ptr;
	long int ir;
	
	ptr = malloc( NR*sizeof(float *) + (iend-istart+1)*NPHI*NLAT*sizeof(float) );
	*msk = (float **) ptr;
	ptr = (*msk + NR);
	(*msk)[istart] = (float *) ptr;
	for (ir = istart+1; ir <= iend; ir++)
		(*msk)[ir] = (*msk)[ir-1] + NPHI*NLAT;
	for (ir = 0; ir < istart; ir++) (*msk)[ir] = NULL;
	for (ir = iend+1; ir<NR; ir++) (*msk)[ir] = NULL;
}

void mask_cyl(float ***msk, double ar)
{
	double L, d, C;
	double rr,s,z,phi,theta;
	long int irs, i,j,k;
	
	L = 0.99/sqrt(1. + ar*ar);	// 0.99 : to avoid touching the outer sphere.
	d = L*ar;
	
	s = 3./NR;	z = pi/NLAT;
	if (z > s) s = z;
	C = nu/(4.*s*s);		// Penalization cannot exceed this value, for Boundary Layers must be resolved...
	if (C > 1.0/dtU) C = 1.0/dtU;	// cannot exceed 1/dt eiter !

	if (L>d) { rr = d; } else { rr = L; }
	irs = r_to_idx(rr)-1;	// ir_start for mask.
	alloc_float(msk, irs, NR-1);
	printf("[mask] Cylinder in Sphere : L=%.4f, d=%.4f, C=%g, ir_start=%d\n",L,d,C,irs);

	for (i = irs; i<NR; i++) {
		rr = r[i];
		for (j=0; j<NPHI; j++) {
			phi = j * (2.*pi/(MRES*NPHI));
			for (k=0; k<NLAT; k++) {
				(*msk)[i][j*NLAT + k] = 0.0;
				s = r[i]*st[k];	z = r[i]*ct[k];	theta = acos( ct[k] );

				if ((s > d)||(z > L)||(z < -L))
					(*msk)[i][j*NLAT + k] = C;
			}
		}
	}
}
#endif

#ifdef COIL
double Icoil = 0.0;
float **Jcoil;

void init_Jcoil(double s0, double z0, double h)
{
	double rr,s,z,phi,theta, j0;
	long int i,j,k, irs,ire;
	
	j0 = 1.0/(pi*h*h);	// convert intensity 1.0 to current density
	j0 *= eta;		// eta = sigma = convert current density to emf

	
	irs = r_to_idx(sqrt(s0*s0 + z0*z0) - h) -1;	// ir_start for coil.
	ire = r_to_idx(sqrt(s0*s0 + z0*z0) + h) +1;	// ir_end for coil.
	if (ire >= NR) runerr("coil does not fit in sphere...\n");
	printf("[coil] irs=%d, ire=%d\n",irs,ire);

	alloc_float(&Jcoil, irs, ire);
	DEB;

	for (i=irs; i<=ire; i++) {
		rr = r[i];
		for (j=0; j<NPHI; j++) {
			phi = j * (2.0*pi/(MRES*NPHI));
			for (k=0; k<NLAT; k++) {
				s = rr*st[k];	z = rr*ct[k];	theta = acos( ct[k] );
				Jcoil[i][j*NLAT + k] = 0.0;

				if ( (s-s0)*(s-s0) + (z-z0)*(z-z0) < h*h ) {
					Jcoil[i][j*NLAT + k] = j0;
				}
			}
		}
	}
}
#endif



#include "xshells_fields.c"
#include "xshells_io.c"


struct VectField B, U, W, J;
#ifdef _MHD_CURRENT_FREE_SMALL_RM_
	struct VectField B0;
#endif
struct PolTor Blm, Ulm, NLb1, NLb2, NLu1, NLu2;
struct StatSpecVect B0lm, J0lm;


int MAKE_MOVIE = 0;	// no output by default.
int SAVE_QUIT = 0;	// for signal catching.



	/* MATRIX INITIALIZATION */

/// Initialize semi-implicit part of Induction equation
///      matrix   MB = (1/dt + 0.5*eta*Lap)
/// and  matrix MB_1 = (1/dt - 0.5*eta*Lap)^(-1)
/// Same matrix MB and MB_1 is used for Poloidal and Toroidal components. (Polidal up to NR-1, Toroidal up to NR-2)
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

/// Initialize semi-implicit part of Induction equation
///      matrix   MUt = (1/dt + 0.5*nu*Lap)	  	MUp = -(1/dt + 0.5*nu*Lap)*Lap
/// and  matrix MUt_1 = (1/dt - 0.5*nu*Lap)^(-1)	MUp_1 = -((1/dt - 0.5*nu*Lap)*Lap)^(-1)
/// BC : boundary condition, can be BC_NO_SLIP or BC_FREE_SLIP
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
/// (d/dt - nu.Lap)(-Lap.Up) = Toroidal[curl_NL]
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
/// (d/dt - nu.Lap) Ut = Poloidal(curl_NL)
	i = NG;		// BC toroidal (free-slip) : d(T/r)/dr = 0 =>  dT/dr = T/r
		dx_1 = 1.0/(r[i+1]-r[i]);	dx_2 = dx_1*dx_1;
					MUt[i-NG].l = 0.0;
		for (l=0; l<=LMAX; l++) MUt[i-NG].d[l] = 1.0/dtU + nu*( -dx_2 - r_1[i]*dx_1 + (1.0 - 0.5*l*(l+1))*r_2[i] );
					MUt[i-NG].u =              nu * dx_2;
	for(i=NG+1; i<NR-1; i++) {
					MUt[i-NG].l =              0.5*nu*Lr[i].l;
		for (l=0; l<=LMAX; l++)	MUt[i-NG].d[l] = 1.0/dtU + 0.5*nu*(Lr[i].d - r_2[i]*l*(l+1));
					MUt[i-NG].u =              0.5*nu*Lr[i].u;
	}
	i = NR-1;	// BC toroidal (free-slip) : d(T/r)/dr = 0 =>  dT/dr = T/r
		dx_1 = 1.0/(r[i]-r[i-1]);	dx_2 = dx_1*dx_1;
					MUt[i-NG].l =              nu * dx_2;
		for (l=0; l<=LMAX; l++) MUt[i-NG].d[l] = 1.0/dtU + nu*( -dx_2 + r_1[i]*dx_1 + (1.0 - 0.5*l*(l+1))*r_2[i] );
					MUt[i-NG].u = 0.0;

	for(i=0; i<NU; i++) {
					MUt_1[i].l =            - MUt[i].l;
		for (l=0; l<=LMAX; l++) MUt_1[i].d[l] = 2.0/dtU - MUt[i].d[l];
					MUt_1[i].u =            - MUt[i].u;
	}
	if (BC == BC_FREE_SLIP) {
		TriDec(MUt_1, 0, NU-1);		// for free-slip : 0 to NU-1
	} else  TriDec(MUt_1, 1, NU-2);		// for no-slip : 1 to NU-2
	
// for NO-SLIP, use : cTriSolve(MUt_1, RHS, Tlm, 1, NU-2)
// for FREE-SLIP, use : cTriSolve(MUt_1, RHS, Tlm, 1, NU-1)
}

// **************************
// ***** TIME STEPPING ******
// **************************

/// substepB : substep advance for B, with NL and NLo set as non-linear terms at time t and t-dtB
inline substepB(struct PolTor *NL, struct PolTor * NLo, complex double **Btmp)
{
	long int i,l;

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

/// substepU : substep advance for U, with NL and NLo set as non-linear terms at time t and t-dtU
inline substepU(struct PolTor *NL, struct PolTor *NLo, complex double **Alm)
{
	long int i,l;

	cTriMulBC(MUt, Ulm.T +NG, Alm, 1, NU-2);	// no slip
//	cTriMul(MUt, Ulm.T +NG, Alm, 0, NU-1);		// free slip
	for (i=1; i<NU-1; i++) {
		for (l=1;l<NLM;l++) {
			Alm[i][l] += 1.5*NL->T[i+NG][l] - 0.5*NLo->T[i+NG][l];
		}
/*#ifdef MASK
		if (mask[i+NG] != NULL) {
			SH_to_spat(Alm[i], (complex double *) W.r[NG+1]);
			for (l=0; l<NPHI*NLAT; l++) W.r[NG+1][l] *= mask[i+NG][l];
			spat_to_SH((complex double *) W.r[NG+1], Alm[i]);
		}
#endif*/
	}
	cTriSolveBC(MUt_1, Alm, Ulm.T +NG, 1, NU-2);	// no slip
//	cTriSolve(MUt_1, Alm, Ulm.T +NG, 0, NU-1);	// free slip

	cPentaMul(MUp, Ulm.P +NG, Alm, 1, NU-2);
	for (i=1; i<NU-1; i++) {
		for (l=1;l<NLM;l++) {
			Alm[i][l] += 1.5*NL->P[i+NG][l] - 0.5*NLo->P[i+NG][l];
		}
/*#ifdef MASK
		if (mask[i+NG] != NULL) {
			SH_to_spat(Alm[i], (complex double *) W.r[NG+1]);
			for (l=0; l<NPHI*NLAT; l++) W.r[NG+1][l] *= mask[i+NG][l];
			spat_to_SH((complex double *) W.r[NG+1], Alm[i]);
		}
#endif*/
	}
	cPentaSolve(MUp_1, Alm, Ulm.P +NG, 1, NU-2);
}


/// updates the forcing at each time-step.
inline void calc_Uforcing(double t)
{
	double domega;

	if (t_forcing > 0.0) {			// no forcing if t_forcing = 0.0
		if (t > t_forcing * 10) {
			domega = 0.0;	t_forcing = 0.0;	// no more forcing afterwards.
		} else {
			t = (t - 3*t_forcing)/t_forcing;
			domega = DeltaOmega * exp(-t*t);	//   /(t_forcing*sqrt(pi)); //constant energy forcing.
		}
		Ulm.T[NG][LM(1,0)]  = r[NG]* domega * Y10_ct;	// set differential rotation of inner core.
	}
}

double step_NS(long int nstep, complex double **Alm)
{
	inline step1(struct PolTor *NLu, struct PolTor *NLuo)
	{
		calc_Uforcing(ftime);
		CALC_U(BC_NO_SLIP);
		calc_Vort(&Ulm, Omega0, &W);
		NL_Fluid(&U, &W, NLu);
		substepU(NLu, NLuo, Alm);
		ftime += dtU;
	}

	while(nstep > 0) {
		nstep--;
		step1(&NLu1, &NLu2);
		step1(&NLu2, &NLu1);
	}
	return 2*nstep*dtU;
}

/// perform nstep steps of Navier-Stokes equation, with temporary array Alm.
double step_MHD(long int nstep, complex double **Alm)
{
	inline step1(struct PolTor *NLu, struct PolTor *NLb, struct PolTor *NLuo, struct PolTor *NLbo)
	{
		
		calc_Uforcing(ftime);
		CALC_U(BC_NO_SLIP);
		calc_Vort(&Ulm, Omega0, &W);
#ifdef _MHD_CURRENT_FREE_SMALL_RM_
		calc_J(&Blm, &J, NULL);
		NL_MHD(&U, &W, &B0, &J, NLu);
		NL_Induction(&U, &B0, &B, NLb);
#else
		calc_B(&Blm, &B, &B0lm);
		calc_J(&Blm, &J, &J0lm);
		NL_MHD(&U, &W, &B, &J, NLu);
		NL_Induction(&U, &B, &J, NLb);
#endif
		substepU(NLu, NLuo, Alm);
		substepB(NLb, NLbo, Alm);
		ftime += dtU;
	}

	while(nstep > 0) {
		nstep--;
		step1(&NLu1, &NLb1, &NLu2, &NLb2);
		step1(&NLu2, &NLb2, &NLu1, &NLb1);
	}
	return 2*nstep*dtU;
}

/// perform nstep steps of full coupled dynamo equation, with temporary array Alm.
double step_Dynamo(long int nstep, complex double **Alm)
{
	inline step1(struct PolTor *NLu, struct PolTor *NLb, struct PolTor *NLuo, struct PolTor *NLbo)
	{
		calc_Uforcing(ftime);
		CALC_U(BC_NO_SLIP);
		calc_Vort(&Ulm, Omega0, &W);
		calc_B(&Blm, &B, NULL);
		calc_J(&Blm, &J, NULL);
		NL_MHD(&U, &W, &B, &J, NLu);
		NL_Induction(&U, &B, &B, NLb);
		substepU(NLu, NLuo, Alm);
		substepB(NLb, NLbo, Alm);
		ftime += dtU;
	}

	while(nstep > 0) {
		nstep--;
		step1(&NLu1, &NLb1, &NLu2, &NLb2);
		step1(&NLu2, &NLb2, &NLu1, &NLb1);
	}
	return 2*nstep*dtU;
}



/// perform nstep steps of Induction equation, with temporary array Alm.
double step_Dyncin(long int nstep, complex double **Alm)
{
	inline step1(struct PolTor *NLb, struct PolTor *NLbo)
	{
		calc_B(&Blm, &B, NULL);
		NL_Induction(&U, &B, &B, NLb);
		substepB(NLb, NLbo, Alm);
	}

	while(nstep > 0) {
		nstep--;
		step1(&NLb1, &NLb2);
		step1(&NLb2, &NLb1);
	}
	return 2*nstep*dtB;
}


/*
void save_n_quit(char* fprefix, char* job, long int it, double t)
{
	char fn[100];

	sprintf(fn,"%s.%s",fprefix,job);
	printf("  => Signal received : saving before quitting...\n");	fflush(stdout);
	save_all(fn, it, t);
	printf("     current state successfully saved to file %s\n",fn);	fflush(stdout);
	exit(0);
}
*/

void read_Par(char *fname, char *job, long int *iter_max, long int *modulo, double *polar_opt_max, double *Ric, double *B0, double *B1)
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
	DeltaOmega = 0.0;	// no differential rotation.
	t_forcing = 0.0;	// permanent forcing.
	eta = 0;	nu = 0;		// if eta or nu is not set, copy the other value. (Pm=1)

	while (!feof(fp))
	{
		fgets(str, 256, fp);
		fputs(str, fpw);		// save file...
//		printf(str);
		if ((str[0] != '#' )&&(strlen(str) > 3)) {	// ignore les lignes commentees (#) et les lignes vides.
			sscanf(str,"%s = %lf",name,&tmp);
//			printf("[%s] = %f\n",name,tmp);
			// JOB
			if (strcmp(name,"job") == 0) 		sscanf(str,"%s = %s",name,job);
			if (strcmp(name,"movie") == 0)		MAKE_MOVIE = tmp;
			// PHYSICAL PARAMS
			if (strcmp(name,"Ric") == 0)		*Ric = tmp;
			if (strcmp(name,"Omega0") == 0)		Omega0 = tmp;
			if (strcmp(name,"nu") == 0)		nu = tmp;
			if (strcmp(name,"eta") == 0)		eta = tmp;
			if (strcmp(name,"B0") == 0)		*B0 = tmp;
			if (strcmp(name,"B1") == 0)		*B1 = tmp;
			// FORCING
			if (strcmp(name,"DeltaOmega") == 0)	DeltaOmega = tmp;
			if (strcmp(name,"t_forcing") == 0)	t_forcing = tmp;
			// NUMERICAL SCHEME
			if (strcmp(name,"NR") == 0)		NR = tmp;
			if (strcmp(name,"NH") == 0)		NH = tmp;
			if (strcmp(name,"dtU") == 0)		dtU = tmp;
			if (strcmp(name,"dtB") == 0)		dtB = tmp;
			// TIME DOMAIN
			if (strcmp(name,"iter_max") == 0)	*iter_max = tmp;
			if (strcmp(name,"modulo") == 0)		*modulo = tmp;
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
	if (eta == 0) eta = nu;
	if (nu == 0) nu = eta;

	fclose(fp);	fclose(fpw);
	sprintf(str, "%s.%s",fname,job);	rename("tmp.par", str);		// rename file.	
	fflush(stdout);		// when writing to a file.
}

/*double j1(double x)
{
	if (x==0.0) return(0.0);
	return (sin(x)/x - cos(x))/x;
}*/

void sig_handler(int sig_num)
{
	switch(sig_num) {
		case 15 : SAVE_QUIT = -2; break;	// TERM signal : save state as soon as possible !
		case 30 : SAVE_QUIT = 1; break;		// Save snapshot
		case 31 : SAVE_QUIT = -1;		// Save & quit
	}
}

int main (int argc, char *argv[])
{
	double (*ptr_step)(long int, complex double **);
	double Ric = 0.0;		// default ic radius
	double polar_opt_max = 0.0;	// default SHT optimization.
	double b0 = 0.0;
	double b1 = 0.0;
	double rr, t0,t1,Rm,Rm2;
	long int i,j,k, im,m,l, it, lmtest;
	long int iter_max, modulo;
	complex double **Alm;		// temp scalar spectral (pol or tor) representation
	FILE *fp;
	struct JobInfo jinfo;
	char command[100] = "xshells.par";
	char job[40];
	clock_t tcpu;

	printf("[XSHELLS] eXtendable Spherical Harmonic Earth-Like Liquid Simulator\n          by Nathanael Schaeffer / LGIT, build %s, %s\n",__DATE__,__TIME__);
#ifdef _NTH_
	printf("          ++ OMP Parallel version, %d threads.\n",_NTH_);
#endif

	read_Par(command, job, &iter_max, &modulo, &polar_opt_max, &Ric, &b0, &b1);
	init_SH(polar_opt_max);
	init_rad_sph(0.0, Ric, 1.0);

	if (b0 != 0.0) init_Bmatrix();
	init_Umatrix(BC_NO_SLIP);

	alloc_DynamicField(&Ulm, &U, &NLu1, &NLu2, NG, NR-1);
	alloc_VectField(&W,NG,NR-1);
	if (b0 != 0.0) {
		alloc_DynamicField(&Blm, &B, &NLb1, &NLb2, 0, NR-1);
		alloc_VectField(&J,0,NR-1);
#ifdef _MHD_CURRENT_FREE_SMALL_RM_
		alloc_VectField(&B0,0,NR-1);
#endif
	}

	printf("[Params] job name : %s\n",job);
	printf("         Ek=%.2e, Ro=%.2e, Pm=%.2e, Re=%.2e, S=%.2e, N=%.2e, M=%.2e, Elsasser=%.2e, Lehnert=%.2e\n",nu/Omega0, Ric*DeltaOmega/Omega0, nu/eta, Ric*DeltaOmega/nu, b0/eta, b0*b0/(Ric*DeltaOmega*Ric*DeltaOmega) , b0*0.5/sqrt(nu*eta), b0*b0/(eta*Omega0), b0/Omega0);
	printf("         dtU.Omega=%.2e, dtU.nu.R^2=%.2e, dtB.eta.R^2=%.2e, dtB/dtU=%.2e\n",dtU*Omega0, dtU*nu, dtB*eta, dtB/dtU);
	fflush(stdout);

	Alm = (complex double **) malloc( NR * sizeof(complex double *));
	for (i=0;i<NR;i++)
		Alm[i] = (complex double *) malloc( NLM * sizeof(complex double));

	if (b0 != 0.0) {
		ptr_step = &step_MHD;
		// init B fields.
		zero_out_field(&Blm);
		#define Set_Poloidal( cval ) Blm.P[i][LM(l,m)] = b0*(cval);
		#define Set_Toroidal( cval ) Blm.T[i][LM(l,m)] = b0*(cval);
		for (i=0; i<NR; i++) {
			rr = r[i];
			#include "inc_B0ini.c"
		}
		sprintf(command,"poltorB0.%s",job);	save_PolTor(command, &Blm, ftime, BC_MAGNETIC);
		#ifdef _MHD_CURRENT_FREE_SMALL_RM_
			PolTor_to_spat(&Blm, &B0, 0, NR-1, BC_MAGNETIC);	// background magnetic field.
			printf("=> Imposed magnetic field B0 : using small Rm MHD integration (assuming current-free B0)\n");
		#else
			printf("=> Imposed magnetic field B0 : using full MHD integration\n");
		#endif
		i = Make_StatSpecVect(&Blm, &B0lm, &J0lm, 1, NR-2, BC_NONE);
		#ifdef INIT_FIELD_NAME
			printf("    B0 set to : \"" INIT_FIELD_NAME "\" (%d modes)\n",i);
		#endif

		zero_out_field(&Blm);
		if (argc > 2) {
			load_PolTor(argv[2], &Blm, &jinfo);
			ftime = jinfo.t;
			printf("    B read from file \"%s\" (t=%f).\n",argv[2],ftime);
		}
	} else {
		printf("=> No imposed magnetic field : using Navier-Stokes integration\n");
		ptr_step = &step_NS;
	}

	// init U fields
	if (t_forcing > 0.0) {
		printf("=> Forcing on time-scale t_forcing=%f\n",t_forcing);
	}
	if (argc > 1) {
		load_PolTor(argv[1], &Ulm, &jinfo);		// load from file.
		ftime = jinfo.t;
		printf("    U field read from file \"%s\" (t=%f).\n",argv[1],ftime);
	} else {
		zero_out_field(&Ulm);
		for (i=NG+1; i<NR-1; i++) {
			for (im=1;im<=MMAX;im++)
				Ulm.T[i][LM(MRES*im+1,im)] = DeltaOmega*1e-8;	// non-zero initial value for m>0
		}
		Ulm.T[NG][LM(1,0)]  = r[NG]*DeltaOmega * Y10_ct;	// differential rotation of inner core.
		sprintf(command,"poltorU0.%s",job);	save_PolTor(command, &Ulm, ftime, BC_NO_SLIP);
	}

#ifdef MASK
	mask_cyl(&mask, 1.0);	// define cylindric mask
#endif
#ifdef COIL
	init_Jcoil(0.4, -0.77, 0.06);
	Icoil = 1.0;
#endif

// Initializing signal handler.
	signal(30,&sig_handler);	signal(31,&sig_handler);	signal(15,&sig_handler);

	printf("let's go for %d iterations !\n",iter_max);
	it =0;		// some pre-requisites are computed here.
		calc_Uforcing(ftime);		// update forcing
		CALC_U(BC_NO_SLIP);
		calc_Vort(&Ulm, Omega0, &W);
		if (b0 == 0.0) {
			NL_Fluid(&U, &W, &NLu2);
		} else {
#ifdef _MHD_CURRENT_FREE_SMALL_RM_
			calc_J(&Blm, &J, NULL);
			NL_MHD(&U, &W, &B0, &J, &NLu2);
			NL_Induction(&U, &B0, &B, &NLb2);
#else
			calc_B(&Blm, &B, &B0lm);
			calc_J(&Blm, &J, &J0lm);
			NL_MHD(&U, &W, &B, &J, &NLu2);
			NL_Induction(&U, &B, &B, &NLb2);
#endif
		}
		t0 = creal(Ulm.T[NG][1]);	t1 = creal(Ulm.P[NG+1][2]);
		printf("[it %d] t=%g, P0=%g, T0=%g\n",it,ftime,t1, t0);	fflush(stdout);

/* MAIN LOOP */
	while (it < iter_max) {
		tcpu = clock();
		(*ptr_step)(modulo, Alm);	// go for time integration !!
		tcpu = clock() - tcpu;

		it++;
		t0 = creal(Ulm.T[NG][1]);	t1 = creal(Ulm.P[NG+1][2]);
		printf("[it %d] t=%g, P0=%g, T0=%g  (cpu=%d)\n",it,ftime,t1, t0, (int) tcpu);	fflush(stdout);
		if (isnan(t0)) runerr("NaN encountered");

//		t0 = t1;
//		printf("%g :: tx=%g\n",t1,log(t1/t0)/(2*dtB));

		if (SAVE_QUIT < 0) {
			printf("  > signal received : save & quit.\n");	fflush(stdout);
			break;
		}
		if (SAVE_QUIT > 0) {
			sprintf(command,"poltorU.%s",job);	save_PolTor(command, &Ulm, ftime, BC_NO_SLIP);
			if (b0 != 0.0) { sprintf(command,"poltorB.%s",job);	save_PolTor(command, &Blm, ftime, BC_MAGNETIC); }
			printf(" > signal received, snapshot dumped.\n");	fflush(stdout);
			SAVE_QUIT = 0;
		}
		if (MAKE_MOVIE != 0) {		// save single-precision data snapshot.
			sprintf(command,"poltorU_%04d.%s",it,job);	save_PolTor_single(command, &Ulm, ftime, BC_NO_SLIP);
			if (b0 != 0.0) { sprintf(command,"poltorB_%04d.%s",it,job);	save_PolTor_single(command, &Blm, ftime, BC_MAGNETIC); }
		}
	}

	// save full double-precision data at the end (for restart)
	sprintf(command,"poltorU.%s",job);	save_PolTor(command, &Ulm, ftime, BC_NO_SLIP);
	if (b0 != 0.0) { sprintf(command,"poltorB.%s",job);	save_PolTor(command, &Blm, ftime, BC_MAGNETIC); }
	printf("[XSHELLS] all done.\n");
}
