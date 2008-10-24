// xshells_fields :: 3D Fields definition and operations, in spectral (PolTor) and physical space (VectField).

// boundary condition constants.
#define BC_NONE 0
#define BC_NO_SLIP 1
#define BC_FREE_SLIP 2
#define BC_MAGNETIC 3

// field definitions
struct VectField {
	double **r,**t,**p;
};
struct PolTor {
	complex double **P,**T;
};


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
	complex double T[NLM];	//  [ l(l+1)/r^2 - 1/r d2/dr2(r .) ] P = -Lap P
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
				T[lm] = (r_2[ir]*l2[lm] - Lr[ir].d)*PT->P[ir][lm] - Lr[ir].l*PT->P[ir-1][lm] - Lr[ir].u*PT->P[ir+1][lm];
			}
			istart++;
		}
	}
	for (ir=istart; ir < iend; ir++) {
		for (lm=0; lm<NLM; lm++) {		// Solenoidal deduced from radial derivative of Poloidal
			Q[lm] = r_1[ir]*l2[lm] * PT->T[ir][lm];
			S[lm] = Wr[ir].l*PT->T[ir-1][lm] + Wr[ir].d*PT->T[ir][lm] + Wr[ir].u*PT->T[ir+1][lm];
			T[lm] = (r_2[ir]*l2[lm] - Lr[ir].d)*PT->P[ir][lm] - Lr[ir].l*PT->P[ir-1][lm] - Lr[ir].u*PT->P[ir+1][lm];
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
			T[lm] = (r_2[ir]*l2[lm] - Lr[ir].d)*PT->P[ir][lm] - Lr[ir].l*PT->P[ir-1][lm] - Lr[ir].u*PT->P[ir+1][lm];
			Q[lm] = r_1[ir]*l2[lm] * PT->T[ir][lm];
			S[lm] = Wr[ir].l*PT->T[ir-1][lm] + Wr[ir].d*PT->T[ir][lm] + Wr[ir].u*PT->T[ir+1][lm];
		}
		(double) Q[1] += Om0;		// add Background Vorticity for Coriolis Force (l=1, m=0)
		(double) S[1] += Om0;		// add Background Vorticity for Coriolis Force (l=1, m=0)
		SH_to_spat(Q,(complex double *) W->r[ir]);
		SHsphtor_to_spat(S, T, (complex double *) W->t[ir], (complex double *) W->p[ir]);
	}
}

/// compute current density field from poloidal/toroidal components of magnetic field [ie poltor_to_rot_spat applied to B]
/// IN: PT : pol/tor components of magnetic field
/// OUT: J : r,theta,phi components of current density field
void calc_J(struct PolTor *PT, struct VectField *J)
{
	complex double Q[NLM];
	complex double S[NLM];
	complex double T[NLM];
	long int ir,lm;

//	Plm -= NG;	Tlm -= NG;	Vr -= NG;	Vt -= NG;	Vp -= NG;	// adjust pointers.
	for (ir=NG+1; ir <= NR-2; ir++) {
		for (lm=0; lm<NLM; lm++) {
			T[lm] = (r_2[ir]*l2[lm] - Lr[ir].d)*PT->P[ir][lm] - Lr[ir].l*PT->P[ir-1][lm] - Lr[ir].u*PT->P[ir+1][lm];
			Q[lm] = r_1[ir]*l2[lm] * PT->T[ir][lm];
			S[lm] = Wr[ir].l*PT->T[ir-1][lm] + Wr[ir].d*PT->T[ir][lm] + Wr[ir].u*PT->T[ir+1][lm];
		}
		SH_to_spat(Q,(complex double *) J->r[ir]);
		SHsphtor_to_spat(S, T, (complex double *) J->t[ir], (complex double *) J->p[ir]);
	}
}

// *****************************
// ***** NON-LINEAR TERMS ******
// *****************************

/// compute rot(VxW + JxB) and its pol-tor components.
/// output : V, B, J unchanged, W = VxW + JxB
///          NL = PolTor[rot(VxW + JxB)]
void NL_MHD(struct VectField *V, struct VectField *W, struct VectField *B, struct VectField *J, struct PolTor *NL)
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
void NL_Fluid(struct VectField *V, struct VectField *W, struct PolTor *NL)
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
void NL_Induction(struct VectField *V, struct VectField *B, struct VectField *VxB, struct PolTor *NL)
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



/// Allocate memory for a dynamic vector field (like U or B) and its pol/tor representation + non-linear term storage.
/// (V, NL1 and NL2 are optional)
void alloc_DynamicField(struct PolTor *PT, struct VectField *V, struct PolTor *NL1, struct PolTor *NL2, long int istart, long int iend)
{
	long int ir, k;

	// how much space is required for NL terms ?	
	k = 6;
	if (NL2 == NULL)  k = 4;
	if (NL1 == NULL)  { k = 2; NL2 = NULL; };	// does not allow to set NL2 if NL1 not set.

	// alloc radial structure.
	PT->P = (complex double **) malloc( k*NR * sizeof(complex double *) );
	PT->T = PT->P + NR;
	if (NL1 != NULL) { NL1->P = PT->P + 2*NR;	NL1->T = PT->P + 3*NR; }
	if (NL2 != NULL) { NL2->P = PT->P + 4*NR;	NL2->T = PT->P + 5*NR; }
	if (V != NULL) {
		V->r = (double **) malloc( 3*NR * sizeof(double *) );
		V->t = V->r + NR;	V->p = V->r + 2*NR;
	}

	// not defined here, set as NULL pointer.
	for (ir = 0; ir < istart; ir++) {
		PT->P[ir] = NULL;	PT->T[ir] = NULL;
		if (V != NULL) { V->r[ir] = NULL;	V->t[ir] = NULL;	V->p[ir] = NULL; }
		if (NL1 != NULL) { NL1->P[ir] = NULL;	NL1->T[ir] = NULL; }
		if (NL2 != NULL) { NL2->P[ir] = NULL;	NL2->T[ir] = NULL; }
	}
	for (ir = istart; ir <= iend; ir++) {		// shell by shell allocation.
		if (V != NULL) {
			V->r[ir] = (double *) fftw_malloc( 3*(NPHI/2+1)*NLAT * sizeof(complex double));	// 3 components
			V->t[ir] = V->r[ir] + 2*(NPHI/2+1)*NLAT;
			V->p[ir] = V->t[ir] + 2*(NPHI/2+1)*NLAT;
		}
		PT->P[ir] = (complex double *) malloc( k*NLM * sizeof(complex double));	// Pol/Tor
		PT->T[ir] = PT->P[ir] + NLM;
		if (NL1 != NULL) {
			NL1->P[ir] = PT->P[ir] + 2*NLM;
			NL1->T[ir] = PT->P[ir] + 3*NLM;
		}
		if (NL2 != NULL) {
			NL2->P[ir] = PT->P[ir] + 4*NLM;
			NL2->T[ir] = PT->P[ir] + 5*NLM;
		}
	}
	for (ir = iend+1; ir<NR; ir++) {
		PT->P[ir] = NULL;	PT->T[ir] = NULL;
		if (V != NULL) { V->r[ir] = NULL;	V->t[ir] = NULL;	V->p[ir] = NULL; }
		if (NL1 != NULL) { NL1->P[ir] = NULL;	NL1->T[ir] = NULL; }
		if (NL2 != NULL) { NL2->P[ir] = NULL;	NL2->T[ir] = NULL; }
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
