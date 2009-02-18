# This file is meta-code for SHT.c (spherical harmonic transform).
# it is intended for "make" to generate C code for similar SHT functions,
# from one generic function + tags.
# > See Makefile and SHT.c
# Basically, there are tags at the beginning of lines that are information
# to keep or remove the line depending on the function to build.
# tags :
# Q : line for scalar transform
# V : line for vector transform (both spheroidal and toroidal)
# S : line for vector transfrom, spheroidal component
# T : line for vector transform, toroidal component.

/////////////////////////////////////////////////////
//   Inverse Spherical Harmonic Transform
// input  : Qlm,Slm,Tlm = spherical harmonics coefficients of Scalar, Spheroidal and Toroidal scalars : 
//          complex double array of size NLM [unmodified]
// output : BrF, BtF, BpF = theta, and phi vector components, spatial/fourrier data : 
//          complex double array of size NLAT*(NPHI/2+1) or double array of size NLAT*(NPHI/2+1)*2
#Q void SH_to_spat(complex double *Qlm, complex double *BrF)
#V void SHsphtor_to_spat(complex double *Slm, complex double *Tlm, complex double *BtF, complex double *BpF)
# {
Q	complex double fe, fo;		// even and odd parts
S	complex double se, so, dse, dso;	// spheroidal even and odd parts
T	complex double te, to, dte, dto;	// toroidal ...
  #if NPHI > 1
Q	#define BR0 BrF
V	#define BT0 BtF
V	#define BP0 BpF
  #else
Q	#define BR0 ((double *)BrF)
V	#define BT0 ((double *)BtF)
V	#define BP0 ((double *)BpF)
  #endif
Q	complex double *Ql;
S	complex double *Sl;
T	complex double *Tl;
Q	double *yl;
V	double *dyl0;	// for m=0
V	struct DtDp *dyl;
	long int i,im,m,l;

	im = 0;		// zonal part : d/dphi = 0;
		m = im*MRES;
	/*#ifdef SHT_DCT
Q		Ql = &Qlm[LiM(0,im)];		// virtual pointer for l=0 and im
Q		yl = ylm_dct[im];
		for (it=0; it<NLAT; it++)
Q			BrF[it] = 0.0;		// zero out array (includes DCT padding)
		for (l=m; l<LTR; l+=2) {
			for (it=0; it<=l; it+=2) {
Q				(double) BrF[it]   += yl[it] *   (double) Ql[l];
Q				(double) BrF[it+1] += yl[it+1] * (double) Ql[l+1];
			}
Q			yl += (l+2 - (m&1));
		}
		if (l==LTR) {
			for (it=0; it<=l; it+=2) {
Q				(double) BrF[it] += yl[it] * (double) Ql[l];
			}
		}
		fftw_execute_r2r(idctm0,(double *) BrF, (double *) BrF);		// iDCT
	#else*/
Q		Ql = &Qlm[LiM(0,im)];	// virtual pointer for l=0 and im
S		Sl = &Slm[LiM(0,im)];	// virtual pointer for l=0 and im
T		Tl = &Tlm[LiM(0,im)];
		i=0;
Q		yl  = ylm[im];
V		dyl0 = (double *) dylm[im];	// only theta derivative (d/dphi = 0 for m=0)
		while (i < NLAT_2) {	// ops : NLAT_2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
Q			fe = 0.0;	fo = 0.0;
S			dse = 0.0;	dso = 0.0;
T			dte = 0.0;	dto = 0.0;
			while (l<LTR) {	// compute even and odd parts
QE				(double) fe += yl[0] * (double) Ql[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
QO				(double) fo += yl[1] * (double) Ql[l+1];	// fo += ylm[im][i*(LMAX-m+1) + (l+1-m)] * Qlm[LiM(l+1,im)];
TO				(double) dto += dyl0[0] * (double) Tl[l];	// m=0 : everything is real.
SE				(double) dso += dyl0[0] * (double) Sl[l];
TE				(double) dte += dyl0[1] * (double) Tl[l+1];
SO				(double) dse += dyl0[1] * (double) Sl[l+1];
				l+=2;
Q				yl+=2;
V				dyl0+=2;
			}
			if (l==LTR) {
QE				(double) fe += yl[0] * (double) Ql[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
TO				(double) dto += dyl0[0] * Tl[l];
SE				(double) dso += dyl0[0] * Sl[l];
Q				yl++;
V				dyl0++;
			}
Q			BR0[i] = fe + fo;
V			BT0[i] = 0.0
S				+ (dse+dso)			// Bt = dS/dt
V				;
V			BP0[i] = 0.0
T		 		- (dte+dto)			// Bp = - dT/dt
V				;
			i++;
QB			BR0[NLAT-i] = fe - fo;
VB			BT0[NLAT-i] = 0.0
SB		 		+ (dse-dso)
VB				;
VB			BP0[NLAT-i] = 0.0
TB				- (dte-dto)
VB				;
Q			yl  += (LMAX-LTR);
V			dyl0 += (LMAX-LTR);
		}
	//#endif
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
	for (im=1; im<=MTR; im++) {
		m = im*MRES;
Q		Ql = &Qlm[LiM(0,im)];	// virtual pointer for l=0 and im
S		Sl = &Slm[LiM(0,im)];	// virtual pointer for l=0 and im
T		Tl = &Tlm[LiM(0,im)];
		i=0;
		while (i<tm[im]) {	// polar optimization
Q			BrF[i] = 0.0;
QB			BrF[NLAT-tm[im] + i] = 0.0;	// south pole zeroes <=> BrF[im*NLAT + NLAT-(i+1)] = 0.0;
V			BtF[i] = 0.0;
VB			BtF[NLAT-tm[im] + i] = 0.0;	// south pole zeroes
V			BpF[i] = 0.0;
VB			BpF[NLAT-tm[im] + i] = 0.0;	// south pole zeroes
			i++;
		}
Q		yl  = ylm[im] + i*(LMAX-m+1);
V		dyl = dylm[im] + i*(LMAX-m+1);
		while (i < NLAT_2) {	// ops : NLAT_2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
Q			fe = 0.0;	fo = 0.0;
S			dse = 0.0;	dso = 0.0;	se = 0.0;	so = 0.0;
T			dte = 0.0;	dto = 0.0;	te = 0.0;	to = 0.0;
			while (l<LTR) {	// compute even and odd parts
QE				fe  += yl[0] * Ql[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
QO				fo  += yl[1] * Ql[l+1];	// fo += ylm[im][i*(LMAX-m+1) + (l+1-m)] * Qlm[LiM(l+1,im)];
TO				dto += dyl[0].t * Tl[l];
TO				te  += dyl[0].p * Tl[l];
SE				dso += dyl[0].t * Sl[l];
SE				se  += dyl[0].p * Sl[l];
TE				dte += dyl[1].t * Tl[l+1];
TE				to  += dyl[1].p * Tl[l+1];
SO				dse += dyl[1].t * Sl[l+1];
SO				so  += dyl[1].p * Sl[l+1];
				l+=2;
Q				yl+=2;
V				dyl+=2;
			}
			if (l==LTR) {
QE				fe  += yl[0] * Ql[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
TO				dto += dyl[0].t * Tl[l];
TO				te  += dyl[0].p * Tl[l];
SE				dso += dyl[0].t * Sl[l];
SE				se  += dyl[0].p * Sl[l];
Q				yl++;
V				dyl++;
			}
Q			BrF[i] = fe + fo;
V			BtF[i] = 		// Bt = dS/dt       + I.m/sint *T
S					(dse+dso)
T					+ I*(te+to)
V					;
V			BpF[i] = 		// Bp = I.m/sint * S - dT/dt
S					I*(se+so)
T					- (dte+dto)
V					;
			i++;
QB			BrF[NLAT-i] = fe - fo;
VB			BtF[NLAT-i] =
SB					(dse-dso)
TB					+ I*(te-to)
VB					;
VB			BpF[NLAT-i] =
SB					I*(se-so)
TB					- (dte-dto)
VB					;
Q			yl  += (LMAX-LTR);
V			dyl += (LMAX-LTR);
		}
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
	}
	for (i=0; i < NLAT*(NPHI/2 -MTR); i++) {	// padding for high m's
Q			BrF[i] = 0.0;
V			BtF[i] = 0.0;	BpF[i] = 0.0;
	}

Q	BrF -= NLAT*(MTR+1);		// restore original pointer
V	BtF -= NLAT*(MTR+1);	BpF -= NLAT*(MTR+1);		// restore original pointers
 #if NPHI>1
Q	fftw_execute_dft_c2r(ifft, BrF, (double *) BrF);
V	fftw_execute_dft_c2r(ifft, BtF, (double *) BtF);
V	fftw_execute_dft_c2r(ifft, BpF, (double *) BpF);
 #endif
# }
