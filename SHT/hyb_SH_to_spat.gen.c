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
//   Inverse Spherical Harmonic Transform (Synthesis) using DCT on a regular theta grid.
// input  : Qlm,Slm,Tlm = spherical harmonics coefficients of Scalar, Spheroidal and Toroidal scalars : 
//          complex double array of size NLM [unmodified]
// output : BrF, BtF, BpF = r, theta, phi vector components, spatial/fourrier data : 
//          complex double array of size NLAT*(NPHI/2+1) or double array of size NLAT*(NPHI/2+1)*2
// MTR_DCT : -1 => no dct
//            0 => dct for m=0 only
//            m => dct up to m, (!!! MTR_DCT <= MTR !!!)
#Q void spat_to_SH_dct(complex double *BrF, complex double *Qlm)
#V void SHsphtor_to_spat_dct(complex double *Slm, complex double *Tlm, complex double *BtF, complex double *BpF)
# {
Q	complex double *Ql;
S	complex double *Sl;
T	complex double *Tl;
Q	double *yl;
V	double *dyl0;
  #if NPHI > 1
V	struct DtDp *dyl;
Q	complex double re,ro;
V	complex double te,to, pe,po;
S	complex double dte, dto;	// spheroidal even and odd parts
T	complex double dpe, dpo;	// toroidal ...
Q	#define BR0 BrF
V	#define BT0 BtF
V	#define BP0 BpF
  #else
Q	double re,ro;
V	double te,to, pe,po;
Q	#define BR0 ((double *)BrF)
V	#define BT0 ((double *)BtF)
V	#define BP0 ((double *)BpF)
  #endif
	long int k,im,m,l;

	im=0;	m=0;
Q		Ql = &Qlm[LiM(0,im)];	// virtual pointer for l=0 and im
S		Sl = &Slm[LiM(0,im)];	// virtual pointer for l=0 and im
T		Tl = &Tlm[LiM(0,im)];
	if (MTR_DCT >= 0) {	// dct for m=0.
Q		yl = ykm_dct[im];
V		dyl0 = (double *) dykm_dct[im];		// only theta derivative (d/dphi = 0 for m=0)
		k=0;
		while (k<LTR) {
			l = k;
Q			re = 0.0;	ro = 0.0;
V			te = 0.0;	to = 0.0;	pe = 0.0;	po = 0.0;
			while(l<LTR) {
QE				re += yl[0]  * (double) Ql[l];
QO				ro += yl[1]  * (double) Ql[l+1];
SE				to += dyl0[0] * (double) Sl[l];
SO				te += dyl0[1] * (double) Sl[l+1];
TO				po -= dyl0[0] * (double) Tl[l];
TE				pe -= dyl0[1] * (double) Tl[l+1];
				l+=2;
Q				yl+=2;
V				dyl0+=2;
			}
			if (l==LTR) {
QE				re += yl[0]  * (double) Ql[l];
SE				to += dyl0[0] * (double) Sl[l];
TO				po -= dyl0[0] * (double) Tl[l];
Q				yl++;
V				dyl0++;
			}
Q			BR0[k] = re;	BR0[k+1] = ro;
V			BT0[k] = te;	BT0[k+1] = to;
V			BP0[k] = pe;	BP0[k+1] = po;
			k+=2;
Q			yl+= (LMAX-LTR);
V			dyl0+= (LMAX-LTR);
		}
QE		if ((LTR&1)==0) {	// k=LTR
QE			BR0[k] = yl[0] * Ql[k];
QE			k++;
QE		}
		while (k<NLAT) {	// dct padding
Q			BR0[k] = 0.0;
V			BT0[k] = 0.0;	BP0[k] = 0.0;
			k++;
		}
	} else {
		k=0;
Q		yl  = ylm[im];
V		dyl0 = (double *) dylm[im];	// only theta derivative (d/dphi = 0 for m=0)
		while (k < NLAT_2) {	// ops : NLAT_2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=0;
Q			re = 0.0;	ro = 0.0;
S			te = 0.0;	to = 0.0;
T			pe = 0.0;	po = 0.0;
			while (l<LTR) {	// compute even and odd parts
QE				re += yl[0] * (double) Ql[l];		// re += ylm[im][k*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
QO				ro += yl[1] * (double) Ql[l+1];	// ro += ylm[im][k*(LMAX-m+1) + (l+1-m)] * Qlm[LiM(l+1,im)];
TO				po += dyl0[0] * (double) Tl[l];	// m=0 : everything is real.
SE				to += dyl0[0] * (double) Sl[l];
TE				pe += dyl0[1] * (double) Tl[l+1];
SO				te += dyl0[1] * (double) Sl[l+1];
				l+=2;
Q				yl+=2;
V				dyl0+=2;
			}
			if (l==LTR) {
QE				re += yl[0] * (double) Ql[l];		// re += ylm[im][k*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
TO				po += dyl0[0] * (double) Tl[l];
SE				to += dyl0[0] * (double) Sl[l];
Q				yl++;
V				dyl0++;
			}
Q			BR0[k] = re + ro;
V			BT0[k] = 0.0
S				+ (te+to)			// Bt = dS/dt
V				;
V			BP0[k] = 0.0
T		 		- (pe+po)			// Bp = - dT/dt
V				;
			k++;
QB			BR0[NLAT-k] = re - ro;
VB			BT0[NLAT-k] = 0.0
SB		 		+ (te-to)
VB				;
VB			BP0[NLAT-k] = 0.0
TB				- (pe-po)
VB				;
Q			yl  += (LMAX-LTR);
V			dyl0 += (LMAX-LTR);
		}
	}
  #if NPHI > 1
	im=1;
Q	BrF += NLAT;
V	BtF += NLAT;	BpF += NLAT;
	while(im<=MTR_DCT) {		// dct for im <= MTR_DCT
		m=im*MRES;
Q		Ql = &Qlm[LiM(0,im)];		// virtual pointer for l=0 and im
S		Sl = &Slm[LiM(0,im)];
T		Tl = &Tlm[LiM(0,im)];
Q		yl = ykm_dct[im];
V		dyl = dykm_dct[im];
		k=0;
		while (k<LTR) {
Q			l = (k < m) ? m : k+(m&1);
V			l = (k < m) ? m : k-(m&1);
Q			re = 0.0;	ro = 0.0;
V			te = 0.0;	to = 0.0;	pe = 0.0;	po = 0.0;
			while(l<LTR) {
QE				re += yl[0] * Ql[l];
QO				ro += yl[1] * Ql[l+1];
VO				te +=
TO					  dyl[0].p * (I*Tl[l])
SO					+ dyl[1].t * Sl[l+1]
VO					;
VO				po +=
SO					  dyl[1].p * (I*Sl[l+1])
TO					- dyl[0].t * Tl[l]
VO					;
VE				pe += 
SE					  dyl[0].p * (I*Sl[l])
TE					- dyl[1].t * Tl[l+1]
VE					;
VE				to += 
TE					  dyl[1].p * (I*Tl[l+1])
SE					+ dyl[0].t * Sl[l]
VE					;
				l+=2;
Q				yl+=2;
V				dyl+=2;
			}
			if (l==LTR) {
QE				re += yl[0]   * Ql[l];
TO				te += dyl[0].p * (I*Tl[l]);
SE				pe += dyl[0].p * (I*Sl[l]);
SE				to += dyl[0].t * Sl[l];
TO				po -= dyl[0].t * Tl[l];
Q				yl++;
V				dyl++;
			}
Q			BrF[k] = re;	BrF[k+1] = ro;
V			BtF[k] = te;	BtF[k+1] = to;
V			BpF[k] = pe;	BpF[k+1] = po;
			k+=2;
Q			yl+= (LMAX-LTR);
V			dyl+= (LMAX-LTR);
		}
QE		if ( ((LTR&1)==0) && ((m&1)==0) ) {	// k=LTR, m even
QE			BrF[k] = yl[0] * Ql[k];
QE			k++;
QE		}
V		if ( ((LTR&1)==0) && ((m&1)==1) ) {	// k=LTR, m odd
SO			BtF[k] =   dyl[1].t * Sl[k];
TE			BpF[k] = - dyl[1].t * Tl[k];
V			k++;
V		}
		while (k<NLAT) {
Q			BrF[k] = 0.0;
V			BtF[k] = 0.0;	BpF[k] = 0.0;
			k++;
		}
		im++;
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
	}
	while(im<=MTR) {	// regular for MTR_DCT < im <= MTR
		m = im*MRES;
Q		Ql = &Qlm[LiM(0,im)];	// virtual pointer for l=0 and im
S		Sl = &Slm[LiM(0,im)];	// virtual pointer for l=0 and im
T		Tl = &Tlm[LiM(0,im)];
		k=0;
		while (k<tm[im]) {	// polar optimization
Q			BrF[k] = 0.0;
QB			BrF[NLAT-tm[im] + k] = 0.0;	// south pole zeroes <=> BrF[im*NLAT + NLAT-(k+1)] = 0.0;
V			BtF[k] = 0.0;
VB			BtF[NLAT-tm[im] + k] = 0.0;	// south pole zeroes
V			BpF[k] = 0.0;
VB			BpF[NLAT-tm[im] + k] = 0.0;	// south pole zeroes
			k++;
		}
Q		yl  = ylm[im] + k*(LMAX-m+1);
V		dyl = dylm[im] + k*(LMAX-m+1);
		while (k < NLAT_2) {	// ops : NLAT_2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
Q			re = 0.0;	ro = 0.0;
S			dte = 0.0;	dto = 0.0;	pe = 0.0;	po = 0.0;
T			dpe = 0.0;	dpo = 0.0;	te = 0.0;	to = 0.0;
			while (l<LTR) {	// compute even and odd parts
QE				re  += yl[0] * Ql[l];		// re += ylm[im][k*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
QO				ro  += yl[1] * Ql[l+1];	// ro += ylm[im][k*(LMAX-m+1) + (l+1-m)] * Qlm[LiM(l+1,im)];
TO				dpo += dyl[0].t * Tl[l];
TO				te  += dyl[0].p * Tl[l];
SE				dto += dyl[0].t * Sl[l];
SE				pe  += dyl[0].p * Sl[l];
TE				dpe += dyl[1].t * Tl[l+1];
TE				to  += dyl[1].p * Tl[l+1];
SO				dte += dyl[1].t * Sl[l+1];
SO				po  += dyl[1].p * Sl[l+1];
				l+=2;
Q				yl+=2;
V				dyl+=2;
			}
			if (l==LTR) {
QE				re  += yl[0] * Ql[l];		// re += ylm[im][k*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
TO				dpo += dyl[0].t * Tl[l];
TO				te  += dyl[0].p * Tl[l];
SE				dto += dyl[0].t * Sl[l];
SE				pe  += dyl[0].p * Sl[l];
Q				yl++;
V				dyl++;
			}
Q			BrF[k] = re + ro;
V			BtF[k] = 		// Bt = dS/dt       + I.m/sint *T
S					(dte+dto)
T					+ I*(te+to)
V					;
V			BpF[k] = 		// Bp = I.m/sint * S - dT/dt
S					I*(pe+po)
T					- (dpe+dpo)
V					;
			k++;
QB			BrF[NLAT-k] = re - ro;
VB			BtF[NLAT-k] =
SB					(dte-dto)
TB					+ I*(te-to)
VB					;
VB			BpF[NLAT-k] =
SB					I*(pe-po)
TB					- (dpe-dpo)
VB					;
Q			yl  += (LMAX-LTR);
V			dyl += (LMAX-LTR);
		}
		im++;
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
	}
	for (k=0; k < NLAT*(NPHI/2 -MTR); k++) {	// padding for high m's
Q			BrF[k] = 0.0;
V			BtF[k] = 0.0;	BpF[k] = 0.0;
	}
Q	BrF -= NLAT*(MTR+1);		// restore original pointer
V	BtF -= NLAT*(MTR+1);	BpF -= NLAT*(MTR+1);	// restore original pointer
  #endif
	if (MTR_DCT >= 0) {
Q		fftw_execute_r2r(idct,(double *) BrF, (double *) BrF);		// iDCT
V		fftw_execute_r2r(idct,(double *) BtF, (double *) BtF);		// iDCT
V		fftw_execute_r2r(idct,(double *) BpF, (double *) BpF);		// iDCT
  #if NPHI>1
		if (MRES & 1) {		// odd m's must be multiplied by sin(theta) which was removed from ylm's
Q			for (im=1; im<=MTR_DCT; im+=2) {	// odd m's
Q				for (k=0; k<NLAT; k++) BrF[im*NLAT + k] *= st[k];
Q			}
V			for (im=0; im<=MTR_DCT; im+=2) {	//even m's
V				for (k=0; k<NLAT; k++) {
V					BtF[im*NLAT + k] *= st[k];	BpF[im*NLAT + k] *= st[k];
V				}
V			}
		} else {	// only even m's
V			for (im=0; im<=MTR_DCT; im++) {
V				for (k=0; k<NLAT; k++) {
V					BtF[im*NLAT + k] *= st[k];	BpF[im*NLAT + k] *= st[k];
V				}
V			}
		}
  #else
V		for (k=0; k<NLAT; k++) {
V			BT0[k] *= st[k];	BP0[k] *= st[k];
V		}
  #endif
	}
  #if NPHI>1
Q	fftw_execute_dft_c2r(ifft, BrF, (double *) BrF);
V	fftw_execute_dft_c2r(ifft, BtF, (double *) BtF);
V	fftw_execute_dft_c2r(ifft, BpF, (double *) BpF);
  #endif
# }