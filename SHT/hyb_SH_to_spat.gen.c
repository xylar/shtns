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
#Q void SH_to_spat(complex double *Qlm, double *Vr)
#V void SHsphtor_to_spat_dct(complex double *Slm, complex double *Tlm, double *Vt, double *Vp)
# {
Q	complex double *Ql;
S	complex double *Sl;
T	complex double *Tl;
Q	complex double *BrF;
  #ifndef SHT_AXISYM
V	complex double *BtF, *BpF;
V	struct DtDp *dyl;
Q	complex double re,ro;
V	complex double te,to, pe,po;
S	complex double dte, dto;	// spheroidal even and odd parts
T	complex double dpe, dpo;	// toroidal ...
Q	#define BR0 BrF
V	#define BT0 BtF
V	#define BP0 BpF
  #else
S	complex double *BtF;
T	complex double *BpF;
Q	double re,ro;
S	double te,to;
T	double pe,po;
Q	#define BR0 ((double *)BrF)
S	#define BT0 ((double *)BtF)
T	#define BP0 ((double *)BpF)
  #endif
Q	double *yl;
V	double *dyl0;
	long int llim;
	long int k,im,m,l;

	llim = LTR;	// copy LTR to a local variable for faster access (inner loop limit)
  #ifndef SHT_AXISYM
Q	BrF = (complex double *) Vr;
V	BtF = (complex double *) Vt;	BpF = (complex double *) Vp;
	if (SHT_FFT > 1) {		// alloc memory for the FFT
Q		BrF = fftw_malloc( (NPHI/2+1)*NLAT * sizeof(complex double) );
V		BtF = fftw_malloc( 2* (NPHI/2+1)*NLAT * sizeof(complex double) );
V		BpF = BtF + (NPHI/2+1)*NLAT;
	}
  #else
Q	BrF = (complex double *) Vr;
S	BtF = (complex double *) Vt;
T	BpF = (complex double *) Vp;
  #endif

	im=0;	m=0;
Q		Ql = Qlm;
S		Sl = Slm;
T		Tl = Tlm;
  #ifndef SHT_NO_DCT
	if (MTR_DCT >= 0) {	// dct for m=0.
Q		yl = ykm_dct[im];
V		dyl0 = (double *) dykm_dct[im];		// only theta derivative (d/dphi = 0 for m=0)
		k=0;
Q			re = 0.0;	ro = 0.0;
S			te = 0.0;	to = 0.0;
T			pe = 0.0;	po = 0.0;
		do {
			l = k;
			do {
QE				re += yl[0]  * (double) Ql[l];
QO				ro += yl[1]  * (double) Ql[l+1];
SE				to += dyl0[0] * (double) Sl[l];
SO				te += dyl0[1] * (double) Sl[l+1];
TO				po -= dyl0[0] * (double) Tl[l];
TE				pe -= dyl0[1] * (double) Tl[l+1];
				l+=2;
Q				yl+=2;
V				dyl0+=2;
			} while(l<llim);
			if (l==llim) {
QE				re += yl[0]  * (double) Ql[l];
SE				to += dyl0[0] * (double) Sl[l];
TO				po -= dyl0[0] * (double) Tl[l];
Q				yl++;
V				dyl0++;
			}
Q			BR0[k] = re;	BR0[k+1] = ro;
S			BT0[k] = te;	BT0[k+1] = to;
T			BP0[k] = pe;	BP0[k+1] = po;
			k+=2;
Q			yl+= (LMAX-LTR);
V			dyl0+= (LMAX-LTR);
Q			re = 0.0;	ro = 0.0;
S			to = 0.0;
T			po = 0.0;
SO			te = dyl0[1] * (double) Sl[k-1];
TE			pe = -dyl0[1] * (double) Tl[k-1];
V			dyl0+=2;
		} while (k<llim);
		if (k==llim) {
QE			re = yl[0] * Ql[k];
SE			to = dyl0[0] * (double) Sl[k];
TO			po = -dyl0[0] * (double) Tl[k];
		}
Q		BR0[k] = re;	BR0[k+1] = ro;
S		BT0[k] = te;	BT0[k+1] = to;
T		BP0[k] = pe;	BP0[k+1] = po;
		k+=2;
		while (k<NLAT) {	// dct padding (NLAT is even)
Q			BR0[k] = 0.0;	BR0[k+1] = 0.0;
S			BT0[k] = 0.0;	BT0[k+1] = 0.0;
T			BP0[k] = 0.0;	BP0[k+1] = 0.0;
			k+=2;
		}
    #ifdef SHT_AXISYM
Q		fftw_execute_r2r(idct_r1,Vr, Vr);		// iDCT m=0
S		fftw_execute_r2r(idct_r1,Vt, Vt);		// iDCT m=0
T		fftw_execute_r2r(idct_r1,Vp, Vp);		// iDCT m=0
V		k=0;	do {
S			Vt[k] *= st_1[k];
T		 	Vp[k] *= st_1[k];
V			k++;
V		} while (k<NLAT);
    #endif
	} else {
  #endif
		k=0;
Q		yl  = ylm[im];
V		dyl0 = (double *) dylm[im];	// only theta derivative (d/dphi = 0 for m=0)
		do {	// ops : NLAT_2 * [ (lmax-m+1) + 4]	: almost twice as fast.
			l=0;
Q			re = 0.0;	ro = 0.0;
S			te = 0.0;	to = 0.0;
T			pe = 0.0;	po = 0.0;
			do {	// compute even and odd parts
Q				re += yl[0] * (double) Ql[l];		// re += ylm[im][k*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
QB				ro += yl[1] * (double) Ql[l+1];	// ro += ylm[im][k*(LMAX-m+1) + (l+1-m)] * Qlm[LiM(l+1,im)];
TB				po += dyl0[0] * (double) Tl[l];	// m=0 : everything is real.
S				to += dyl0[0] * (double) Sl[l];
T				pe += dyl0[1] * (double) Tl[l+1];
SB				te += dyl0[1] * (double) Sl[l+1];
				l+=2;
Q				yl+=2;
V				dyl0+=2;
			} while (l<llim);
			if (l==llim) {
Q				re += yl[0] * (double) Ql[l];		// re += ylm[im][k*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
TB				po += dyl0[0] * (double) Tl[l];
S				to += dyl0[0] * (double) Sl[l];
Q				yl++;
V				dyl0++;
			}
Q			BR0[k] = re + ro;
S			BT0[k] = te + to;			// Bt = dS/dt
T			BP0[k] = -(pe+po);			// Bp = - dT/dt
			k++;
QB			BR0[NLAT-k] = re - ro;
SB			BT0[NLAT-k] = te - to;
TB			BP0[NLAT-k] = -(pe-po);
Q			yl  += (LMAX-LTR);
V			dyl0 += (LMAX-LTR);
		} while (k < NLAT_2);
  #ifndef SHT_NO_DCT
	}
  #endif

  #ifndef SHT_AXISYM
	im=1;
Q	BrF += NLAT;
V	BtF += NLAT;	BpF += NLAT;
    #ifndef SHT_NO_DCT
	while(im<=MTR_DCT) {		// dct for im <= MTR_DCT
		m=im*MRES;
Q		Ql = &Qlm[LiM(0,im)];		// virtual pointer for l=0 and im
S		Sl = &Slm[LiM(0,im)];
T		Tl = &Tlm[LiM(0,im)];
Q		yl = ykm_dct[im];
V		dyl = dykm_dct[im];
		k=0;	l=m;
		do {
Q			re = 0.0;	ro = 0.0;
V			te = 0.0;	to = 0.0;	pe = 0.0;	po = 0.0;
			while(l<llim) {
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
			if (l==llim) {
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
V			l = (k < m) ? m : k+(m&1);
			k+=2;
Q			yl+= (LMAX-LTR);
V			dyl+= (LMAX-LTR);
Q			l = (k < m) ? m : k-(m&1);
V		} while (k<=llim+1-(m&1));
Q		} while (k<=llim);
QE		if (l==llim) {
QE			BrF[k] = yl[0]   * Ql[l];
QE			BrF[k+1] = 0.0;
QE			k+=2;
QE		}
		while (k<NLAT) {		// NLAT even
Q			BrF[k] = 0.0;	BrF[k+1] = 0.0;
V			BtF[k] = 0.0;	BpF[k] = 0.0;
V			BtF[k+1] = 0.0;	BpF[k+1] = 0.0;
			k+=2;
		}
		im++;
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
	}
    #endif
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
		do {	// ops : NLAT_2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
Q			re = 0.0;	ro = 0.0;
S			dte = 0.0;	dto = 0.0;	pe = 0.0;	po = 0.0;
T			dpe = 0.0;	dpo = 0.0;	te = 0.0;	to = 0.0;
			while (l<llim) {	// compute even and odd parts
Q				re  += yl[0] * Ql[l];		// re += ylm[im][k*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
QB				ro  += yl[1] * Ql[l+1];	// ro += ylm[im][k*(LMAX-m+1) + (l+1-m)] * Qlm[LiM(l+1,im)];
TB				dpo += dyl[0].t * Tl[l];
TB				te  += dyl[0].p * Tl[l];
S				dto += dyl[0].t * Sl[l];
S				pe  += dyl[0].p * Sl[l];
T				dpe += dyl[1].t * Tl[l+1];
T				to  += dyl[1].p * Tl[l+1];
SB				dte += dyl[1].t * Sl[l+1];
SB				po  += dyl[1].p * Sl[l+1];
				l+=2;
Q				yl+=2;
V				dyl+=2;
			}
			if (l==llim) {
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
		} while (k < NLAT_2);
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

    if (NPHI>1) {
    #ifndef SHT_NO_DCT
		if (MTR_DCT >= 0) {
Q			fftw_execute_r2r(idct,(double *) BrF, (double *) BrF);		// iDCT
V			fftw_execute_r2r(idct,(double *) BtF, (double *) BtF);		// iDCT
V			fftw_execute_r2r(idct,(double *) BpF, (double *) BpF);		// iDCT
V			k=0;	do {		// m=0
V				((double *)BtF)[2*k] *= st_1[k];		((double *)BpF)[2*k] *= st_1[k];
V				k++;
V			} while(k<NLAT);
			if (MRES & 1) {		// odd m's must be multiplied by sin(theta) which was removed from ylm's
Q				for (im=1; im<=MTR_DCT; im+=2) {	// odd m's
Q					k=0;	do {	BrF[im*NLAT + k] *= st_1[k];	k++;	} while(k<NLAT);
Q				}
V				for (im=2; im<=MTR_DCT; im+=2) {	//even m's
V					k=0;	do {
V						BtF[im*NLAT + k] *= st_1[k];	BpF[im*NLAT + k] *= st_1[k];
V						k++;
V					} while(k<NLAT);
V				}
V			} else {	// only even m's
V				for (im=1; im<=MTR_DCT; im++) {
V					k=0;	do {
V						BtF[im*NLAT + k] *= st_1[k];	BpF[im*NLAT + k] *= st_1[k];
V						k++;
V					} while(k<NLAT);
V				}
			}
		}
    #endif
Q		fftw_execute_dft_c2r(ifft, BrF, Vr);
V		fftw_execute_dft_c2r(ifft, BtF, Vt);
V		fftw_execute_dft_c2r(ifft, BpF, Vp);
		if (SHT_FFT > 1) {		// free memory
Q		    fftw_free(BrF);
V		    fftw_free(BtF);	// this frees also BpF.
		}
    } else {
		k=1;	do {	// compress complex to real
Q			Vr[k] = (double) BrF[k];
V			Vt[k] = (double) BtF[k];
V			Vp[k] = (double) BpF[k];
			k++;
		} while(k<NLAT);
    #ifndef SHT_NO_DCT
		if (MTR_DCT >= 0) {
Q			fftw_execute_r2r(idct_r1,Vr, Vr);		// iDCT m=0
V			fftw_execute_r2r(idct_r1,Vt, Vt);		// iDCT m=0
V			fftw_execute_r2r(idct_r1,Vp, Vp);		// iDCT m=0
V			k=0;	do {
V				Vt[k] *= st_1[k]; 	Vp[k] *= st_1[k];
V				k++;
V			} while (k<NLAT);
		}
    #endif
    }
  #endif

Q	#undef BR0
V	#undef BT0
V	#undef BP0
# }
