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
#Q void spat_to_SH_dct(complex double *BrF, complex double *Qlm)
#V void SHsphtor_to_spat_dct(complex double *Slm, complex double *Tlm, complex double *BtF, complex double *BpF)
# {
Q	complex double *Ql;
S	complex double *Sl;
T	complex double *Tl;
Q	double *yl;
V	struct DtDp *dyl;
  #if NPHI > 1
Q	#define BR0 BrF
V	#define BT0 BtF
V	#define BP0 BpF
  #else
Q	#define BR0 ((double *)BrF)
V	#define BT0 ((double *)BtF)
V	#define BP0 ((double *)BpF)
  #endif
	long int k,im,m,l;

	im=0;	m=0;
Q		Ql = &Qlm[LiM(0,im)];		// virtual pointer for l=0 and im
S		Sl = &Slm[LiM(0,im)];
T		Tl = &Tlm[LiM(0,im)];

Q		yl = ylm_dct[im];
V		dyl = dylm_dct[im];
		for (k=0; k<NLAT; k++) {
Q			BR0[k] = 0.0;		// zero out array (includes DCT padding)
V			BT0[k] = 0.0;	BP0[k] = 0.0;	// zero out array (includes DCT padding)
		}
		for (l=m; l<LTR; l+=2) {
			for (k=0; k<=l; k+=2) {
Q				BR0[k]   += yl[k]   * (double) Ql[l];
Q				BR0[k+1] += yl[k+1] * (double) Ql[l+1];
S				BT0[k]   += dyl[k].t   * (double) Sl[l+1];
T				BP0[k]   -= dyl[k].t   * (double) Tl[l+1];
S				BT0[k+1] += dyl[k+1].t * (double) Sl[l];
T				BP0[k+1] -= dyl[k+1].t * (double) Tl[l];
			}
Q			yl += l+2;	//(l+2 - (m&1));
V			dyl += l+2;	//(l+1 + (m&1));
		}
		if (l==LTR) {
			for (k=0; k<=l; k+=2) {
Q				BR0[k]   += yl[k]   * (double) Ql[l];
S				BT0[k+1] += dyl[k+1].t * (double) Sl[l];
T				BP0[k+1] -= dyl[k+1].t * (double) Tl[l];
			}
		}
  #if NPHI > 1
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
	for (im=1;im<=MTR;im++) {
		m=im*MRES;
Q		Ql = &Qlm[LiM(0,im)];		// virtual pointer for l=0 and im
S		Sl = &Slm[LiM(0,im)];
T		Tl = &Tlm[LiM(0,im)];

Q		yl = ylm_dct[im];
V		dyl = dylm_dct[im];
		for (k=0; k<NLAT; k++) {
Q			BrF[k] = 0.0;		// zero out array (includes DCT padding)
V			BtF[k] = 0.0;	BpF[k] = 0.0;	// zero out array (includes DCT padding)
		}
		for (l=m; l<LTR; l+=2) {
			for(k=0; k<=l; k+=2) {
Q				BrF[k]   += yl[k]   * Ql[l];
Q				BrF[k+1] += yl[k+1] * Ql[l+1];
V				BtF[k]   += 
T					  dyl[k].p   * (I*Tl[l])
S					+ dyl[k].t   * Sl[l+1]
V					;
V				BpF[k]   += 
S					  dyl[k].p   * (I*Sl[l])
T					- dyl[k].t   * Tl[l+1]
V					;
V				BtF[k+1] += 
T					  dyl[k+1].p * (I*Tl[l+1])
S					+ dyl[k+1].t * Sl[l]
V					;
V				BpF[k+1] +=
S					  dyl[k+1].p * (I*Sl[l+1])
T					- dyl[k+1].t * Tl[l]
V					;
			}
V			if (k<=l+1) {
S				BtF[k]   += dyl[k].t   * Sl[l+1];
T				BpF[k]   -= dyl[k].t   * Tl[l+1];
V			}
Q			yl += (l+2 - (m&1));
V			dyl += l+2;	//(l+1 + (m&1));
		}
		if (l==LTR) {
			for (k=0; k<=l; k+=2) {
Q				BrF[k]   += yl[k]   * Ql[l];
T				BtF[k]   += dyl[k].p   * (I*Tl[l]);
S				BpF[k]   += dyl[k].p   * (I*Sl[l]);
S				BtF[k+1] += dyl[k+1].t * Sl[l];
T				BpF[k+1] -= dyl[k+1].t * Tl[l];
			}
		}
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
	}
	for (k=0; k < NLAT*(NPHI/2 -MTR); k++) {	// FFT padding for high m's
Q		BrF[k] = 0.0;
V		BtF[k] = 0.0;	BpF[k] = 0.0;
	}
Q	BrF -= NLAT*(MTR+1);		// restore original pointer
V	BtF -= NLAT*(MTR+1);	BpF -= NLAT*(MTR+1);	// restore original pointer
  #endif
Q	fftw_execute_r2r(idct,(double *) BrF, (double *) BrF);		// iDCT
V	fftw_execute_r2r(idct,(double *) BtF, (double *) BtF);		// iDCT
V	fftw_execute_r2r(idct,(double *) BpF, (double *) BpF);		// iDCT
  #if NPHI>1
	if (MRES & 1) {		// odd m's must be multiplied by sin(theta) which was removed from ylm's
Q		for (im=1; im<=MTR; im+=2) {	// odd m's
Q			for (k=0; k<NLAT; k++) BrF[im*NLAT + k] *= st[k];
Q		}
V		for (im=0; im<=MTR; im+=2) {	//even m's
V			for (k=0; k<NLAT; k++) {
V				BtF[im*NLAT + k] *= st[k];	BpF[im*NLAT + k] *= st[k];
V			}
V		}
	} else {	// only even m's
V		for (im=0; im<=MTR; im++) {
V			for (k=0; k<NLAT; k++) {
V				BtF[im*NLAT + k] *= st[k];	BpF[im*NLAT + k] *= st[k];
V			}
V		}
	}
Q	fftw_execute_dft_c2r(ifft, BrF, (double *) BrF);
V	fftw_execute_dft_c2r(ifft, BtF, (double *) BtF);
V	fftw_execute_dft_c2r(ifft, BpF, (double *) BpF);
  #else
V	for (k=0; k<NLAT; k++) {
V		BT0[k] *= st[k];	BP0[k] *= st[k];
V	}
  #endif
# }
