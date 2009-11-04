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
//          complex double array of size NLAT_2*(NPHI/2+1) or double array of size NLAT_2*(NPHI/2+1)*2
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
V	complex double *BtF, *BpF;
Q	double *yl;
V	double *dyl0;
  #ifndef SHT_AXISYM
V	struct DtDp *dyl;
Q	complex double re;
V	complex double to, pe, dto, dpe;
Q	#define BR0 BrF
V	#define BT0 BtF
V	#define BP0 BpF
  #else
Q	double re;
V	double to, pe;
Q	#define BR0 ((double *)BrF)
V	#define BT0 ((double *)BtF)
V	#define BP0 ((double *)BpF)
  #endif
V  	double sgn_sph, sgn_tor;
	long int llim;
	long int k,im,m,l;

	llim = LTR;	// copy LTR to a local variable for faster access (inner loop limit)
Q	BrF = (complex double *) Vr;
V	BtF = (complex double *) Vt;	BpF = (complex double *) Vp;
V	sgn_tor = -1.0;		sgn_sph = 1.0;

  #ifndef SHT_AXISYM
	if (SHT_FFT > 1) {		// alloc memory for the FFT
Q		BrF = fftw_malloc( (NPHI/2+1)*NLAT_2 * sizeof(complex double) );
V		BtF = fftw_malloc( 2* (NPHI/2+1)*NLAT_2 * sizeof(complex double) );
V		BpF = BtF + (NPHI/2+1)*NLAT_2;
	}
  #endif

V	if (parity) {	// antisymmetric scalar or vector.
V		Sl = Slm;		Slm = Tlm;		Tlm = Sl;		// exchange Slm and Tlm pointers.
V		Sl = BtF;		BtF = BpF;		BpF = Sl;		// exchange BtF and BpF pointers.
V		sgn_tor = 1.0;		sgn_sph = -1.0;
V	}

	im=0;	m=0;
Q		Ql = Qlm + parity;
S		Sl = Slm;
T		Tl = Tlm;
		k=0;
Q		yl  = ylm[im] + parity;
V		dyl0 = (double *) dylm[im];
		do {	// ops : 
			l = 0;
Q			re = 0.0;
V			to = 0.0;		pe = 0.0;
			do {	// compute even and odd parts
Q				re += yl[0] * (double) Ql[l];		// re += ylm[im][k*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
S				to += dyl0[0] * (double) Sl[l];
T				pe += dyl0[1] * (double) Tl[l+1];
				l+=2;
Q				yl+=2;
V				dyl0+=2;
			} while (l<llim);
			if (l==llim) {
Q				if (parity == 0)
Q					re += yl[0] * (double) Ql[l];		// re += ylm[im][k*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
S				to += dyl0[0] * (double) Sl[l];
Q				yl++;
V				dyl0++;
			}
Q			BR0[k] = re;
V			BT0[k] = to * sgn_sph;	// Bt = dS/dt
V			BP0[k] = pe * sgn_tor;	// Bp = - dT/dt
			k++;
Q			yl  += (LMAX-LTR);
V			dyl0 += (LMAX-LTR);
		} while (k < NLAT_2);

  #ifndef SHT_AXISYM
	im=1;
Q	BrF += NLAT_2;
V	BtF += NLAT_2;	BpF += NLAT_2;
	while(im<=MTR) {	// regular for MTR_DCT < im <= MTR
		m = im*MRES;
Q		Ql = &Qlm[LiM(0,im)] + parity;	// virtual pointer for l=0 and im
V		Sl = &Slm[LiM(0,im)];	// virtual pointer for l=0 and im
V		Tl = &Tlm[LiM(0,im)];
		k=0;
		while (k<tm[im]) {	// polar optimization
Q			BrF[k] = 0.0;
V			BtF[k] = 0.0;
V			BpF[k] = 0.0;
			k++;
		}
Q		yl  = ylm[im + parity*(MMAX+1)] + k*(LMAX-m+1) + parity;
V		dyl = dylm[im] + k*(LMAX-m+1);
		do {	// ops : NLAT_2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
Q			re = 0.0;
V			dto = 0.0;	pe = 0.0;	dpe = 0.0;	to = 0.0;
			while (l<llim) {	// compute even and odd parts
Q				re  += yl[0] * Ql[l];		// re += ylm[im][k*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
V				dto += dyl[0].t * Sl[l];
V				pe  += dyl[0].p * Sl[l];
V				dpe += dyl[1].t * Tl[l+1];
V				to  += dyl[1].p * Tl[l+1];
				l+=2;
Q				yl+=2;
V				dyl+=2;
			}
			if (l==llim) {
Q				if (parity == 0)
Q					re  += yl[0] * Ql[l];		// re += ylm[im][k*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
V				dto += dyl[0].t * Sl[l];
V				pe  += dyl[0].p * Sl[l];
Q				yl++;
V				dyl++;
			}
Q			BrF[k] = re;
V			BtF[k] = I*to + dto*sgn_sph;		// Bt = I.m/sint *T + dS/dt
V			BpF[k] = I*pe + dpe*sgn_tor;		// Bp = I.m/sint *S - dT/dt
			k++;
Q			yl  += (LMAX-LTR);
V			dyl += (LMAX-LTR);
		} while (k < NLAT_2);
		im++;
Q		BrF += NLAT_2;
V		BtF += NLAT_2;	BpF += NLAT_2;
	}
	for (k=0; k < NLAT_2*(NPHI/2 -MTR); k++) {	// padding for high m's
Q			BrF[k] = 0.0;
V			BtF[k] = 0.0;	BpF[k] = 0.0;
	}
Q	BrF -= NLAT_2*(MTR+1);		// restore original pointer
V	BtF -= NLAT_2*(MTR+1);	BpF -= NLAT_2*(MTR+1);	// restore original pointer
V	if (parity) {
V		Sl = BtF;		BtF = BpF;		BpF = Sl;		// restore original BtF and BpF pointers for FFTW and compression.
V	}

    if (NPHI>1) {
Q		fftw_execute_dft_c2r(ifft_eo, BrF, Vr);
V		fftw_execute_dft_c2r(ifft_eo, BtF, Vt);
V		fftw_execute_dft_c2r(ifft_eo, BpF, Vp);
		if (SHT_FFT > 1) {		// free memory
Q		    fftw_free(BrF);
V	    	fftw_free(BtF);	// this frees also BpF.
		}
    } else {
		k=1;	do {	// compress complex to real
Q			Vr[k] = (double) BrF[k];
V			Vt[k] = (double) BtF[k];
V			Vp[k] = (double) BpF[k];
			k++;
		} while(k<NLAT_2);
  #endif
  #ifndef SHT_AXISYM
    }
  #endif

Q	#undef BR0
V	#undef BT0
V	#undef BP0
# }
