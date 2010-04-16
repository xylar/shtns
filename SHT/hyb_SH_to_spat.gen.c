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
  #ifdef SHT_AXISYM
/// Only the \b axisymmetric (m=0) component is transformed, and the resulting spatial fields have size NLAT.
  #endif

/// Truncation and spatial discretization are defined by \ref shtns_set_size and \ref shtns_precompute.
Q/// \param[in] Qlm = spherical harmonics coefficients :
Q/// complex double arrays of size NLM [unmodified].
S/// \param[in] Slm = spherical harmonics coefficients of \b Spheroidal scalar :
S/// complex double array of size NLM [unmodified].
T/// \param[in] Tlm = spherical harmonics coefficients of \b Toroidal scalar :
T/// complex double array of size NLM [unmodified].
Q/// \param[out] Vr = spatial scalar field : double array.
  #ifndef SHT_AXISYM
V/// \param[out] Vt, Vp = (theta,phi)-components of spatial vector : double arrays.
  #else
S/// \param[out] Vt = theta-component of spatial vector : double array.
T/// \param[out] Vp = phi-component of spatial vector : double array.  
  #endif
  #ifdef SHT_VAR_LTR
/// \param[in] ltr = specify maximum degree of spherical harmonic. ltr must be at most LMAX, and all spherical harmonic degree higher than ltr are ignored. 
  #endif

// MTR_DCT : -1 => no dct
//            0 => dct for m=0 only
//            m => dct up to m, (!!! MTR_DCT <= MTR !!!)
#Q void SH_to_spat(complex double *Qlm, double *Vr)
#V void SHsphtor_to_spat(complex double *Slm, complex double *Tlm, double *Vt, double *Vp)
# {
Q	complex double *Ql;
S	complex double *Sl;
T	complex double *Tl;
Q	v2d *BrF;
  #ifndef SHT_AXISYM
V	v2d *BtF, *BpF;
V	struct DtDp *dyl;
Q	complex double re,ro;
V	complex double te,to, pe,po;
V	complex double dte,dto, dpe,dpo;
Q	#define BR0 ((complex double *)BrF)
V	#define BT0 ((complex double *)BtF)
V	#define BP0 ((complex double *)BpF)
  #else
S	v2d *BtF;
T	v2d *BpF;
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
Q	BrF = (v2d *) Vr;
V	BtF = (v2d *) Vt;	BpF = (v2d *) Vp;
	if (SHT_FFT > 1) {		// alloc memory for the FFT
Q		BrF = fftw_malloc( (NPHI/2+1)*NLAT * sizeof(complex double) );
V		BtF = fftw_malloc( 2* (NPHI/2+1)*NLAT * sizeof(complex double) );
V		BpF = BtF + (NPHI/2+1)*NLAT;
	}
  #else
Q	BrF = (v2d*) Vr;
S	BtF = (v2d*) Vt;
T	BpF = (v2d*) Vp;
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
QE			re = yl[0] * (double) Ql[k];
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
V		#ifdef _GCC_VEC_
V			v2d sin_1 = ((v2d *)st_1)[k];
S			((v2d *)Vt)[k] *= sin_1;
T		 	((v2d *)Vp)[k] *= sin_1;
V		#else
V			double sin_1 = st_1[2*k]; 	double sin_2 = st_1[2*k+1];
S			Vt[2*k] *= sin_1;		Vt[2*k+1] *= sin_2;
T			Vp[2*k] *= sin_1;		Vp[2*k+1] *= sin_2;
V		#endif
V			k++;
V		} while (k<NLAT_2);
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
Q		v2d* Ql = (v2d*) &Qlm[LiM(0,im)];		// virtual pointer for l=0 and im
S		v2d* Sl = (v2d*) &Slm[LiM(0,im)];
T		v2d* Tl = (v2d*) &Tlm[LiM(0,im)];
Q		yl = ykm_dct[im];
V		dyl = dykm_dct[im];
		k=0;	l=m;
		do {
Q			v2d re = vdup(0.0);		v2d ro = vdup(0.0);
V			v2d te = vdup(0.0);		v2d to = vdup(0.0);		v2d pe = vdup(0.0);		v2d po = vdup(0.0);
V			v2d dte = vdup(0.0);	v2d dto = vdup(0.0);	v2d dpe = vdup(0.0);	v2d dpo = vdup(0.0);
			while(l<llim) {
QE				re  += vdup(yl[0]) * Ql[l];
QO				ro  += vdup(yl[1]) * Ql[l+1];
TO				dte += vdup(dyl[0].p) * Tl[l];
TO				po  -= vdup(dyl[0].t) * Tl[l];
SE				dpe += vdup(dyl[0].p) * Sl[l];
SE				to  += vdup(dyl[0].t) * Sl[l];
SO				dpo += vdup(dyl[1].p) * Sl[l+1];
SO				te  += vdup(dyl[1].t) * Sl[l+1];
TE				dto += vdup(dyl[1].p) * Tl[l+1];
TE				pe  -= vdup(dyl[1].t) * Tl[l+1];
				l+=2;
Q				yl+=2;
V				dyl+=2;
			}
			if (l==llim) {
QE				re  += vdup(yl[0])   * Ql[l];
TO				dte += vdup(dyl[0].p) * Tl[l];
TO				po  -= vdup(dyl[0].t) * Tl[l];
SE				dpe += vdup(dyl[0].p) * Sl[l];
SE				to  += vdup(dyl[0].t) * Sl[l];
Q				yl++;
V				dyl++;
			}
Q			BrF[k] = re;	BrF[k+1] = ro;
V			BtF[k] = addi(te, dte);		BtF[k+1] = addi(to, dto);
V			BpF[k] = addi(pe, dpe);		BpF[k+1] = addi(po, dpo);
V			l = (k < m) ? m : k+(m&1);
			k+=2;
Q			yl+= (LMAX-LTR);
V			dyl+= (LMAX-LTR);
Q			l = (k < m) ? m : k-(m&1);
V		} while (k<=llim+1-(m&1));
Q		} while (k<=llim);
QE		if (l==llim) {
QE			BrF[k] = vdup(yl[0])   * Ql[l];
QE			BrF[k+1] = vdup(0.0);
QE			k+=2;
QE		}
		while (k<NLAT) {		// NLAT even
Q			BrF[k] = vdup(0.0); 	BrF[k+1] = vdup(0.0);
V			BtF[k] = vdup(0.0);		BpF[k] = vdup(0.0);
V			BtF[k+1] = vdup(0.0);	BpF[k+1] = vdup(0.0);
			k+=2;
		}
		im++;
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
	}
    #endif
	while(im<=MTR) {	// regular for MTR_DCT < im <= MTR
		m = im*MRES;
Q		v2d* Ql = (v2d*) &Qlm[LiM(0,im)];	// virtual pointer for l=0 and im
S		v2d* Sl = (v2d*) &Slm[LiM(0,im)];	// virtual pointer for l=0 and im
T		v2d* Tl = (v2d*) &Tlm[LiM(0,im)];
		k=0;
		while (k<tm[im]) {	// polar optimization
Q			BrF[k] = vdup(0.0);
QB			BrF[NLAT-tm[im] + k] = vdup(0.0);	// south pole zeroes <=> BrF[im*NLAT + NLAT-(k+1)] = 0.0;
V			BtF[k] = vdup(0.0);		BpF[k] = vdup(0.0);
VB			BtF[NLAT-tm[im] + k] = vdup(0.0);		BpF[NLAT-tm[im] + k] = vdup(0.0);	// south pole zeroes
			k++;
		}
Q		yl  = ylm[im];
V		dyl = dylm[im];
		do {	// ops : NLAT_2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
Q			v2d re = vdup(0.0); 	v2d ro = vdup(0.0);
V			v2d dte = vdup(0.0); 	v2d dto = vdup(0.0); 	v2d pe = vdup(0.0); 	v2d po = vdup(0.0);
V			v2d dpe = vdup(0.0); 	v2d dpo = vdup(0.0); 	v2d te = vdup(0.0); 	v2d to = vdup(0.0);
			while (l<llim) {	// compute even and odd parts
Q				re  += vdup(yl[0]) * Ql[l];		// re += ylm[im][k*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
QB				ro  += vdup(yl[1]) * Ql[l+1];	// ro += ylm[im][k*(LMAX-m+1) + (l+1-m)] * Qlm[LiM(l+1,im)];
TB				dpo -= vdup(dyl[0].t) * Tl[l];
TB				te  += vdup(dyl[0].p) * Tl[l];
S				dto += vdup(dyl[0].t) * Sl[l];
S				pe  += vdup(dyl[0].p) * Sl[l];
T				dpe -= vdup(dyl[1].t) * Tl[l+1];
T				to  += vdup(dyl[1].p) * Tl[l+1];
SB				dte += vdup(dyl[1].t) * Sl[l+1];
SB				po  += vdup(dyl[1].p) * Sl[l+1];
				l+=2;
Q				yl+=2;
V				dyl+=2;
			}
			if (l==llim) {
QE				re  += vdup(yl[0]) * Ql[l];		// re += ylm[im][k*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
TO				dpo -= vdup(dyl[0].t) * Tl[l];
TO				te  += vdup(dyl[0].p) * Tl[l];
SE				dto += vdup(dyl[0].t) * Sl[l];
SE				pe  += vdup(dyl[0].p) * Sl[l];
Q				yl++;
V				dyl++;
			}
Q			BrF[k] = re + ro;
V			BtF[k] = addi(dte+dto, te+to);		// Bt = dS/dt       + I.m/sint *T
V			BpF[k] = addi(dpe+dpo, pe+po);		// Bp = I.m/sint * S - dT/dt
			k++;
QB			BrF[NLAT-k] = re - ro;
VB			BtF[NLAT-k] = addi(dte-dto, te-to);
VB			BpF[NLAT-k] = addi(dpe-dpo, pe-po);
Q			yl  += (LMAX-LTR);
V			dyl += (LMAX-LTR);
		} while (k < NLAT_2);
		im++;
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
	}
	for (k=0; k < NLAT*(NPHI/2 -MTR); k++) {	// padding for high m's
Q			BrF[k] = vdup(0.0);
V			BtF[k] = vdup(0.0);	BpF[k] = vdup(0.0);
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
V				double sin_1 = st_1[k];		double sin_2 = st_1[k+1];
V				((double *)BtF)[2*k] *= sin_1;		((double *)BpF)[2*k] *= sin_1;
V				((double *)BtF)[2*k+2] *= sin_2;		((double *)BpF)[2*k+2] *= sin_2;
V				k+=2;
V			} while(k<NLAT);
Q			if (MRES & 1) {		// odd m's must be divided by sin(theta)
Q				for (im=1; im<=MTR_DCT; im+=2) {	// odd m's
Q					k=0;	do {
Q						((v2d *)BrF)[im*NLAT + k] *= vdup(st_1[k]);		((v2d *)BrF)[im*NLAT + k+1] *= vdup(st_1[k+1]);
Q						k+=2;
Q					} while(k<NLAT);
Q				}
Q			}
V			l = (MRES & 1) + 1;		// im-stride (l=1 if MRES even, l=2 if MRES odd)
V			for (im=l; im<=MTR_DCT; im+=l) {	//even m's must be divided by sin(theta)
V				k=0;	do {
V					s2d sin_1 = vdup(st_1[k]);		s2d sin_2 = vdup(st_1[k+1]);
V					((v2d *)BtF)[im*NLAT + k] *= sin_1;		((v2d *)BpF)[im*NLAT + k] *= sin_1;
V					((v2d *)BtF)[im*NLAT + k+1] *= sin_2;	((v2d *)BpF)[im*NLAT + k+1] *= sin_2;
V					k+=2;
V				} while(k<NLAT);
V			}
		}
    #endif
Q		fftw_execute_dft_c2r(ifft, (complex double *) BrF, Vr);
V		fftw_execute_dft_c2r(ifft, (complex double *) BtF, Vt);
V		fftw_execute_dft_c2r(ifft, (complex double *) BpF, Vp);
		if (SHT_FFT > 1) {		// free memory
Q		    fftw_free(BrF);
V		    fftw_free(BtF);	// this frees also BpF.
		}
    } else {
		k=1;	do {	// compress complex to real
Q			Vr[k] = ((double *)BrF)[2*k];
V			Vt[k] = ((double *)BtF)[2*k];
V			Vp[k] = ((double *)BpF)[2*k];
			k++;
		} while(k<NLAT);
    #ifndef SHT_NO_DCT
		if (MTR_DCT >= 0) {
Q			fftw_execute_r2r(idct_r1,Vr, Vr);		// iDCT m=0
V			fftw_execute_r2r(idct_r1,Vt, Vt);		// iDCT m=0
V			fftw_execute_r2r(idct_r1,Vp, Vp);		// iDCT m=0
V			k=0;	do {
V		#ifdef _GCC_VEC_
V				v2d sin_1 = ((v2d *)st_1)[k];
V				((v2d *)Vt)[k] *= sin_1; 	((v2d *)Vp)[k] *= sin_1;
V		#else
V			double sin_1 = st_1[2*k]; 	double sin_2 = st_1[2*k+1];
V			Vt[2*k] *= sin_1;		Vt[2*k+1] *= sin_2;
V			Vp[2*k] *= sin_1;		Vp[2*k+1] *= sin_2;
V		#endif
V				k++;
V			} while (k<NLAT_2);
		}
    #endif
    }
  #endif

Q	#undef BR0
V	#undef BT0
V	#undef BP0
# }
