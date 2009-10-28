# This file is meta-code for SHT.c (spherical harmonic transform).
# it is intended for "make" to generate C code for 3 similar SHT functions,
# (namely spat_to_SH [Q tag]), spat_to_SHsphtor [V tag], spat_to_SH3 [both Q&V tags])
# from one generic function + tags.
# Basically, there are tags at the beginning of lines (Q,V) that are information
# to keep or remove the line depending on the function to build. (Q for scalar, V for vector, # for comment)
#
//////////////////////////////////////////////////
//  Spherical Harmonics Transform
// input  : BrF, BtF, BpF = spatial/fourrier data : complex double array of size NLAT*(NPHI/2+1) or double array of size NLAT*(NPHI/2+1)*2
// output : Qlm, Slm, Tlm = spherical harmonics coefficients : complex double array of size NLM
#Q void spat_to_SH(double *Vr, complex double *Qlm)
#V void spat_to_SHsphtor(double *Vt, double *Vp, complex double *Slm, complex double *Tlm)
# {
Q	complex double *Ql;		// virtual pointers for given im
V	complex double *Sl, *Tl;	// virtual pointers for given im
Q	complex double *BrF;		// contains the Fourier transformed data
V	complex double *BtF, *BpF;	// contains the Fourier transformed data
Q	double *zl;
V	double *dzl0;
V	struct DtDp *dzl;
	long int ni;
	long int i,i0, im,l;
Q	complex double q0,q1;
V	complex double s0,t0,s1,t1;
  #ifndef SHT_AXISYM
QB	complex double reo[2*NLAT_2] SSE;	// symmetric (even) and anti-symmetric (odd) parts, interleaved.
VB	complex double tpeo[4*NLAT_2] SSE;	// theta and phi even and odd parts
Q	#define reo0 ((double*)reo)
V	#define tpeo0 ((double*)tpeo)
  #else
QB	double reo0[2*NLAT_2] SSE;	// symmetric (even) and anti-symmetric (odd) parts, interleaved.
VB	double tpeo0[4*NLAT_2] SSE;	// theta and phi even and odd parts
  #endif

// defines how to access even and odd parts of data
QB	#define re	reo[2*i]
QB	#define ro	reo[2*i+1]
VB	#define te	tpeo[4*i]
VB	#define to	tpeo[4*i+1]
VB	#define pe	tpeo[4*i+2]
VB	#define po	tpeo[4*i+3]
QB	#define re0	reo0[2*i]
QB	#define ro0	reo0[2*i+1]
VB	#define te0	tpeo0[4*i]
VB	#define to0	tpeo0[4*i+1]
VB	#define pe0	tpeo0[4*i+2]
VB	#define po0	tpeo0[4*i+3]
  #ifndef SHT_AXISYM
Q1	#define re	BrF[i]
Q1	#define ro	BrF[i]
V1	#define te	BtF[i]
V1	#define to	BtF[i]
V1	#define pe	BpF[i]
V1	#define po	BpF[i]
Q1	#define re0	BrF[i]
Q1	#define ro0	BrF[i]
V1	#define te0	BtF[i]
V1	#define to0	BtF[i]
V1	#define pe0	BpF[i]
V1	#define po0	BpF[i]
  #else
Q1	#define re	((double *)BrF)[i]
Q1	#define ro	((double *)BrF)[i]
V1	#define te	((double *)BtF)[i]
V1	#define to	((double *)BtF)[i]
V1	#define pe	((double *)BpF)[i]
V1	#define po	((double *)BpF)[i]
Q1	#define re0	((double *)BrF)[i]
Q1	#define ro0	((double *)BrF)[i]
V1	#define te0	((double *)BtF)[i]
V1	#define to0	((double *)BtF)[i]
V1	#define pe0	((double *)BpF)[i]
V1	#define po0	((double *)BpF)[i]
  #endif

	ni = NLAT_2;	// copy NLAT_2 to a local variable for faster access (inner loop limit)
Q    	BrF = (complex double *) Vr;
V    	BtF = (complex double *) Vt;	BpF = (complex double *) Vp;

  #ifndef SHT_AXISYM
	if (SHT_FFT > 0) {
	    if (SHT_FFT > 1) {		// alloc memory for the FFT
Q	    	BrF = fftw_malloc( (NPHI/2+1)*NLAT * sizeof(complex double) );
V	    	BtF = fftw_malloc( 2* (NPHI/2+1)*NLAT * sizeof(complex double) );
V	    	BpF = BtF + (NPHI/2+1)*NLAT;
	    }
Q	    fftw_execute_dft_r2c(fft,Vr, BrF);
V	    fftw_execute_dft_r2c(fft,Vt, BtF);
V	    fftw_execute_dft_r2c(fft,Vp, BpF);
	}
  #endif

	im = 0;		// dzl.p = 0.0 : and evrything is REAL
  #ifndef SHT_NO_DCT
	if (MTR_DCT >= 0) {
	#ifndef SHT_AXISYM
Q		#define BR0	((double *)reo)
V		#define BT0	((double *)tpeo)
V		#define BP0	((double *)tpeo + NLAT)
V		l = (NPHI==1) ? 1 : 2;		// stride of source data.
V		i=0;	i0=0;	do {
V			((double *)BtF)[i0] *= st_1[i]; 	((double *)BpF)[i0] *= st_1[i];
V			i++;	i0+=l;
V		} while (i<NLAT);
Q		fftw_execute_r2r(dct_m0,(double *) BrF, BR0);		// DCT out-of-place.
V		fftw_execute_r2r(dct_m0,(double *) BtF, BT0);		// DCT out-of-place.
V		fftw_execute_r2r(dct_m0,(double *) BpF, BP0);		// DCT out-of-place.
	#else
Q		#define BR0	((double *)BrF)
V		#define BT0	((double *)BtF)
V		#define BP0	((double *)BpF)
V		i=0;	do {
V			BT0[i] *= st_1[i]; 	BP0[i] *= st_1[i];
V			i++;
V		} while (i<NLAT);
Q		fftw_execute_r2r(dct_r1,(double *) BrF, (double *) BrF);	// DCT in-place.
V		fftw_execute_r2r(dct_r1,(double *) BtF, (double *) BtF);	// DCT in-place.
V		fftw_execute_r2r(dct_r1,(double *) BpF, (double *) BpF);	// DCT in-place.
	#endif
Q		l=0;
V		l=1;
Q		Ql = Qlm;		// virtual pointer for l=0 and im
Q		zl = zlm_dct0;
V		Sl = Slm;	Tl = Tlm;
V		dzl0 = dzlm_dct0;
V		Sl[0] = 0.0;	Tl[0] = 0.0;
#		qs0 = 0.0;	qs1 = 0.0;			// sum of first Ql's
		do {	// l has parity of m
Q			q0 = 0.0;	q1 = 0.0;
V			s0 = 0.0;	t1 = 0.0;	t0 = 0.0;	s1 = 0.0;
Q			i=l;	// l < NLAT
V			i=l-1;	// l > 0
#			qs0 += BR0[l];	qs1 += BR0[l+1];
#			q0 = qs0*zl[i];	q1 = qs1*zl[i+1];
			do {
Q				q0 += BR0[i]   * zl[0];
Q				q1 += BR0[i+1] * zl[1];
V				s0 += BT0[i]   * dzl0[0];
V				t0 -= BP0[i]   * dzl0[0];
V				s1 += BT0[i+1] * dzl0[1];
V				t1 -= BP0[i+1] * dzl0[1];
Q				zl+=2;
V				dzl0+=2;
				i+=2;
			} while(i<NLAT);
Q			Ql[l] = q0;	Ql[l+1] = q1;
V			Sl[l] = s0;	Sl[l+1] = s1;
V			Tl[l] = t0;	Tl[l+1] = t1;
			l+=2;
		} while(l<LTR);
		if (l == LTR) {
Q			q0 = 0.0;
V			s0 = 0.0;	t0 = 0.0;
Q			i=l;	// l < NLAT
V			i=l-1;
			do {
Q				q0 += BR0[i] * zl[0];
V				s0 += BT0[i] * dzl0[0];
V				t0 -= BP0[i] * dzl0[0];
Q				zl+=2;
V				dzl0+=2;
				i+=2;
			} while(i<NLAT);
Q			Ql[l] = q0;
V			Sl[l] = s0;	Tl[l] = t0;
			l++;
		}
Q	#undef BR0
V	#undef BT0
V	#undef BP0
  #ifdef SHT_VAR_LTR
		while( l<=LMAX ) {
Q			Ql[l] = 0.0;
V			Sl[l] = 0.0;	Tl[l] = 0.0;
			l++;
		}
  #endif
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
	} else {
  #endif
  #ifndef SHT_AXISYM
	    if (NPHI>1) {
		i=0;	l=NLAT-1;
 B		do {	// compute symmetric and antisymmetric parts.
QB			re0 = (double) BrF[i] + (double) BrF[l];
QB			ro0 = (double) BrF[i] - (double) BrF[l];
VB			te0 = (double) BtF[i] + (double) BtF[l];
VB			to0 = (double) BtF[i] - (double) BtF[l];
VB			pe0 = (double) BpF[i] + (double) BpF[l];
VB			po0 = (double) BpF[i] - (double) BpF[l];
 B			i++;	l--;
 B		} while(i<l);
    #ifndef SHT_NLAT_EVEN
 B		if (i == l) {		// NLAT is odd : special equator handling
QB			re0 = (double) BrF[i];	ro0 = 0.0;
VB			te0 = (double) BtF[i];	to0 = 0.0;
VB			pe0 = (double) BpF[i];	po0 = 0.0;
 B		}
    #endif
	    } else {
  #endif
		i=0;	l=NLAT-1;
 B		do {	// compute symmetric and antisymmetric parts.
QB			re0 = ((double*)BrF)[i] + ((double*)BrF)[l];
QB			ro0 = ((double*)BrF)[i] - ((double*)BrF)[l];
VB			te0 = ((double*)BtF)[i] + ((double*)BtF)[l];
VB			to0 = ((double*)BtF)[i] - ((double*)BtF)[l];
VB			pe0 = ((double*)BpF)[i] + ((double*)BpF)[l];
VB			po0 = ((double*)BpF)[i] - ((double*)BpF)[l];
 B			i++;	l--;
 B		} while(i<l);
    #ifndef SHT_NLAT_EVEN
 B		if (i == l) {		// NLAT is odd : special equator handling
QB			re0 = ((double*)BrF)[i];	ro0 = 0.0;
VB			te0 = ((double*)BtF)[i];	to0 = 0.0;
VB			pe0 = ((double*)BpF)[i];	po0 = 0.0;
 B		}
    #endif
  #ifndef SHT_AXISYM
	    }
  #endif
Q		l=0;
V		l=1;		// l=0 is zero for the vector transform.
Q		Ql = Qlm;		// virtual pointer for l=0 and im
V		Sl = Slm;	Tl = Tlm;		// virtual pointer for l=0 and im
Q		zl = zlm[im];
V		dzl0 = (double *) dzlm[im];		// only theta derivative (d/dphi = 0 for m=0)
QB		BrF += NLAT;
VB		BtF += NLAT;	BpF += NLAT;
V		Sl[0] = 0.0;	Tl[0] = 0.0;
		do {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
			i=0;
QE			q0 = 0.0;
QO			q1 = 0.0;
VE			s0 = 0.0;	t1 = 0.0;
VO			t0 = 0.0;	s1 = 0.0;
			do {
QE				q0 += zl[0] * re0;	// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
QO				q1 += zl[1] * ro0;	// Qlm[LiM(l+1,im)] += zlm[im][(l+1-m)*NLAT/2 + i] * fm[i];
VE				s0 += dzl0[0] * to0;
VO				t0 -= dzl0[0] * po0;
VO				s1 += dzl0[1] * te0;
VE				t1 -= dzl0[1] * pe0;
Q				zl +=2;
V				dzl0 +=2;
				i++;
			} while(i < ni);
QE			Ql[l] = q0;
QO			Ql[l+1] = q1;
VE			Sl[l] = s1;	Tl[l+1] = t0;
VO			Tl[l] = t1;	Sl[l+1] = s0;
			l+=2;
		} while (l<LTR);
		if (l==LMAX) {
QE			q0 = 0.0;
VE			s1 = 0.0;
VO			t1 = 0.0;
			i=0;	do {
QE				q0 += zl[0] * re0;		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
VE				s1 += dzl0[0] * te0;
VO				t1 -= dzl0[0] * pe0;
Q				zl ++;
V				dzl0 ++;
				i++;
			} while(i<ni);
QE			Ql[l] = q0;
VE			Sl[l] = s1;
VO			Tl[l] = t1;
  #ifdef SHT_VAR_LTR
		} else {
		    if (l==LTR) {
QE			q0 = 0.0;
VE			s1 = 0.0;
VO			t1 = 0.0;
			i=0;	do {
QE				q0 += zl[0] * re0;		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
VO				s1 += dzl0[1] * te0;
VE				t1 -= dzl0[1] * pe0;
Q				zl +=2;
V				dzl0 +=2;
				i++;
			} while(i<ni);
QE			Ql[l] = q0;
VE			Sl[l] = s1;
VO			Tl[l] = t1;
			l++;
		    }
		    while( l<=LMAX ) {
Q			Ql[l] = 0.0;
V			Sl[l] = 0.0;	Tl[l] = 0.0;
			l++;
		    }
  #endif
		}
  #ifndef SHT_NO_DCT
	}
  #endif
  #ifndef SHT_AXISYM
	for (im=1;im<=MTR;im++) {
		i0 = tm[im];
 B		i=i0;	l=NLAT-1-i;
 B		do {	// compute symmetric and antisymmetric parts.
QB			re = BrF[i] + BrF[l];	ro = BrF[i] - BrF[l];
VB			te = BtF[i] + BtF[l];	to = BtF[i] - BtF[l];
VB			pe = BpF[i] + BpF[l];	po = BpF[i] - BpF[l];
 B			i++;	l--;
 B		} while (i < l);
    #ifndef SHT_NLAT_EVEN
 B		if (i == l) {		// NLAT is odd : special equator handling
QB			re = BrF[i];	ro = 0.0;
VB			te = BtF[i];	to = 0.0;
VB			pe = BpF[i];	po = 0.0;
 B		}
    #endif
		l=im*MRES;
Q		Ql = &Qlm[LiM(0,im)];	// virtual pointer for l=0 and im
V		Sl = &Slm[LiM(0,im)];	Tl = &Tlm[LiM(0,im)];
Q		zl = zlm[im];
V		dzl = dzlm[im];
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
		while (l<LTR) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
Q			zl += 2*i0;
V			dzl += 2*i0;
QE			q0 = 0.0;
QO			q1 = 0.0;
VE			s0 = 0.0;	t1 = 0.0;		// Slm[LiM(l,im)] = 0.0;	Slm[LiM(l+1,im)] = 0.0;
VO			t0 = 0.0;	s1 = 0.0;
			i=i0;	do {		// tm[im] : polar optimization
QE				q0 += re * zl[0];		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
QO				q1 += ro * zl[1];	// Qlm[LiM(l+1,im)] += zlm[im][(l+1-m)*NLAT/2 + i] * fm[i];
VE				s0 += dzl[0].t *to - dzl[0].p *pe*I;		// ref: these E. Dormy p 72.
VO				t0 -= dzl[0].t *po + dzl[0].p *te*I;
VO				s1 += dzl[1].t *te - dzl[1].p *po*I;
VE				t1 -= dzl[1].t *pe + dzl[1].p *to*I;
Q				zl +=2;
V				dzl +=2;
				i++;
			} while (i < ni);
QE			Ql[l] = q0;
QO			Ql[l+1] = q1;
VE			Sl[l] = s0;	Tl[l+1] = t1;
VO			Tl[l] = t0;	Sl[l+1] = s1; 
			l+=2;
		}
		if (l==LMAX) {
Q			zl += i0;
V			dzl += i0;
QE			q0 = 0.0;	// Qlm[LiM(l,im)] = 0.0;
VE			s0 = 0.0;
VO			t0 = 0.0;
			i=i0;	do {		// tm[im] : polar optimization
QE				q0 += zl[0] * re;	// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
VE				s0 += dzl[0].t *to - dzl[0].p *pe*I;
VO				t0 -= dzl[0].t *po + dzl[0].p *te*I;
Q				zl++;
V				dzl++;
				i++;
			} while(i<ni);
QE			Ql[l] = q0;
VE			Sl[l] = s0;
VO			Tl[l] = t0;
    #ifdef SHT_VAR_LTR
		} else {
		    if (l==LTR) {
Q			zl += 2*i0;
V			dzl += 2*i0;
QE			q0 = 0.0;
VE			s0 = 0.0;
VO			t0 = 0.0;
			i=i0;	do {		// tm[im] : polar optimization
QE				q0 += re * zl[0];		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
VE				s0 += dzl[0].t *to - dzl[0].p *pe*I;		// ref: these E. Dormy p 72.
VO				t0 -= dzl[0].t *po + dzl[0].p *te*I;
Q				zl +=2;
V				dzl +=2;
				i++;
			} while (i<ni);
QE			Ql[l] = q0;
VE			Sl[l] = s0;
VO			Tl[l] = t0;
			l++;
		    }
		    while( l<=LMAX ) {
Q			Ql[l] = 0.0;
V			Sl[l] = 0.0;	Tl[l] = 0.0;
			l++;
		    }
    #endif
		}
	}

  	if (SHT_FFT > 1) {		// free memory
Q	    fftw_free(BrF - NLAT*(MTR+1));
V	    fftw_free(BtF - NLAT*(MTR+1));	// this frees also BpF.
	}
  #endif

Q	#undef re
Q	#undef ro
V	#undef te
V	#undef to
V	#undef pe
V	#undef po
Q	#undef re0
Q	#undef ro0
V	#undef te0
V	#undef to0
V	#undef pe0
V	#undef po0
Q	#undef reo0
V	#undef teo0
V	#undef peo0
Q	#undef reo0
V	#undef tpeo0
# }
