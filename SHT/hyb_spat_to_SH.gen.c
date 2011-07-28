/*
 * Copyright (c) 2010-2011 Centre National de la Recherche Scientifique.
 * written by Nathanael Schaeffer (CNRS, ISTerre, Grenoble, France).
 * 
 * nathanael.schaeffer@ujf-grenoble.fr
 * 
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software. You can use,
 * modify and/or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 * 
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * 
 */

# This file is meta-code for SHT.c (spherical harmonic transform).
# it is intended for "make" to generate C code for 3 similar SHT functions,
# (namely spat_to_SH [Q tag]), spat_to_SHsphtor [V tag], spat_to_SH3 [both Q&V tags])
# from one generic function + tags.
# Basically, there are tags at the beginning of lines (Q,V) that are information
# to keep or remove the line depending on the function to build. (Q for scalar, V for vector, # for comment)
#
//////////////////////////////////////////////////
  #ifdef SHT_AXISYM
/// The spatial field is assumed to be \b axisymmetric (spatial size NLAT), and only the m=0 harmonics are written to output.
  #endif

/// Truncation and spatial discretization are defined by \ref shtns_set_size and \ref shtns_precompute.
Q/// \param[in] Vr = spatial scalar field : double array.
V/// \param[in] Vt, Vp = spatial (theta, phi) vector components : double arrays.
Q/// \param[out] Qlm = spherical harmonics coefficients :
Q/// complex double arrays of size NLM.
V/// \param[out] Slm,Tlm = spherical harmonics coefficients of \b Spheroidal and \b Toroidal scalars :
V/// complex double arrays of size NLM.
  #ifdef SHT_VAR_LTR
/// \param[in] ltr = specify maximum degree of spherical harmonic. ltr must be at most LMAX, and all spherical harmonic degree higher than ltr are set to zero. 
  #endif

QX	void GEN3(spat_to_SH_,ID_NME,SUFFIX)(shtns_cfg shtns, double *Vr, complex double *Qlm SUPARG) {
VX	void GEN3(spat_to_SHsphtor_,ID_NME,SUFFIX)(shtns_cfg shtns, double *Vt, double *Vp, complex double *Slm, complex double *Tlm SUPARG) {
3	void GEN3(spat_to_SHqst_,ID_NME,SUFFIX)(shtns_cfg shtns, double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm SUPARG) {

Q	complex double *BrF;		// contains the Fourier transformed data
V	complex double *BtF, *BpF;	// contains the Fourier transformed data
Q	double *zl;
V	double *dzl0;
V	struct DtDp *dzl;
	long int ni, llim, imlim;
	long int i,i0, im,l;
Q	complex double q0,q1;
V	complex double s0,t0,s1,t1;
  #ifndef SHT_AXISYM
QB	v2d reo[2*NLAT_2];	// symmetric (even) and anti-symmetric (odd) parts, interleaved.
VB	v2d tpeo[4*NLAT_2];	// theta and phi even and odd parts
Q	#define reo0 ((double*)reo)
V	#define tpeo0 ((double*)tpeo)
  #else
QB	double reo0[2*NLAT_2] SSE;	// symmetric (even) and anti-symmetric (odd) parts, interleaved.
VB	double tpeo0[4*NLAT_2] SSE;	// theta and phi even and odd parts
  #endif

// defines how to access even and odd parts of data
QB	#define re(i)	reo[2*(i)]
QB	#define ro(i)	reo[2*(i)+1]
VB	#define te(i)	tpeo[4*(i)]
VB	#define to(i)	tpeo[4*(i)+1]
VB	#define pe(i)	tpeo[4*(i)+2]
VB	#define po(i)	tpeo[4*(i)+3]
QB	#define re0(i)	reo0[2*(i)+1]
QB	#define ro0(i)	reo0[2*(i)]
VB	#define te0(i)	tpeo0[4*(i)]
VB	#define to0(i)	tpeo0[4*(i)+1]
VB	#define pe0(i)	tpeo0[4*(i)+2]
VB	#define po0(i)	tpeo0[4*(i)+3]

	ni = NLAT_2;	// copy NLAT_2 to a local variable for faster access (inner loop limit)
	llim = LTR;		imlim = MTR;
	#ifdef SHT_VAR_LTR
		if (imlim*MRES > llim) imlim = llim/MRES;
	#endif
Q		BrF = (complex double *) Vr;
V		BtF = (complex double *) Vt;	BpF = (complex double *) Vp;

  #ifndef SHT_AXISYM
	if (SHT_FFT > 0) {
	    if (SHT_FFT > 1) {		// alloc memory for the FFT
	    	long int nspat = ((NPHI>>1) +1)*NLAT;
QX	    	BrF = fftw_malloc( nspat * sizeof(complex double) );
VX	    	BtF = fftw_malloc( 2* nspat * sizeof(complex double) );
VX	    	BpF = BtF + nspat;
3	    	BrF = fftw_malloc( 3* nspat * sizeof(complex double) );
3	    	BtF = BrF + nspat;		BpF = BtF + nspat;
	    }
Q	    fftw_execute_dft_r2c(shtns->fft,Vr, BrF);
V	    fftw_execute_dft_r2c(shtns->fft,Vt, BtF);
V	    fftw_execute_dft_r2c(shtns->fft,Vp, BpF);
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
V		double* st_1 = shtns->st_1;
V		i=0;	i0=0;	do {
V			double sin_1 = st_1[i];
V			((double *)BtF)[i0] *= sin_1; 	((double *)BpF)[i0] *= sin_1;
V			i++;	i0+=l;
V		} while (i<NLAT);
Q		fftw_execute_r2r(shtns->dct_m0,(double *) BrF, BR0);		// DCT out-of-place.
V		fftw_execute_r2r(shtns->dct_m0,(double *) BtF, BT0);		// DCT out-of-place.
V		fftw_execute_r2r(shtns->dct_m0,(double *) BpF, BP0);		// DCT out-of-place.
	#else
Q		#define BR0	((double *)BrF)
V		#define BT0	((double *)BtF)
V		#define BP0	((double *)BpF)
Q		if (NPHI > 1) {
Q			s2d np = vdup(NPHI);
Q			i=0;	do {
Q				((v2d*) BrF)[i] *= np;
Q				i++;
Q			} while (i<ni);
Q		}
V		i=0;	do {
V			s2d np = vdup(NPHI);
V		#ifdef _GCC_VEC_
V			v2d sin_1 = ((v2d *)st_1)[i] * np;
V			((v2d*) BtF)[i] *= sin_1; 	((v2d*) BpF)[i] *= sin_1;
V		#else
V			double sin_1 = st_1[2*i] * np;		double sin_2 = st_1[2*i+1] * np;
V			BT0[2*i] *= sin_1;		BT0[2*i+1] *= sin_2;
V			BP0[2*i] *= sin_1;		BP0[2*i+1] *= sin_2;
V		#endif
V			i++;
V		} while (i<ni);
Q		fftw_execute_r2r(shtns->dct_r1,(double *) BrF, (double *) BrF);	// DCT in-place.
V		fftw_execute_r2r(shtns->dct_r1,(double *) BtF, (double *) BtF);	// DCT in-place.
V		fftw_execute_r2r(shtns->dct_r1,(double *) BpF, (double *) BpF);	// DCT in-place.
	#endif
		long int klim = shtns->klim;
		l=0;
Q		v2d* Ql = (v2d*) Qlm;
V		v2d* Sl = (v2d*) Slm;	v2d* Tl = (v2d*) Tlm;
Q		zl = shtns->zlm_dct0;
V		dzl0 = shtns->dzlm_dct0;
V	#ifndef _GCC_VEC_
V		s1 = 0.0;	t1 = 0.0;		// l=0 : Sl = Tl = 0
V	#else
V		v2d s = vdup(0.0);		v2d t = vdup(0.0);		// l=0 : Sl = Tl = 0
V	#endif
	#ifdef SHT_VAR_LTR
		i = (llim * SHT_NL_ORDER) + 2;		// sum truncation
		if (i < klim) klim = i;
		while(l < llim) {
	#else
		do {		// l < LMAX
	#endif
			i=l;	// l < klim
	  #ifndef _GCC_VEC_
V			Sl[l] = s1;		Tl[l] = t1;
Q			q0 = 0.0;	q1 = 0.0;
V			s0 = 0.0;	t1 = 0.0;	t0 = 0.0;	s1 = 0.0;
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
			} while(i<klim);
		#ifdef SHT_VAR_LTR
Q			zl += (shtns->klim-i);
V			dzl0 += (shtns->klim-i);
		#endif
Q			Ql[l] = q0;		Ql[l+1] = q1;
V			Sl[l+1] = s0;	Tl[l+1] = t0;
	  #else
V			Sl[l] = vhi_to_cplx(s);		Tl[l] = vhi_to_cplx(t);
Q			v2d q = vdup(0.0);
V			s = vdup(0.0);		t = vdup(0.0);
			i >>= 1;	// i = i/2
			do {
Q				q += ((v2d*) zl)[0] * ((v2d*) BR0)[i];
V				s += ((v2d*) dzl0)[0] * ((v2d*) BT0)[i];
V				t -= ((v2d*) dzl0)[0] * ((v2d*) BP0)[i];
Q				zl +=2;
V				dzl0 +=2;
				i++;
			} while(2*i < klim);
		#ifdef SHT_VAR_LTR
Q			zl += (shtns->klim-2*i);
V			dzl0 += (shtns->klim-2*i);
		#endif
Q			Ql[l]   = vlo_to_cplx(q);		Ql[l+1] = vhi_to_cplx(q);
V			Sl[l+1] = vlo_to_cplx(s);		Tl[l+1] = vlo_to_cplx(t);
	  #endif
			l+=2;
	#ifndef SHT_VAR_LTR
		} while(l<llim);
	#else
		}
	#endif
		if (l == llim) {
V	#ifndef _GCC_VEC_
V			Sl[l] = s1;		Tl[l] = t1;
V	#else
V			((v2d*) Sl)[l] = vhi_to_cplx(s);		((v2d*) Tl)[l] = vhi_to_cplx(t);
V	#endif
Q			q0 = 0.0;
Q			i=l;	// l < klim
Q			do {
Q				q0 += BR0[i] * zl[0];
Q				zl+=2;
Q				i+=2;
Q			} while(i<klim);
Q			((complex double *) Ql)[l] = q0;
			l++;
		}
Q	#undef BR0
V	#undef BT0
V	#undef BP0
  #ifdef SHT_VAR_LTR
		while( l<=LMAX ) {
Q			Ql[l] = vdup(0.0);
V			Sl[l] = vdup(0.0);		Tl[l] = vdup(0.0);
			l++;
		}
  #endif
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
	} else {
  #endif
		i=0;
QE		double r0 = 0.0;
Q		zl = shtns->zlm[0];
  #ifndef SHT_AXISYM
		i0 = (NPHI==1) ? 1 : 2;		// stride of source data.
 B		do {	// compute symmetric and antisymmetric parts.
QB			double a = ((double*)BrF)[i*i0];		double b = ((double*)BrF)[(NLAT-1)*i0 -i*i0];
QB			ro0(i) = a-b;		re0(i) = a+b;
QB			r0 += zl[i] * (a+b);
VB			double c = ((double*)BtF)[i*i0];		double d = ((double*)BtF)[(NLAT-1)*i0 -i*i0];
VB			te0(i) = c+d;		to0(i) = c-d;
VB			double e = ((double*)BpF)[i*i0];		double f = ((double*)BpF)[(NLAT-1)*i0 -i*i0];
VB			pe0(i) = e+f;		po0(i) = e-f;
 B			i++;
 B		} while(i<ni);
  #else
 B		do {	// compute symmetric and antisymmetric parts.
 B			double np = NPHI;
QB			double a = ((double*)BrF)[i];		double b = ((double*)BrF)[NLAT-1-i];
QB			ro0(i) = (a-b)*np;		re0(i) = (a+b)*np;
QB			r0 += zl[i] * ((a+b)*np);
VB			double c = ((double*)BtF)[i];		double d = ((double*)BtF)[NLAT-1-i];
VB			te0(i) = (c+d)*np;		to0(i) = (c-d)*np;
VB			double e = ((double*)BpF)[i];		double f = ((double*)BpF)[NLAT-1-i];
VB			pe0(i) = (e+f)*np;		po0(i) = (e-f)*np;
 B			i++;
 B		} while(i<ni);
  #endif
Q		zl += ni + (ni&1);		// SSE alignement
		l=1;			// l=0 is zero for the vector transform.
Q		v2d* Ql = (v2d*) Qlm;		// virtual pointer for l=0 and im
V		v2d* Sl = (v2d*) Slm;	v2d* Tl = (v2d*) Tlm;		// virtual pointer for l=0 and im
V		dzl0 = (double *) shtns->dzlm[0];		// only theta derivative (d/dphi = 0 for m=0)
QB		BrF += NLAT;
VB		BtF += NLAT;	BpF += NLAT;
Q		((complex double *)Ql)[0] = r0;
V		Sl[0] = vdup(0.0);	Tl[0] = vdup(0.0);	// l=0 is zero for the vector transform.
	#ifdef SHT_VAR_LTR
		while (l<llim) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
	#else
		do {
	#endif
			i=0;
  #ifndef _GCC_VEC_
QE			double q0 = 0.0;
QO			double q1 = 0.0;
VE			double s0 = 0.0;	double t1 = 0.0;
VO			double t0 = 0.0;	double s1 = 0.0;
			do {
QE				q0 += zl[0] * ro0(i);	// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
QO				q1 += zl[1] * re0(i);	// Qlm[LiM(l+1,im)] += zlm[im][(l+1-m)*NLAT/2 + i] * fm[i];
VO				s0 += dzl0[0] * te0(i);
VE				t0 -= dzl0[0] * pe0(i);
VE				s1 += dzl0[1] * to0(i);
VO				t1 -= dzl0[1] * po0(i);
Q				zl +=2;
V				dzl0 +=2;
				i++;
			} while(i < ni);
QE			Ql[l] = q0;
QO			Ql[l+1] = q1;
VE			Sl[l] = s0;		Tl[l+1] = t1;
VO			Tl[l] = t0;		Sl[l+1] = s1;
  #else
Q			v2d q = vdup(0.0);
V			v2d s = vdup(0.0);		v2d t = vdup(0.0);
			do {
Q				q += ((v2d*) zl)[0] * ((v2d*) reo0)[i];
V				s += ((v2d*) dzl0)[0] * ((v2d*) tpeo0)[2*i];
V				t -= ((v2d*) dzl0)[0] * ((v2d*) tpeo0)[2*i+1];
Q				zl +=2;
V				dzl0 +=2;
				i++;
			} while(i < ni);
Q			Ql[l] = vlo_to_cplx(q);		Ql[l+1] = vhi_to_cplx(q);
V			Sl[l] = vlo_to_cplx(s);		Sl[l+1] = vhi_to_cplx(s);
V			Tl[l] = vlo_to_cplx(t);		Tl[l+1] = vhi_to_cplx(t);
  #endif
			l+=2;
	#ifndef SHT_VAR_LTR
		} while (l<llim);
	#else
		}
	#endif
		if (l==llim) {
			long int lstride=1;
	  #ifdef SHT_VAR_LTR
			if (l != LMAX) lstride=2;
	  #endif
QE			double q0 = 0.0;
VE			double s0 = 0.0;
VO			double t0 = 0.0;
			i=0;	do {
QE				q0 += zl[0] * ro0(i);		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
VE				s0 += dzl0[0] * te0(i);
VO				t0 -= dzl0[0] * pe0(i);
Q				zl += lstride;
V				dzl0 += lstride;
				i++;
			} while(i<ni);
QE			((complex double *)Ql)[l] = q0;
VE			((complex double *)Sl)[l] = s0;
VO			((complex double *)Tl)[l] = t0;
	  #ifdef SHT_VAR_LTR
	  		l++;
		}
	    while( l<=LMAX ) {
Q			Ql[l] = vdup(0.0);
V			Sl[l] = vdup(0.0);	Tl[l] = vdup(0.0);
			l++;
      #endif
		}
  #ifndef SHT_NO_DCT
	}
  #endif
  #ifndef SHT_AXISYM
	for (im=1;im<=imlim;im++) {
		i0 = shtns->tm[im];
 B		i=i0;
 B		do {	// compute symmetric and antisymmetric parts.
		  #ifdef SHT_3COMP
QB			s2d sin = vdup(shtns->st[i]);
QB			v2d q0 = ((v2d *)BrF)[i];	v2d q1 = ((v2d *)BrF)[NLAT-1-i];		re(i) = (q0+q1)*sin;	ro(i) = (q0-q1)*sin;
		  #else
QB			v2d q0 = ((v2d *)BrF)[i];	v2d q1 = ((v2d *)BrF)[NLAT-1-i];		re(i) = q0+q1;	ro(i) = q0-q1;		  
		  #endif
VB			v2d t0 = ((v2d *)BtF)[i];	v2d t1 = ((v2d *)BtF)[NLAT-1-i];		te(i) = t0+t1;	to(i) = t0-t1;
VB			v2d s0 = ((v2d *)BpF)[i];	v2d s1 = ((v2d *)BpF)[NLAT-1-i];		pe(i) = s0+s1;	po(i) = s0-s1;
 B			i++;
 B		} while (i<ni);
		l = LiM(shtns, 0,im);
Q		v2d* Ql = (v2d*) &Qlm[l];	// virtual pointer for l=0 and im
V		v2d* Sl = (v2d*) &Slm[l];		v2d* Tl = (v2d*) &Tlm[l];
		l=im*MRES;
3		double m_1 = 1.0/l;
Q		zl = shtns->zlm[im];
V		dzl = shtns->dzlm[im];
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
		while (l<llim) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
QE			v2d q0 = vdup(0.0);
QO			v2d q1 = vdup(0.0);
VE			v2d s0 = vdup(0.0);	v2d t1 = vdup(0.0);		v2d s0i = vdup(0.0);	v2d t1i = vdup(0.0);
VO			v2d t0 = vdup(0.0);	v2d s1 = vdup(0.0);		v2d t0i = vdup(0.0);	v2d s1i = vdup(0.0);
			i=i0;	do {		// tm[im] : polar optimization
VE				s0  += vdup(dzl[0].t) *to(i);		// ref: these E. Dormy p 72.
VO				t0  -= vdup(dzl[0].t) *po(i);
VO				s1  += vdup(dzl[1].t) *te(i);
VE				t1  -= vdup(dzl[1].t) *pe(i);
VE				s0i -= vdup(dzl[0].p) *pe(i);
VO				t0i -= vdup(dzl[0].p) *te(i);
VO				s1i -= vdup(dzl[1].p) *po(i);
VE				t1i -= vdup(dzl[1].p) *to(i);
			  #ifndef SHT_3COMP
QE				q0  += re(i) * vdup(zl[0]);		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
QO				q1  += ro(i) * vdup(zl[1]);	// Qlm[LiM(l+1,im)] += zlm[im][(l+1-m)*NLAT/2 + i] * fm[i];
			  #else
QE				q0  += re(i) * vdup(dzl[0].p);		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
QO				q1  += ro(i) * vdup(dzl[1].p);	// Qlm[LiM(l+1,im)] += zlm[im][(l+1-m)*NLAT/2 + i] * fm[i];
			  #endif
Q				zl +=2;
V				dzl +=2;
				i++;
			} while (i < ni);
3			q0 *= vdup((l*(l+1))*m_1);
3			q1 *= vdup(((l+1)*(l+2))*m_1);
VE			Sl[l] = addi(s0,s0i);	Tl[l+1] = addi(t1,t1i);
VO			Tl[l] = addi(t0,t0i);	Sl[l+1] = addi(s1,s1i);
QE			Ql[l] = q0;
QO			Ql[l+1] = q1;
			l+=2;
		}
		if (l==llim) {
			long int lstride=1;
	  #ifdef SHT_VAR_LTR
			if (l != LMAX) lstride=2;
	  #endif
QE			v2d q0 = vdup(0.0);	// Qlm[LiM(l,im)] = 0.0;
VE			v2d s0 = vdup(0.0);	v2d s0i = vdup(0.0);
VO			v2d t0 = vdup(0.0);	v2d t0i = vdup(0.0);
			i=i0;	do {		// tm[im] : polar optimization
VE				s0  += vdup(dzl[0].t) *to(i);
VO				t0  -= vdup(dzl[0].t) *po(i);
VE				s0i -= vdup(dzl[0].p) *pe(i);
VO				t0i -= vdup(dzl[0].p) *te(i);
			  #ifndef SHT_3COMP
QE				q0  += re(i) * vdup(zl[0]);		// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
			  #else
QE				q0  += vdup(dzl[0].p) * re(i);	// Qlm[LiM(l,im)] += zlm[im][(l-m)*NLAT/2 + i] * fp[i];
			  #endif
Q				zl  += lstride;
V				dzl += lstride;
				i++;
			} while(i<ni);
3			q0 *= vdup((l*(l+1))*m_1);
VE			Sl[l] = addi(s0,s0i);
VO			Tl[l] = addi(t0,t0i);
QE			Ql[l] = q0;
	  #ifdef SHT_VAR_LTR
	  		l++;
		}
	    while( l<=LMAX ) {
Q			Ql[l] = vdup(0.0);
V			Sl[l] = vdup(0.0);	Tl[l] = vdup(0.0);
			l++;
      #endif
		}
	}
	#ifdef SHT_VAR_LTR
	if (imlim < MMAX) {
		im = imlim+1;
		l = LiM(shtns, im*MRES, im);
		do {
Q			((v2d*)Qlm)[l] = vdup(0.0);
V			((v2d*)Slm)[l] = vdup(0.0);		((v2d*)Tlm)[l] = vdup(0.0);
		} while(++l < NLM);
	}
	#endif

  	if (SHT_FFT > 1) {		// free memory
Q	    fftw_free(BrF - NLAT*(imlim+1));
VX	    fftw_free(BtF - NLAT*(imlim+1));	// this frees also BpF.
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
  }
