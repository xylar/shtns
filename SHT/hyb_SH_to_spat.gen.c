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

/// Truncation and spatial discretization are defined by \ref shtns_create and \ref shtns_set_grid_*
/// \param[in] shtns = a configuration created by \ref shtns_create with a grid set by shtns_set_grid_*
Q/// \param[in] Qlm = spherical harmonics coefficients :
Q/// complex double arrays of size shtns->nlm [unmodified].
S/// \param[in] Slm = spherical harmonics coefficients of \b Spheroidal scalar :
S/// complex double array of size shtns->nlm [unmodified].
T/// \param[in] Tlm = spherical harmonics coefficients of \b Toroidal scalar :
T/// complex double array of size shtns->nlm [unmodified].
Q/// \param[out] Vr = spatial scalar field : double array.
  #ifndef SHT_AXISYM
V/// \param[out] Vt, Vp = (theta,phi)-components of spatial vector : double arrays.
  #else
S/// \param[out] Vt = theta-component of spatial vector : double array.
T/// \param[out] Vp = phi-component of spatial vector : double array.  
  #endif
  #ifdef SHT_VAR_LTR
/// \param[in] llim = specify maximum degree of spherical harmonic. llim must be at most LMAX, and all spherical harmonic degree higher than llim are ignored. 
  #else
/// \param[in] llim MUST be shtns->lmax.
  #endif

// 3 components not possible with DCT acceleration.
3	#define SHT_NO_DCT
	#define NWAY 1

// MTR_DCT : -1 => no dct
//            0 => dct for m=0 only
//            m => dct up to m, (!!! MTR_DCT <= MTR !!!)

3	void GEN3(SHqst_to_spat_,ID_NME,SUFFIX)(shtns_cfg shtns, complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp, long int llim) {
QX	void GEN3(SH_to_spat_,ID_NME,SUFFIX)(shtns_cfg shtns, complex double *Qlm, double *Vr, long int llim) {
  #ifndef SHT_GRAD
VX	void GEN3(SHsphtor_to_spat_,ID_NME,SUFFIX)(shtns_cfg shtns, complex double *Slm, complex double *Tlm, double *Vt, double *Vp, long int llim) {
  #else
	#ifndef SHT_AXISYM
S	void GEN3(SHsph_to_spat_,ID_NME,SUFFIX)(shtns_cfg shtns, complex double *Slm, double *Vt, double *Vp, long int llim) {
T	void GEN3(SHtor_to_spat_,ID_NME,SUFFIX)(shtns_cfg shtns, complex double *Tlm, double *Vt, double *Vp, long int llim) {
	#else
S	void GEN3(SHsph_to_spat_,ID_NME,SUFFIX)(shtns_cfg shtns, complex double *Slm, double *Vt, long int llim) {
T	void GEN3(SHtor_to_spat_,ID_NME,SUFFIX)(shtns_cfg shtns, complex double *Tlm, double *Vp, long int llim) {
	#endif
  #endif

Q	v2d *BrF;
  #ifndef SHT_AXISYM
V	v2d *BtF, *BpF;
V	struct DtDp *dyl;
Q	complex double re,ro;
V	complex double te,to, pe,po;
V	complex double dte,dto, dpe,dpo;
Q	#define BR0(i) ((double *)BrF)[2*(i)]
V	#define BT0(i) ((double *)BtF)[2*(i)]
V	#define BP0(i) ((double *)BpF)[2*(i)]
  #else
S	v2d *BtF;
T	v2d *BpF;
Q	double re,ro;
S	double te,to;
T	double pe,po;
Q	#define BR0(i) ((double *)BrF)[i]
S	#define BT0(i) ((double *)BtF)[i]
T	#define BP0(i) ((double *)BpF)[i]
  #endif
Q	double *yl;
V	double *dyl0;
	long int k,m,l;
	long int im, imlim, imlim_dct;
  #ifdef _GCC_VEC_
Q	v2d Ql0[(llim+3)>>1];		// we need some zero-padding.
S	v2d Sl0[(llim+2)>>1];
T	v2d Tl0[(llim+2)>>1];
  #endif

  #ifndef SHT_AXISYM
Q	BrF = (v2d *) Vr;
V	BtF = (v2d *) Vt;	BpF = (v2d *) Vp;
	if (SHT_FFT > 1) {		// alloc memory for the FFT
		unsigned long nspat = ((NPHI>>1) +1)*NLAT;
QX		BrF = fftw_malloc( nspat * sizeof(complex double) );
VX		BtF = fftw_malloc( 2* nspat * sizeof(complex double) );
VX		BpF = BtF + nspat;
3		BrF = fftw_malloc( 3* nspat * sizeof(complex double) );
3		BtF = BrF + nspat;		BpF = BtF + nspat;
	}
  #else
Q	BrF = (v2d*) Vr;
S	BtF = (v2d*) Vt;
T	BpF = (v2d*) Vp;
  #endif

	imlim = MTR;
	#ifdef SHT_VAR_LTR
		if (MTR*MRES > (int) llim) imlim = ((int) llim)/MRES;		// 32bit mul and div should be faster
	#endif
	im=0;	m=0;
  #ifdef _GCC_VEC_
  	{	// store the m=0 coefficients in an efficient & vectorizable way.
		l=1;
Q		double* Ql = (double*) &Ql0;
S		double* Sl = (double*) &Sl0;
T		double* Tl = (double*) &Tl0;
Q		Ql[0] = (double) Qlm[0];		// l=0
		do {		// for m=0, compress the complex Q,S,T to double
Q			Ql[l] = (double) Qlm[l];	//	Ql[l+1] = (double) Qlm[l+1];
S			Sl[l-1] = (double) Slm[l];	//	Sl[l] = (double) Slm[l+1];
T			Tl[l-1] = (double) Tlm[l];	//	Tl[l] = (double) Tlm[l+1];
			l++;
		} while(l<=llim);
Q		Ql[l] = 0.0;
S		Sl[l-1] = 0.0;
T		Tl[l-1] = 0.0;
	}
  #endif
  #ifndef SHT_NO_DCT
	imlim_dct = MTR_DCT;
	#ifdef SHT_VAR_LTR
		if (MTR_DCT*MRES > (int) llim) imlim_dct = ((int) llim)/MRES;		// 32bit mul and div should be faster
	#endif
	#ifndef _GCC_VEC_
Q		double* Ql = (double*) Qlm;
S		double* Sl = (double*) Slm;
T		double* Tl = (double*) Tlm;
	#endif
Q		yl = shtns->ykm_dct[im];
V		dyl0 = (double *) shtns->dykm_dct[im];		// only theta derivative (d/dphi = 0 for m=0)
		k=0;
V		l = 1;
		do {
Q			l = k;
	#ifndef _GCC_VEC_
Q			re = 0.0;	ro = 0.0;
S			te = 0.0;	to = 0.0;
T			pe = 0.0;	po = 0.0;
			while(l<llim) {
QE				re += yl[0]  * Ql[2*l];
QO				ro += yl[1]  * Ql[2*l+2];
SO				te += dyl0[0] * Sl[2*l];
TE				pe -= dyl0[0] * Tl[2*l];
SE				to += dyl0[1] * Sl[2*l+2];
TO				po -= dyl0[1] * Tl[2*l+2];
				l+=2;
Q				yl+=2;
V				dyl0+=2;
			}
			if (l==llim) {
QE				re += yl[0]  * Ql[2*l];
SO				te += dyl0[0] * Sl[2*l];
TE				pe -= dyl0[0] * Tl[2*l];
Q				yl+=2;
V				dyl0+=2;
			}
Q			BR0(k) = re;	BR0(k+1) = ro;
VX		#ifndef SHT_AXISYM
VX			BT0(k) = 0.0;	BT0(k+1) = 0.0;			// required for tor or sph only transform
VX			BP0(k) = 0.0;	BP0(k+1) = 0.0;
VX		#endif
S			BT0(k) = te;	BT0(k+1) = to;
T			BP0(k) = pe;	BP0(k+1) = po;
	#else
Q			v2d r[NWAY];
S			v2d t[NWAY];
T			v2d p[NWAY];
			for (int j=0; j<NWAY; j++) {
Q				r[j] = vdup(0.0);
S				t[j] = vdup(0.0);
T				p[j] = vdup(0.0);
			}
			l >>= 1;	// l = l/2;
			do {
				for (int j=0; j<NWAY; j++) {
Q					r[j] += ((v2d*) yl)[j]   * Ql0[l];		// { re, ro }
S					t[j] += ((v2d*) dyl0)[j] * Sl0[l];		// { te, to }
T					p[j] -= ((v2d*) dyl0)[j] * Tl0[l];		// { pe, po }
				}
				l++;
Q				yl+=2*NWAY;
V				dyl0+=2*NWAY;
Q			} while(2*l <= llim);
V			} while(2*l < llim);
			for (int j=0; j<NWAY; j++) {
		#ifndef SHT_AXISYM
Q				BR0(k+2*j) = vlo_to_dbl(r[j]);		BR0(k+1+2*j) = vhi_to_dbl(r[j]);
VX				BT0(k+2*j) = 0.0;					BT0(k+1+2*j) = 0.0;
S				BT0(k+2*j) = vlo_to_dbl(t[j]);		BT0(k+1+2*j) = vhi_to_dbl(t[j]);
VX				BP0(k+2*j) = 0.0;					BP0(k+1+2*j) = 0.0;
T				BP0(k+2*j) = vlo_to_dbl(p[j]);		BP0(k+1+2*j) = vhi_to_dbl(p[j]);
		#else
Q				*((v2d*)(((double*)BrF)+k+2*j)) = r[j];
S				*((v2d*)(((double*)BtF)+k+2*j)) = t[j];
T				*((v2d*)(((double*)BpF)+k+2*j)) = p[j];
		#endif
			}
	#endif
			k+=2*NWAY;
		#ifdef SHT_VAR_LTR
Q			yl  += ((LMAX>>1) - (llim>>1))*2;
V			dyl0 += (((LMAX+1)>>1) - ((llim+1)>>1))*2;
		#endif
V			l = k-1;
Q		} while (k<=llim);
V		} while (k<=llim+1);
		while (k<NLAT) {	// dct padding (NLAT is even)
Q			BR0(k) = 0.0;	BR0(k+1) = 0.0;
		#ifndef SHT_AXISYM
V			BT0(k) = 0.0;	BT0(k+1) = 0.0;			// required for tor or sph only transform
V			BP0(k) = 0.0;	BP0(k+1) = 0.0;
		#else
S			BT0(k) = 0.0;	BT0(k+1) = 0.0;
T			BP0(k) = 0.0;	BP0(k+1) = 0.0;
		#endif
			k+=2;
		}
    #ifdef SHT_AXISYM
Q		fftw_execute_r2r(shtns->idct_r1,Vr, Vr);		// iDCT m=0
S		fftw_execute_r2r(shtns->idct_r1,Vt, Vt);		// iDCT m=0
T		fftw_execute_r2r(shtns->idct_r1,Vp, Vp);		// iDCT m=0
V		double* st_1 = shtns->st_1;
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
  #else		// ifndef SHT_NO_DCT
	#ifndef _GCC_VEC_
Q		double* Ql = (double*) Qlm;
S		double* Sl = (double*) Slm;
T		double* Tl = (double*) Tlm;
	#endif
		k=0;
Q		yl  = shtns->ylm[im];
V		dyl0 = (double *) shtns->dylm[im];	// only theta derivative (d/dphi = 0 for m=0)
		do {	// ops : NLAT_2 * [ (lmax-m+1) + 4]	: almost twice as fast.
	#ifndef _GCC_VEC_
			l=1;
Q			re = yl[0] * Ql[0];		// re += ylm[im][k*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
Q			yl++;
Q			ro = 0.0;
S			te = 0.0;	to = 0.0;
T			pe = 0.0;	po = 0.0;
			while (l<llim) {	// compute even and odd parts
Q				ro += yl[0] * Ql[2*l];		// re += ylm[im][k*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
QB				re += yl[1] * Ql[2*l+2];	// ro += ylm[im][k*(LMAX-m+1) + (l+1-m)] * Qlm[LiM(l+1,im)];
T				pe -= dyl0[0] * Tl[2*l];	// m=0 : everything is real.
SB				te += dyl0[0] * Sl[2*l];
TB				po -= dyl0[1] * Tl[2*l+2];
S				to += dyl0[1] * Sl[2*l+2];
				l+=2;
Q				yl+=2;
V				dyl0+=2;
			}
			if (l==llim) {
Q				ro += yl[0] * Ql[2*l];
TB				pe -= dyl0[0] * Tl[2*l];
S				te += dyl0[0] * Sl[2*l];
V				dyl0+=2;
			}
Q			yl++;
Q			BR0(k) = re + ro;
V		#ifndef SHT_AXISYM
V			BT0(k) = 0.0;		BP0(k) = 0.0;		// required for partial tor or sph transform
V		#endif
S			BT0(k) = te + to;			// Bt = dS/dt
T			BP0(k) = pe + po;			// Bp = - dT/dt
			k++;
QB			BR0(NLAT-k) = re - ro;
VB		#ifndef SHT_AXISYM
VB			BT0(NLAT-k) = 0.0;		BP0(NLAT-k) = 0.0;		// required for partial tor or sph transform
VB		#endif
SB			BT0(NLAT-k) = te - to;
TB			BP0(NLAT-k) = pe - po;
	#else
			l=0;
Q			v2d r[NWAY];
S			v2d t[NWAY];
T			v2d p[NWAY];
			for (int j=0; j<NWAY; j++) {
Q				r[j] = vdup(0.0);
S				t[j] = vdup(0.0);
T				p[j] = vdup(0.0);
			}
			do {
				for (int j=0; j<NWAY; j++) {
Q					r[j] += ((v2d*) yl)[j]   * Ql0[l];		// { re, ro }
S					t[j] += ((v2d*) dyl0)[j] * Sl0[l];		// { te, to }
T					p[j] -= ((v2d*) dyl0)[j] * Tl0[l];		// { pe, po }
				}
				l++;
Q				yl+=2*NWAY;
V				dyl0+=2*NWAY;
QX			} while (2*l <= llim);
V			} while (2*l < llim);
3			if (2*l == llim) {
3				for (int j=0; j<NWAY; j++) {
3					r[j] += ((v2d*) yl)[j]   * Ql0[l];		// { re, ro }
3				}
3				yl+=2*NWAY;
3			}
		#if __SSE3__
	/*	alternate code, which may be faster (slightly) on SSE3.	*/
		  for (int j=0; j<NWAY; j++) {
Q			r[j] = addi(r[j],r[j]);		// { re-ro , re+ro }
S			t[j] = addi(t[j],t[j]);		// { te-to , te+to }
T			p[j] = addi(p[j],p[j]);		// { pe-po , pe+po }
		  }
		  for (int j=0; j<NWAY; j++) {
Q			BR0(k+j) = vhi_to_dbl(r[j]);
V		#ifndef SHT_AXISYM
V			BT0(k+j) = 0.0;		BP0(k+j) = 0.0;		// required for partial tor or sph transform
V		#endif
S			BT0(k+j) = vhi_to_dbl(t[j]);	// Bt = dS/dt
T			BP0(k+j) = vhi_to_dbl(p[j]);	// Bp = - dT/dt
QB			BR0(NLAT-NWAY-k+j) = vlo_to_dbl(r[NWAY-1-j]);
VB		#ifndef SHT_AXISYM
VB			BT0(NLAT-NWAY-k+j) = 0.0;		BP0(NLAT-NWAY-k+j) = 0.0;		// required for partial tor or sph transform
VB		#endif
SB			BT0(NLAT-NWAY-k+j) = vlo_to_dbl(t[NWAY-1-j]);
TB			BP0(NLAT-NWAY-k+j) = vlo_to_dbl(p[NWAY-1-j]);
		  }
			k+=NWAY;
		#else
		  for (int j=0; j<NWAY; j++) {
Q			BR0(k+j) = vhi_to_dbl(r[j]) + vlo_to_dbl(r[j]);
V		#ifndef SHT_AXISYM
V			BT0(k+j) = 0.0;		BP0(k+j) = 0.0;		// required for partial tor or sph transform
V		#endif
S			BT0(k+j) = vhi_to_dbl(t[j]) + vlo_to_dbl(t[j]);	// Bt = dS/dt
T			BP0(k+j) = vhi_to_dbl(p[j]) + vlo_to_dbl(p[j]);	// Bp = - dT/dt
QB			BR0(NLAT-NWAY-k+j) = vlo_to_dbl(r[j]) - vhi_to_dbl(r[j]);
VB		#ifndef SHT_AXISYM
VB			BT0(NLAT-NWAY-k+j) = 0.0;		BP0(NLAT-NWAY-k+j) = 0.0;		// required for partial tor or sph transform
VB		#endif
SB			BT0(NLAT-NWAY-k+j) = vlo_to_dbl(t[j]) - vhi_to_dbl(t[j]);
TB			BP0(NLAT-NWAY-k+j) = vlo_to_dbl(p[j]) - vhi_to_dbl(p[j]);
		  }
		  k+=NWAY;
		#endif
	#endif
		#ifdef SHT_VAR_LTR
Q			yl  += ((LMAX>>1) - (llim>>1))*2;
V			dyl0 += (((LMAX+1)>>1) - ((llim+1)>>1))*2;
		#endif
		} while (k < NLAT_2);
  #endif		// ifndef SHT_NO_DCT

  #ifndef SHT_AXISYM
	im=1;
Q	BrF += NLAT;
V	BtF += NLAT;	BpF += NLAT;
    #ifndef SHT_NO_DCT
	while(im<=imlim_dct) {		// dct for im <= MTR_DCT
		m=im*MRES;
		l = LiM(shtns, 0,im);
Q		v2d* Ql = (v2d*) &Qlm[l];		// virtual pointer for l=0 and im
S		v2d* Sl = (v2d*) &Slm[l];
T		v2d* Tl = (v2d*) &Tlm[l];
Q		yl = shtns->ykm_dct[im];
V		dyl = shtns->dykm_dct[im];
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
		#ifdef SHT_VAR_LTR
Q			yl+= (LMAX-llim);
V			dyl+= (LMAX-llim);
		#endif
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

	while(im<=imlim) {	// regular for MTR_DCT < im <= MTR
		m = im*MRES;
3		double m_1 = 1.0/m;
		l = LiM(shtns, 0,im);
Q		v2d* Ql = (v2d*) &Qlm[l];	// virtual pointer for l=0 and im
S		v2d* Sl = (v2d*) &Slm[l];	// virtual pointer for l=0 and im
T		v2d* Tl = (v2d*) &Tlm[l];
		k=0;	l=shtns->tm[im];
		while (k<l) {	// polar optimization
Q			BrF[k] = vdup(0.0);
QB			BrF[NLAT-l + k] = vdup(0.0);	// south pole zeroes <=> BrF[im*NLAT + NLAT-(k+1)] = 0.0;
V			BtF[k] = vdup(0.0);		BpF[k] = vdup(0.0);
VB			BtF[NLAT-l + k] = vdup(0.0);		BpF[NLAT-l + k] = vdup(0.0);	// south pole zeroes
			k++;
		}
Q		yl  = shtns->ylm[im];
V		dyl = shtns->dylm[im];
		do {	// ops : NLAT_2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
Q			v2d re = vdup(0.0); 	v2d ro = vdup(0.0);
V			v2d dte = vdup(0.0); 	v2d dto = vdup(0.0); 	v2d pe = vdup(0.0); 	v2d po = vdup(0.0);
V			v2d dpe = vdup(0.0); 	v2d dpo = vdup(0.0); 	v2d te = vdup(0.0); 	v2d to = vdup(0.0);
			while (l<llim) {	// compute even and odd parts
QX				re  += vdup(yl[0]) * Ql[l];		// re += ylm[im][k*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
QX				ro  += vdup(yl[1]) * Ql[l+1];	// ro += ylm[im][k*(LMAX-m+1) + (l+1-m)] * Qlm[LiM(l+1,im)];
TB				dpo -= vdup(dyl[0].t) * Tl[l];
S				dto += vdup(dyl[0].t) * Sl[l];
TB				te  += vdup(dyl[0].p) * Tl[l];
S				pe  += vdup(dyl[0].p) * Sl[l];
3				re  += vdup(dyl[0].p) * Ql[l];
T				dpe -= vdup(dyl[1].t) * Tl[l+1];
SB				dte += vdup(dyl[1].t) * Sl[l+1];
T				to  += vdup(dyl[1].p) * Tl[l+1];
SB				po  += vdup(dyl[1].p) * Sl[l+1];
3				ro  += vdup(dyl[1].p) * Ql[l+1];
				l+=2;
Q				yl+=2;
V				dyl+=2;
			}
			if (l==llim) {
QX				re  += vdup(yl[0]) * Ql[l];		// re += ylm[im][k*(LMAX-m+1) + (l-m)] * Qlm[LiM(l,im)];
TO				dpo -= vdup(dyl[0].t) * Tl[l];
SE				dto += vdup(dyl[0].t) * Sl[l];
TO				te  += vdup(dyl[0].p) * Tl[l];
SE				pe  += vdup(dyl[0].p) * Sl[l];
3				re  += vdup(dyl[0].p) * Ql[l];
Q				yl++;
V				dyl++;
			}
3			s2d qv = vdup(shtns->st[k]*m_1);
V			BtF[k] = addi(dte+dto, te+to);		// Bt = dS/dt       + I.m/sint *T
VB			BtF[NLAT-1-k] = addi(dte-dto, te-to);
3			re *= qv;	ro *= qv;
V			BpF[k] = addi(dpe+dpo, pe+po);		// Bp = I.m/sint * S - dT/dt
VB			BpF[NLAT-1-k] = addi(dpe-dpo, pe-po);
Q			BrF[k] = re + ro;
QB			BrF[NLAT-1-k] = re - ro;
			k++;
		#ifdef SHT_VAR_LTR
Q			yl  += (LMAX-llim);
V			dyl += (LMAX-llim);
		#endif
		} while (k < NLAT_2);
		im++;
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
	}
	for (k=0; k < NLAT*((NPHI>>1) -imlim); k++) {	// padding for high m's
Q			BrF[k] = vdup(0.0);
V			BtF[k] = vdup(0.0);	BpF[k] = vdup(0.0);
	}
Q	BrF -= NLAT*(imlim+1);		// restore original pointer
V	BtF -= NLAT*(imlim+1);	BpF -= NLAT*(imlim+1);	// restore original pointer

    if (NPHI>1) {
    #ifndef SHT_NO_DCT
		if (MTR_DCT >= 0) {
Q			fftw_execute_r2r(shtns->idct,(double *) BrF, (double *) BrF);		// iDCT
V			fftw_execute_r2r(shtns->idct,(double *) BtF, (double *) BtF);		// iDCT
V			fftw_execute_r2r(shtns->idct,(double *) BpF, (double *) BpF);		// iDCT
V			double* st_1 = shtns->st_1;
V			k=0;	do {		// m=0
V				double sin_1 = st_1[k];		double sin_2 = st_1[k+1];
V				((double *)BtF)[2*k] *= sin_1;   	((double *)BpF)[2*k] *= sin_1;
V				((double *)BtF)[2*k+2] *= sin_2; 	((double *)BpF)[2*k+2] *= sin_2;
V				k+=2;
V			} while(k<NLAT);
Q			if (MRES & 1) {		// odd m's must be divided by sin(theta)
Q				double* st_1 = shtns->st_1;
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
Q		fftw_execute_dft_c2r(shtns->ifft, (complex double *) BrF, Vr);
V		fftw_execute_dft_c2r(shtns->ifft, (complex double *) BtF, Vt);
V		fftw_execute_dft_c2r(shtns->ifft, (complex double *) BpF, Vp);
		if (SHT_FFT > 1) {		// free memory
Q			fftw_free(BrF);
VX			fftw_free(BtF);	// this frees also BpF.
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
Q			fftw_execute_r2r(shtns->idct_r1,Vr, Vr);		// iDCT m=0
S			fftw_execute_r2r(shtns->idct_r1,Vt, Vt);		// iDCT m=0
T			fftw_execute_r2r(shtns->idct_r1,Vp, Vp);		// iDCT m=0
V			double* st_1 = shtns->st_1;
V			k=0;	do {
V		#ifdef _GCC_VEC_
V				v2d sin_1 = ((v2d *)st_1)[k];
S				((v2d *)Vt)[k] *= sin_1;
T				((v2d *)Vp)[k] *= sin_1;
V		#else
V			double sin_1 = st_1[2*k]; 	double sin_2 = st_1[2*k+1];
S			Vt[2*k] *= sin_1;		Vt[2*k+1] *= sin_2;
T			Vp[2*k] *= sin_1;		Vp[2*k+1] *= sin_2;
V		#endif
V				k++;
V			} while (k<NLAT_2);
		}
    #endif
    }
  #endif

	#undef NWAY
Q	#undef BR0
V	#undef BT0
V	#undef BP0
  }
