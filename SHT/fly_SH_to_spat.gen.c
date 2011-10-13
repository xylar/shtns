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

3	void GEN3(SHqst_to_spat_fly,NWAY,SUFFIX)(shtns_cfg shtns, complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp, long int llim) {
QX	void GEN3(SH_to_spat_fly,NWAY,SUFFIX)(shtns_cfg shtns, complex double *Qlm, double *Vr, long int llim) {
  #ifndef SHT_GRAD
VX	void GEN3(SHsphtor_to_spat_fly,NWAY,SUFFIX)(shtns_cfg shtns, complex double *Slm, complex double *Tlm, double *Vt, double *Vp, long int llim) {
  #else
S	void GEN3(SHsph_to_spat_fly,NWAY,SUFFIX)(shtns_cfg shtns, complex double *Slm, double *Vt, double *Vp, long int llim) {
T	void GEN3(SHtor_to_spat_fly,NWAY,SUFFIX)(shtns_cfg shtns, complex double *Tlm, double *Vt, double *Vp, long int llim) {
  #endif

Q	v2d *BrF;
  #ifndef SHT_AXISYM
V	v2d *BtF, *BpF;
Q	#define BR0(i) ((double *)BrF)[2*(i)]
V	#define BT0(i) ((double *)BtF)[2*(i)]
V	#define BP0(i) ((double *)BpF)[2*(i)]
Q	#define qr vdup(creal(Ql[l]));
Q	#define qi vdup(cimag(Ql[l]));
S	#define sr vdup(creal(Sl[l]));
S	#define si vdup(cimag(Sl[l]));
T	#define tr vdup(creal(Tl[l]));
T	#define ti vdup(cimag(Tl[l]));
  #else
S	v2d *BtF;
T	v2d *BpF;
Q	#define BR0(i) ((double *)BrF)[i]
S	#define BT0(i) ((double *)BtF)[i]
T	#define BP0(i) ((double *)BpF)[i]
  #endif
	long int nk, k,m,l;
	long int im, imlim;
	double *alm, *al, *ct, *st;
Q	double Ql0[llim+1];
S	double Sl0[llim+1];
T	double Tl0[llim];

  #ifndef SHT_AXISYM
Q	BrF = (v2d *) Vr;
V	BtF = (v2d *) Vt;	BpF = (v2d *) Vp;
	if (shtns->ncplx_fft > 0) {		// alloc memory for the FFT
QX		BrF = VMALLOC( shtns->ncplx_fft * sizeof(complex double) );
VX		BtF = VMALLOC( 2* shtns->ncplx_fft * sizeof(complex double) );
VX		BpF = BtF + shtns->ncplx_fft;
3		BrF = VMALLOC( 3* shtns->ncplx_fft * sizeof(complex double) );
3		BtF = BrF + shtns->ncplx_fft;		BpF = BtF + shtns->ncplx_fft;
	}
	imlim = MTR;
	#ifdef SHT_VAR_LTR
		if (MTR*MRES > (int) llim) imlim = ((int) llim)/MRES;		// 32bit mul and div should be faster
	#endif
  #else
	#ifdef SHT_GRAD
S		if (Vp != NULL) { k=0; do { ((s2d*)Vp)[k]=vdup(0.0); } while(++k<NLAT/2); Vp[NLAT-1]=0.0; }
T		if (Vt != NULL) { k=0; do { ((s2d*)Vt)[k]=vdup(0.0); } while(++k<NLAT/2); Vt[NLAT-1]=0.0; }
	#endif
Q	BrF = (v2d*) Vr;
S	BtF = (v2d*) Vt;
T	BpF = (v2d*) Vp;
  #endif

	ct = shtns->ct;		st = shtns->st;
	im=0;
 		l=1;
		alm = shtns->alm[0];
Q		Ql0[0] = (double) Qlm[0];		// l=0
		do {		// for m=0, compress the complex Q,S,T to double
Q			Ql0[l] = (double) Qlm[l];	//	Ql[l+1] = (double) Qlm[l+1];
S			Sl0[l-1] = (double) Slm[l];	//	Sl[l] = (double) Slm[l+1];
T			Tl0[l-1] = (double) Tlm[l];	//	Tl[l] = (double) Tlm[l+1];
			l++;
		} while(l<=llim);
		k=0;	nk = NLAT_2;
	#if _GCC_VEC_
		nk = (nk+1)>>1;
	#endif
		do {
			l=0;	al = alm;
			s2d cost[NWAY], y0[NWAY], y1[NWAY];
V			s2d sint[NWAY], dy0[NWAY], dy1[NWAY];
Q			s2d re[NWAY], ro[NWAY];
S			s2d te[NWAY], to[NWAY];
T			s2d pe[NWAY], po[NWAY];
			for (int j=0; j<NWAY; j++) {
				cost[j] = ((s2d*)(ct))[j+k];
V				sint[j] = ((s2d*)(st))[j+k];
				y0[j] = vdup(al[0]);
V				dy0[j] = vdup(0.0);
Q				re[j] = y0[j] * vdup(Ql0[0]);
S				to[j] = dy0[j];
T				po[j] = dy0[j];
			}
			for (int j=0; j<NWAY; j++) {
				y1[j] = vdup(al[0]*al[1]) * cost[j];
V				dy1[j] = -vdup(al[0]*al[1]) * sint[j];
			}
			for (int j=0; j<NWAY; j++) {
Q				ro[j] = y1[j] * vdup(Ql0[1]);
S				te[j] = dy1[j] * vdup(Sl0[0]);
T				pe[j] = -dy1[j] * vdup(Tl0[0]);
			}
			al+=2;	l+=2;
			while(l<llim) {
				for (int j=0; j<NWAY; j++) {
					y0[j]  = vdup(al[1])*cost[j]*y1[j] + vdup(al[0])*y0[j];
V					dy0[j] = vdup(al[1])*(cost[j]*dy1[j] - y1[j]*sint[j]) + vdup(al[0])*dy0[j];
				}
				for (int j=0; j<NWAY; j++) {
Q					re[j] += y0[j] * vdup(Ql0[l]);
S					to[j] += dy0[j] * vdup(Sl0[l-1]);
T					po[j] -= dy0[j] * vdup(Tl0[l-1]);
				}
				for (int j=0; j<NWAY; j++) {
					y1[j]  = vdup(al[3])*cost[j]*y0[j] + vdup(al[2])*y1[j];
V					dy1[j] = vdup(al[3])*(cost[j]*dy0[j] - y0[j]*sint[j]) + vdup(al[2])*dy1[j];
				}
				for (int j=0; j<NWAY; j++) {
Q					ro[j] += y1[j] * vdup(Ql0[l+1]);
S					te[j] += dy1[j] * vdup(Sl0[l]);
T					pe[j] -= dy1[j] * vdup(Tl0[l]);
				}
				al+=4;	l+=2;
			}
			if (l==llim) {
				for (int j=0; j<NWAY; j++) {
					y0[j]  = vdup(al[1])*cost[j]*y1[j] + vdup(al[0])*y0[j];
V					dy0[j] = vdup(al[1])*(cost[j]*dy1[j] - y1[j]*sint[j]) + vdup(al[0])*dy0[j];
				}
				for (int j=0; j<NWAY; j++) {
Q					re[j] += y0[j] * vdup(Ql0[l]);
S					to[j] += dy0[j] * vdup(Sl0[l-1]);
T					po[j] -= dy0[j] * vdup(Tl0[l-1]);
				}
			}
		#if _GCC_VEC_
			for (int j=0; j<NWAY; j++) {
Q				BrF[j+k] = re[j]+ro[j];
Q				BrF[NLAT_2-1-(j+k)] = vxchg(re[j]-ro[j]);
S				BtF[j+k] = te[j]+to[j];
S				BtF[NLAT_2-1-(j+k)] = vxchg(te[j]-to[j]);
T				BpF[j+k] = pe[j]+po[j];
T				BpF[NLAT_2-1-(j+k)] = vxchg(pe[j]-po[j]);
			}
		#else
			for (int j=0; j<NWAY; j++) {
Q				BR0(k+j) = (re[j]+ro[j]);
Q				BR0(NLAT-k-1-j) = (re[j]-ro[j]);
S				BT0(k+j) = (te[j]+to[j]);
S				BT0(NLAT-k-1-j) = (te[j]-to[j]);
T				BP0(k+j) = (pe[j]+po[j]);
T				BP0(NLAT-k-1-j) = (pe[j]-po[j]);
			}
		#endif
			k+=NWAY;
		} while (k < nk);

  #ifndef SHT_AXISYM
//	#undef NWAY
//V	#define NWAY 1
//QX	#define NWAY 2
	#if _GCC_VEC_
Q		BrF += NLAT_2;
V		BtF += NLAT_2;	BpF += NLAT_2;
	#else
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
	#endif
  if (llim <= SHT_L_RESCALE_FLY) {
	for(im=1; im<=imlim; im++) {
		m = im*MRES;
		l = LiM(shtns, 0,im);
		alm = shtns->alm[im];
Q		complex double* Ql = &Qlm[l];	// virtual pointer for l=0 and im
S		complex double* Sl = &Slm[l];	// virtual pointer for l=0 and im
T		complex double* Tl = &Tlm[l];
		k=0;	l=shtns->tm[im];
	#if _GCC_VEC_
		l>>=1;		// stay on a 16 byte boundary
		while (k<l) {	// polar optimization
Q			BrF[k] = vdup(0.0);				BrF[(NPHI-2*im)*NLAT_2 + k] = vdup(0.0);
Q			BrF[NLAT_2-l+k] = vdup(0.0);	BrF[(NPHI+1-2*im)*NLAT_2 -l+k] = vdup(0.0);
V			BtF[k] = vdup(0.0);				BtF[(NPHI-2*im)*NLAT_2 + k] = vdup(0.0);
V			BtF[NLAT_2-l+k] = vdup(0.0);	BtF[(NPHI+1-2*im)*NLAT_2 -l+k] = vdup(0.0);
V			BpF[k] = vdup(0.0);				BpF[(NPHI-2*im)*NLAT_2 + k] = vdup(0.0);
V			BpF[NLAT_2-l+k] = vdup(0.0);	BpF[(NPHI+1-2*im)*NLAT_2 -l+k] = vdup(0.0);
			k++;
		}
	#else
		while (k<l) {	// polar optimization
Q			BrF[k] = 0.0;		BrF[NLAT-l+k] = vdup(0.0);
V			BtF[k] = 0.0;		BtF[NLAT-l+k] = vdup(0.0);
V			BpF[k] = 0.0;		BpF[NLAT-l+k] = vdup(0.0);
			k++;
		}
	#endif
		do {
			al = alm;
			s2d cost[NWAY], y0[NWAY], y1[NWAY];
V			s2d st2[NWAY], dy0[NWAY], dy1[NWAY];
Q			s2d rer[NWAY], rei[NWAY], ror[NWAY], roi[NWAY];
V			s2d ter[NWAY], tei[NWAY], tor[NWAY], toi[NWAY];
V			s2d per[NWAY], pei[NWAY], por[NWAY], poi[NWAY];
			for (int j=0; j<NWAY; j++) {
				cost[j] = ((s2d*)(st))[k+j];
				y0[j] = vdup(al[0]);
V				st2[j] = cost[j]*cost[j]*vdup(1.0/m);
V				y0[j] *= vdup(m);		// for the vector transform, compute ylm*m/sint
			}
Q			l=m;
V			l=m-1;
			do {		// sin(theta)^m
				if (l&1) for (int j=0; j<NWAY; j++) y0[j] *= cost[j];
				for (int j=0; j<NWAY; j++) cost[j] *= cost[j];
			} while(l >>= 1);
			for (int j=0; j<NWAY; j++) {
				cost[j] = ((s2d*)(ct))[k+j];
V				dy0[j] = cost[j]*y0[j];
Q				ror[j] = vdup(0.0);		roi[j] = vdup(0.0);
Q				rer[j] = vdup(0.0);		rei[j] = vdup(0.0);
			}
			for (int j=0; j<NWAY; j++) {
				y1[j]  = (vdup(al[1])*y0[j]) *cost[j];		//	y1[j] = vdup(al[1])*cost[j]*y0[j];
S				por[j] = vdup(0.0);		tei[j] = vdup(0.0);
T				tor[j] = vdup(0.0);		pei[j] = vdup(0.0);
V				dy1[j] = (vdup(al[1])*y0[j]) *(cost[j]*cost[j] - st2[j]);		//	dy1[j] = vdup(al[1])*(cost[j]*dy0[j] - y0[j]*st2[j]);
S				poi[j] = vdup(0.0);		ter[j] = vdup(0.0);
T				toi[j] = vdup(0.0);		per[j] = vdup(0.0);
			}
			l=m;		al+=2;
			#define Y0 y0[j]
			#define Y1 y1[j]
V			#define DY0 dy0[j]
V			#define DY1 dy1[j]
			while (l<llim) {	// compute even and odd parts
Q				for (int j=0; j<NWAY; j++) {	rer[j] += Y0 * qr;		rei[j] += Y0 * qi;	}
S				for (int j=0; j<NWAY; j++) {	tor[j] += DY0 * sr;		pei[j] += Y0 * sr;	}
S				for (int j=0; j<NWAY; j++) {	toi[j] += DY0 * si;		per[j] -= Y0 * si;	}
T				for (int j=0; j<NWAY; j++) {	por[j] -= DY0 * tr;		tei[j] += Y0 * tr;	}
T				for (int j=0; j<NWAY; j++) {	poi[j] -= DY0 * ti;		ter[j] -= Y0 * ti;	}
				l++;
				for (int j=0; j<NWAY; j++) {
					y0[j] = vdup(al[1])*cost[j]*y1[j] + vdup(al[0])*y0[j];
V					dy0[j] = vdup(al[1])*(cost[j]*dy1[j] - y1[j]*st2[j]) + vdup(al[0])*dy0[j];
				}
Q				for (int j=0; j<NWAY; j++) {	ror[j] += Y1 * qr;		roi[j] += Y1 * qi;	}
S				for (int j=0; j<NWAY; j++) {	ter[j] += DY1 * sr;		poi[j] += Y1 * sr;	}
S				for (int j=0; j<NWAY; j++) {	tei[j] += DY1 * si;		por[j] -= Y1 * si;	}
T				for (int j=0; j<NWAY; j++) {	per[j] -= DY1 * tr;		toi[j] += Y1 * tr;	}
T				for (int j=0; j<NWAY; j++) {	pei[j] -= DY1 * ti;		tor[j] -= Y1 * ti;	}
				l++;
				for (int j=0; j<NWAY; j++) {
					y1[j] = vdup(al[3])*cost[j]*y0[j] + vdup(al[2])*y1[j];
V					dy1[j] = vdup(al[3])*(cost[j]*dy0[j] - y0[j]*st2[j]) + vdup(al[2])*dy1[j];
				}
				al+=4;
			}
			if (l==llim) {
Q				for (int j=0; j<NWAY; j++) {	rer[j] += Y0 * qr;		rei[j] += Y0 * qi;	}
S				for (int j=0; j<NWAY; j++) {	tor[j] += DY0 * sr;		pei[j] += Y0 * sr;	}
S				for (int j=0; j<NWAY; j++) {	toi[j] += DY0 * si;		per[j] -= Y0 * si;	}
T				for (int j=0; j<NWAY; j++) {	por[j] -= DY0 * tr;		tei[j] += Y0 * tr;	}
T				for (int j=0; j<NWAY; j++) {	poi[j] -= DY0 * ti;		ter[j] -= Y0 * ti;	}
			}
			#undef Y0
			#undef Y1
V			#undef DY0
V			#undef DY1
3			for (int j=0; j<NWAY; j++) {
3				cost[j]  = ((s2d*)(st))[k+j];
3				cost[j] *= vdup(1.0/m);
3			}
3			for (int j=0; j<NWAY; j++) {  rer[j] *= cost[j];  ror[j] *= cost[j];	rei[j] *= cost[j];  roi[j] *= cost[j];  }
		#if _GCC_VEC_
			for (int j=0; j<NWAY; j++) {
Q				y0[j] = rer[j] + ror[j];	y1[j] = rei[j] + roi[j];	// north
Q				((complex double *)BrF)[k+j] =                      vlo_to_dbl(y0[j]) - vhi_to_dbl(y1[j]) + I*(vhi_to_dbl(y0[j]) + vlo_to_dbl(y1[j]));
Q				((complex double *)BrF)[(NPHI-2*im)*NLAT_2 + k+j] = vlo_to_dbl(y0[j]) + vhi_to_dbl(y1[j]) + I*(vhi_to_dbl(y0[j]) - vlo_to_dbl(y1[j]));
Q				y0[j] = rer[j] - ror[j];	y1[j] = rei[j] - roi[j];	// south
Q				((complex double *)BrF)[NLAT_2-1 -k-j] =                vhi_to_dbl(y0[j]) - vlo_to_dbl(y1[j]) + I*(vlo_to_dbl(y0[j]) + vhi_to_dbl(y1[j]));
Q				((complex double *)BrF)[(NPHI+1-2*im)*NLAT_2 -1- k-j] = vhi_to_dbl(y0[j]) + vlo_to_dbl(y1[j]) + I*(vlo_to_dbl(y0[j]) - vhi_to_dbl(y1[j]));
V				y0[j] = ter[j] + tor[j];	y1[j] = tei[j] + toi[j];	// north
V				((complex double *)BtF)[k+j] =                      vlo_to_dbl(y0[j]) - vhi_to_dbl(y1[j]) + I*(vhi_to_dbl(y0[j]) + vlo_to_dbl(y1[j]));
V				((complex double *)BtF)[(NPHI-2*im)*NLAT_2 + k+j] = vlo_to_dbl(y0[j]) + vhi_to_dbl(y1[j]) + I*(vhi_to_dbl(y0[j]) - vlo_to_dbl(y1[j]));
V				y0[j] = ter[j] - tor[j];	y1[j] = tei[j] - toi[j];	// south
V				((complex double *)BtF)[NLAT_2-1 -k-j] =                vhi_to_dbl(y0[j]) - vlo_to_dbl(y1[j]) + I*(vlo_to_dbl(y0[j]) + vhi_to_dbl(y1[j]));
V				((complex double *)BtF)[(NPHI+1-2*im)*NLAT_2 -1- k-j] = vhi_to_dbl(y0[j]) + vlo_to_dbl(y1[j]) + I*(vlo_to_dbl(y0[j]) - vhi_to_dbl(y1[j]));
V				y0[j] = per[j] + por[j];	y1[j] = pei[j] + poi[j];	// north
V				((complex double *)BpF)[k+j] =                      vlo_to_dbl(y0[j]) - vhi_to_dbl(y1[j]) + I*(vhi_to_dbl(y0[j]) + vlo_to_dbl(y1[j]));
V				((complex double *)BpF)[(NPHI-2*im)*NLAT_2 + k+j] = vlo_to_dbl(y0[j]) + vhi_to_dbl(y1[j]) + I*(vhi_to_dbl(y0[j]) - vlo_to_dbl(y1[j]));
V				y0[j] = per[j] - por[j];	y1[j] = pei[j] - poi[j];	// south
V				((complex double *)BpF)[NLAT_2-1 -k-j] =                vhi_to_dbl(y0[j]) - vlo_to_dbl(y1[j]) + I*(vlo_to_dbl(y0[j]) + vhi_to_dbl(y1[j]));
V				((complex double *)BpF)[(NPHI+1-2*im)*NLAT_2 -1- k-j] = vhi_to_dbl(y0[j]) + vlo_to_dbl(y1[j]) + I*(vlo_to_dbl(y0[j]) - vhi_to_dbl(y1[j]));
			}
		#else
			for (int j=0; j<NWAY; j++) {
Q				BrF[k+j] = (rer[j]+ror[j]) + I*(rei[j]+roi[j]);
Q				BrF[NLAT-k-1-j] = (rer[j]-ror[j]) + I*(rei[j]-roi[j]);
V				BtF[k+j] = (ter[j]+tor[j]) + I*(tei[j]+toi[j]);
V				BtF[NLAT-1-k-j] = (ter[j]-tor[j]) + I*(tei[j]-toi[j]);
V				BpF[k+j] = (per[j]+por[j]) + I*(pei[j]+poi[j]);
V				BpF[NLAT-1-k-j] = (per[j]-por[j]) + I*(pei[j]-poi[j]);
			}
		#endif
			k+=NWAY;
		} while (k < nk);
	#if _GCC_VEC_
Q		BrF += NLAT_2;
V		BtF += NLAT_2;	BpF += NLAT_2;
	#else
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
	#endif
	}
  } else {		// llim > SHT_L_RESCALE_FLY
	for(im=1; im<=imlim; im++) {
		m = im*MRES;
		l = LiM(shtns, 0,im);
		alm = shtns->alm[im];
Q		complex double* Ql = &Qlm[l];	// virtual pointer for l=0 and im
S		complex double* Sl = &Slm[l];	// virtual pointer for l=0 and im
T		complex double* Tl = &Tlm[l];
		k=0;	l=shtns->tm[im];
		l>>=1;		// stay on a 16 byte boundary
		while (k<l) {	// polar optimization
Q			BrF[2*k] = vdup(0.0);		BrF[2*k+1] = vdup(0.0);
Q			BrF[NLAT+2*(k-l)] = vdup(0.0);		BrF[NLAT+2*(k-l)+1] = vdup(0.0);
V			BtF[2*k] = vdup(0.0);		BtF[2*k+1] = vdup(0.0);
V			BpF[2*k] = vdup(0.0);		BpF[2*k+1] = vdup(0.0);
V			BtF[NLAT+2*(k-l)] = vdup(0.0);		BtF[NLAT+2*(k-l)+1] = vdup(0.0);
V			BpF[NLAT+2*(k-l)] = vdup(0.0);		BpF[NLAT+2*(k-l)+1] = vdup(0.0);
			k++;
		}
	#if _GCC_VEC_
		k *= 2;
	#endif
		do {
			al = alm;
			s2d cost[NWAY], y0[NWAY], y1[NWAY], scale[NWAY];
V			s2d st2[NWAY], dy0[NWAY], dy1[NWAY];
Q			s2d rer[NWAY], rei[NWAY], ror[NWAY], roi[NWAY];
V			s2d ter[NWAY], tei[NWAY], tor[NWAY], toi[NWAY];
V			s2d per[NWAY], pei[NWAY], por[NWAY], poi[NWAY];
			for (int j=0; j<NWAY; j++) {
				cost[j] = ((s2d*)(st))[k+j];
				y0[j] = vdup(al[0]);
V				st2[j] = cost[j]*cost[j]*vdup(1.0/m);
V				y0[j] *= vdup(m);		// for the vector transform, compute ylm*m/sint
			}
Q			l=m;
V			l=m-1;
			for (int j=0; j<NWAY; j++) {
				y0[j] *= vdup(SHT_LEG_SCALEF);
				scale[j] = vdup(1.0/SHT_LEG_SCALEF);
			}
			long int ll = l >> 8;
			do {		// sin(theta)^m
				if (l&1) for (int j=0; j<NWAY; j++) scale[j] *= cost[j];
				l >>= 1;
				for (int j=0; j<NWAY; j++) cost[j] *= cost[j];
			} while(l > ll);
			while(l > 0) {
				l--;
				for (int j=0; j<NWAY; j++) scale[j] *= cost[j];
			}
			for (int j=0; j<NWAY; j++) {
				cost[j] = ((s2d*)(ct))[k+j];
V				dy0[j] = cost[j]*y0[j];
Q				ror[j] = vdup(0.0);		roi[j] = vdup(0.0);
Q				rer[j] = vdup(0.0);		rei[j] = vdup(0.0);
			}
			for (int j=0; j<NWAY; j++) {
				y1[j]  = (vdup(al[1])*y0[j]) *cost[j];		//	y1[j] = vdup(al[1])*cost[j]*y0[j];
S				por[j] = vdup(0.0);		tei[j] = vdup(0.0);
T				tor[j] = vdup(0.0);		pei[j] = vdup(0.0);
V				dy1[j] = (vdup(al[1])*y0[j]) *(cost[j]*cost[j] - st2[j]);		//	dy1[j] = vdup(al[1])*(cost[j]*dy0[j] - y0[j]*st2[j]);
S				poi[j] = vdup(0.0);		ter[j] = vdup(0.0);
T				toi[j] = vdup(0.0);		per[j] = vdup(0.0);
			}
			l=m;		al+=2;
			if (SHT_ACCURACY > 0) {		// this saves a lot of work, as rescale is not needed here.
				while ((vlo_to_dbl(y1[NWAY-1]) < SHT_ACCURACY*SHT_LEG_SCALEF)&&(vlo_to_dbl(y1[NWAY-1]) > -SHT_ACCURACY*SHT_LEG_SCALEF)&&(l<llim)) {
					for (int j=0; j<NWAY; j++) {
						y0[j] = vdup(al[1])*cost[j]*y1[j] + vdup(al[0])*y0[j];
V						dy0[j] = vdup(al[1])*(cost[j]*dy1[j] - y1[j]*st2[j]) + vdup(al[0])*dy0[j];
					}
					for (int j=0; j<NWAY; j++) {					
						y1[j] = vdup(al[3])*cost[j]*y0[j] + vdup(al[2])*y1[j];
V						dy1[j] = vdup(al[3])*(cost[j]*dy0[j] - y0[j]*st2[j]) + vdup(al[2])*dy1[j];				
					}
					l+=2;	al+=4;
				}
			}
			#define Y0 (y0[j]*scale[j])
			#define Y1 (y1[j]*scale[j])
V			#define DY0 (dy0[j]*scale[j])
V			#define DY1 (dy1[j]*scale[j])
			while (l<llim) {	// compute even and odd parts
				for (int j=0; j<NWAY; j++) {
Q					rer[j] += Y0 * qr;		rei[j] += Y0 * qi;
S					pei[j] += Y0 * sr;		per[j] -= Y0 * si;
T					tei[j] += Y0 * tr;		ter[j] -= Y0 * ti;
				}
V				for (int j=0; j<NWAY; j++) {
S					tor[j] += DY0 * sr;		toi[j] += DY0 * si;
T					por[j] -= DY0 * tr;		poi[j] -= DY0 * ti;
V				}
				l++;
				for (int j=0; j<NWAY; j++) {
					y0[j] = vdup(al[1])*cost[j]*y1[j] + vdup(al[0])*y0[j];
V					dy0[j] = vdup(al[1])*(cost[j]*dy1[j] - y1[j]*st2[j]) + vdup(al[0])*dy0[j];
				}
				for (int j=0; j<NWAY; j++) {
Q					ror[j] += Y1 * qr;		roi[j] += Y1 * qi;
S					poi[j] += Y1 * sr;		por[j] -= Y1 * si;
T					toi[j] += Y1 * tr;		tor[j] -= Y1 * ti;
				}
V				for (int j=0; j<NWAY; j++) {
S					ter[j] += DY1 * sr;		tei[j] += DY1 * si;
T					per[j] -= DY1 * tr;		pei[j] -= DY1 * ti;
V				}
				l++;
				for (int j=0; j<NWAY; j++) {
					y1[j] = vdup(al[3])*cost[j]*y0[j] + vdup(al[2])*y1[j];
V					dy1[j] = vdup(al[3])*(cost[j]*dy0[j] - y0[j]*st2[j]) + vdup(al[2])*dy1[j];
				}
				al+=4;
			}
			if (l==llim) {
				for (int j=0; j<NWAY; j++) {
Q					rer[j] += Y0 * qr;		rei[j] += Y0 * qi;
S					pei[j] += Y0 * sr;		per[j] -= Y0 * si;
T					tei[j] += Y0 * tr;		ter[j] -= Y0 * ti;
				}
V				for (int j=0; j<NWAY; j++) {
S					tor[j] += DY0 * sr;		toi[j] += DY0 * si;
T					por[j] -= DY0 * tr;		poi[j] -= DY0 * ti;
V				}
			}
			#undef Y0
			#undef Y1
V			#undef DY0
V			#undef DY1
3			for (int j=0; j<NWAY; j++) {
3				cost[j]  = ((s2d*)(st))[k+j];
3				cost[j] *= vdup(1.0/m);
3			}
3			for (int j=0; j<NWAY; j++) {  rer[j] *= cost[j];  ror[j] *= cost[j];	rei[j] *= cost[j];  roi[j] *= cost[j];  }
		#if _GCC_VEC_
			for (int j=0; j<NWAY; j++) {
Q				((complex double *)BrF)[2*(k+j)] = vlo_to_dbl(rer[j]+ror[j]) + I*vlo_to_dbl(rei[j]+roi[j]);
Q				((complex double *)BrF)[2*(k+j)+1] = vhi_to_dbl(rer[j]+ror[j]) + I*vhi_to_dbl(rei[j]+roi[j]);
Q				((complex double *)BrF)[NLAT-2-2*(k+j)] = vhi_to_dbl(rer[j]-ror[j]) + I*vhi_to_dbl(rei[j]-roi[j]);
Q				((complex double *)BrF)[NLAT-1-2*(k+j)] = vlo_to_dbl(rer[j]-ror[j]) + I*vlo_to_dbl(rei[j]-roi[j]);
V				((complex double *)BtF)[2*(k+j)] = vlo_to_dbl(ter[j]+tor[j]) + I*vlo_to_dbl(tei[j]+toi[j]);
V				((complex double *)BtF)[2*(k+j)+1] = vhi_to_dbl(ter[j]+tor[j]) + I*vhi_to_dbl(tei[j]+toi[j]);
V				((complex double *)BtF)[NLAT-2-2*(k+j)] = vhi_to_dbl(ter[j]-tor[j]) + I*vhi_to_dbl(tei[j]-toi[j]);
V				((complex double *)BtF)[NLAT-1-2*(k+j)] = vlo_to_dbl(ter[j]-tor[j]) + I*vlo_to_dbl(tei[j]-toi[j]);
V				((complex double *)BpF)[2*(k+j)] = vlo_to_dbl(per[j]+por[j]) + I*vlo_to_dbl(pei[j]+poi[j]);
V				((complex double *)BpF)[2*(k+j)+1] = vhi_to_dbl(per[j]+por[j]) + I*vhi_to_dbl(pei[j]+poi[j]);
V				((complex double *)BpF)[NLAT-2-2*(k+j)] = vhi_to_dbl(per[j]-por[j]) + I*vhi_to_dbl(pei[j]-poi[j]);
V				((complex double *)BpF)[NLAT-1-2*(k+j)] = vlo_to_dbl(per[j]-por[j]) + I*vlo_to_dbl(pei[j]-poi[j]);
			}
		#else
			for (int j=0; j<NWAY; j++) {
Q				BrF[k+j] = (rer[j]+ror[j]) + I*(rei[j]+roi[j]);
Q				BrF[NLAT-k-1-j] = (rer[j]-ror[j]) + I*(rei[j]-roi[j]);
V				BtF[k+j] = (ter[j]+tor[j]) + I*(tei[j]+toi[j]);
V				BtF[NLAT-1-k-j] = (ter[j]-tor[j]) + I*(tei[j]-toi[j]);
V				BpF[k+j] = (per[j]+por[j]) + I*(pei[j]+poi[j]);
V				BpF[NLAT-1-k-j] = (per[j]-por[j]) + I*(pei[j]-poi[j]);
			}
		#endif
			k+=NWAY;
		} while (k < nk);
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
	}
  }
  #if _GCC_VEC_
	for (k=0; k < NLAT_2*(NPHI-1-2*imlim); k++) {	// padding for high m's
Q		BrF[k] = vdup(0.0);
V		BtF[k] = vdup(0.0);	BpF[k] = vdup(0.0);
	}
Q	BrF -= NLAT_2*(imlim+1);		// restore original pointer
V	BtF -= NLAT_2*(imlim+1);	BpF -= NLAT_2*(imlim+1);
  #else
	for (k=0; k < NLAT*((NPHI>>1) -imlim); k++) {	// padding for high m's
Q			BrF[k] = vdup(0.0);
V			BtF[k] = vdup(0.0);	BpF[k] = vdup(0.0);
	}
Q	BrF -= NLAT*(imlim+1);		// restore original pointer
V	BtF -= NLAT*(imlim+1);	BpF -= NLAT*(imlim+1);	// restore original pointer
  #endif
    // NPHI > 1 as SHT_AXISYM is not defined.
	if (shtns->ncplx_fft >= 0) {
	#if _GCC_VEC_
Q		fftw_execute_dft(shtns->ifftc, (complex double *) BrF, (complex double *) Vr);
V		fftw_execute_dft(shtns->ifftc, (complex double *) BtF, (complex double *) Vt);
V		fftw_execute_dft(shtns->ifftc, (complex double *) BpF, (complex double *) Vp);
	#else
Q		fftw_execute_dft_c2r(shtns->ifft, (complex double *) BrF, Vr);
V		fftw_execute_dft_c2r(shtns->ifft, (complex double *) BtF, Vt);
V		fftw_execute_dft_c2r(shtns->ifft, (complex double *) BpF, Vp);
	#endif
		if (shtns->ncplx_fft > 0) {		// free memory
Q			VFREE(BrF);
VX			VFREE(BtF);		// this frees also BpF.
		}
	}
  #endif

Q	#undef BR0
V	#undef BT0
V	#undef BP0
  }
