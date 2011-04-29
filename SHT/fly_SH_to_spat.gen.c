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

#void SH_to_spat_fly(complex double *Qlm, double *Vr)
#{
Q	v2d *BrF;
  #ifndef SHT_AXISYM
V	v2d *BtF, *BpF;
Q	#define BR0(i) ((double *)BrF)[2*(i)]
V	#define BT0(i) ((double *)BtF)[2*(i)]
V	#define BP0(i) ((double *)BpF)[2*(i)]
  #else
S	v2d *BtF;
T	v2d *BpF;
Q	#define BR0(i) ((double *)BrF)[i]
S	#define BT0(i) ((double *)BtF)[i]
T	#define BP0(i) ((double *)BpF)[i]
  #endif
	long int llim;
	long int k,im,m,l;
	double *al;
Q	double Ql0[LTR+1];
S	double Sl0[LTR+1];
T	double Tl0[LTR];

  #ifndef SHT_AXISYM
Q	BrF = (v2d *) Vr;
V	BtF = (v2d *) Vt;	BpF = (v2d *) Vp;
	if (SHT_FFT > 1) {		// alloc memory for the FFT
		long int nspat = ((NPHI>>1) +1)*NLAT;
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

	llim = LTR;		// copy LTR to a local variable for faster access (inner loop limit)
	im=0;	m=0;
 		l=1;
Q		Ql0[0] = (double) Qlm[0];		// l=0
		do {		// for m=0, compress the complex Q,S,T to double
Q			Ql0[l] = (double) Qlm[l];	//	Ql[l+1] = (double) Qlm[l+1];
S			Sl0[l-1] = (double) Slm[l];	//	Sl[l] = (double) Slm[l+1];
T			Tl0[l-1] = (double) Tlm[l];	//	Tl[l] = (double) Tlm[l+1];
			l++;
		} while(l<=llim);
		k=0;
		do {
			l=0;	al = al0;
			s2d cost[NWAY], y0[NWAY], y1[NWAY];
V			s2d sint[NWAY], dy0[NWAY], dy1[NWAY];
Q			s2d re[NWAY], ro[NWAY];
S			s2d te[NWAY], to[NWAY];
T			s2d pe[NWAY], po[NWAY];
			for (int j=0; j<NWAY; j++) {
				cost[j] = ((s2d*)(ct+k))[j];
V				sint[j] = ((s2d*)(st+k))[j];
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
		  #ifndef SHT_AXISYM
			for (int j=0; j<NWAY; j++) {
Q				BR0(k+2*j) = vlo_to_dbl(re[j]+ro[j]);
Q				BR0(k+2*j+1) = vhi_to_dbl(re[j]+ro[j]);
Q				BR0(NLAT-k-2-2*j) = vhi_to_dbl(re[j]-ro[j]);
Q				BR0(NLAT-k-1-2*j) = vlo_to_dbl(re[j]-ro[j]);
S				BT0(k+2*j) = vlo_to_dbl(te[j]+to[j]);
S				BT0(k+2*j+1) = vhi_to_dbl(te[j]+to[j]);
S				BT0(NLAT-k-2-2*j) = vhi_to_dbl(te[j]-to[j]);
S				BT0(NLAT-k-1-2*j) = vlo_to_dbl(te[j]-to[j]);
T				BP0(k+2*j) = vlo_to_dbl(pe[j]+po[j]);
T				BP0(k+2*j+1) = vhi_to_dbl(pe[j]+po[j]);
T				BP0(NLAT-k-2-2*j) = vhi_to_dbl(pe[j]-po[j]);
T				BP0(NLAT-k-1-2*j) = vlo_to_dbl(pe[j]-po[j]);
			}
		  #else
			for (int j=0; j<NWAY; j++) {
Q				((v2d*)(((double*)BrF)+k))[j] = re[j]+ro[j];
Q				*((v2d*)(((double*)BrF)+NLAT-k-2-2*j)) = vxchg(re[j]-ro[j]);
S				((v2d*)(((double*)BtF)+k))[j] = te[j]+to[j];
S				*((v2d*)(((double*)BtF)+NLAT-k-2-2*j)) = vxchg(te[j]-to[j]);
T				((v2d*)(((double*)BpF)+k))[j] = pe[j]+po[j];
T				*((v2d*)(((double*)BpF)+NLAT-k-2-2*j)) = vxchg(pe[j]-po[j]);
			}
		  #endif
			k+=2*NWAY;
		#else
			for (int j=0; j<NWAY; j++) {
Q				BR0(k+j) = (re[j]+ro[j]);
Q				BR0(NLAT-k-1-j) = (re[j]-ro[j]);
S				BT0(k+j) = (te[j]+to[j]);
S				BT0(NLAT-k-1-j) = (te[j]-to[j]);
T				BP0(k+j) = (pe[j]+po[j]);
T				BP0(NLAT-k-1-j) = (pe[j]-po[j]);
			}
			k+=NWAY;
		#endif
		} while (k < NLAT_2);

  #ifndef SHT_AXISYM
	im=1;
//	#undef NWAY
//V	#define NWAY 1
//QX	#define NWAY 2
Q	BrF += NLAT;
V	BtF += NLAT;	BpF += NLAT;
	while(im<=MTR) {	// regular for MTR_DCT < im <= MTR
		m = im*MRES;
Q		v2d* Ql = (v2d*) &Qlm[LiM(0,im)];	// virtual pointer for l=0 and im
S		v2d* Sl = (v2d*) &Slm[LiM(0,im)];	// virtual pointer for l=0 and im
T		v2d* Tl = (v2d*) &Tlm[LiM(0,im)];
		k=0;	l=tm[im];
	#if _GCC_VEC_
		l=(l>>1)*2;		// stay on a 16 byte boundary
	#endif
		while (k<l) {	// polar optimization
Q			BrF[k] = vdup(0.0);
Q			BrF[NLAT-l + k] = vdup(0.0);	// south pole zeroes <=> BrF[im*NLAT + NLAT-(k+1)] = 0.0;
V			BtF[k] = vdup(0.0);		BpF[k] = vdup(0.0);
V			BtF[NLAT-l + k] = vdup(0.0);		BpF[NLAT-l + k] = vdup(0.0);	// south pole zeroes
			k++;
		}
		do {
			al = alm[im];
			s2d cost[NWAY], y0[NWAY], y1[NWAY];
V			s2d st2[NWAY], dy0[NWAY], dy1[NWAY];
Q			v2d re[NWAY], ro[NWAY];
		#if _GCC_VEC_
Q			v2d rex[NWAY], rox[NWAY];
		#endif
V			v2d te[NWAY], to[NWAY], dpe[NWAY], dpo[NWAY];
V			v2d pe[NWAY], po[NWAY], dte[NWAY], dto[NWAY];
			for (int j=0; j<NWAY; j++) {
				cost[j] = ((s2d*)(st+k))[j];
				y0[j] = vdup(al[0]);
V				st2[j] = cost[j]*cost[j]*vdup(1.0/m);
V				y0[j] *= vdup(m);		// for the vector transform, compute ylm*m/sint
			}
Q			l=m;
V			l=m-1;
		#ifndef SHT_LARGE_L
			do {		// sin(theta)^m
				if (l&1) for (int j=0; j<NWAY; j++) y0[j] *= cost[j];
				for (int j=0; j<NWAY; j++) cost[j] *= cost[j];
			} while(l >>= 1);
			#define Y0 y0[j]
			#define Y1 y1[j]
V			#define DY0 dy0[j]
V			#define DY1 dy1[j]
		#else
			s2d scale[NWAY];
			for (int j=0; j<NWAY; j++) {
				y0[j] *= vdup(SHT_LEG_SCALEF);
				scale[j] = vdup(1.0/SHT_LEG_SCALEF);
			}

			int ll = l >> 8;
			do {		// sin(theta)^m
				if (l&1) for (int j=0; j<NWAY; j++) scale[j] *= cost[j];
				l >>= 1;
				for (int j=0; j<NWAY; j++) cost[j] *= cost[j];
			} while(l > ll);
			while(l > 0) {
				l--;
				for (int j=0; j<NWAY; j++) scale[j] *= cost[j];
			}
			#define Y0 (y0[j]*scale[j])
			#define Y1 (y1[j]*scale[j])
V			#define DY0 (dy0[j]*scale[j])
V			#define DY1 (dy1[j]*scale[j])
		#endif
			for (int j=0; j<NWAY; j++) {
				cost[j] = ((s2d*)(ct+k))[j];
V				dy0[j] = cost[j]*y0[j];
			}
			l=m;
Q			v2d q = ((v2d*)Ql)[l];
S			v2d s = ((v2d*)Sl)[l];
T			v2d t = ((v2d*)Tl)[l];
			for (int j=0; j<NWAY; j++) {
Q				re[j]  = Y0 * q;
S				pe[j]  = Y0 * s;
S				dto[j] = DY0 * s;
T				dpo[j] = -DY0 * t;
T				te[j]  = -Y0 * t;
T				dpe[j] = vdup(0.0);		to[j] = vdup(0.0);
S				po[j] = vdup(0.0);		dte[j] = vdup(0.0);
Q				ro[j] = vdup(0.0);
			}
		#if _GCC_VEC_
			for (int j=0; j<NWAY; j++) {
Q				rox[j] = vdup(0.0);
Q				rex[j] = Y0 * vxchg(q);
S				dpe[j] = subadd(dpe[j], Y0 * vxchg(s));
S				to[j]  = subadd(to[j], DY0 * vxchg(s));
T				po[j]  = subadd(po[j], DY0 * vxchg(t));
T				dte[j] = subadd(dte[j], Y0 * vxchg(t));
			}
		#endif
			l++;
			for (int j=0; j<NWAY; j++) {
				y1[j]  = (vdup(al[1])*y0[j]) *cost[j];		//	y1[j] = vdup(al[1])*cost[j]*y0[j];
V				dy1[j] = (vdup(al[1])*y0[j]) *(cost[j]*cost[j] - st2[j]);		//	dy1[j] = vdup(al[1])*(cost[j]*dy0[j] - y0[j]*st2[j]);
			}
			al+=2;
			while (l<llim) {	// compute even and odd parts
Q				q = ((v2d*)Ql)[l];
S				s = ((v2d*)Sl)[l];
T				t = ((v2d*)Tl)[l];
				for (int j=0; j<NWAY; j++) {
Q					ro[j]  += Y1 * q;
S					po[j]  += Y1 * s;
S					dte[j] += DY1 * s;
T					dpe[j] -= DY1 * t;
T					to[j]  -= Y1 * t;
				}
		#if _GCC_VEC_
				for (int j=0; j<NWAY; j++) {
Q					rox[j] += Y1 * vxchg(q);
S					dpo[j] = subadd(dpo[j], Y1 * vxchg(s));
S					te[j]  = subadd(te[j], DY1 * vxchg(s));
T					pe[j]  = subadd(pe[j], DY1 * vxchg(t));
T					dto[j] = subadd(dto[j], Y1 * vxchg(t));
				}
		#endif
				for (int j=0; j<NWAY; j++) {
					y0[j] = vdup(al[1])*cost[j]*y1[j] + vdup(al[0])*y0[j];
V					dy0[j] = vdup(al[1])*(cost[j]*dy1[j] - y1[j]*st2[j]) + vdup(al[0])*dy0[j];
				}
Q				q = ((v2d*)Ql)[l+1];
S				s = ((v2d*)Sl)[l+1];
T				t = ((v2d*)Tl)[l+1];
				for (int j=0; j<NWAY; j++) {
Q					re[j]  += Y0 * q;
S					pe[j]  += Y0 * s;
S					dto[j] += DY0 * s;
T					dpo[j] -= DY0 * t;
T					te[j]  -= Y0 * t;
				}
		#if _GCC_VEC_
				for (int j=0; j<NWAY; j++) {
Q					rex[j] += Y0 * vxchg(q);
S					dpe[j] = subadd(dpe[j], Y0 * vxchg(s));
S					to[j]  = subadd(to[j], DY0 * vxchg(s));
T					po[j]  = subadd(po[j], DY0 * vxchg(t));
T					dte[j] = subadd(dte[j], Y0 * vxchg(t));
				}
		#endif
				l+=2;
				for (int j=0; j<NWAY; j++) {
					y1[j] = vdup(al[3])*cost[j]*y0[j] + vdup(al[2])*y1[j];
V					dy1[j] = vdup(al[3])*(cost[j]*dy0[j] - y0[j]*st2[j]) + vdup(al[2])*dy1[j];
				}
				al+=4;
			}
			if (l==llim) {
Q				q = ((v2d*)Ql)[l];
S				s = ((v2d*)Sl)[l];
T				t = ((v2d*)Tl)[l];
				for (int j=0; j<NWAY; j++) {
Q					ro[j]  += Y1 * q;
S					po[j]  += Y1 * s;
S					dte[j] += DY1 * s;
T					dpe[j] -= DY1 * t;
T					to[j]  -= Y1 * t;
				}
		#if _GCC_VEC_
				for (int j=0; j<NWAY; j++) {
Q					rox[j] += Y1 * vxchg(q);
S					dpo[j] = subadd(dpo[j], Y1 * vxchg(s));
S					te[j]  = subadd(te[j], DY1 * vxchg(s));
T					pe[j]  = subadd(pe[j], DY1 * vxchg(t));
T					dto[j] = subadd(dto[j], Y1 * vxchg(t));
				}
		#endif
			}
3			for (int j=0; j<NWAY; j++) {
3				cost[j]  = ((s2d*)(st+k))[j];
3				cost[j] *= vdup(1.0/m);
3			}
3			for (int j=0; j<NWAY; j++) {  re[j] *= cost[j];  ro[j] *= cost[j];  }
		#if _GCC_VEC_
			for (int j=0; j<NWAY; j++) {
3				rex[j] *= cost[j];  rox[j] *= cost[j];
Q				((complex double *)BrF)[k+2*j] = vlo_to_dbl(re[j]+ro[j]) +I*vlo_to_dbl(rex[j]+rox[j]);
Q				((complex double *)BrF)[k+1+2*j] = vhi_to_dbl(rex[j]+rox[j]) +I*vhi_to_dbl(re[j]+ro[j]);
Q				((complex double *)BrF)[NLAT-k-2-2*j] = vhi_to_dbl(rex[j]-rox[j]) +I*vhi_to_dbl(re[j]-ro[j]);
Q				((complex double *)BrF)[NLAT-k-1-2*j] = vlo_to_dbl(re[j]-ro[j]) +I*vlo_to_dbl(rex[j]-rox[j]);
V				((complex double *)BtF)[k+2*j] = -I*vlo_to_dbl(te[j]+to[j]) + vlo_to_dbl(dte[j]+dto[j]);
V				((complex double *)BtF)[k+1+2*j] = vhi_to_dbl(te[j]+to[j]) + I*vhi_to_dbl(dte[j]+dto[j]);
V				((complex double *)BtF)[NLAT-k-2-2*j] = vhi_to_dbl(te[j]-to[j]) + I*vhi_to_dbl(dte[j]-dto[j]);
V				((complex double *)BtF)[NLAT-k-1-2*j] = -I*vlo_to_dbl(te[j]-to[j]) + vlo_to_dbl(dte[j]-dto[j]);
V				((complex double *)BpF)[k+2*j] = I*vlo_to_dbl(pe[j]+po[j]) + vlo_to_dbl(dpe[j]+dpo[j]);
V				((complex double *)BpF)[k+1+2*j] = -vhi_to_dbl(pe[j]+po[j]) + I*vhi_to_dbl(dpe[j]+dpo[j]);
V				((complex double *)BpF)[NLAT-k-2-2*j] = -vhi_to_dbl(pe[j]-po[j]) + I*vhi_to_dbl(dpe[j]-dpo[j]);
V				((complex double *)BpF)[NLAT-k-1-2*j] = I*vlo_to_dbl(pe[j]-po[j]) + vlo_to_dbl(dpe[j]-dpo[j]);
			}
			k += 2*NWAY;
		#else
			for (int j=0; j<NWAY; j++) {
Q				BrF[k+j] = (re[j]+ro[j]);
Q				BrF[NLAT-k-1-j] = (re[j]-ro[j]);
V				BtF[k+j] = addi(dte[j]+dto[j], -(te[j]+to[j]));		// Bt = dS/dt       + I.m/sint *T
VB				BtF[NLAT-1-k-j] = addi(dte[j]-dto[j], -(te[j]-to[j]));
V				BpF[k+j] = addi(dpe[j]+dpo[j], pe[j]+po[j]);		// Bp = I.m/sint * S - dT/dt
VB				BpF[NLAT-1-k-j] = addi(dpe[j]-dpo[j], pe[j]-po[j]);
			}
			k+=NWAY;
		#endif
		} while (k < NLAT_2);
		im++;
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
	}
	for (k=0; k < NLAT*((NPHI>>1) -MTR); k++) {	// padding for high m's
Q			BrF[k] = vdup(0.0);
V			BtF[k] = vdup(0.0);	BpF[k] = vdup(0.0);
	}
Q	BrF -= NLAT*(MTR+1);		// restore original pointer
V	BtF -= NLAT*(MTR+1);	BpF -= NLAT*(MTR+1);	// restore original pointer

    if (NPHI>1) {
Q		fftw_execute_dft_c2r(ifft, (complex double *) BrF, Vr);
V		fftw_execute_dft_c2r(ifft, (complex double *) BtF, Vt);
V		fftw_execute_dft_c2r(ifft, (complex double *) BpF, Vp);
		if (SHT_FFT > 1) {		// free memory
Q			fftw_free(BrF);
VX			fftw_free(BtF);		// this frees also BpF.
		}
    } else {
		k=1;	do {	// compress complex to real
Q			Vr[k] = ((double *)BrF)[2*k];
V			Vt[k] = ((double *)BtF)[2*k];
V			Vp[k] = ((double *)BpF)[2*k];
			k++;
		} while(k<NLAT);
    }
  #endif

	#undef Y0
	#undef Y1
V	#undef DY0
V	#undef DY1
Q	#undef BR0
V	#undef BT0
V	#undef BP0
# }
