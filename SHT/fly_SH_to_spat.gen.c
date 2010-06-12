/*
 * Copyright (c) 2010 Centre National de la Recherche Scientifique.
 * written by Nathanael Schaeffer (CNRS, LGIT, Grenoble, France).
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
			s2d cost = *((s2d*)(ct+k));
QX		#ifdef _GCC_VEC_
QX			s2d ctb = *((s2d*)(ct+k+2));	// 4 latitudes/iteration with SSE2
QX		#else
QX			s2d ctb = *((s2d*)(ct+k+1));	// 2 latitudes/iteration without SSE2
QX		#endif
V			s2d sint = *((s2d*)(st+k));
			s2d y0 = vdup(al[0]);
QX			s2d y0b = y0;
V			s2d dy0 = vdup(0.0);
Q			s2d re = y0 * vdup(Ql0[0]);
QX			s2d reb = y0b * vdup(Ql0[0]);
S			s2d to = dy0;
T			s2d po = dy0;
			s2d y1 = vdup(al[1]) * cost * y0;
QX			s2d y1b = vdup(al[1]) * ctb * y0b;
V			s2d dy1 = -vdup(al[1]) * sint * y0;
Q			s2d ro = y1 * vdup(Ql0[1]);
QX			s2d rob = y1b * vdup(Ql0[1]);
S			s2d te = dy1 * vdup(Sl0[0]);
T			s2d pe = -dy1 * vdup(Tl0[0]);
			al+=2;	l+=2;
			while(l<llim) {
				y0  = vdup(al[1])*cost*y1 + vdup(al[0])*y0;
QX				y0b  = vdup(al[1])*ctb*y1b + vdup(al[0])*y0b;
V				dy0 = vdup(al[1])*(cost*dy1 - y1*sint) + vdup(al[0])*dy0;
Q				re += y0 * vdup(Ql0[l]);
QX				reb += y0b * vdup(Ql0[l]);
S				to += dy0 * vdup(Sl0[l-1]);
T				po -= dy0 * vdup(Tl0[l-1]);
				y1  = vdup(al[3])*cost*y0 + vdup(al[2])*y1;
QX				y1b  = vdup(al[3])*ctb*y0b + vdup(al[2])*y1b;
V				dy1 = vdup(al[3])*(cost*dy0 - y0*sint) + vdup(al[2])*dy1;
Q				ro += y1 * vdup(Ql0[l+1]);
QX				rob += y1b * vdup(Ql0[l+1]);
S				te += dy1 * vdup(Sl0[l]);
T				pe -= dy1 * vdup(Tl0[l]);
				al+=4;	l+=2;
			}
			if (l==llim) {
				y0  = vdup(al[1])*cost*y1 + vdup(al[0])*y0;
QX				y0b  = vdup(al[1])*ctb*y1b + vdup(al[0])*y0b;
V				dy0 = vdup(al[1])*(cost*dy1 - y1*sint) + vdup(al[0])*dy0;
Q				re += y0 * vdup(Ql0[l]);
QX				reb += y0b * vdup(Ql0[l]);
S				to += dy0 * vdup(Sl0[l-1]);
T				po -= dy0 * vdup(Tl0[l-1]);
			}
		#if _GCC_VEC_
		  #ifndef SHT_AXISYM
Q			BR0(k) = vlo_to_dbl(re+ro);
Q			BR0(k+1) = vhi_to_dbl(re+ro);
QX			BR0(k+2) = vlo_to_dbl(reb+rob);
QX			BR0(k+3) = vhi_to_dbl(reb+rob);
QX			BR0(NLAT-k-4) = vhi_to_dbl(reb-rob);
QX			BR0(NLAT-k-3) = vlo_to_dbl(reb-rob);
Q			BR0(NLAT-k-2) = vhi_to_dbl(re-ro);
Q			BR0(NLAT-k-1) = vlo_to_dbl(re-ro);
S			BT0(k) = vlo_to_dbl(te+to);
S			BT0(k+1) = vhi_to_dbl(te+to);
S			BT0(NLAT-k-2) = vhi_to_dbl(te-to);
S			BT0(NLAT-k-1) = vlo_to_dbl(te-to);
T			BP0(k) = vlo_to_dbl(pe+po);
T			BP0(k+1) = vhi_to_dbl(pe+po);
T			BP0(NLAT-k-2) = vhi_to_dbl(pe-po);
T			BP0(NLAT-k-1) = vlo_to_dbl(pe-po);
		  #else
Q			*((v2d*)(((double*)BrF)+k)) = re+ro;
QX			*((v2d*)(((double*)BrF)+k+2)) = reb+rob;
QX			*((v2d*)(((double*)BrF)+NLAT-k-4)) = vxchg(reb-rob);
Q			*((v2d*)(((double*)BrF)+NLAT-k-2)) = vxchg(re-ro);
S			*((v2d*)(((double*)BtF)+k)) = te+to;
S			*((v2d*)(((double*)BtF)+NLAT-k-2)) = vxchg(te-to);
S			*((v2d*)(((double*)BpF)+k)) = pe+po;
S			*((v2d*)(((double*)BpF)+NLAT-k-2)) = vxchg(pe-po);
		  #endif
			k+=2;
QX			k+=2;
		#else
Q			BR0(k) = (re+ro);
QX			BR0(k+1) = (reb+rob);
QX			BR0(NLAT-k-2) = (reb-rob);
Q			BR0(NLAT-k-1) = (re-ro);
S			BT0(k) = (te+to);
S			BT0(NLAT-k-1) = (te-to);
T			BP0(k) = (pe+po);
T			BP0(NLAT-k-1) = (pe-po);
			k++;
QX			k++;
		#endif
		} while (k < NLAT_2);

  #ifndef SHT_AXISYM
	im=1;
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
QX			al = alm[im];
V			al = dlm[im];
			s2d cost  = *((s2d*)(st+k));
			s2d y0 = vdup(al[0]);
Q			l=m;
V			l=m-1;	y0 *= vdup(m);		// for the vector transform, compute ylm*m/sint
			do {
				if (l&1) y0 *= cost;
				cost *= cost;
			} while(l >>= 1);
			cost = *((s2d*)(ct+k));
			l=m;
V			s2d dy = cost*y0;
Q			v2d q = ((v2d*)Ql)[l];
S			v2d s = ((v2d*)Sl)[l];
T			v2d t = ((v2d*)Tl)[l];
Q			v2d re  = y0 * q;
S			v2d pe  = y0 * s;
S			v2d dto = dy * s;
T			v2d dpo = -dy * t;
T			v2d te  = -y0 * t;
T			v2d dpe = vdup(0.0);	v2d to = vdup(0.0);
S			v2d po = vdup(0.0);		v2d dte = vdup(0.0);
Q			v2d ro = vdup(0.0);
		#if _GCC_VEC_
Q			v2d rox = vdup(0.0);
Q			v2d rex = y0 * vxchg(q);
S			dpe = subadd(dpe, y0 * vxchg(s));
S			to  = subadd(to, dy * vxchg(s));
T			po  = subadd(po, dy * vxchg(t));
T			dte = subadd(dte, y0 * vxchg(t));
		#endif
			l++;
			s2d y1 = vdup(0.0);
			while (l<llim) {	// compute even and odd parts
				y1 = vdup(al[1])*cost*y0 + vdup(al[0])*y1;
V				dy = vdup(al[3])*cost*y1 + vdup(al[2])*y0;
Q				q = ((v2d*)Ql)[l];
S				s = ((v2d*)Sl)[l];
T				t = ((v2d*)Tl)[l];
Q				ro  += y1 * q;
S				po  += y1 * s;
S				dte += dy * s;
T				dpe -= dy * t;
T				to  -= y1 * t;
		#if _GCC_VEC_
Q				rox += y1 * vxchg(q);
S				dpo = subadd(dpo, y1 * vxchg(s));
S				te  = subadd(te, dy * vxchg(s));
T				pe  = subadd(pe, dy * vxchg(t));
T				dto = subadd(dto, y1 * vxchg(t));
		#endif
QX				al+=2;
V				al+=4;
				l++;
				y0   = vdup(al[1])*cost*y1 + vdup(al[0])*y0;
V				dy = vdup(al[3])*cost*y0 + vdup(al[2])*y1;
Q				q = ((v2d*)Ql)[l];
S				s = ((v2d*)Sl)[l];
T				t = ((v2d*)Tl)[l];
Q				re  += y0 * q;
S				pe  += y0 * s;
S				dto += dy * s;
T				dpo -= dy * t;
T				te  -= y0 * t;
		#if _GCC_VEC_
Q				rex += y0 * vxchg(q);
S				dpe = subadd(dpe, y0 * vxchg(s));
S				to  = subadd(to, dy * vxchg(s));
T				po  = subadd(po, dy * vxchg(t));
T				dte = subadd(dte, y0 * vxchg(t));
		#endif
QX				al+=2;
V				al+=4;
				l++;
			}
			if (l==llim) {
				y1 = vdup(al[1])*cost*y0 + vdup(al[0])*y1;
V				dy = vdup(al[3])*cost*y1 + vdup(al[2])*y0;
Q				q = ((v2d*)Ql)[l];
S				s = ((v2d*)Sl)[l];
T				t = ((v2d*)Tl)[l];
Q				ro  += y1 * q;
S				po  += y1 * s;
S				dte += dy * s;
T				dpe -= dy * t;
T				to  -= y1 * t;
		#if _GCC_VEC_
Q				rox += y1 * vxchg(q);
S				dpo = subadd(dpo, y1 * vxchg(s));
S				te  = subadd(te, dy * vxchg(s));
T				pe  = subadd(pe, dy * vxchg(t));
T				dto = subadd(dto, y1 * vxchg(t));
		#endif
			}
3			cost  = *((s2d*)(st+k));
3			cost *= vdup(1.0/m);
3			re *= cost;		ro *= cost;
		#if _GCC_VEC_
3			rex *= cost;		rox *= cost;
Q			((complex double *)BrF)[k] = vlo_to_dbl(re+ro) +I*vlo_to_dbl(rex+rox);
Q			((complex double *)BrF)[k+1] = vhi_to_dbl(rex+rox) +I*vhi_to_dbl(re+ro);
Q			((complex double *)BrF)[NLAT-k-2] = vhi_to_dbl(rex-rox) +I*vhi_to_dbl(re-ro);
Q			((complex double *)BrF)[NLAT-k-1] = vlo_to_dbl(re-ro) +I*vlo_to_dbl(rex-rox);

V			((complex double *)BtF)[k] = -I*vlo_to_dbl(te+to) + vlo_to_dbl(dte+dto);
V			((complex double *)BtF)[k+1] = vhi_to_dbl(te+to) + I*vhi_to_dbl(dte+dto);
V			((complex double *)BtF)[NLAT-k-2] = vhi_to_dbl(te-to) + I*vhi_to_dbl(dte-dto);
V			((complex double *)BtF)[NLAT-k-1] = -I*vlo_to_dbl(te-to) + vlo_to_dbl(dte-dto);

V			((complex double *)BpF)[k] = I*vlo_to_dbl(pe+po) + vlo_to_dbl(dpe+dpo);
V			((complex double *)BpF)[k+1] = -vhi_to_dbl(pe+po) + I*vhi_to_dbl(dpe+dpo);
V			((complex double *)BpF)[NLAT-k-2] = -vhi_to_dbl(pe-po) + I*vhi_to_dbl(dpe-dpo);
V			((complex double *)BpF)[NLAT-k-1] = I*vlo_to_dbl(pe-po) + vlo_to_dbl(dpe-dpo);
			k += 2;
		#else
Q			BrF[k] = (re+ro);
Q			BrF[NLAT-k-1] = (re-ro);
V			BtF[k] = addi(dte+dto, -(te+to));		// Bt = dS/dt       + I.m/sint *T
VB			BtF[NLAT-1-k] = addi(dte-dto, -(te-to));
V			BpF[k] = addi(dpe+dpo, pe+po);		// Bp = I.m/sint * S - dT/dt
VB			BpF[NLAT-1-k] = addi(dpe-dpo, pe-po);
			k++;
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

Q	#undef BR0
V	#undef BT0
V	#undef BP0
# }
