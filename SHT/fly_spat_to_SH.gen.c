/*
 * Copyright (c) 2010-2013 Centre National de la Recherche Scientifique.
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

QX	static void GEN3(spat_to_SH_fly,NWAY,SUFFIX)(shtns_cfg shtns, double *Vr, complex double *Qlm, long int llim) {
VX	static void GEN3(spat_to_SHsphtor_fly,NWAY,SUFFIX)(shtns_cfg shtns, double *Vt, double *Vp, complex double *Slm, complex double *Tlm, long int llim) {
3	static void GEN3(spat_to_SHqst_fly,NWAY,SUFFIX)(shtns_cfg shtns, double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm, long int llim) {

Q	s2d *BrF;		// contains the Fourier transformed data
V	s2d *BtF, *BpF;	// contains the Fourier transformed data
	double *alm, *al;
	s2d *wg, *ct, *st;
V	double *l_2;
	long int nk, k, vnlat, l,m;
  #ifndef SHT_AXISYM
	s2d sgn_mask;
V	double m_1;
	unsigned imlim, im;
  #endif
  #if _GCC_VEC_
Q	rnd qq[2*llim];
V	rnd ss[2*llim];
V	rnd tt[2*llim];
  #else
Q	double qq[llim];
V	double ss[llim];
V	double tt[llim];
  #endif

Q	s2d rer[(NLAT_2+VSIZE-1)/VSIZE+NWAY*(VSIZE2/VSIZE)-1];
Q	s2d ror[(NLAT_2+VSIZE-1)/VSIZE+NWAY*(VSIZE2/VSIZE)-1];
V	s2d ter[(NLAT_2+VSIZE-1)/VSIZE+NWAY*(VSIZE2/VSIZE)-1];
V	s2d tor[(NLAT_2+VSIZE-1)/VSIZE+NWAY*(VSIZE2/VSIZE)-1];
V	s2d per[(NLAT_2+VSIZE-1)/VSIZE+NWAY*(VSIZE2/VSIZE)-1];
V	s2d por[(NLAT_2+VSIZE-1)/VSIZE+NWAY*(VSIZE2/VSIZE)-1];
  #ifndef SHT_AXISYM
Q	s2d rei[(NLAT_2+VSIZE-1)/VSIZE+NWAY*(VSIZE2/VSIZE)-1];
Q	s2d roi[(NLAT_2+VSIZE-1)/VSIZE+NWAY*(VSIZE2/VSIZE)-1];
V	s2d tei[(NLAT_2+VSIZE-1)/VSIZE+NWAY*(VSIZE2/VSIZE)-1];
V	s2d toi[(NLAT_2+VSIZE-1)/VSIZE+NWAY*(VSIZE2/VSIZE)-1];
V	s2d pei[(NLAT_2+VSIZE-1)/VSIZE+NWAY*(VSIZE2/VSIZE)-1];
V	s2d poi[(NLAT_2+VSIZE-1)/VSIZE+NWAY*(VSIZE2/VSIZE)-1];
  #endif

Q	BrF = (s2d *) Vr;
V	BtF = (s2d *) Vt;	BpF = (s2d *) Vp;
  #ifndef SHT_AXISYM
	if (shtns->fftc_mode >= 0) {
	    if (shtns->fftc_mode == 0) {	// in-place
Q			fftw_execute_dft(shtns->fftc,(complex double*)BrF, (complex double*)BrF);
V			fftw_execute_dft(shtns->fftc,(complex double*)BtF, (complex double*)BtF);
V			fftw_execute_dft(shtns->fftc,(complex double*)BpF, (complex double*)BpF);
		} else {	// alloc memory for the transpose FFT
			unsigned long nv = shtns->nspat /VSIZE;
QX			BrF = (s2d*) VMALLOC( nv * sizeof(s2d) );
VX			BtF = (s2d*) VMALLOC( 2*nv * sizeof(s2d) );
VX			BpF = BtF + nv;
3			BrF = (s2d*) VMALLOC( 3*nv * sizeof(s2d) );
3			BtF = BrF + nv;		BpF = BtF + nv;
Q			fftw_execute_split_dft(shtns->fftc, Vr+NPHI, Vr, ((double*)BrF)+1, ((double*)BrF));
V			fftw_execute_split_dft(shtns->fftc, Vt+NPHI, Vt, ((double*)BtF)+1, ((double*)BtF));
V			fftw_execute_split_dft(shtns->fftc, Vp+NPHI, Vp, ((double*)BpF)+1, ((double*)BpF));
	    }
	}
	imlim = MTR;
	#ifdef SHT_VAR_LTR
		if (imlim*MRES > (unsigned) llim) imlim = ((unsigned) llim)/MRES;		// 32bit mul and div should be faster
	#endif
  #endif

	vnlat = ((unsigned) NLAT)/VSIZE;		// vector size.
	nk = NLAT_2;	// copy NLAT_2 to a local variable for faster access (inner loop limit)
	wg = (s2d*) shtns->wg;		ct = (s2d*) shtns->ct;		st = (s2d*) shtns->st;
	#if _GCC_VEC_
	  nk = ((unsigned) nk+(VSIZE2-1))/VSIZE2;
	#endif
V	l_2 = shtns->l_2;
	//	im = 0;		// dzl.p = 0.0 : and evrything is REAL
		k=0;
		alm = shtns->blm[0];
Q		s2d r0 = vdup(0.0);
		do {	// compute symmetric and antisymmetric parts. (do not weight here, it is cheaper to weight y0)
Q			s2d a = BrF[k];		s2d b = vxchg(BrF[vnlat-1-k]);
Q			rer[k] = a+b;		ror[k] = a-b;
Q			r0 += (a+b)*wg[k];
V			s2d c = BtF[k];		s2d d = vxchg(BtF[vnlat-1-k]);
V			ter[k] = c+d;		tor[k] = c-d;
V			s2d e = BpF[k];		s2d f = vxchg(BpF[vnlat-1-k]);
V			per[k] = e+f;		por[k] = e-f;
		} while(++k < nk*(VSIZE2/VSIZE));
		for (k=nk*(VSIZE2/VSIZE); k<(nk+NWAY)*(VSIZE2/VSIZE)-1; ++k) {
Q			rer[k] = vdup(0.0);		ror[k] = vdup(0.0);
V			ter[k] = vdup(0.0);		tor[k] = vdup(0.0);
V			per[k] = vdup(0.0);		por[k] = vdup(0.0);
		}
		#if _GCC_VEC_
Q			Qlm[0] = (vlo_to_dbl(r0) + vhi_to_dbl(r0)) * alm[0];					// l=0 is done.
		#else
Q			Qlm[0] = r0 * alm[0];				// l=0 is done.
		#endif
V		Slm[0] = 0.0;		Tlm[0] = 0.0;		// l=0 is zero for the vector transform.
		k = 0;
		for (l=0;l<llim;++l) {
Q			qq[l] = vall(0.0);
V			ss[l] = vall(0.0);		tt[l] = vall(0.0);
		}
		do {
			al = alm;
			rnd cost[NWAY], y0[NWAY], y1[NWAY];
V			rnd sint[NWAY], dy0[NWAY], dy1[NWAY];
Q			rnd rerk[NWAY], rork[NWAY];		// help the compiler to cache into registers.
V			rnd terk[NWAY], tork[NWAY], perk[NWAY], pork[NWAY];
			for (int j=0; j<NWAY; ++j) {
				cost[j] = vread(ct, k+j);
				y0[j] = vall(al[0]) * vread(wg, k+j);		// weight of Gauss quadrature appears here
V				dy0[j] = vall(0.0);
V				sint[j] = -vread(st, k+j);
				y1[j] =  (vall(al[1])*y0[j]) * cost[j];
V				dy1[j] = (vall(al[1])*y0[j]) * sint[j];
Q				rerk[j] = vread(rer, k+j);		rork[j] = vread(ror, k+j);		// cache into registers.
V				terk[j] = vread(ter, k+j);		tork[j] = vread(tor, k+j);
V				perk[j] = vread(per, k+j);		pork[j] = vread(por, k+j);
			}
			al+=2;	l=1;
			while(l<llim) {
				for (int j=0; j<NWAY; ++j) {
V					dy0[j] = vall(al[1])*(cost[j]*dy1[j] + y1[j]*sint[j]) + vall(al[0])*dy0[j];
					y0[j]  = vall(al[1])*(cost[j]*y1[j]) + vall(al[0])*y0[j];
				}
				for (int j=0; j<NWAY; ++j) {
Q					qq[l-1] += y1[j] * rork[j];
V					ss[l-1] += dy1[j] * terk[j];
V					tt[l-1] -= dy1[j] * perk[j];
				}
				for (int j=0; j<NWAY; ++j) {
V					dy1[j] = vall(al[3])*(cost[j]*dy0[j] + y0[j]*sint[j]) + vall(al[2])*dy1[j];
					y1[j]  = vall(al[3])*(cost[j]*y0[j]) + vall(al[2])*y1[j];
				}
				for (int j=0; j<NWAY; ++j) {
Q					qq[l] += y0[j] * rerk[j];
V					ss[l] += dy0[j] * tork[j];
V					tt[l] -= dy0[j] * pork[j];
				}
				al+=4;	l+=2;
			}
			if (l==llim) {
				for (int j=0; j<NWAY; ++j) {
Q					qq[l-1] += y1[j] * rork[j];
V					ss[l-1] += dy1[j] * terk[j];
V					tt[l-1] -= dy1[j] * perk[j];
				}
			}
			k+=NWAY;
		} while (k < nk);
		for (l=1; l<=llim; ++l) {
			#if _GCC_VEC_
			  #ifdef __AVX__
Q			    s2d qa = _mm256_castpd256_pd128(qq[l-1]) + _mm256_extractf128_pd(qq[l-1],1);
Q				Qlm[l] = vlo_to_dbl(qa) + vhi_to_dbl(qa);
V				rnd a = _mm256_hadd_pd(ss[l-1], tt[l-1]);
V				s2d b = (_mm256_castpd256_pd128(a) + _mm256_extractf128_pd(a,1)) * vdup(l_2[l]);
V				Slm[l] = vlo_to_dbl(b);		Tlm[l] = vhi_to_dbl(b);
			  #elif defined __SSE3__
Q				Qlm[l] = vlo_to_dbl(qq[l-1]) + vhi_to_dbl(qq[l-1]);
V				s2d a = _mm_hadd_pd(ss[l-1], tt[l-1]) * vdup(l_2[l]);
V				Slm[l] = vlo_to_dbl(a);		Tlm[l] = vhi_to_dbl(a);
			  #else
Q				Qlm[l] = vlo_to_dbl(qq[l-1]) + vhi_to_dbl(qq[l-1]);
V				Slm[l] = (vlo_to_dbl(ss[l-1]) + vhi_to_dbl(ss[l-1]))*l_2[l];
V				Tlm[l] = (vlo_to_dbl(tt[l-1]) + vhi_to_dbl(tt[l-1]))*l_2[l];
			  #endif
			#else
Q				Qlm[l] = qq[l-1];
V				Slm[l] = ss[l-1]*l_2[l];		Tlm[l] = tt[l-1]*l_2[l];
			#endif
		}
		#ifdef SHT_VAR_LTR
			for (l=llim+1; l<= LMAX; ++l) {
Q				Qlm[l] = 0.0;
V				Slm[l] = 0.0;		Tlm[l] = 0.0;
			}
		#endif

  #ifndef SHT_AXISYM
	#if _GCC_VEC_
		sgn_mask = SIGN_MASK_HI;		// build mask to change sign of hi value using xorpd (used in CFFT_TO_2REAL)
	#endif
	for (k=nk*(VSIZE2/VSIZE); k<(nk+NWAY)*(VSIZE2/VSIZE)-1; ++k) {		// never written, so this is now done for all m's (real parts already zero)
Q		rei[k] = vdup(0.0);		roi[k] = vdup(0.0);
V		tei[k] = vdup(0.0);		toi[k] = vdup(0.0);
V		pei[k] = vdup(0.0);		poi[k] = vdup(0.0);
	}
	for (im=1;im<=imlim;++im) {
Q		BrF += vnlat;
V		BtF += vnlat;	BpF += vnlat;
		l = shtns->tm[im] / VSIZE2;
		alm = shtns->blm[im];
		m = im*MRES;
 	#if _GCC_VEC_
 		k=(VSIZE2/VSIZE)*l;
		do {	// compute symmetric and antisymmetric parts, and reorganize data.
			s2d nr, ni, sr, si, tn, ts;
3			s2d sin = st[k];
Q			ni = BrF[k];				nr = BrF[(NPHI-2*im)*vnlat + k];
Q			si = BrF[vnlat-1 - k];		sr = BrF[(NPHI-2*im)*vnlat +vnlat-1-k];
Q			CFFT_TO_2REAL(nr, ni, sr, si, sgn_mask)
QX			rer[k] = nr+sr;		ror[k] = nr-sr;		rei[k] = ni+si;		roi[k] = ni-si;
3			rer[k] = (nr+sr)*sin;		ror[k] = (nr-sr)*sin;		rei[k] = (ni+si)*sin;		roi[k] = (ni-si)*sin;

V			ni = BtF[k];				nr = BtF[(NPHI-2*im)*vnlat + k];
V			si = BtF[vnlat-1 - k];		sr = BtF[(NPHI-2*im)*vnlat +vnlat-1-k];
V			CFFT_TO_2REAL(nr, ni, sr, si, sgn_mask)
V			ter[k] = nr+sr;		tor[k] = nr-sr;		tei[k] = ni+si;		toi[k] = ni-si;

V			ni = BpF[k];				nr = BpF[(NPHI-2*im)*vnlat + k];
V			si = BpF[vnlat-1 - k];		sr = BpF[(NPHI-2*im)*vnlat +vnlat-1-k];
V			CFFT_TO_2REAL(nr, ni, sr, si, sgn_mask)
V			per[k] = nr+sr;		por[k] = nr-sr;		pei[k] = ni+si;		poi[k] = ni-si;
		} while(++k < nk*(VSIZE2/VSIZE));
	#else
		k = (l>>1)*2;		// k must be even here.
		do {	// compute symmetric and antisymmetric parts, and reorganize data.
			double an, bn, ani, bni, bs, as, bsi, asi, t;
3			double sina = st[k];	double sinb = st[k+1];
Q			ani = BrF[k];		bni = BrF[k+1];		// north
Q			an = BrF[(NPHI-2*im)*vnlat + k];	bn = BrF[(NPHI-2*im)*vnlat + k+1];
Q			t = ani-an;	an += ani;		ani = bn-bni;		bn += bni;		bni = t;
3			an *= sina;		ani*= sina;		bn *= sinb;		bni *= sinb;
Q			bsi = BrF[vnlat-2 -k];	asi = BrF[vnlat-1 -k];	// south
Q			bs = BrF[(NPHI-2*im)*vnlat +vnlat-2-k];		as = BrF[(NPHI-2*im)*vnlat +vnlat-1-k];
Q			t = bsi-bs;		bs += bsi;		bsi = as-asi;		as += asi;		asi = t;
3			as *= sina;		asi*= sina;		bs *= sinb;		bsi *= sinb;
Q			rer[k] = an+as;		rei[k] = ani+asi;		rer[k+1] = bn+bs;		rei[k+1] = bni+bsi;
Q			ror[k] = an-as;		roi[k] = ani-asi;		ror[k+1] = bn-bs;		roi[k+1] = bni-bsi;

V			ani = BtF[k];		bni = BtF[k+1];		// north
V			an = BtF[(NPHI-2*im)*vnlat + k];	bn = BtF[(NPHI-2*im)*vnlat + k+1];
V			t = ani-an;	an += ani;		ani = bn-bni;		bn += bni;		bni = t;
V			bsi = BtF[vnlat-2 -k];	asi = BtF[vnlat-1 -k];	// south
V			bs = BtF[(NPHI-2*im)*vnlat +vnlat-2-k];		as = BtF[(NPHI-2*im)*vnlat +vnlat-1-k];
V			t = bsi-bs;		bs += bsi;		bsi = as-asi;		as += asi;		asi = t;
V			ter[k] = an+as;		tei[k] = ani+asi;		ter[k+1] = bn+bs;		tei[k+1] = bni+bsi;
V			tor[k] = an-as;		toi[k] = ani-asi;		tor[k+1] = bn-bs;		toi[k+1] = bni-bsi;

V			ani = BpF[k];		bni = BpF[k+1];		// north
V			an = BpF[(NPHI-2*im)*vnlat + k];	bn = BpF[(NPHI-2*im)*vnlat + k+1];
V			t = ani-an;	an += ani;		ani = bn-bni;		bn += bni;		bni = t;
V			bsi = BpF[vnlat-2 -k];	asi = BpF[vnlat-1 -k];	// south
V			bs = BpF[(NPHI-2*im)*vnlat +vnlat-2-k];		as = BpF[(NPHI-2*im)*vnlat +vnlat-1-k];
V			t = bsi-bs;		bs += bsi;		bsi = as-asi;		as += asi;		asi = t;
V			per[k] = an+as;		pei[k] = ani+asi;		per[k+1] = bn+bs;		pei[k+1] = bni+bsi;
V			por[k] = an-as;		poi[k] = ani-asi;		por[k+1] = bn-bs;		poi[k+1] = bni-bsi;
			k+=2;
 		} while (k<nk);
	#endif
V		m_1 = 1.0/m;
		k=l;
		#if _GCC_VEC_
Q			rnd* q = qq;
V			rnd* s = ss;		rnd* t = tt;
		#else
			l = LiM(shtns, m, im);
Q			double* q = (double *) &Qlm[l];
V			double* s = (double *) &Slm[l];
V			double* t = (double *) &Tlm[l];
		#endif
		for (l=llim-m; l>=0; l--) {
Q			q[0] = vall(0.0);		q[1] = vall(0.0);		q+=2;
V			s[0] = vall(0.0);		s[1] = vall(0.0);		s+=2;
V			t[0] = vall(0.0);		t[1] = vall(0.0);		t+=2;
		}
		do {
		#if _GCC_VEC_
Q			rnd* q = qq;
V			rnd* s = ss;		rnd* t = tt;
		#else
			l = LiM(shtns, m, im);
Q			double* q = (double *) &Qlm[l];
V			double* s = (double *) &Slm[l];
V			double* t = (double *) &Tlm[l];
		#endif
			al = alm;
			rnd cost[NWAY], y0[NWAY], y1[NWAY];
V			rnd st2[NWAY], dy0[NWAY], dy1[NWAY];
Q			rnd rerk[NWAY], reik[NWAY], rork[NWAY], roik[NWAY];		// help the compiler to cache into registers.
V			rnd terk[NWAY], teik[NWAY], tork[NWAY], toik[NWAY];
V			rnd perk[NWAY], peik[NWAY], pork[NWAY], poik[NWAY];
			for (int j=0; j<NWAY; ++j) {
				cost[j] = vread(st, k+j);
				y0[j] = vall(0.5);
V				st2[j] = cost[j]*cost[j]*vall(-m_1);
V				y0[j] *= vall(m);		// for the vector transform, compute ylm*m/sint
			}
Q			l=m;
V			l=m-1;
			long int ny = 0;	// exponent to extend double precision range.
		if ((int)llim <= SHT_L_RESCALE_FLY) {
			do {		// sin(theta)^m
				if (l&1) for (int j=0; j<NWAY; ++j) y0[j] *= cost[j];
				for (int j=0; j<NWAY; ++j) cost[j] *= cost[j];
			} while(l >>= 1);
		} else {
			long int nsint = 0;
			do {		// sin(theta)^m		(use rescaling to avoid underflow)
				if (l&1) {
					for (int j=0; j<NWAY; ++j) y0[j] *= cost[j];
					ny += nsint;
					if (vlo(y0[0]) < (SHT_ACCURACY+1.0/SHT_SCALE_FACTOR)) {
						ny--;
						for (int j=0; j<NWAY; ++j) y0[j] *= vall(SHT_SCALE_FACTOR);
					}
				}
				for (int j=0; j<NWAY; ++j) cost[j] *= cost[j];
				nsint += nsint;
				if (vlo(cost[0]) < 1.0/SHT_SCALE_FACTOR) {
					nsint--;
					for (int j=0; j<NWAY; ++j) cost[j] *= vall(SHT_SCALE_FACTOR);
				}
			} while(l >>= 1);
		}
			for (int j=0; j<NWAY; ++j) {
				y0[j] *= vall(al[0]);
				cost[j] = vread(ct, k+j);
V				dy0[j] = cost[j]*y0[j];
				y1[j]  = (vall(al[1])*y0[j]) *cost[j];
V				dy1[j] = (vall(al[1])*y0[j]) *(cost[j]*cost[j] + st2[j]);
			}
			l=m;	al+=2;
			while ((ny<0) && (l<llim)) {		// ylm treated as zero and ignored if ny < 0
				for (int j=0; j<NWAY; ++j) {
V					dy0[j] = vall(al[1])*(cost[j]*dy1[j] + y1[j]*st2[j]) + vall(al[0])*dy0[j];
					y0[j] = vall(al[1])*(cost[j]*y1[j]) + vall(al[0])*y0[j];
				}
				for (int j=0; j<NWAY; ++j) {
V					dy1[j] = vall(al[3])*(cost[j]*dy0[j] + y0[j]*st2[j]) + vall(al[2])*dy1[j];				
					y1[j] = vall(al[3])*(cost[j]*y0[j]) + vall(al[2])*y1[j];
				}
				l+=2;	al+=4;
				if (fabs(vlo(y0[NWAY-1])) > SHT_ACCURACY*SHT_SCALE_FACTOR + 1.0) {		// rescale when value is significant
					++ny;
					for (int j=0; j<NWAY; ++j) {
						y0[j] *= vall(1.0/SHT_SCALE_FACTOR);		y1[j] *= vall(1.0/SHT_SCALE_FACTOR);
V						dy0[j] *= vall(1.0/SHT_SCALE_FACTOR);		dy1[j] *= vall(1.0/SHT_SCALE_FACTOR);
					}
				}
			}
		  if (ny == 0) {
Q			q+=2*(l-m);
V			s+=2*(l-m);		t+=2*(l-m);
			for (int j=0; j<NWAY; ++j) {	// prefetch
				y0[j] *= vread(wg, k+j);		y1[j] *= vread(wg, k+j);		// weight appears here (must be after the previous accuracy loop).
V				dy0[j] *= vread(wg, k+j);		dy1[j] *= vread(wg, k+j);
Q				rerk[j] = vread( rer, k+j);		reik[j] = vread( rei, k+j);		rork[j] = vread( ror, k+j);		roik[j] = vread( roi, k+j);
V				terk[j] = vread( ter, k+j);		teik[j] = vread( tei, k+j);		tork[j] = vread( tor, k+j);		toik[j] = vread( toi, k+j);
V				perk[j] = vread( per, k+j);		peik[j] = vread( pei, k+j);		pork[j] = vread( por, k+j);		poik[j] = vread( poi, k+j);
			}
			while (l<llim) {	// compute even and odd parts
Q				for (int j=0; j<NWAY; ++j)	q[0] += y0[j] * rerk[j];		// real even
Q				for (int j=0; j<NWAY; ++j)	q[1] += y0[j] * reik[j];		// imag even
V				for (int j=0; j<NWAY; ++j)	s[0] += dy0[j] * tork[j]  + y0[j] * peik[j];
V				for (int j=0; j<NWAY; ++j)	s[1] += dy0[j] * toik[j]  - y0[j] * perk[j];
V				for (int j=0; j<NWAY; ++j)	t[0] -= dy0[j] * pork[j]  - y0[j] * teik[j];
V				for (int j=0; j<NWAY; ++j)	t[1] -= dy0[j] * poik[j]  + y0[j] * terk[j];
				for (int j=0; j<NWAY; ++j) {
V					dy0[j] = vall(al[1])*(cost[j]*dy1[j] + y1[j]*st2[j]) + vall(al[0])*dy0[j];
					y0[j] = vall(al[1])*(cost[j]*y1[j]) + vall(al[0])*y0[j];
				}
Q				for (int j=0; j<NWAY; ++j)	q[2] += y1[j] * rork[j];		// real odd
Q				for (int j=0; j<NWAY; ++j)	q[3] += y1[j] * roik[j];		// imag odd
V				for (int j=0; j<NWAY; ++j)	s[2] += dy1[j] * terk[j]  + y1[j] * poik[j];
V				for (int j=0; j<NWAY; ++j)	s[3] += dy1[j] * teik[j]  - y1[j] * pork[j];
V				for (int j=0; j<NWAY; ++j)	t[2] -= dy1[j] * perk[j]  - y1[j] * toik[j];
V				for (int j=0; j<NWAY; ++j)	t[3] -= dy1[j] * peik[j]  + y1[j] * tork[j];
Q				q+=4;
V				s+=4;	t+=4;
				for (int j=0; j<NWAY; ++j) {
V					dy1[j] = vall(al[3])*(cost[j]*dy0[j] + y0[j]*st2[j]) + vall(al[2])*dy1[j];
					y1[j] = vall(al[3])*(cost[j]*y0[j]) + vall(al[2])*y1[j];
				}
				l+=2;	al+=4;
			}
			if (l==llim) {
Q				for (int j=0; j<NWAY; ++j)	q[0] += y0[j] * rerk[j];		// real even
Q				for (int j=0; j<NWAY; ++j)	q[1] += y0[j] * reik[j];		// imag even
V				for (int j=0; j<NWAY; ++j)	s[0] += dy0[j] * tork[j]  + y0[j] * peik[j];
V				for (int j=0; j<NWAY; ++j)	s[1] += dy0[j] * toik[j]  - y0[j] * perk[j];
V				for (int j=0; j<NWAY; ++j)	t[0] -= dy0[j] * pork[j]  - y0[j] * teik[j];
V				for (int j=0; j<NWAY; ++j)	t[1] -= dy0[j] * poik[j]  + y0[j] * terk[j];
			}
		  }
			k+=NWAY;
		} while (k < nk);
		l = LiM(shtns, m, im);
Q		v2d *Ql = (v2d*) &Qlm[l];
V		v2d *Sl = (v2d*) &Slm[l];
V		v2d *Tl = (v2d*) &Tlm[l];
		#if _GCC_VEC_
			for (l=0; l<=llim-m; ++l) {
			#ifdef __AVX__
Q				rnd qa = _mm256_hadd_pd(qq[2*l], qq[2*l+1]);
QX				Ql[l] = _mm256_castpd256_pd128(qa) + _mm256_extractf128_pd(qa,1);
3				Ql[l] = (_mm256_castpd256_pd128(qa) + _mm256_extractf128_pd(qa,1)) * vdup(m_1);
V				rnd sa = _mm256_hadd_pd(ss[2*l], ss[2*l+1]);		rnd ta = _mm256_hadd_pd(tt[2*l], tt[2*l+1]);
V				sa = (_mm256_permute2f128_pd(sa, ta, 0x02) + _mm256_permute2f128_pd(sa, ta, 0x13)) * vall(l_2[l+m]);
V				Sl[l] = _mm256_extractf128_pd(sa,1);	Tl[l] = _mm256_castpd256_pd128(sa);
			#elif defined __SSE3__
QX				Ql[l] = _mm_hadd_pd(qq[2*l], qq[2*l+1]);
3				Ql[l] = _mm_hadd_pd(qq[2*l], qq[2*l+1]) * vdup(m_1);
V				Sl[l] = _mm_hadd_pd(ss[2*l], ss[2*l+1]) * vdup(l_2[l+m]);
V				Tl[l] = _mm_hadd_pd(tt[2*l], tt[2*l+1]) * vdup(l_2[l+m]);
			#else
Q				s2d qa = _mm_unpacklo_pd(qq[2*l], qq[2*l+1]);		s2d qb = _mm_unpackhi_pd(qq[2*l], qq[2*l+1]);
QX				Ql[l] = qa + qb;
3				Ql[l] = (qa + qb)*vdup(m_1);
V				s2d sa = _mm_unpacklo_pd(ss[2*l], ss[2*l+1]);		s2d sb = _mm_unpackhi_pd(ss[2*l], ss[2*l+1]);
V				s2d ta = _mm_unpacklo_pd(tt[2*l], tt[2*l+1]);		s2d tb = _mm_unpackhi_pd(tt[2*l], tt[2*l+1]);
V				Sl[l] = (sa + sb)*vdup(l_2[l+m]);
V				Tl[l] = (ta + tb)*vdup(l_2[l+m]);
			#endif
			}
		#else
V			for (l=0; l<=llim-m; ++l) {
3				Ql[l] *= m_1;
V				Sl[l] *= l_2[l+m];
V				Tl[l] *= l_2[l+m];
V			}
		#endif
		#ifdef SHT_VAR_LTR
			for (l=llim+1-m; l<=LMAX-m; ++l) {
Q				Ql[l] = vdup(0.0);
V				Sl[l] = vdup(0.0);		Tl[l] = vdup(0.0);
			}
		#endif
	}
	#ifdef SHT_VAR_LTR
	if (imlim < MMAX) {
		im = imlim+1;
		l = LiM(shtns, im*MRES, im);
		do {
Q			((v2d*)Qlm)[l] = vdup(0.0);
V			((v2d*)Slm)[l] = vdup(0.0);		((v2d*)Tlm)[l] = vdup(0.0);
		} while(++l < shtns->nlm);
	}
	#endif

  	if (shtns->fftc_mode > 0) {		// free memory
Q	    VFREE(BrF - vnlat*imlim);
VX	    VFREE(BtF - vnlat*imlim);	// this frees also BpF.
	}
  #endif

  }
