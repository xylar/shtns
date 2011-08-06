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

/// Truncation and spatial discretization are defined by \ref shtns_create and \ref shtns_set_grid_*
/// \param[in] shtns = a configuration created by \ref shtns_create with a grid set by shtns_set_grid_*
Q/// \param[in] Vr = spatial scalar field : double array.
V/// \param[in] Vt, Vp = spatial (theta, phi) vector components : double arrays.
Q/// \param[out] Qlm = spherical harmonics coefficients :
Q/// complex double arrays of size shtns->nlm.
V/// \param[out] Slm,Tlm = spherical harmonics coefficients of \b Spheroidal and \b Toroidal scalars :
V/// complex double arrays of size shtns->nlm.
  #ifdef SHT_VAR_LTR
/// \param[in] llim = specify maximum degree of spherical harmonic. ltr must be at most LMAX, and all spherical harmonic degree higher than ltr are set to zero. 
  #endif

QX	void GEN3(spat_to_SH_fly,NWAY,SUFFIX)(shtns_cfg shtns, double *Vr, complex double *Qlm, long int llim) {
VX	void GEN3(spat_to_SHsphtor_fly,NWAY,SUFFIX)(shtns_cfg shtns, double *Vt, double *Vp, complex double *Slm, complex double *Tlm, long int llim) {
3	void GEN3(spat_to_SHqst_fly,NWAY,SUFFIX)(shtns_cfg shtns, double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm, long int llim) {

Q	complex double *BrF;		// contains the Fourier transformed data
V	complex double *BtF, *BpF;	// contains the Fourier transformed data
	double *al, *wg, *ct, *st;
V	double *l_2;
	long int ni, k;
	long int i,i0, m,l;
	long int imlim, im;
  #if _GCC_VEC_
Q	s2d qq[2*llim];
V	s2d ss[2*llim];
V	s2d tt[2*llim];
  #else
Q	double qq[llim+1];
V	double ss[llim+1];
V	double tt[llim+1];
  #endif

  #ifndef SHT_AXISYM
Q	double rei[NLAT_2+2*NWAY] SSE;
Q	double rer[NLAT_2+2*NWAY] SSE;
Q	double ror[NLAT_2+2*NWAY] SSE;
Q	double roi[NLAT_2+2*NWAY] SSE;
V	double ter[NLAT_2+2*NWAY] SSE;
V	double tor[NLAT_2+2*NWAY] SSE;
V	double per[NLAT_2+2*NWAY] SSE;
V	double por[NLAT_2+2*NWAY] SSE;
V	double tei[NLAT_2+2*NWAY] SSE;
V	double toi[NLAT_2+2*NWAY] SSE;
V	double pei[NLAT_2+2*NWAY] SSE;
V	double poi[NLAT_2+2*NWAY] SSE;
  #else
Q	double rer[NLAT_2+2*NWAY] SSE;
Q	double ror[NLAT_2+2*NWAY] SSE;
V	double ter[NLAT_2+2*NWAY] SSE;
V	double tor[NLAT_2+2*NWAY] SSE;
V	double per[NLAT_2+2*NWAY] SSE;
V	double por[NLAT_2+2*NWAY] SSE;
  #endif

	ni = NLAT_2;	// copy NLAT_2 to a local variable for faster access (inner loop limit)
	imlim = MTR;
	#ifdef SHT_VAR_LTR
		if (MTR*MRES > (int) llim) imlim = ((int) llim)/MRES;		// 32bit mul and div should be faster
	#endif
Q	BrF = (complex double *) Vr;
V	BtF = (complex double *) Vt;	BpF = (complex double *) Vp;

  #ifndef SHT_AXISYM
	if (SHT_FFT > 0) {
	    if (SHT_FFT > 1) {		// alloc memory for the FFT
	    	unsigned long nspat = ((NPHI>>1) +1)*NLAT;
QX			BrF = fftw_malloc( nspat * sizeof(complex double) );
VX			BtF = fftw_malloc( 2* nspat * sizeof(complex double) );
VX			BpF = BtF + nspat;
3			BrF = fftw_malloc( 3* nspat * sizeof(complex double) );
3			BtF = BrF + nspat;		BpF = BtF + nspat;
	    }
Q		fftw_execute_dft_r2c(shtns->fft,Vr, BrF);
V		fftw_execute_dft_r2c(shtns->fft,Vt, BtF);
V		fftw_execute_dft_r2c(shtns->fft,Vp, BpF);
	}
  #endif

	wg = shtns->wg;		ct = shtns->ct;		st = shtns->st;
V	l_2 = shtns->l_2;
	im = 0;		// dzl.p = 0.0 : and evrything is REAL
		i=0;
Q		double r0 = 0.0;
  #ifndef SHT_AXISYM
		i0 = (NPHI==1) ? 1 : 2;		// stride of source data.
 		do {	// compute symmetric and antisymmetric parts.
 			double w = wg[i];
Q			double a = ((double*)BrF)[i*i0];		double b = ((double*)BrF)[(NLAT-1)*i0 -i*i0];
Q			ror[i] = (a-b)*w;		rer[i] = (a+b)*w;
Q			r0 += ((a+b)*w);
V			double c = ((double*)BtF)[i*i0];		double d = ((double*)BtF)[(NLAT-1)*i0 -i*i0];
V			ter[i] = (c+d)*w;		tor[i] = (c-d)*w;
V			double e = ((double*)BpF)[i*i0];		double f = ((double*)BpF)[(NLAT-1)*i0 -i*i0];
V			per[i] = (e+f)*w;		por[i] = (e-f)*w;
 			i++;
 		} while(i<ni);
  #else
 		do {	// compute symmetric and antisymmetric parts.
 			double np = wg[i]*NPHI;
Q			double a = ((double*)BrF)[i];		double b = ((double*)BrF)[NLAT-1-i];
Q			ror[i] = (a-b)*np;		rer[i] = (a+b)*np;
Q			r0 += ((a+b)*np);
V			double c = ((double*)BtF)[i];		double d = ((double*)BtF)[NLAT-1-i];
V			ter[i] = (c+d)*np;		tor[i] = (c-d)*np;
V			double e = ((double*)BpF)[i];		double f = ((double*)BpF)[NLAT-1-i];
V			per[i] = (e+f)*np;		por[i] = (e-f)*np;
			i++;
 		} while(i<ni);
  #endif
		do {
Q			rer[i] = 0.0;		ror[i] = 0.0;
V			ter[i] = 0.0;		tor[i] = 0.0;
V			per[i] = 0.0;		por[i] = 0.0;
			i++;
	#if _GCC_VEC_
		} while(i<ni+2*NWAY-1);
	#else
		} while(i<ni+NWAY-1);
	#endif
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
Q		Qlm[0] = r0 * shtns->bl0[0];					// l=0 is done.
V		Slm[0] = 0.0;		Tlm[0] = 0.0;		// l=0 is zero for the vector transform.
		k = 0;
		for (l=1;l<=llim;l++) {
Q			qq[l] = vdup(0.0);
V			ss[l] = vdup(0.0);		tt[l] = vdup(0.0);
		}
		do {
			al = shtns->bl0;
			s2d cost[NWAY], y0[NWAY], y1[NWAY];
V			s2d sint[NWAY], dy0[NWAY], dy1[NWAY];
Q			s2d rerk[NWAY], rork[NWAY];		// help the compiler to cache into registers.
V			s2d terk[NWAY], tork[NWAY], perk[NWAY], pork[NWAY];
			for (int j=0; j<NWAY; j++) {
				cost[j] = ((s2d*)(ct+k))[j];
V				sint[j] = ((s2d*)(st+k))[j];
				y0[j] = vdup(al[0]);
V				dy0[j] = vdup(0.0);
				y1[j] = vdup(al[0]*al[1]) * cost[j];
V				dy1[j] = -vdup(al[0]*al[1]) * sint[j];
Q				rerk[j] = ((s2d*)(rer+k))[j];		rork[j] = ((s2d*)(ror+k))[j];		// cache into registers.
V				terk[j] = ((s2d*)(ter+k))[j];		tork[j] = ((s2d*)(tor+k))[j];
V				perk[j] = ((s2d*)(per+k))[j];		pork[j] = ((s2d*)(por+k))[j];
			}
			al+=2;	l=1;
			while(l<llim) {
				for (int j=0; j<NWAY; j++) {
					y0[j]  = vdup(al[1])*cost[j]*y1[j] + vdup(al[0])*y0[j];
V					dy0[j] = vdup(al[1])*(cost[j]*dy1[j] - y1[j]*sint[j]) + vdup(al[0])*dy0[j];
				}
				for (int j=0; j<NWAY; j++) {
Q					qq[l] += y1[j] * rork[j];
Q					qq[l+1]   += y0[j] * rerk[j];
V					ss[l] += dy1[j] * terk[j];
V					tt[l] -= dy1[j] * perk[j];
V					ss[l+1] += dy0[j] * tork[j];
V					tt[l+1] -= dy0[j] * pork[j];
				}
				for (int j=0; j<NWAY; j++) {
					y1[j]  = vdup(al[3])*cost[j]*y0[j] + vdup(al[2])*y1[j];
V					dy1[j] = vdup(al[3])*(cost[j]*dy0[j] - y0[j]*sint[j]) + vdup(al[2])*dy1[j];
				}
				al+=4;	l+=2;
			}
			if (l==llim) {
				for (int j=0; j<NWAY; j++) {
Q					qq[l] += y1[j] * rork[j];
V					ss[l] += dy1[j] * terk[j];
V					tt[l] -= dy1[j] * perk[j];
				}
			}
		#if _GCC_VEC_
			k+=2*NWAY;
		#else
			k+=NWAY;
		#endif
		} while (k < ni);
		for (l=1; l<=llim; l++) {
			#if _GCC_VEC_
Q				Qlm[l] = vlo_to_dbl(qq[l]) + vhi_to_dbl(qq[l]);
V				Slm[l] = (vlo_to_dbl(ss[l]) + vhi_to_dbl(ss[l]))*l_2[l];
V				Tlm[l] = (vlo_to_dbl(tt[l]) + vhi_to_dbl(tt[l]))*l_2[l];
			#else
Q				Qlm[l] = qq[l];
V				Slm[l] = ss[l]*l_2[l];		Tlm[l] = tt[l]*l_2[l];
			#endif
		}
		#ifdef SHT_VAR_LTR
			for (l=llim+1; l<= LMAX; l++) {
Q				Qlm[l] = 0.0;
V				Slm[l] = 0.0;		Tlm[l] = 0.0;
			}
		#endif

  #ifndef SHT_AXISYM
	for (im=1;im<=imlim;im++) {
		i0 = shtns->tm[im];
		m = im*MRES;
	#if _GCC_VEC_
		i0=(i0>>1)*2;		// stay on a 16 byte boundary
	#endif
 		i=i0;
 		do {	// compute symmetric and antisymmetric parts, and reorganize data.
3			s2d sin = vdup(st[i]);
Q			v2d r0 = ((v2d *)BrF)[i];	v2d r1 = ((v2d *)BrF)[NLAT-1-i];
V			v2d t0 = ((v2d *)BtF)[i];	v2d t1 = ((v2d *)BtF)[NLAT-1-i];
V			v2d p0 = ((v2d *)BpF)[i];	v2d p1 = ((v2d *)BpF)[NLAT-1-i];
3			r0 *= sin;		r1 *= sin;
			#if _GCC_VEC_
V				ter[i] = vlo_to_dbl(t0+t1);		tei[i] = vhi_to_dbl(t0+t1);
V				tor[i] = vlo_to_dbl(t0-t1);		toi[i] = vhi_to_dbl(t0-t1);
V				per[i] = vlo_to_dbl(p0+p1);		pei[i] = vhi_to_dbl(p0+p1);
V				por[i] = vlo_to_dbl(p0-p1);		poi[i] = vhi_to_dbl(p0-p1);
Q				rer[i] = vlo_to_dbl(r0+r1);		rei[i] = vhi_to_dbl(r0+r1);
Q				ror[i] = vlo_to_dbl(r0-r1);		roi[i] = vhi_to_dbl(r0-r1);
			#else
V				ter[i] = creal(t0+t1);		tei[i] = cimag(t0+t1);
V				tor[i] = creal(t0-t1);		toi[i] = cimag(t0-t1);
V				per[i] = creal(p0+p1);		pei[i] = cimag(p0+p1);
V				por[i] = creal(p0-p1);		poi[i] = cimag(p0-p1);
Q				rer[i] = creal(r0+r1);		rei[i] = cimag(r0+r1);
Q				ror[i] = creal(r0-r1);		roi[i] = cimag(r0-r1);
			#endif
 			i++;
 		} while (i<ni);
		do {
Q			rer[i] = 0.0;		rei[i] = 0.0;		ror[i] = 0.0;		roi[i] = 0.0;
V			ter[i] = 0.0;		tei[i] = 0.0;		tor[i] = 0.0;		toi[i] = 0.0;
V			per[i] = 0.0;		pei[i] = 0.0;		por[i] = 0.0;		poi[i] = 0.0;
			i++;
	#if _GCC_VEC_
		} while (i < ni+2*NWAY-1);
	#else
		} while (i < ni+NWAY-1);
	#endif
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;

		k=i0;		// i0 must be even.
		#if _GCC_VEC_
Q			s2d* q = qq;
V			s2d* s = ss;		s2d* t = tt;
		#else
			l = LiM(shtns, m, im);
Q			double* q = (double *) &Qlm[l];
V			double* s = (double *) &Slm[l];
V			double* t = (double *) &Tlm[l];
		#endif
		for (l=llim-m; l>=0; l--) {
Q			q[0] = vdup(0.0);		q[1] = vdup(0.0);		q+=2;
V			s[0] = vdup(0.0);		s[1] = vdup(0.0);		s+=2;
V			t[0] = vdup(0.0);		t[1] = vdup(0.0);		t+=2;
		}
	  if (llim <= SHT_L_RESCALE_FLY) {
		#define Y0 y0[j]
		#define Y1 y1[j]
V		#define DY0 dy0[j]
V		#define DY1 dy1[j]
		do {
		#if _GCC_VEC_
Q			s2d* q = qq;
V			s2d* s = ss;		s2d* t = tt;
		#else
			l = LiM(shtns, m, im);
Q			double* q = (double *) &Qlm[l];
V			double* s = (double *) &Slm[l];
V			double* t = (double *) &Tlm[l];
		#endif
			al = shtns->blm[im];
			s2d cost[NWAY], y0[NWAY], y1[NWAY];
V			s2d st2[NWAY], dy0[NWAY], dy1[NWAY];
Q			s2d rerk[NWAY], reik[NWAY], rork[NWAY], roik[NWAY];		// help the compiler to cache into registers.
V			s2d terk[NWAY], teik[NWAY], tork[NWAY], toik[NWAY];
V			s2d perk[NWAY], peik[NWAY], pork[NWAY], poik[NWAY];
			for (int j=0; j<NWAY; j++) {
				cost[j] = ((s2d*)(st+k))[j];
				y0[j] = vdup(al[0]) * ((s2d*)(wg+k))[j];		// weight appears here.
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
				cost[j] = ((s2d*)(ct+k))[j];
V				dy0[j] = cost[j]*y0[j];
			}
			for (int j=0; j<NWAY; j++) {		// help the compiler to cache spatial data into registers.
Q				rerk[j] = ((s2d*)(rer+k))[j];		reik[j] = ((s2d*)(rei+k))[j];		rork[j] = ((s2d*)(ror+k))[j];		roik[j] = ((s2d*)(roi+k))[j];
V				terk[j] = ((s2d*)(ter+k))[j];		teik[j] = ((s2d*)(tei+k))[j];		tork[j] = ((s2d*)(tor+k))[j];		toik[j] = ((s2d*)(toi+k))[j];
V				perk[j] = ((s2d*)(per+k))[j];		peik[j] = ((s2d*)(pei+k))[j];		pork[j] = ((s2d*)(por+k))[j];		poik[j] = ((s2d*)(poi+k))[j];
			}
			l=m;
Q			for (int j=0; j<NWAY; j++)	q[0] += Y0 * rerk[j];		// real even
Q			for (int j=0; j<NWAY; j++)	q[1] += Y0 * reik[j];		// imag even
V			for (int j=0; j<NWAY; j++)	s[0] += DY0 * tork[j]  + Y0 * peik[j];
V			for (int j=0; j<NWAY; j++)	s[1] += DY0 * toik[j]  - Y0 * perk[j];
V			for (int j=0; j<NWAY; j++)	t[0] -= DY0 * pork[j]  - Y0 * teik[j];
V			for (int j=0; j<NWAY; j++)	t[1] -= DY0 * poik[j]  + Y0 * terk[j];
Q			q+=2;
V			s+=2;	t+=2;
			l++;
			for (int j=0; j<NWAY; j++) {
				y1[j]  = (vdup(al[1])*y0[j]) *cost[j];		//	y1[j] = vdup(al[1])*cost[j]*y0[j];
V				dy1[j] = (vdup(al[1])*y0[j]) *(cost[j]*cost[j] - st2[j]);		//	dy1[j] = vdup(al[1])*(cost[j]*dy0[j] - y0[j]*st2[j]);
			}
			al+=2;
			while (l<llim) {	// compute even and odd parts
Q				for (int j=0; j<NWAY; j++)	q[0] += Y1 * rork[j];		// real odd
Q				for (int j=0; j<NWAY; j++)	q[1] += Y1 * roik[j];		// imag odd
V				for (int j=0; j<NWAY; j++)	s[0] += DY1 * terk[j]  + Y1 * poik[j];
V				for (int j=0; j<NWAY; j++)	s[1] += DY1 * teik[j]  - Y1 * pork[j];
V				for (int j=0; j<NWAY; j++)	t[0] -= DY1 * perk[j]  - Y1 * toik[j];
V				for (int j=0; j<NWAY; j++)	t[1] -= DY1 * peik[j]  + Y1 * tork[j];
				for (int j=0; j<NWAY; j++) {
					y0[j] = vdup(al[1])*cost[j]*y1[j] + vdup(al[0])*y0[j];
V					dy0[j] = vdup(al[1])*(cost[j]*dy1[j] - y1[j]*st2[j]) + vdup(al[0])*dy0[j];
				}
Q				for (int j=0; j<NWAY; j++)	q[2] += Y0 * rerk[j];		// real even
Q				for (int j=0; j<NWAY; j++)	q[3] += Y0 * reik[j];		// imag even
V				for (int j=0; j<NWAY; j++)	s[2] += DY0 * tork[j]  + Y0 * peik[j];
V				for (int j=0; j<NWAY; j++)	s[3] += DY0 * toik[j]  - Y0 * perk[j];
V				for (int j=0; j<NWAY; j++)	t[2] -= DY0 * pork[j]  - Y0 * teik[j];
V				for (int j=0; j<NWAY; j++)	t[3] -= DY0 * poik[j]  + Y0 * terk[j];
Q				q+=4;
V				s+=4;	t+=4;
				l+=2;
				for (int j=0; j<NWAY; j++) {
					y1[j] = vdup(al[3])*cost[j]*y0[j] + vdup(al[2])*y1[j];
V					dy1[j] = vdup(al[3])*(cost[j]*dy0[j] - y0[j]*st2[j]) + vdup(al[2])*dy1[j];
				}
				al+=4;
			}
			if (l==llim) {
Q				for (int j=0; j<NWAY; j++)	q[0] += Y1 * rork[j];		// real odd
Q				for (int j=0; j<NWAY; j++)	q[1] += Y1 * roik[j];		// imag odd
V				for (int j=0; j<NWAY; j++)	s[0] += DY1 * terk[j]  + Y1 * poik[j];
V				for (int j=0; j<NWAY; j++)	s[1] += DY1 * teik[j]  - Y1 * pork[j];
V				for (int j=0; j<NWAY; j++)	t[0] -= DY1 * perk[j]  - Y1 * toik[j];
V				for (int j=0; j<NWAY; j++)	t[1] -= DY1 * peik[j]  + Y1 * tork[j];
			}
		#if _GCC_VEC_
			k += 2*NWAY;
		#else
			k+=NWAY;
		#endif
		} while (k < ni);
		#undef Y0
		#undef Y1
V		#undef DY0
V		#undef DY1
	  } else {		// llim > SHT_L_RESCALE_FLY
		#define Y0 (y0[j]*scale[j])
		#define Y1 (y1[j]*scale[j])
V		#define DY0 (dy0[j]*scale[j])
V		#define DY1 (dy1[j]*scale[j])
		do {
		#if _GCC_VEC_
Q			s2d* q = qq;
V			s2d* s = ss;		s2d* t = tt;
		#else
			l = LiM(shtns, m, im);
Q			double* q = (double *) &Qlm[l];
V			double* s = (double *) &Slm[l];
V			double* t = (double *) &Tlm[l];
		#endif
			al = shtns->blm[im];
			s2d cost[NWAY], y0[NWAY], y1[NWAY], scale[NWAY];
V			s2d st2[NWAY], dy0[NWAY], dy1[NWAY];
Q			s2d rerk[NWAY], reik[NWAY], rork[NWAY], roik[NWAY];		// help the compiler to cache into registers.
V			s2d terk[NWAY], teik[NWAY], tork[NWAY], toik[NWAY];
V			s2d perk[NWAY], peik[NWAY], pork[NWAY], poik[NWAY];
			for (int j=0; j<NWAY; j++) {
				cost[j] = ((s2d*)(st+k))[j];
				y0[j] = vdup(al[0]) * ((s2d*)(wg+k))[j];		// weight appears here.
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
				cost[j] = ((s2d*)(ct+k))[j];
V				dy0[j] = cost[j]*y0[j];
			}
			l=m;
			for (int j=0; j<NWAY; j++) {
Q				rerk[j] = ((s2d*)(rer+k))[j];		reik[j] = ((s2d*)(rei+k))[j];		rork[j] = ((s2d*)(ror+k))[j];		roik[j] = ((s2d*)(roi+k))[j];
Q				q[0] += Y0 * rerk[j];		// real even
Q				q[1] += Y0 * reik[j];		// imag even
V				terk[j] = ((s2d*)(ter+k))[j];		teik[j] = ((s2d*)(tei+k))[j];		tork[j] = ((s2d*)(tor+k))[j];		toik[j] = ((s2d*)(toi+k))[j];
V				perk[j] = ((s2d*)(per+k))[j];		peik[j] = ((s2d*)(pei+k))[j];		pork[j] = ((s2d*)(por+k))[j];		poik[j] = ((s2d*)(poi+k))[j];
V				s[0] += DY0 * tork[j]  + Y0 * peik[j];
V				s[1] += DY0 * toik[j]  - Y0 * perk[j];
V				t[0] -= DY0 * pork[j]  - Y0 * teik[j];
V				t[1] -= DY0 * poik[j]  + Y0 * terk[j];
			}
Q			q+=2;
V			s+=2;	t+=2;
			l++;
			for (int j=0; j<NWAY; j++) {
				y1[j]  = (vdup(al[1])*y0[j]) *cost[j];		//	y1[j] = vdup(al[1])*cost[j]*y0[j];
V				dy1[j] = (vdup(al[1])*y0[j]) *(cost[j]*cost[j] - st2[j]);		//	dy1[j] = vdup(al[1])*(cost[j]*dy0[j] - y0[j]*st2[j]);
			}
			al+=2;
			while (l<llim) {	// compute even and odd parts
				for (int j=0; j<NWAY; j++) {
Q					q[0] += Y1 * rork[j];		// real odd
Q					q[1] += Y1 * roik[j];		// imag odd
V					s[0] += DY1 * terk[j]  + Y1 * poik[j];
V					s[1] += DY1 * teik[j]  - Y1 * pork[j];
V					t[0] -= DY1 * perk[j]  - Y1 * toik[j];
V					t[1] -= DY1 * peik[j]  + Y1 * tork[j];
				}
				for (int j=0; j<NWAY; j++) {
					y0[j] = vdup(al[1])*cost[j]*y1[j] + vdup(al[0])*y0[j];
V					dy0[j] = vdup(al[1])*(cost[j]*dy1[j] - y1[j]*st2[j]) + vdup(al[0])*dy0[j];
				}
				for (int j=0; j<NWAY; j++) {
Q					q[2] += Y0 * rerk[j];		// real even
Q					q[3] += Y0 * reik[j];		// imag even
V					s[2] += DY0 * tork[j]  + Y0 * peik[j];
V					s[3] += DY0 * toik[j]  - Y0 * perk[j];
V					t[2] -= DY0 * pork[j]  - Y0 * teik[j];
V					t[3] -= DY0 * poik[j]  + Y0 * terk[j];
				}
Q				q+=4;
V				s+=4;	t+=4;
				l+=2;
				for (int j=0; j<NWAY; j++) {
					y1[j] = vdup(al[3])*cost[j]*y0[j] + vdup(al[2])*y1[j];
V					dy1[j] = vdup(al[3])*(cost[j]*dy0[j] - y0[j]*st2[j]) + vdup(al[2])*dy1[j];
				}
				al+=4;
			}
			if (l==llim) {
				for (int j=0; j<NWAY; j++) {
Q					q[0] += Y1 * rork[j];		// real odd
Q					q[1] += Y1 * roik[j];		// imag odd
V					s[0] += DY1 * terk[j]  + Y1 * poik[j];
V					s[1] += DY1 * teik[j]  - Y1 * pork[j];
V					t[0] -= DY1 * perk[j]  - Y1 * toik[j];
V					t[1] -= DY1 * peik[j]  + Y1 * tork[j];
				}
			}
		#if _GCC_VEC_
			k += 2*NWAY;
		#else
			k+=NWAY;
		#endif
		} while (k < ni);
		#undef Y0
		#undef Y1
V		#undef DY0
V		#undef DY1
	  }
		l = LiM(shtns, m, im);
Q		complex double *Ql = &Qlm[l];
V		complex double *Sl = &Slm[l];
V		complex double *Tl = &Tlm[l];
3		double m_1 = 1.0/m;
		#if _GCC_VEC_
			for (l=0; l<=llim-m; l++) {
Q				Ql[l] = vlo_to_dbl(qq[2*l]) + vhi_to_dbl(qq[2*l]) + I*(vlo_to_dbl(qq[2*l+1]) + vhi_to_dbl(qq[2*l+1]));
3				Ql[l] *= m_1;
V				Sl[l] = (vlo_to_dbl(ss[2*l]) + vhi_to_dbl(ss[2*l]) + I*(vlo_to_dbl(ss[2*l+1]) + vhi_to_dbl(ss[2*l+1])))*l_2[l+m];
V				Tl[l] = (vlo_to_dbl(tt[2*l]) + vhi_to_dbl(tt[2*l]) + I*(vlo_to_dbl(tt[2*l+1]) + vhi_to_dbl(tt[2*l+1])))*l_2[l+m];
			}
		#else
V			for (l=0; l<=llim-m; l++) {
3				Ql[l] *= m_1;
V				Sl[l] *= l_2[l+m];
V				Tl[l] *= l_2[l+m];
V			}
		#endif
		#ifdef SHT_VAR_LTR
			for (l=llim+1-m; l<=LMAX-m; l++) {
Q				Ql[l] = 0.0;
V				Sl[l] = 0.0;		Tl[l] = 0.0;
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

  	if (SHT_FFT > 1) {		// free memory
Q	    fftw_free(BrF - NLAT*(imlim+1));
VX	    fftw_free(BtF - NLAT*(imlim+1));	// this frees also BpF.
	}
  #endif

  }
