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

#Q void spat_to_SH(double *Vr, complex double *Qlm)
#V void spat_to_SHsphtor(double *Vt, double *Vp, complex double *Slm, complex double *Tlm)
# {
Q	complex double *BrF;		// contains the Fourier transformed data
V	complex double *BtF, *BpF;	// contains the Fourier transformed data
	double *al;
	long int llim, m;
	long int ni, k;
	long int i,i0, im,l;
  #if _GCC_VEC_
Q	s2d qq[2*LMAX];
V	s2d ss[2*LMAX];
V	s2d tt[2*LMAX];
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

	llim = LTR;
	ni = NLAT_2;	// copy NLAT_2 to a local variable for faster access (inner loop limit)
Q    	BrF = (complex double *) Vr;
V    	BtF = (complex double *) Vt;	BpF = (complex double *) Vp;

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
Q	    fftw_execute_dft_r2c(fft,Vr, BrF);
V	    fftw_execute_dft_r2c(fft,Vt, BtF);
V	    fftw_execute_dft_r2c(fft,Vp, BpF);
	}
  #endif

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
		} while(i<ni+2*(NWAY-1));
	#else
		} while(i<ni+NWAY-1);
	#endif
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
Q		Qlm[0] = r0 * al0[0];					// l=0 is done.
V		Slm[0] = 0.0;		Tlm[0] = 0.0;		// l=0 is zero for the vector transform.
		k = 0;
		for (l=1;l<=llim;l++) {
Q			qq[l] = vdup(0.0);
V			ss[l] = vdup(0.0);		tt[l] = vdup(0.0);
		}
		do {
			al = al0;
			s2d cost[NWAY], y0[NWAY], y1[NWAY];
V			s2d sint[NWAY], dy0[NWAY], dy1[NWAY];
			for (int j=0; j<NWAY; j++) {
				cost[j] = ((s2d*)(ct+k))[j];
V				sint[j] = ((s2d*)(st+k))[j];
				y0[j] = vdup(al[0]);
V				dy0[j] = vdup(0.0);
				y1[j] = vdup(al[0]*al[1]) * cost[j];
V				dy1[j] = -vdup(al[0]*al[1]) * sint[j];
			}
			al+=2;	l=1;
			while(l<llim) {
				for (int j=0; j<NWAY; j++) {
					y0[j]  = vdup(al[1])*cost[j]*y1[j] + vdup(al[0])*y0[j];
V					dy0[j] = vdup(al[1])*(cost[j]*dy1[j] - y1[j]*sint[j]) + vdup(al[0])*dy0[j];
				}
				for (int j=0; j<NWAY; j++) {
Q					qq[l] += y1[j] * ((s2d*)(ror+k))[j];
Q					qq[l+1]   += y0[j] * ((s2d*)(rer+k))[j];
V					ss[l] += dy1[j] * ((s2d*)(ter+k))[j];
V					tt[l] -= dy1[j] * ((s2d*)(per+k))[j];
V					ss[l+1] += dy0[j] * ((s2d*)(tor+k))[j];
V					tt[l+1] -= dy0[j] * ((s2d*)(por+k))[j];
				}
				for (int j=0; j<NWAY; j++) {
					y1[j]  = vdup(al[3])*cost[j]*y0[j] + vdup(al[2])*y1[j];
V					dy1[j] = vdup(al[3])*(cost[j]*dy0[j] - y0[j]*sint[j]) + vdup(al[2])*dy1[j];
				}
				al+=4;	l+=2;
			}
			if (l==llim) {
				for (int j=0; j<NWAY; j++) {
Q					qq[l] += y1[j] * ((s2d*)(ror+k))[j];
V					ss[l] += dy1[j] * ((s2d*)(ter+k))[j];
V					tt[l] -= dy1[j] * ((s2d*)(per+k))[j];
				}
			}
		#if _GCC_VEC_
			k+=2*NWAY;
		#else
			k+=NWAY;
		#endif
		} while (k < NLAT_2);
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
	for (im=1;im<=MTR;im++) {
		i0 = tm[im];
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
		} while (i < ni+2*(NWAY-1));
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
Q			double* q = (double *) &Qlm[LiM(m,im)];
V			double* s = (double *) &Slm[LiM(m,im)];
V			double* t = (double *) &Tlm[LiM(m,im)];
		#endif
		for (l=llim-m; l>=0; l--) {
Q			q[0] = vdup(0.0);		q[1] = vdup(0.0);		q+=2;
V			s[0] = vdup(0.0);		s[1] = vdup(0.0);		s+=2;
V			t[0] = vdup(0.0);		t[1] = vdup(0.0);		t+=2;
		}
		do {
		#if _GCC_VEC_
Q			s2d* q = qq;
V			s2d* s = ss;		s2d* t = tt;
		#else
Q			double* q = (double *) &Qlm[LiM(m,im)];
V			double* s = (double *) &Slm[LiM(m,im)];
V			double* t = (double *) &Tlm[LiM(m,im)];
		#endif
			al = alm[im];
			s2d cost[NWAY], y0[NWAY], y1[NWAY];
V			s2d st2[NWAY], dy0[NWAY], dy1[NWAY];
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
			l=m;
			for (int j=0; j<NWAY; j++) {
Q				q[0] += y0[j] * ((s2d*)(rer+k))[j];		// real even
Q				q[1] += y0[j] * ((s2d*)(rei+k))[j];		// imag even
V				s[0] += dy0[j] * ((s2d*)(tor+k))[j]  + y0[j] * ((s2d*)(pei+k))[j];
V				s[1] += dy0[j] * ((s2d*)(toi+k))[j]  - y0[j] * ((s2d*)(per+k))[j];
V				t[0] -= dy0[j] * ((s2d*)(por+k))[j]  - y0[j] * ((s2d*)(tei+k))[j];
V				t[1] -= dy0[j] * ((s2d*)(poi+k))[j]  + y0[j] * ((s2d*)(ter+k))[j];
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
Q					q[0] += y1[j] * ((s2d*)(ror+k))[j];		// real odd
Q					q[1] += y1[j] * ((s2d*)(roi+k))[j];		// imag odd
V					s[0] += dy1[j] * ((s2d*)(ter+k))[j]  + y1[j] * ((s2d*)(poi+k))[j];
V					s[1] += dy1[j] * ((s2d*)(tei+k))[j]  - y1[j] * ((s2d*)(por+k))[j];
V					t[0] -= dy1[j] * ((s2d*)(per+k))[j]  - y1[j] * ((s2d*)(toi+k))[j];
V					t[1] -= dy1[j] * ((s2d*)(pei+k))[j]  + y1[j] * ((s2d*)(tor+k))[j];
				}
				for (int j=0; j<NWAY; j++) {
					y0[j] = vdup(al[1])*cost[j]*y1[j] + vdup(al[0])*y0[j];
V					dy0[j] = vdup(al[1])*(cost[j]*dy1[j] - y1[j]*st2[j]) + vdup(al[0])*dy0[j];
				}
				for (int j=0; j<NWAY; j++) {
Q					q[2] += y0[j] * ((s2d*)(rer+k))[j];		// real even
Q					q[3] += y0[j] * ((s2d*)(rei+k))[j];		// imag even
V					s[2] += dy0[j] * ((s2d*)(tor+k))[j]  + y0[j] * ((s2d*)(pei+k))[j];
V					s[3] += dy0[j] * ((s2d*)(toi+k))[j]  - y0[j] * ((s2d*)(per+k))[j];
V					t[2] -= dy0[j] * ((s2d*)(por+k))[j]  - y0[j] * ((s2d*)(tei+k))[j];
V					t[3] -= dy0[j] * ((s2d*)(poi+k))[j]  + y0[j] * ((s2d*)(ter+k))[j];
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
Q					q[0] += y1[j] * ((s2d*)(ror+k))[j];		// real odd
Q					q[1] += y1[j] * ((s2d*)(roi+k))[j];		// imag odd
V					s[0] += dy1[j] * ((s2d*)(ter+k))[j]  + y1[j] * ((s2d*)(poi+k))[j];
V					s[1] += dy1[j] * ((s2d*)(tei+k))[j]  - y1[j] * ((s2d*)(por+k))[j];
V					t[0] -= dy1[j] * ((s2d*)(per+k))[j]  - y1[j] * ((s2d*)(toi+k))[j];
V					t[1] -= dy1[j] * ((s2d*)(pei+k))[j]  + y1[j] * ((s2d*)(tor+k))[j];
				}
			}
		#if _GCC_VEC_
			k += 2*NWAY;
		#else
			k+=NWAY;
		#endif
		} while (k < NLAT_2);
Q		complex double *Ql = &Qlm[LiM(m,im)];
V		complex double *Sl = &Slm[LiM(m,im)];
V		complex double *Tl = &Tlm[LiM(m,im)];
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

  	if (SHT_FFT > 1) {		// free memory
Q	    fftw_free(BrF - NLAT*(MTR+1));
VX	    fftw_free(BtF - NLAT*(MTR+1));	// this frees also BpF.
	}
  #endif

# }
