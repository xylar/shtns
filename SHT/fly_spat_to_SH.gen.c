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
Q	complex double q0,q1;
V	complex double s0,t0,s1,t1;
  #if _GCC_VEC_
	s2d qt[2*LMAX];
  #endif
  
  #ifndef SHT_AXISYM
Q	double rei[NLAT_2+2*NWAY] SSE;
Q	double rer[NLAT_2+2*NWAY] SSE;
Q	double ror[NLAT_2+2*NWAY] SSE;
Q	double roi[NLAT_2+2*NWAY] SSE;
  #else
Q	double rer[NLAT_2+2*NWAY] SSE;
Q	double ror[NLAT_2+2*NWAY] SSE;
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
QE		double r0 = 0.0;
  #ifndef SHT_AXISYM
		i0 = (NPHI==1) ? 1 : 2;		// stride of source data.
 B		do {	// compute symmetric and antisymmetric parts.
 B			double w = wg[i];
QB			double a = ((double*)BrF)[i*i0];		double b = ((double*)BrF)[(NLAT-1)*i0 -i*i0];
QB			ror[i] = (a-b)*w;		rer[i] = (a+b)*w;
QB			r0 += ((a+b)*w);
 B			i++;
 B		} while(i<ni);
  #else
 B		do {	// compute symmetric and antisymmetric parts.
 B			double np = wg[i]*NPHI;
QB			double a = ((double*)BrF)[i];		double b = ((double*)BrF)[NLAT-1-i];
QB			ror[i] = (a-b)*np;		rer[i] = (a+b)*np;
QB			r0 += ((a+b)*np);
 B			i++;
 B		} while(i<ni);
  #endif
		do {
			rer[i] = 0.0;		ror[i] = 0.0;
			i++;
	#if _GCC_VEC_
		} while(i<ni+2*(NWAY-1));
	#else
		} while(i<ni+NWAY-1);
	#endif
Q		BrF += NLAT;
V		BtF += NLAT;	BpF += NLAT;
Q		Qlm[0] = r0 * al0[0];		// l=0 is done.
		k = 0;
		for (l=1;l<=llim;l++) qt[l] = vdup(0.0);
		do {
			al = al0;
			s2d cost[NWAY], y0[NWAY], y1[NWAY];
			for (int j=0; j<NWAY; j++) {
				cost[j] = ((s2d*)(ct+k))[j];
				y0[j] = vdup(al[0]);
				y1[j] = vdup(al[0]*al[1]) * cost[j];
			}
			al+=2;	l=1;
			while(l<llim) {
				for (int j=0; j<NWAY; j++) {
					y0[j]  = vdup(al[1])*cost[j]*y1[j] + vdup(al[0])*y0[j];
V					dy0[j] = vdup(al[1])*(cost[j]*dy1[j] - y1[j]*sint[j]) + vdup(al[0])*dy0[j];
				}
				for (int j=0; j<NWAY; j++) {
Q					qt[l] += y1[j] * ((s2d*)(ror+k))[j];
Q					qt[l+1]   += y0[j] * ((s2d*)(rer+k))[j];
				}
				for (int j=0; j<NWAY; j++) {
					y1[j]  = vdup(al[3])*cost[j]*y0[j] + vdup(al[2])*y1[j];
V					dy1[j] = vdup(al[3])*(cost[j]*dy0[j] - y0[j]*sint[j]) + vdup(al[2])*dy1[j];
				}
				al+=4;	l+=2;
			}
			if (l==llim) {
				for (int j=0; j<NWAY; j++) {
Q					qt[l] += y1[j] * ((s2d*)(ror+k))[j];
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
				Qlm[l] = vlo_to_dbl(qt[l]) + vhi_to_dbl(qt[l]);
			#else
				Qlm[l] = qt[l];
			#endif
		}
		#ifdef SHT_VAR_LTR
			for (l=llim+1; l<= LMAX; l++) {
Q				Qlm[l] = 0.0;
			}
		#endif

  #ifndef SHT_AXISYM
	for (im=1;im<=MTR;im++) {
		i0 = tm[im];
		m = im*MRES;
	#if _GCC_VEC_
		i0=(i0>>1)*2;		// stay on a 16 byte boundary
	#endif
 B		i=i0;
 B		do {	// compute symmetric and antisymmetric parts, and reorganize data.
QB			v2d q0 = ((v2d *)BrF)[i];	v2d q1 = ((v2d *)BrF)[NLAT-1-i];
			#if _GCC_VEC_
				rer[i] = vlo_to_dbl(q0+q1);		rei[i] = vhi_to_dbl(q0+q1);
				ror[i] = vlo_to_dbl(q0-q1);		roi[i] = vhi_to_dbl(q0-q1);
			#else
				rer[i] = creal(q0+q1);		rei[i] = cimag(q0+q1);
				ror[i] = creal(q0-q1);		roi[i] = cimag(q0-q1);
			#endif
 B			i++;
 B		} while (i<ni);
		do {
			rer[i] = 0.0;		rei[i] = 0.0;		ror[i] = 0.0;		roi[i] = 0.0;
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
			s2d* q = qt;
		#else
			double* q = (double *) &Qlm[LiM(m,im)];
		#endif
		for (l=llim-m; l>=0; l--) {
			q[0] = vdup(0.0);		q[1] = vdup(0.0);
			q+=2;
		}
		do {
		#if _GCC_VEC_
			s2d* q = qt;
		#else
			double* q = (double *) &Qlm[LiM(m,im)];
		#endif
			al = alm[im];
			s2d cost[NWAY], y0[NWAY], y1[NWAY];
			for (int j=0; j<NWAY; j++) {
				cost[j] = ((s2d*)(st+k))[j];
				y0[j] = vdup(al[0]) * ((s2d*)(wg+k))[j];		// weight appears here.
			}
			l=m;
			do {		// sin(theta)^m
				if (l&1) for (int j=0; j<NWAY; j++) y0[j] *= cost[j];
				for (int j=0; j<NWAY; j++) cost[j] *= cost[j];
			} while(l >>= 1);
			for (int j=0; j<NWAY; j++) {
				cost[j] = ((s2d*)(ct+k))[j];
			}
			l=m;
			for (int j=0; j<NWAY; j++) {
				q[0] += y0[j] * ((s2d*)(rer+k))[j];		// real even
				q[1] += y0[j] * ((s2d*)(rei+k))[j];		// imag even
			}
			q+=2;	l++;
			for (int j=0; j<NWAY; j++) {
				y1[j]  = (vdup(al[1])*y0[j]) *cost[j];		//	y1[j] = vdup(al[1])*cost[j]*y0[j];
			}
			al+=2;
			while (l<llim) {	// compute even and odd parts
				for (int j=0; j<NWAY; j++) {
					q[0] += y1[j] * ((s2d*)(ror+k))[j];		// real odd
					q[1] += y1[j] * ((s2d*)(roi+k))[j];		// imag odd
				}
				for (int j=0; j<NWAY; j++) {
					y0[j] = vdup(al[1])*cost[j]*y1[j] + vdup(al[0])*y0[j];
				}
				for (int j=0; j<NWAY; j++) {
					q[2] += y0[j] * ((s2d*)(rer+k))[j];		// real even
					q[3] += y0[j] * ((s2d*)(rei+k))[j];		// imag even
				}
				q+=4;	l+=2;
				for (int j=0; j<NWAY; j++) {
					y1[j] = vdup(al[3])*cost[j]*y0[j] + vdup(al[2])*y1[j];
				}
				al+=4;
			}
			if (l==llim) {
				for (int j=0; j<NWAY; j++) {
					q[0] += y1[j] * ((s2d*)(ror+k))[j];		// real odd
					q[1] += y1[j] * ((s2d*)(roi+k))[j];		// imag odd
				}
			}
		#if _GCC_VEC_
			k += 2*NWAY;
		#else
			k+=NWAY;
		#endif
		} while (k < NLAT_2);
		complex double *Ql = &Qlm[LiM(m,im)];
		#if _GCC_VEC_
			for (l=0; l<=llim-m; l++) {
				Ql[l] = vlo_to_dbl(qt[2*l]) + vhi_to_dbl(qt[2*l]) + I*(vlo_to_dbl(qt[2*l+1]) + vhi_to_dbl(qt[2*l+1]));
			}
		#endif
		#ifdef SHT_VAR_LTR
			for (l=llim+1-m; l<=LMAX-m; l++) {
				Ql[l] = 0.0;
			}
		#endif
	}

  	if (SHT_FFT > 1) {		// free memory
Q	    fftw_free(BrF - NLAT*(MTR+1));
VX	    fftw_free(BtF - NLAT*(MTR+1));	// this frees also BpF.
	}
  #endif

# }
