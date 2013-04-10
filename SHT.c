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

/********************************************************************
 * SHTns : Spherical Harmonic Transform for numerical simulations.  *
 *    written by Nathanael Schaeffer / CNRS                         *
 ********************************************************************/

/** \internal \file SHT.c
 * \brief main source file for SHTns.
 * This files contains initialization code and also some partial transforms (point or latitudinal evaluations)
 */

#include <stdio.h>
#include <string.h>
// global variables definitions
#include "sht_private.h"

// cycle counter from FFTW
#include "cycle.h"

// chained list of sht_setup : start with NULL
shtns_cfg sht_data = NULL;
#ifdef _OPENMP
  int omp_threads = 1;	// multi-thread disabled by default.
  #if HAVE_LIBFFTW3_OMP
	#define OMP_FFTW
  #endif
#else
  #define omp_threads 1
#endif

/// \internal Abort program with error message.
static void shtns_runerr(const char * error_text)
{
	printf("*** [" PACKAGE_NAME "] Run-time error : %s\n",error_text);
	exit(1);
}

/* PUBLIC useful functions */

/// returns the l=0, m=0 SH coefficient corresponding to a uniform value of 1.
double sh00_1(shtns_cfg shtns) {
	return shtns->Y00_1;
}
/// returns the l=1, m=0 SH coefficient corresponding to cos(theta).
double sh10_ct(shtns_cfg shtns) {
	return shtns->Y10_ct;
}
/// returns the l=1, m=1 SH coefficient corresponding to sin(theta).cos(phi).
double sh11_st(shtns_cfg shtns) {
	return shtns->Y11_st;
}
/// returns the l,m SH coefficient corresponding to unit energy.
double shlm_e1(shtns_cfg shtns, int l, int m) {
	double x = shtns->Y00_1/sqrt(4.*M_PI);
	if (SHT_NORM == sht_schmidt) x *= sqrt(2*l+1);
	if ((m!=0)&&((shtns->norm & SHT_REAL_NORM)==0)) x *= sqrt(0.5);
	return(x);
}

/// \code return (mmax+1)*(lmax+1) - mres*(mmax*(mmax+1))/2; \endcode */
/// \ingroup init
long nlm_calc(long lmax, long mmax, long mres)
{
	if (mmax*mres > lmax) mmax = lmax/mres;
	return (mmax+1)*(lmax+1) - ((mmax*mres)*(mmax+1))/2;	// this is wrong if lmax < mmax*mres
}


/*  LEGENDRE FUNCTIONS  */
#include "sht_legendre.c"


/*
	SHT FUNCTIONS
*/

/** \addtogroup rotation Rotation of SH fields.
Rotation around axis other than Z should be considered of beta quality (they have been tested but may still contain bugs).
They also require \c mmax = \c lmax. They use an Algorithm inspired by the pseudospectral rotation described in
Gimbutas Z. and Greengard L. 2009 "A fast and stable method for rotating spherical harmonic expansions" <i>Journal of Computational Physics</i>.
doi:<a href="http://dx.doi.org/10.1016/j.jcp.2009.05.014">10.1016/j.jcp.2009.05.014</a>

These functions do only require a call to \ref shtns_create, but not to \ref shtns_set_grid.
*/
//@{

/// Rotate a SH representation Qlm around the z-axis by angle alpha (in radians),
/// which is the same as rotating the reference frame by angle -alpha.
/// Result is stored in Rlm (which can be the same array as Qlm).
void SH_Zrotate(shtns_cfg shtns, complex double *Qlm, double alpha, complex double *Rlm)
{
	int im, l, lmax, mmax, mres;

	lmax = shtns->lmax;		mmax = shtns->mmax;		mres = shtns->mres;

	if (Rlm != Qlm) {		// copy m=0 which does not change.
		l=0;	do { Rlm[l] = Qlm[l]; } while(++l <= lmax);
	}
	if (mmax > 0) {
		im=1; do {
			complex double eima = cos(im*mres*alpha) - I*sin(im*mres*alpha);		// rotate reference frame by angle -alpha
			for (l=im*mres; l<=lmax; ++l)	Rlm[LiM(shtns, l, im)] = Qlm[LiM(shtns, l, im)] * eima;
		} while(++im <= mmax);
	}
}

//@}

/** \internal rotation kernel used by SH_Yrotate90(), SH_Xrotate90() and SH_rotate().
 Algorithm based on the pseudospectral rotation[1] :
 - rotate around Z by angle dphi0.
 - synthetize for each l the spatial description for phi=0 and phi=pi on an equispaced latitudinal grid.
 - Fourier ananlyze as data on the equator to recover the m in the 90 degrees rotated frame.
 - rotate around new Z by angle dphi1.
 [1] Gimbutas Z. and Greengard L. 2009 "A fast and stable method for rotating spherical harmonic expansions" Journal of Computational Physics. **/
static void SH_rotK90(shtns_cfg shtns, complex double *Qlm, complex double *Rlm, double dphi0, double dphi1)
{
	fftw_plan fft;
	complex double *q;
	double *q0;
	long int k, m, l;
	int lmax, ntheta;
	int nrembed, ncembed;

	lmax = shtns->lmax;
	ntheta = ((lmax+2)>>1)*2;
	m = 2* sizeof(double)*(2*ntheta+2)*lmax;
	q0 = fftw_malloc(m);		memset(q0, 0, m);		// alloc & zero out.

	// rotate around Z by dphi0
	if (dphi0 != 0.0) {
		SH_Zrotate(shtns, Qlm, dphi0, Rlm);
		Qlm = Rlm;
	} else {
		Rlm[0] = Qlm[0];		// l=0 is rotation invariant.
	}

  #pragma omp parallel private(k,m,l) num_threads(shtns->nthreads)
  {
	double yl[lmax+1];
	// compute q(l) on the meridian phi=0 and phi=pi. (rotate around X)
	#pragma omp for schedule(static)
	for (k=0; k<ntheta/2; ++k) {
		double cost= cos(((0.5*M_PI)*(2*k+1))/ntheta);
		double sint_1 = 1.0/sqrt((1.0-cost)*(1.0+cost));
		m=0;
			legendre_sphPlm_array(shtns, lmax, m, cost, yl+m);
			double sgnt = -1.0;
			for (l=1; l<=lmax; ++l) {
				double qr = creal(Qlm[LiM(shtns, l, m)]) * yl[l];
				q0[k*2*lmax +2*(l-1)] = qr;
				q0[(ntheta-1-k)*2*lmax +2*(l-1)] = sgnt*qr;
				q0[(ntheta+k)*2*lmax +2*(l-1)] = sgnt*qr;
				q0[(2*ntheta-1-k)*2*lmax +2*(l-1)] = qr;
				sgnt *= -1.0;
			}
		#ifndef _GCC_VEC_
		double sgnm = 1.0;
		for (m=1; m<=lmax; ++m) {
			legendre_sphPlm_array(shtns, lmax, m, cost, yl+m);
			double sgnt = 1.0;
			sgnm *= -1.0;
			for (l=m; l<=lmax; ++l) {
				double qr = creal(Qlm[LiM(shtns, l, m)]) * yl[l];
				double qi = cimag(Qlm[LiM(shtns, l, m)]) * m*yl[l]*sint_1;
				qr += qr;	qi += qi;		// x2 for m>0
				q0[k*2*lmax +2*(l-1)] += qr;						// q0
				q0[k*2*lmax +2*(l-1)+1] -= qi;						// dq0
				q0[(ntheta-1-k)*2*lmax +2*(l-1)] += sgnt*qr;
				q0[(ntheta-1-k)*2*lmax +2*(l-1)+1] -= sgnt*qi;
				q0[(ntheta+k)*2*lmax +2*(l-1)] += (sgnm*sgnt)*qr;
				q0[(ntheta+k)*2*lmax +2*(l-1)+1] += (sgnm*sgnt)*qi;
				q0[(2*ntheta-1-k)*2*lmax +2*(l-1)] += sgnm*qr;
				q0[(2*ntheta-1-k)*2*lmax +2*(l-1)+1] += sgnm*qi;
				sgnt *= -1.0;
			}
		}
		#else
		s2d sgnm = SIGN_MASK_HI;
		s2d sgnflip = SIGN_MASK_2;
		for (m=1; m<=lmax; ++m) {
			legendre_sphPlm_array(shtns, lmax, m, cost, yl+m);
			s2d sgnt = vdup(0.0);
			s2d m_st = vset(2.0, -2*m*sint_1);		// x2 for m>0
			sgnm = _mm_xor_pd(sgnm, sgnflip);	// (-1)^m
			for (l=m; l<=lmax; ++l) {
				v2d qc = ((v2d*)Qlm)[LiM(shtns, l, m)] * vdup(yl[l]) * m_st;	// (q0, dq0)
				((v2d*)q0)[k*lmax +(l-1)] += qc;
				((v2d*)q0)[(ntheta-1-k)*lmax +(l-1)] += _mm_xor_pd(sgnt, qc);
				qc = _mm_xor_pd(sgnm, qc);
				((v2d*)q0)[(ntheta+k)*lmax +(l-1)] += _mm_xor_pd( sgnt, qc );
				((v2d*)q0)[(2*ntheta-1-k)*lmax +(l-1)] += qc;
				sgnt = _mm_xor_pd(sgnt, sgnflip);	// (-1)^(l+m)
			}
		}
		#endif
	}
  }

	// perform FFT
	#ifdef OMP_FFTW
		k = (lmax < 63) ? 1 : shtns->nthreads;
		fftw_plan_with_nthreads(k);
	#endif
	q = (complex double*) q0;
	ntheta*=2;		nrembed = ntheta+2;		ncembed = nrembed/2;
	fft = fftw_plan_many_dft_r2c(1, &ntheta, 2*lmax, q0, &nrembed, 2*lmax, 1, q, &ncembed, 2*lmax, 1, FFTW_ESTIMATE);
	fftw_execute_dft_r2c(fft, q0, q);
	fftw_destroy_plan(fft);

	double yl[lmax+1];		double dyl[lmax+1];
	m=0;
		//legendre_sphPlm_deriv_array(shtns, lmax, m, 0.0, 1.0, yl+m, dyl+m);
		legendre_sphPlm_deriv_array_equ(shtns, lmax, m, yl+m, dyl+m);
		for (l=1; l<lmax; l+=2) {
			Rlm[LiM(shtns, l,m)] =  -creal(q[m*2*lmax +2*(l-1)+1])/(dyl[l]*ntheta);
			Rlm[LiM(shtns, l+1,m)] =  creal(q[m*2*lmax +2*l])/(yl[l+1]*ntheta);
		}
		if (l==lmax) {
			Rlm[LiM(shtns, l,m)] =  -creal(q[m*2*lmax +2*(l-1)+1])/(dyl[l]*ntheta);
		}
	dphi1 += M_PI/ntheta;	// shift rotation angle by angle of first synthesis latitude.
	for (m=1; m<=lmax; ++m) {
		//legendre_sphPlm_deriv_array(shtns, lmax, m, 0.0, 1.0, yl+m, dyl+m);
		legendre_sphPlm_deriv_array_equ(shtns, lmax, m, yl+m, dyl+m);
		complex double eimdp = (cos(m*dphi1) - I*sin(m*dphi1))/(ntheta);
		for (l=m; l<lmax; l+=2) {
			Rlm[LiM(shtns, l,m)] =  eimdp*q[m*2*lmax +2*(l-1)]*(1./yl[l]);
			Rlm[LiM(shtns, l+1,m)] =  eimdp*q[m*2*lmax +2*l+1]*(-1./dyl[l+1]);
		}
		if (l==lmax) {
			Rlm[LiM(shtns, l,m)] =  eimdp*q[m*2*lmax +2*(l-1)]*(1./yl[l]);
		}
	}
	fftw_free(q0);
}


/// \addtogroup rotation
//@{

/// rotate Qlm by 90 degrees around X axis and store the result in Rlm.
/// shtns->mres MUST be 1, and lmax=mmax.
void SH_Xrotate90(shtns_cfg shtns, complex double *Qlm, complex double *Rlm)
{
	int lmax= shtns->lmax;
	if ((shtns->mres != 1) || (shtns->mmax < lmax)) shtns_runerr("truncature makes rotation not closed.");

	if (lmax == 1) {
		Rlm[0] = Qlm[0];	// l=0 is invariant.
		int l=1;													// rotation matrix for rotX(90), l=1 : m=[0, 1r, 1i]
			double q0 = creal(Qlm[LiM(shtns, l, 0)]);
			Rlm[LiM(shtns, l, 0)] = sqrt(2.0) * cimag(Qlm[LiM(shtns, l, 1)]);			//[m=0]     0        0    sqrt(2)
			Rlm[LiM(shtns, l ,1)] = creal(Qlm[LiM(shtns, l, 1)]) - I*(sqrt(0.5)*q0);	//[m=1r]    0        1      0
		return;																			//[m=1i] -sqrt(2)/2  0      0
	}

	SH_rotK90(shtns, Qlm, Rlm, 0.0,  -M_PI/2);
}

/// rotate Qlm by 90 degrees around Y axis and store the result in Rlm.
/// shtns->mres MUST be 1, and lmax=mmax.
void SH_Yrotate90(shtns_cfg shtns, complex double *Qlm, complex double *Rlm)
{
	int lmax= shtns->lmax;
	if ((shtns->mres != 1) || (shtns->mmax < lmax)) shtns_runerr("truncature makes rotation not closed.");

	if (lmax == 1) {
		Rlm[0] = Qlm[0];	// l=0 is invariant.
		int l=1;											// rotation matrix for rotY(90), l=1 : m=[0, 1r, 1i]
			double q0 = creal(Qlm[LiM(shtns, l, 0)]);									//[m=0]       0       0    sqrt(2)
			Rlm[LiM(shtns, l, 0)] = sqrt(2.0) * creal(Qlm[LiM(shtns, l, 1)]);			//[m=1r] -sqrt(2)/2   0      0
			Rlm[LiM(shtns, l ,1)] = I*cimag(Qlm[LiM(shtns, l, 1)]) - sqrt(0.5) * q0;	//[m=1i]      0       0      1
		return;
	}

	SH_rotK90(shtns, Qlm, Rlm, -M_PI/2, 0.0);
}

/// rotate Qlm around Y axis by arbitrary angle, using composition of rotations. Store the result in Rlm.
void SH_Yrotate(shtns_cfg shtns, complex double *Qlm, double alpha, complex double *Rlm)
{
	if ((shtns->mres != 1) || (shtns->mmax < shtns->lmax)) shtns_runerr("truncature makes rotation not closed.");

	SH_rotK90(shtns, Qlm, Rlm, 0.0, M_PI/2 + alpha);	// Zrotate(pi/2) + Yrotate90 + Zrotate(pi+alpha)
	SH_rotK90(shtns, Rlm, Rlm, 0.0, M_PI/2);			// Yrotate90 + Zrotate(pi/2)
}

//@}


/** \addtogroup operators Special operators
 * Apply special operators in spectral space: multiplication by cos(theta), sin(theta).d/dtheta.
*/
//@{

/// fill mx with the coefficients for multiplication by cos(theta)
/// \param mx : an array of 2*NLM double that will be filled with the matrix coefficients.
/// xq[lm] = mx[2*lm] * q[lm-1] + mx[2*lm+1] * q[lm+1];
void mul_ct_matrix(shtns_cfg shtns, double* mx)
{
	long int im,l,lm;
	double a_1;

	if (SHT_NORM == sht_schmidt) {
		lm=0;
		for (im=0; im<=MMAX; im++) {
			double* al = shtns->alm[im];
			long int m=im*MRES;
			mx[2*lm] = 0.0;
			a_1 = 1.0 / al[1];
			l=m;
			while(++l < LMAX) {
				al+=2;				
				mx[2*lm+2] = a_1;
				a_1 = 1.0 / al[1];
				mx[2*lm+1] = -a_1*al[0];        // = -al[2*(lm+1)] / al[2*(lm+1)+1];
				lm++;
			}
			if (l == LMAX) {	// the last one needs to be computed.
				mx[2*lm+2] = a_1;
				mx[2*lm+1] = sqrt((l+m)*(l-m))/(2*l+1);
				lm++;
			}
			mx[2*lm +1] = 0.0;
			lm++;
		}
	} else {
		lm=0;
		for (im=0; im<=MMAX; im++) {
			double* al = shtns->alm[im];
			l=im*MRES;
			mx[2*lm] = 0.0;
			while(++l <= LMAX) {
				a_1 = 1.0 / al[1];
				mx[2*lm+1] = a_1;		// specific to orthonormal.
				mx[2*lm+2] = a_1;
				lm++;	al+=2;
			}
			mx[2*lm +1] = 0.0;
			lm++;
		}
	}
}

/// fill mx with the coefficients of operator sin(theta).d/dtheta
/// \param mx : an array of 2*NLM double that will be filled with the matrix coefficients.
/// stdq[lm] = mx[2*lm] * q[lm-1] + mx[2*lm+1] * q[lm+1];
void st_dt_matrix(shtns_cfg shtns, double* mx)
{
	mul_ct_matrix(shtns, mx);
	for (int lm=0; lm<NLM; lm++) {
		mx[2*lm]   *=   shtns->li[lm] - 1;
		mx[2*lm+1] *= -(shtns->li[lm] + 2);
	}
}

/// Multiplication of Qlm by a matrix involving l+1 and l-1 only.
/// The result is stored in Rlm, which MUST be different from Qlm.
/// mx is an array of 2*NLM values as returned by \ref mul_ct_matrix or \ref st_dt_matrix
/// compute: Rlm[lm] = mx[2*lm] * Qlm[lm-1] + mx[2*lm+1] * Qlm[lm+1];
void SH_mul_mx(shtns_cfg shtns, double* mx, complex double *Qlm, complex double *Rlm)
{
	long int nlmlim, lm;
	v2d* vq = (v2d*) Qlm;
	v2d* vr = (v2d*) Rlm;
	nlmlim = NLM-1;
	lm = 0;
		s2d mxu = vdup(mx[1]);
		vr[0] = mxu*vq[1];
	for (lm=1; lm<nlmlim; lm++) {
		s2d mxl = vdup(mx[2*lm]);		s2d mxu = vdup(mx[2*lm+1]);
		vr[lm] = mxl*vq[lm-1] + mxu*vq[lm+1];
	}
	lm = nlmlim;
		s2d mxl = vdup(mx[2*lm]);
		vr[lm] = mxl*vq[lm-1];
}

//@}

// truncation at LMAX and MMAX
#define LTR LMAX
#define MTR MMAX

/** \addtogroup local Local and partial evaluation of SH fields.
 * These do only require a call to \ref shtns_create, but not to \ref shtns_set_grid.
 * These functions are not optimized and can be relatively slow, but they provide good
 * reference implemenation for the transforms.
*/
//@{

/// Evaluate scalar SH representation \b Qlm at physical point defined by \b cost = cos(theta) and \b phi
double SH_to_point(shtns_cfg shtns, complex double *Qlm, double cost, double phi)
{
	double yl[LMAX+1];
	double vr0, vr1;
	long int l,m,im;

	vr0 = 0.0;		vr1 = 0.0;
	m=0;	im=0;
		legendre_sphPlm_array(shtns, LTR, im, cost, &yl[m]);
		for (l=m; l<LTR; l+=2) {
			vr0 += yl[l] * creal( Qlm[l] );
			vr1 += yl[l+1] * creal( Qlm[l+1] );
		}
		if (l==LTR) {
			vr0 += yl[l] * creal( Qlm[l] );
		}
		vr0 += vr1;
	if (MTR>0) {
		complex double eip, eimp;
		eip = cos(phi*MRES) + I*sin(phi*MRES);	eimp = 2.0;
		im = 1;  do {
			m = im*MRES;
			legendre_sphPlm_array(shtns, LTR, im, cost, &yl[m]);
			v2d* Ql = (v2d*) &Qlm[LiM(shtns, 0,im)];	// virtual pointer for l=0 and im
			v2d vrm0 = vdup(0.0);		v2d vrm1 = vdup(0.0);
			for (l=m; l<LTR; l+=2) {
				vrm0 += vdup(yl[l]) * Ql[l];
				vrm1 += vdup(yl[l+1]) * Ql[l+1];
			}
//			eimp = 2.*(cos(m*phi) + I*sin(m*phi));
			eimp *= eip;			// not so accurate, but it should be enough for rendering uses.
			vrm0 += vrm1;
			if (l==LTR) {
				vrm0 += vdup(yl[l]) * Ql[l];
			}
			vr0 += vcplx_real(vrm0)*creal(eimp) - vcplx_imag(vrm0)*cimag(eimp);
		} while(++im <= MTR);
	}
	return vr0;
}

void SH_to_grad_point(shtns_cfg shtns, complex double *DrSlm, complex double *Slm, double cost, double phi,
					   double *gr, double *gt, double *gp)
{
	double yl[LMAX+1];
	double dtyl[LMAX+1];
	double vtt, vpp, vr0, vrm;
	long int l,m,im;

	const double sint = sqrt((1.-cost)*(1.+cost));
	vtt = 0.;  vpp = 0.;  vr0 = 0.;  vrm = 0.;
	m=0;	im=0;
		legendre_sphPlm_deriv_array(shtns, LTR, im, cost, sint, &yl[m], &dtyl[m]);
		for (l=m; l<=LTR; ++l) {
			vr0 += yl[l] * creal( DrSlm[l] );
			vtt += dtyl[l] * creal( Slm[l] );
		}
	if (MTR>0) {
		complex double eip, eimp, imeimp;
		eip = cos(phi*MRES) + I*sin(phi*MRES);	eimp = 2.0;
		im=1;  do {
			m = im*MRES;
			legendre_sphPlm_deriv_array(shtns, LTR, im, cost, sint, &yl[m], &dtyl[m]);
//			eimp = 2.*(cos(m*phi) + I*sin(m*phi));
			eimp *= eip;		// not so accurate, but it should be enough for rendering uses.
			imeimp = eimp*m*I;
			l = LiM(shtns, 0,im);
			v2d* Ql = (v2d*) &DrSlm[l];		v2d* Sl = (v2d*) &Slm[l];
			v2d qm = vdup(0.0);
			v2d dsdt = vdup(0.0);		v2d dsdp = vdup(0.0);
			for (l=m; l<=LTR; ++l) {
				qm += vdup(yl[l]) * Ql[l];
				dsdt += vdup(dtyl[l]) * Sl[l];
				dsdp += vdup(yl[l]) * Sl[l];
			}
			vrm += vcplx_real(qm)*creal(eimp) - vcplx_imag(qm)*cimag(eimp);			// dS/dr
			vtt += vcplx_real(dsdt)*creal(eimp) - vcplx_imag(dsdt)*cimag(eimp);		// dS/dt
			vpp += vcplx_real(dsdp)*creal(imeimp) - vcplx_imag(dsdp)*cimag(imeimp);	// + I.m/sint *S
		} while (++im <= MTR);
		vr0 += vrm*sint;
	}
	*gr = vr0;	// Gr = dS/dr
	*gt = vtt;	// Gt = dS/dt
	*gp = vpp;	// Gp = I.m/sint *S
}

/// Evaluate vector SH representation \b Qlm at physical point defined by \b cost = cos(theta) and \b phi
void SHqst_to_point(shtns_cfg shtns, complex double *Qlm, complex double *Slm, complex double *Tlm, double cost, double phi,
					   double *vr, double *vt, double *vp)
{
	double yl[LMAX+1];
	double dtyl[LMAX+1];
	double vtt, vpp, vr0, vrm;
	long int l,m,im;

	const double sint = sqrt((1.-cost)*(1.+cost));
	vtt = 0.;  vpp = 0.;  vr0 = 0.;  vrm = 0.;
	m=0;	im=0;
		legendre_sphPlm_deriv_array(shtns, LTR, im, cost, sint, &yl[m], &dtyl[m]);
		for (l=m; l<=LTR; ++l) {
			vr0 += yl[l] * creal( Qlm[l] );
			vtt += dtyl[l] * creal( Slm[l] );
			vpp -= dtyl[l] * creal( Tlm[l] );
		}
	if (MTR>0) {
		complex double eip, eimp, imeimp;
		eip = cos(phi*MRES) + I*sin(phi*MRES);	eimp = 2.0;
		im=1;  do {
			m = im*MRES;
			legendre_sphPlm_deriv_array(shtns, LTR, im, cost, sint, &yl[m], &dtyl[m]);
//			eimp = 2.*(cos(m*phi) + I*sin(m*phi));
			eimp *= eip;		// not so accurate, but it should be enough for rendering uses.
			imeimp = eimp*m*I;
			l = LiM(shtns, 0,im);
			v2d* Ql = (v2d*) &Qlm[l];	v2d* Sl = (v2d*) &Slm[l];	v2d* Tl = (v2d*) &Tlm[l];
			v2d qm = vdup(0.0);
			v2d dsdt = vdup(0.0);		v2d dtdt = vdup(0.0);
			v2d dsdp = vdup(0.0);		v2d dtdp = vdup(0.0);
			for (l=m; l<=LTR; ++l) {
				qm += vdup(yl[l]) * Ql[l];
				dsdt += vdup(dtyl[l]) * Sl[l];
				dtdt += vdup(dtyl[l]) * Tl[l];
				dsdp += vdup(yl[l]) * Sl[l];
				dtdp += vdup(yl[l]) * Tl[l];
			}
			vrm += vcplx_real(qm)*creal(eimp) - vcplx_imag(qm)*cimag(eimp);
			vtt += (vcplx_real(dtdp)*creal(imeimp) - vcplx_imag(dtdp)*cimag(imeimp))	// + I.m/sint *T
					+ (vcplx_real(dsdt)*creal(eimp) - vcplx_imag(dsdt)*cimag(eimp));	// + dS/dt
			vpp += (vcplx_real(dsdp)*creal(imeimp) - vcplx_imag(dsdp)*cimag(imeimp))	// + I.m/sint *S
					- (vcplx_real(dtdt)*creal(eimp) - vcplx_imag(dtdt)*cimag(eimp));	// - dT/dt
		} while (++im <= MTR);
		vr0 += vrm*sint;
	}
	*vr = vr0;
	*vt = vtt;	// Bt = I.m/sint *T  + dS/dt
	*vp = vpp;	// Bp = I.m/sint *S  - dT/dt
}
//@}
	
#undef LTR
#undef MTR


/*
	SYNTHESIS AT A GIVEN LATITUDE
	(does not require a previous call to shtns_set_grid)
*/

fftw_plan ifft_lat = NULL;		///< fftw plan for SHqst_to_lat
long int nphi_lat = 0;			///< nphi of previous SHqst_to_lat
double* ylm_lat = NULL;
double* dylm_lat;
double ct_lat = 2.0;
double st_lat;

/// synthesis at a given latitude, on nphi equispaced longitude points.
/// vr, vt, and vp arrays must have nphi+2 doubles allocated (fftw requirement).
/// It does not require a previous call to shtns_set_grid, but it is NOT thread-safe.
/// \ingroup local
void SHqst_to_lat(shtns_cfg shtns, complex double *Qlm, complex double *Slm, complex double *Tlm, double cost,
					double *vr, double *vt, double *vp, int nphi, int ltr, int mtr)
{
	complex double vst, vtt, vsp, vtp, vrr;
	complex double *vrc, *vtc, *vpc;
	long int m, l, j;

	if (ltr > LMAX) ltr=LMAX;
	if (mtr > MMAX) mtr=MMAX;
	if (mtr*MRES > ltr) mtr=ltr/MRES;
	if (mtr*2*MRES >= nphi) mtr = (nphi-1)/(2*MRES);

	vrc = (complex double *) vr;
	vtc = (complex double *) vt;
	vpc = (complex double *) vp;

	if ((nphi != nphi_lat)||(ifft_lat == NULL)) {
		if (ifft_lat != NULL) fftw_destroy_plan(ifft_lat);
		#ifdef OMP_FFTW
			fftw_plan_with_nthreads(1);
		#endif
		ifft_lat = fftw_plan_dft_c2r_1d(nphi, vrc, vr, FFTW_ESTIMATE);
		nphi_lat = nphi;
	}
	if (ylm_lat == NULL) {
		ylm_lat = (double *) malloc(sizeof(double)* NLM*2);
		dylm_lat = ylm_lat + NLM;
	}
	if (cost != ct_lat) {		// don't recompute if same latitude (ie equatorial disc rendering)
		st_lat = sqrt((1.-cost)*(1.+cost));	// sin(theta)
		for (m=0,j=0; m<=mtr; ++m) {
			legendre_sphPlm_deriv_array(shtns, ltr, m, cost, st_lat, &ylm_lat[j], &dylm_lat[j]);
			j += LMAX -m*MRES +1;
		}
	}

	for (m = 0; m<nphi/2+1; ++m) {	// init with zeros
		vrc[m] = 0.0;	vtc[m] = 0.0;	vpc[m] = 0.0;
	}
	j=0;
	m=0;
		vrr=0;	vtt=0;	vst=0;
		for(l=m; l<=ltr; ++l, ++j) {
			vrr += ylm_lat[j] * creal(Qlm[j]);
			vst += dylm_lat[j] * creal(Slm[j]);
			vtt += dylm_lat[j] * creal(Tlm[j]);
		}
		j += (LMAX-ltr);
		vrc[m] = vrr;
		vtc[m] =  vst;	// Vt =   dS/dt
		vpc[m] = -vtt;	// Vp = - dT/dt
	for (m=MRES; m<=mtr*MRES; m+=MRES) {
		vrr=0;	vtt=0;	vst=0;	vsp=0;	vtp=0;
		for(l=m; l<=ltr; ++l, ++j) {
			vrr += ylm_lat[j] * Qlm[j];
			vst += dylm_lat[j] * Slm[j];
			vtt += dylm_lat[j] * Tlm[j];
			vsp += ylm_lat[j] * Slm[j];
			vtp += ylm_lat[j] * Tlm[j];
		}
		j+=(LMAX-ltr);
		vrc[m] = vrr*st_lat;
		vtc[m] = I*m*vtp + vst;	// Vt = I.m/sint *T  + dS/dt
		vpc[m] = I*m*vsp - vtt;	// Vp = I.m/sint *S  - dT/dt
	}
	fftw_execute_dft_c2r(ifft_lat,vrc,vr);
	fftw_execute_dft_c2r(ifft_lat,vtc,vt);
	fftw_execute_dft_c2r(ifft_lat,vpc,vp);
//	free(ylm_lat);
}

/// synthesis at a given latitude, on nphi equispaced longitude points.
/// vr arrays must have nphi+2 doubles allocated (fftw requirement).
/// It does not require a previous call to shtns_set_grid, but it is NOT thread-safe.
/// \ingroup local
void SH_to_lat(shtns_cfg shtns, complex double *Qlm, double cost,
					double *vr, int nphi, int ltr, int mtr)
{
	complex double vrr;
	complex double *vrc;
	long int m, l, j;

	if (ltr > LMAX) ltr=LMAX;
	if (mtr > MMAX) mtr=MMAX;
	if (mtr*MRES > ltr) mtr=ltr/MRES;
	if (mtr*2*MRES >= nphi) mtr = (nphi-1)/(2*MRES);

	vrc = (complex double *) vr;

	if ((nphi != nphi_lat)||(ifft_lat == NULL)) {
		if (ifft_lat != NULL) fftw_destroy_plan(ifft_lat);
		#ifdef OMP_FFTW
			fftw_plan_with_nthreads(1);
		#endif
		ifft_lat = fftw_plan_dft_c2r_1d(nphi, vrc, vr, FFTW_ESTIMATE);
		nphi_lat = nphi;
	}
	if (ylm_lat == NULL) {
		ylm_lat = (double *) malloc(sizeof(double)* NLM*2);
		dylm_lat = ylm_lat + NLM;
	}
	if (cost != ct_lat) {		// don't recompute if same latitude (ie equatorial disc rendering)
		st_lat = sqrt((1.-cost)*(1.+cost));	// sin(theta)
		for (m=0,j=0; m<=mtr; ++m) {
			legendre_sphPlm_deriv_array(shtns, ltr, m, cost, st_lat, &ylm_lat[j], &dylm_lat[j]);
			j += LMAX -m*MRES +1;
		}
	}

	for (m = 0; m<nphi/2+1; ++m) {	// init with zeros
		vrc[m] = 0.0;
	}
	j=0;
	m=0;
		vrr=0;
		for(l=m; l<=ltr; ++l, ++j) {
			vrr += ylm_lat[j] * creal(Qlm[j]);
		}
		j += (LMAX-ltr);
		vrc[m] = vrr;
	for (m=MRES; m<=mtr*MRES; m+=MRES) {
		vrr=0;
		for(l=m; l<=ltr; ++l, ++j) {
			vrr += ylm_lat[j] * Qlm[j];
		}
		j+=(LMAX-ltr);
		vrc[m] = vrr*st_lat;
	}
	fftw_execute_dft_c2r(ifft_lat,vrc,vr);
//	free(ylm_lat);
}

/*
	INTERNAL INITIALIZATION FUNCTIONS
*/

char* sht_name[SHT_NALG] = {"dct", "mem", "s+v", "fly1", "fly2", "fly3", "fly4", "fly6", "fly8", "omp1", "omp2", "omp3", "omp4", "omp6", "omp8" };
char* sht_var[SHT_NVAR] = {"std", "ltr"};
char *sht_type[SHT_NTYP] = {"syn", "ana", "vsy", "van", "gsp", "gto", "v3s", "v3a" };
int sht_npar[SHT_NTYP] = {2, 2, 4, 4, 3, 3, 6, 6};

extern void* sht_array[SHT_NALG][SHT_NTYP];
extern void* sht_array_l[SHT_NALG][SHT_NTYP];
extern void* sht_array_m0[SHT_NALG][SHT_NTYP];
extern void* sht_array_m0l[SHT_NALG][SHT_NTYP];

// big array holding all sht functions, variants and algorithms
void* sht_func[SHT_NVAR][SHT_NTYP][SHT_NALG];

/// \internal use on-the-fly alogorithm (guess without measuring)
static void set_sht_fly(shtns_cfg shtns, int typ_start)
{
	int algo = SHT_FLY2;
	if (shtns->nthreads > 1) algo = SHT_OMP2;
	for (int it=typ_start; it<SHT_NTYP; it++) {
		for (int v=0; v<SHT_NVAR; v++)
			shtns->fptr[v][it] = sht_func[v][it][algo];
	}
}

/// \internal choose memory algorithm everywhere.
static void set_sht_mem(shtns_cfg shtns) {
	for (int it=0; it<SHT_NTYP; it++) {
		for (int v=0; v<SHT_NVAR; v++)
			shtns->fptr[v][it] = sht_func[v][it][SHT_MEM];
	}
}

/// \internal copy all algos to sht_func array (should be called by set_grid befor choosing variants).
/// if nphi is 1, axisymmetric algorythms are used.
static void init_sht_array_func(shtns_cfg shtns)
{
	int it, j;
	int alg_lim = SHT_NALG-1;

	if (shtns->nlat_2 < 8*VSIZE2) {		// limit available on-the-fly algorithm to avoid overflow (and segfaults).
		it = shtns->nlat_2 / VSIZE2;
		switch(it) {
			case 0 : alg_lim = SHT_FLY1-1; break;
			case 1 : alg_lim = SHT_FLY1; break;
			case 2 : alg_lim = SHT_FLY2; break;
			case 3 : alg_lim = SHT_FLY3; break;
			case 4 : ;
			case 5 : alg_lim = SHT_FLY4; break;
			default : alg_lim = SHT_FLY6;
		}
	}

	for (it=0; it<SHT_NTYP; it++) {
		for (j=0; j<=alg_lim; j++) {
			sht_func[SHT_LTR][it][j] = sht_array_l[j][it];
			if (j >= SHT_FLY1) {
				sht_func[SHT_STD][it][j] = sht_array_l[j][it];		// on-the-fly only exist in LTR version
			} else  sht_func[SHT_STD][it][j] = sht_array[j][it];

			if ((shtns->nphi == 1)&&(j != SHT_SV)) {		// axisymmetric transform requested (but not for special case S+V)
				sht_func[SHT_LTR][it][j] = sht_array_m0l[j][it];
				if (j >= SHT_FLY1) {
					sht_func[SHT_STD][it][j] = sht_array_m0l[j][it];		// on-the-fly only exist in LTR version
				} else  sht_func[SHT_STD][it][j] = sht_array_m0[j][it];
			}
		}
		for (j=alg_lim+1; j<SHT_NALG; j++) {
			sht_func[SHT_LTR][it][j] = NULL;
			sht_func[SHT_STD][it][j] = NULL;
		}
	}
	set_sht_mem(shtns);		// default transform is MEM
}

/// \internal return the smallest power of 2 larger than n.
static int next_power_of_2(int n)
{
	int f = 1;
	if ( (n<=0) || (n>(1<<(sizeof(int)*8-2))) ) return 0;
	while (f<n) f*=2;
	return f;
}

/// \internal find the closest integer that is larger than n and that contains only factors up to fmax.
/// fmax is 7 for optimal FFTW fourier transforms.
/// return only even integers for n>fmax.
static int fft_int(int n, int fmax)
{
	int k,f;

	if (n<=fmax) return n;
	if (fmax<2) return 0;
	if (fmax==2) return next_power_of_2(n);

	n -= 2-(n&1);		// only even n
	do {
		n+=2;	f=2;
		while ((2*f <= n) && ((n&f)==0)) f *= 2;		// no divisions for factor 2.
		k=3;
		while ((k<=fmax) &&  (k*f <= n)) {
			while ((k*f <= n) && (n%(k*f)==0)) f *= k;
			k+=2;
		}
	} while (f != n);

	k = next_power_of_2(n);			// what is the closest power of 2 ?
	if ((k-n)*33 < n) return k;		// rather choose power of 2 if not too far (3%)

	return n;
}


/// \internal returns an aproximation of the memory usage in mega bytes for the scalar matrices.
/// \ingroup init
static double sht_mem_size(int lmax, int mmax, int mres, int nlat)
{
	double s = 1./(1024*1024);
	s *= ((nlat+1)/2) * sizeof(double) * nlm_calc(lmax, mmax, mres);
	return s;
}

/// \internal return the number of shtns config that contain a reference to the memory location *pp.
static int ref_count(shtns_cfg shtns, void* pp)
{
	shtns_cfg s2 = sht_data;
	void* p;
	long int nref, offset;

	if ((pp==NULL) || (shtns==NULL)) return -1;		// error.

	p = *(void**)pp;			// the pointer to memory location that we want to free
	if (p == NULL) return 0;	// nothing to do.

	offset = (char*)pp - (char*)shtns;			// relative location in shtns_info structure.
	nref = 0;		// reference count.
	while (s2 != NULL) {		// scan all configs.
		if ( *(void**)(((char*)s2) + offset) == p ) nref++;
		s2 = s2->next;
	}
	return nref;		// will be >= 1 (as shtns is included in search)
}

/// \internal check if the memory location *pp is referenced by another sht config before freeing it.
/// returns the number of other references, or -1 on error. If >=1, the memory location could not be freed.
/// If the return value is 0, the ressource has been freed.
static int free_unused(shtns_cfg shtns, void* pp)
{
	int n = ref_count(shtns, pp);	// howmany shtns config do reference this memory location ?
	if (n <= 0) return n;		// nothing to free.
	if (n == 1) {				// no other reference found...
		void** ap = (void**) pp;
		free(*ap);		// ...so we can free it...
		*ap = NULL;		// ...and mark as unaloccated.
	}
	return (n-1);
}

/// \internal Free matrices if on-the-fly has been selected.
static void free_unused_matrices(shtns_cfg shtns)
{
	long int marray_size = (MMAX+1)*sizeof(double) + (MIN_ALIGNMENT-1);
	int count[SHT_NTYP];
	int it, iv, ia, iv2;

	for (it=0; it<SHT_NTYP; it++) {
		count[it] = 0;
		for (iv=0; iv<SHT_NVAR; iv++) {
			if (shtns->fptr[iv][it] != NULL) {
				for (ia=0; ia<SHT_SV; ia++)		// count occurences to mem algo
					for (iv2=0; iv2<SHT_NVAR; iv2++)
						if (sht_func[iv2][it][ia] == shtns->fptr[iv][it]) count[it]++;
			}
		}
		#if SHT_VERBOSE > 1
			printf(" %d ",count[it]);
		#endif
	}

	if (shtns->dylm == NULL) {		// no vector transform : do not try to free
		count[SHT_TYP_VAN] = -1;		count[SHT_TYP_VSY] = -1;
	}

	if (count[SHT_TYP_3AN] == 0) {		// analysis may be freed.
		if (count[SHT_TYP_VAN] == 0) {
			PRINT_VERB("freeing vector analysis matrix\n");
			free_unused(shtns, &shtns->dzlm);
		}
		if (count[SHT_TYP_SAN] == 0) {
			PRINT_VERB("freeing scalar analysis matrix\n");
			free_unused(shtns, &shtns->zlm);
		}
	} else if (shtns->mmax > 0) {	// scalar may be reduced to m=0
		if ((count[SHT_TYP_SAN] == 0) && (ref_count(shtns, &shtns->zlm) == 1)) {
			PRINT_VERB("keeping scalar analysis matrix up to m=0\n");
			shtns->zlm = realloc(shtns->zlm, ((LMAX+1)*NLAT_2 +3)*sizeof(double) + marray_size );
		}
	}
	if (count[SHT_TYP_3SY] == 0) {		// synthesis may be freed.
		if (count[SHT_TYP_VSY] + count[SHT_TYP_GSP] + count[SHT_TYP_GTO] == 0) {
			PRINT_VERB("freeing vector synthesis matrix\n");
			free_unused(shtns, &shtns->dylm);
		}
		if (count[SHT_TYP_SSY] == 0) {
			PRINT_VERB("freeing scalar synthesis matrix\n");
			free_unused(shtns, &shtns->ylm);
		}
	} else if (shtns->mmax > 0) {	// scalar may be reduced to m=0
		if ((count[SHT_TYP_SSY] == 0) && (ref_count(shtns, &shtns->ylm) == 1)) {
			PRINT_VERB("keeping scalar synthesis matrix up to m=0\n");
			shtns->ylm = realloc(shtns->ylm, (LMAX+2)*NLAT_2*sizeof(double) + marray_size );
		}
	}
}


/// \internal allocate arrays for SHT related to a given grid.
static void alloc_SHTarrays(shtns_cfg shtns, int on_the_fly, int vect, int analys)
{
	long int im, l0;
	long int size, marray_size, lstride;

	l0 = ((NLAT+1)>>1)*2;		// round up to even
	shtns->ct = (double *) VMALLOC( sizeof(double) * l0*3 );			/// ct[] (including st and st_1)
	shtns->st = shtns->ct + l0;		shtns->st_1 = shtns->ct + 2*l0;
	
	shtns->ylm = NULL;		shtns->dylm = NULL;			// synthesis
	shtns->zlm = NULL;		shtns->dzlm = NULL;			// analysis

	if (on_the_fly == 0) {		// Allocate legendre functions lookup tables.
		marray_size = (MMAX+1)*sizeof(double*) + (MIN_ALIGNMENT-1);		// for sse2 alignement

		/* ylm */
		lstride = (LMAX+1);		lstride += (lstride&1);		// even stride.
		size = sizeof(double) * ((NLM-(LMAX+1)+lstride)*NLAT_2);
		if (MMAX == 0) size += 3*sizeof(double);			// some overflow needed.
		shtns->ylm = (double **) malloc( marray_size + size );
		shtns->ylm[0] = (double *) PTR_ALIGN( shtns->ylm + (MMAX+1) );
		if (MMAX>0) shtns->ylm[1] = shtns->ylm[0] + NLAT_2*lstride;
		for (im=1; im<MMAX; im++) shtns->ylm[im+1] = shtns->ylm[im] + NLAT_2*(LMAX+1-im*MRES);
		/* dylm */		
		if (vect) {
			lstride = LMAX;		lstride += (lstride&1);		// even stride.
			size = sizeof(struct DtDp) * ((NLM-(LMAX+1) +lstride/2)*NLAT_2);
			if (MMAX == 0) size += 2*sizeof(double);	// some overflow needed.
			shtns->dylm = (struct DtDp **) malloc( marray_size + size );
			shtns->dylm[0] = (struct DtDp *) PTR_ALIGN( shtns->dylm + (MMAX+1) );
			if (MMAX>0) shtns->dylm[1] = shtns->dylm[0] + (lstride/2)*NLAT_2;		// phi-derivative is zero for m=0
			for (im=1; im<MMAX; im++) shtns->dylm[im+1] = shtns->dylm[im] + NLAT_2*(LMAX+1-im*MRES);
		}
		if (analys) {
			/* zlm */
			size = sizeof(double) * (NLM*NLAT_2 + (NLAT_2 & 1));
			if (MMAX == 0) size += 2*(NLAT_2 & 1)*((LMAX+1) & 1) * sizeof(double);
			shtns->zlm = (double **) malloc( marray_size + size );
			shtns->zlm[0] = (double *) PTR_ALIGN( shtns->zlm + (MMAX+1) );
			if (MMAX>0) shtns->zlm[1] = shtns->zlm[0] + NLAT_2*(LMAX+1) + (NLAT_2&1);
			for (im=1; im<MMAX; im++) shtns->zlm[im+1] = shtns->zlm[im] + NLAT_2*(LMAX+1-im*MRES);
			/* dzlm */
			if (vect) {
				size = sizeof(struct DtDp)* (NLM-1)*NLAT_2;		// remove l=0
				shtns->dzlm = (struct DtDp **) malloc( marray_size + size );
				shtns->dzlm[0] = (struct DtDp *) PTR_ALIGN( shtns->dzlm + (MMAX+1) );
				if (MMAX>0) shtns->dzlm[1] = shtns->dzlm[0] + NLAT_2*(LMAX);
				for (im=1; im<MMAX; im++) shtns->dzlm[im+1] = shtns->dzlm[im] + NLAT_2*(LMAX+1-im*MRES);
			}
		}
	}
	#if SHT_VERBOSE > 1
		printf("          Memory used for Ylm and Zlm matrices = %.3f Mb x2\n",3.0*sizeof(double)*NLM*NLAT_2/(1024.*1024.));
	#endif
}


/// \internal free arrays allocated by init_SH_dct
static void free_SH_dct(shtns_cfg shtns)
{
	if (shtns->zlm_dct0 == NULL) return;

	if (ref_count(shtns, &shtns->dzlm_dct0) == 1) VFREE(shtns->dzlm_dct0);
	if (ref_count(shtns, &shtns->zlm_dct0) == 1)  VFREE(shtns->zlm_dct0);
	shtns->dzlm_dct0 = NULL;	shtns->zlm_dct0 = NULL;

	free_unused(shtns, &shtns->dykm_dct);
	free_unused(shtns, &shtns->ykm_dct);

	if (ref_count(shtns, &shtns->idct) == 1)    fftw_destroy_plan(shtns->idct);		// free unused dct plans
	if (ref_count(shtns, &shtns->dct_m0) == 1)  fftw_destroy_plan(shtns->dct_m0);
	shtns->idct = NULL;		shtns->dct_m0 = NULL;
}

/// \internal free arrays allocated by alloc_SHTarrays.
static void free_SHTarrays(shtns_cfg shtns)
{
	free_SH_dct(shtns);
	free_unused(shtns, &shtns->ylm);
	free_unused(shtns, &shtns->dylm);
	free_unused(shtns, &shtns->zlm);
	free_unused(shtns, &shtns->dzlm);

	if (ref_count(shtns, &shtns->ct) == 1)	VFREE(shtns->ct);
	shtns->ct = NULL;		shtns->st = NULL;

	if (ref_count(shtns, &shtns->fft) == 1)  fftw_destroy_plan(shtns->fft);
	if (ref_count(shtns, &shtns->ifft) == 1) fftw_destroy_plan(shtns->ifft);
	shtns->fft = NULL;		shtns->ifft = NULL;		shtns->ncplx_fft = -1;	// no fft
}

#ifndef HAVE_FFTW_COST
	// substitute undefined symbol in fftw older than 3.3
	#define fftw_cost(a) 0.0
#endif

/// \internal initialize FFTs using FFTW.
/// \param[in] layout defines the spatial layout (see \ref spat).
/// \param[in] on_the_fly is one, if only on-the-fly transform are considered.
/// returns the number of double to be allocated for a spatial field.
static void planFFT(shtns_cfg shtns, int layout, int on_the_fly)
{
	double cost_fft_ip, cost_fft_oop, cost_ifft_ip, cost_ifft_oop;
	complex double *ShF;
	double *Sh;
	fftw_plan fft2, ifft2, fft, ifft;
	int nfft, ncplx, nreal;
	int theta_inc, phi_inc, phi_embed;
	int in_place = 1;		// try to use in-place fft.

	if (NPHI <= 2*MMAX) shtns_runerr("the sampling condition Nphi > 2*Mmax is not met.");

	#ifdef OMP_FFTW
		if ((shtns->fftw_plan_mode & (FFTW_EXHAUSTIVE | FFTW_PATIENT)) && (omp_threads > 1)) {
			shtns->fftw_plan_mode = FFTW_PATIENT;
			fftw_plan_with_nthreads(omp_threads);
		} else fftw_plan_with_nthreads(shtns->nthreads);
	#endif

	shtns->fft = NULL;		shtns->ifft = NULL;
	if (NPHI==1) 	// no FFT needed.
	{
		shtns->fftc_mode = -1;		// no FFT
		#if SHT_VERBOSE > 0
			printf("        => no fft : Mmax=0, Nphi=1, Nlat=%d\n",NLAT);
		#endif
		shtns->nspat = NLAT;
		shtns->ncplx_fft = -1;	// no fft.
	}
	else	/* NPHI > 1 */
	{
		theta_inc=1;  phi_inc=NLAT;  phi_embed=2*(NPHI/2+1);	// SHT_NATIVE_LAYOUT is the default.
		if (layout & SHT_THETA_CONTIGUOUS) {	theta_inc=1;  phi_inc=NLAT;  phi_embed=NPHI;	}
		if (layout & SHT_PHI_CONTIGUOUS)   {	phi_inc=1;  theta_inc=NPHI;  phi_embed=NPHI;	}
		nfft = NPHI;
		ncplx = NPHI/2 +1;
		nreal = phi_embed;
		if ((theta_inc != 1)||(phi_inc != NLAT)||(nreal < 2*ncplx))  in_place = 0;		// we need to do the fft out-of-place.

		#if SHT_VERBOSE > 0
			printf("        => using FFTW : Mmax=%d, Nphi=%d, Nlat=%d  (data layout : phi_inc=%d, theta_inc=%d, phi_embed=%d)\n",MMAX,NPHI,NLAT,phi_inc,theta_inc,phi_embed);
			if (NPHI <= (SHT_NL_ORDER+1)*MMAX)	printf("     !! Warning : anti-aliasing condition Nphi > %d*Mmax is not met !\n", SHT_NL_ORDER+1);
			if (NPHI != fft_int(NPHI,7))		printf("     !! Warning : Nphi is not optimal for FFTW !\n");
		#endif

	// Allocate dummy Spatial Fields.
		ShF = (complex double *) VMALLOC(ncplx * NLAT * sizeof(complex double));
		Sh = (double *) VMALLOC(ncplx * NLAT * sizeof(complex double));
		fft = NULL;		ifft = NULL;	fft2 = NULL;	ifft2 = NULL;

	// complex fft for fly transform is a bit different.
		if (layout & SHT_PHI_CONTIGUOUS) {		// out-of-place split dft
			fftw_iodim dim, many;
			shtns->fftc_mode = 1;
			dim.n = NPHI;    	dim.os = 1;			dim.is = NLAT;		// complex transpose
			many.n = NLAT/2;	many.os = 2*NPHI;	many.is = 2;
			shtns->ifftc = fftw_plan_guru_split_dft(1, &dim, 1, &many, ((double*)ShF)+1, (double*)ShF, Sh+NPHI, Sh, shtns->fftw_plan_mode);
			dim.n = NPHI;    	dim.is = 1;			dim.os = NLAT;		// complex transpose
			many.n = NLAT/2;	many.is = 2*NPHI;	many.os = 2;
			shtns->fftc = fftw_plan_guru_split_dft(1, &dim, 1, &many,  Sh+NPHI, Sh, ((double*)ShF)+1, (double*)ShF, shtns->fftw_plan_mode);
			#if SHT_VERBOSE > 1
				printf("          fftw cost ifftc=%lg,  fftc=%lg  ",fftw_cost(shtns->ifftc), fftw_cost(shtns->fftc));	fflush(stdout);
			#endif
		} else {		// use only in-place here, supposed to be faster.
			shtns->fftc_mode = 0;
			shtns->ifftc = fftw_plan_many_dft(1, &nfft, NLAT/2, ShF, &nfft, NLAT/2, 1, ShF, &nfft, NLAT/2, 1, FFTW_BACKWARD, shtns->fftw_plan_mode);
			shtns->fftc = shtns->ifftc;		// same thing, with m>0 and m<0 exchanged.
			#if SHT_VERBOSE > 1
				printf("          fftw cost ifftc=%lg  ",fftw_cost(shtns->ifftc));	fflush(stdout);
			#endif
		}

	#if _GCC_VEC_
	  if (on_the_fly == 0)
	#endif
	  {		// the real ffts are required if _GCC_VEC_ == 0 or if on_the_fly is zero.
	// IFFT : unnormalized.  FFT : must be normalized.
		cost_fft_ip = 0.0;	cost_ifft_ip = 0.0;		cost_fft_oop = 0.0;		cost_ifft_oop = 0.0;
		if (in_place) {		// in-place FFT (if allowed)
			ifft2 = fftw_plan_many_dft_c2r(1, &nfft, NLAT, ShF, &ncplx, NLAT, 1, (double*) ShF, &nreal, phi_inc, theta_inc, shtns->fftw_plan_mode);
			#if SHT_VERBOSE > 1
				printf("          in-place cost : ifft=%lg  ",fftw_cost(ifft2));	fflush(stdout);
			#endif
			if (ifft2 != NULL) {
				fft2 = fftw_plan_many_dft_r2c(1, &nfft, NLAT, (double*) ShF, &nreal, phi_inc, theta_inc, ShF, &ncplx, NLAT, 1, shtns->fftw_plan_mode);
				#if SHT_VERBOSE > 1
					printf("fft=%lg\n",fftw_cost(fft2));	fflush(stdout);
				#endif
				if (fft2 != NULL) {
					cost_fft_ip = fftw_cost(fft2);		cost_ifft_ip = fftw_cost(ifft2);
				}
			}
		}
		if ( (in_place == 0) || (cost_fft_ip * cost_ifft_ip > 0.0) )
		{		// out-of-place FFT
			ifft = fftw_plan_many_dft_c2r(1, &nfft, NLAT, ShF, &ncplx, NLAT, 1, Sh, &nreal, phi_inc, theta_inc, shtns->fftw_plan_mode);
			#if SHT_VERBOSE > 1
				printf("          oop cost : ifft=%lg  ",fftw_cost(ifft));	fflush(stdout);
			#endif
			if (ifft == NULL) shtns_runerr("[FFTW] ifft planning failed !");
			fft = fftw_plan_many_dft_r2c(1, &nfft, NLAT, Sh, &nreal, phi_inc, theta_inc, ShF, &ncplx, NLAT, 1, shtns->fftw_plan_mode);
			#if SHT_VERBOSE > 1
				printf("fft=%lg\n",fftw_cost(fft));	fflush(stdout);
			#endif
			if (fft == NULL) shtns_runerr("[FFTW] fft planning failed !");
			cost_fft_oop = fftw_cost(fft);		cost_ifft_oop = fftw_cost(ifft);
		}
		if ( (cost_fft_ip * cost_ifft_ip > 0.0) && (cost_fft_oop * cost_ifft_oop > 0.0) ) {		// both have been succesfully timed.
			if ( cost_fft_oop + SHT_NL_ORDER*cost_ifft_oop < cost_fft_ip + SHT_NL_ORDER*cost_ifft_ip )
				in_place = 0;		// disable in-place, because out-of-place is faster.
		}

		if (in_place) {
			/* IN-PLACE FFT */
			if (fft != NULL)  fftw_destroy_plan(fft);
			if (ifft != NULL) fftw_destroy_plan(ifft);
			fft = fft2;		ifft = ifft2;
			shtns->ncplx_fft = 0;		// fft is done in-place, no allocation needed.
		} else {
			/* OUT-OF-PLACE FFT */
			if (fft2 != NULL)  fftw_destroy_plan(fft2);		// use OUT-OF-PLACE FFT
			if (ifft2 != NULL) fftw_destroy_plan(ifft2);
			shtns->ncplx_fft = ncplx * NLAT;		// fft is done out-of-place, store allocation size.
			phi_embed = NPHI;
			#if SHT_VERBOSE > 1
				printf("        ** out-of-place fft **\n");
			#endif
		}
		shtns->fft = fft;		shtns->ifft = ifft;
	  }
		shtns->nspat = phi_embed * NLAT;
		VFREE(Sh);		VFREE(ShF);

		#if SHT_VERBOSE > 2
			printf(" *** fft plan :\n");
			fftw_print_plan(fft);
			printf("\n *** ifft plan :\n");
			fftw_print_plan(ifft);
			printf("\n");
		#endif
	}

	shtns->dct_m0 = NULL;	shtns->idct = NULL;		// set dct plans to uninitialized.
}

/// \internal initialize DCTs using FFTW. Must be called if MTR_DCT is changed.
static void planDCT(shtns_cfg shtns)
{
	double *Sh;
	int ndct = NLAT;
	fftw_r2r_kind r2r_kind;
	fftw_iodim dims, hdims[2];
	double Sh0[NLAT] SSE;				// temp storage on the stack, aligned.

	#ifdef OMP_FFTW
		if (shtns->fftw_plan_mode & (FFTW_EXHAUSTIVE | FFTW_PATIENT)) {
			fftw_plan_with_nthreads(omp_threads);
		} else fftw_plan_with_nthreads(1);
	#endif
	// Allocate dummy Spatial Fields.
	Sh = (double *) VMALLOC((NPHI/2 +1) * NLAT*2 * sizeof(double));

	if (shtns->dct_m0 == NULL) {		// allocate only once since it does not change.
		int stride = (NPHI>1) ? 2 : 1;
		r2r_kind = FFTW_REDFT10;		// out-of-place.
		shtns->dct_m0 = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, stride, stride*NLAT, Sh0, &ndct, 1, NLAT, &r2r_kind, shtns->fftw_plan_mode);	// out-of-place.
		if (shtns->dct_m0 == NULL)
			shtns_runerr("[FFTW] dct_m0 planning failed !");
	#if SHT_VERBOSE > 2
		printf(" *** dct_m0 plan :\n");		fftw_print_plan(shtns->dct_m0);	printf("\n");
	#endif
	}

	if (NPHI > 1) {		// complex data for NPHI>1, recompute as it does depend on MTR_DCT
		if (shtns->idct != NULL) fftw_destroy_plan(shtns->idct);

		dims.n = NLAT;	dims.is = 2;	dims.os = 2;		// real and imaginary part.
		hdims[0].n = MTR_DCT+1;	hdims[0].is = 2*NLAT; 	hdims[0].os = 2*NLAT;
		hdims[1].n = 2;			hdims[1].is = 1; 	hdims[1].os = 1;
		r2r_kind = FFTW_REDFT01;
		shtns->idct = fftw_plan_guru_r2r(1, &dims, 2, hdims, Sh, Sh, &r2r_kind, shtns->fftw_plan_mode);
	} else {		// NPHI == 1
		if (shtns->idct == NULL) {
			r2r_kind = FFTW_REDFT01;
			shtns->idct = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, 1, NLAT, Sh, &ndct, 1, NLAT, &r2r_kind, shtns->fftw_plan_mode);
		}
	}

	VFREE(Sh);
	if (shtns->idct == NULL)  shtns_runerr("[FFTW] idct planning failed !");
	#if SHT_VERBOSE > 2
		printf(" *** idct plan :\n");	fftw_print_plan(shtns->idct);	printf("\n");
	#endif

}

/// \internal SET MTR_DCT and updates fftw_plan for DCT's
static void Set_MTR_DCT(shtns_cfg shtns, int m)
{
	if ((shtns->zlm_dct0 == NULL)||(m == MTR_DCT)) return;
	if ( m < 0 ) {	// don't use dct
		shtns->mtr_dct = -1;
	} else {
		if (m>MMAX) m=MMAX;
		shtns->mtr_dct = m;
		planDCT(shtns);
	}
}

/// \internal returns the m-truncation of DCT part of synthesis
static int Get_MTR_DCT(shtns_cfg shtns) {
	return MTR_DCT;
}


/// \internal Sets the value tm[im] used for polar optimiation on-the-fly.
static void PolarOptimize(shtns_cfg shtns, double eps)
{
	int im, m, l, it;
	double v;
	double y[LMAX+1];

	for (im=0;im<=MMAX;im++)	shtns->tm[im] = 0;

	if (eps > 0.0) {
		for (im=1;im<=MMAX;im++) {
			m = im*MRES;
			it = shtns->tm[im-1] -1;	// tm[im] is monotonic.
			do {
				it++;
				legendre_sphPlm_array(shtns, LMAX, im, shtns->ct[it], y+m);
				v = 0.0;
				for (l=m; l<=LMAX; l++) {
					double ya = fabs(y[l]);
					if ( v < ya )	v = ya;
				}
			} while (v < eps);
			shtns->tm[im] = it;
		}
	#if SHT_VERBOSE > 0
		printf("        + polar optimization threshold = %.1e\n",eps);
	#endif
	#if SHT_VERBOSE > 1
		printf("          tm[im]=");
		for (im=0;im<=MMAX;im++)
			printf(" %d",shtns->tm[im]);
		printf("\n");
	#endif
	}
}

/// \internal Perform some optimization on the SHT matrices.
static void OptimizeMatrices(shtns_cfg shtns, double eps)
{
	unsigned short *tm;
	double **ylm, **zlm;
	struct DtDp** dylm;
	struct DtDp** dzlm;
	int im,m,l,it;
	int vector = (shtns->dylm != NULL);

	tm = shtns->tm;
	ylm = shtns->ylm;		dylm = shtns->dylm;
	zlm = shtns->zlm;		dzlm = shtns->dzlm;

	for (im=0;im<=MMAX;im++)	tm[im] = 0;
/// POLAR OPTIMIZATION : analyzing coefficients, some can be safely neglected.
	if (eps > 0.0) {
		for (im=1;im<=MMAX;im++) {
			m = im*MRES;
			tm[im] = NLAT_2;
			for (l=m;l<=LMAX;l++) {
				it = tm[im-1];		// tm[im] is monotonic.
				while( fabs(ylm[im][it*(LMAX-m+1) + (l-m)]) < eps ) { it++; }
				if (tm[im] > it) tm[im] = it;
			}
		}
#if SHT_VERBOSE > 0
		printf("        + polar optimization threshold = %.1e\n",eps);
#endif
#if SHT_VERBOSE > 1
		printf("          tm[im]=");
		for (im=0;im<=MMAX;im++)
			printf(" %d",tm[im]);
		printf("\n");
#endif

		for (im=1; im<=MMAX; im++) {	//	im >= 1
		  if (tm[im] > 0) {		// we can remove the data corresponding to polar values.
			m = im*MRES;
			ylm[im]  += tm[im]*(LMAX-m+1);		// shift pointers (still one block for each m)
			if (vector)  dylm[im] += tm[im]*(LMAX-m+1);
			if (zlm[0] != NULL) {
			  for (l=m; l<LMAX; l+=2) {
				for (it=0; it<NLAT_2-tm[im]; it++) {	// copy data to avoid cache misses.
					zlm[im][(l-m)*(NLAT_2-tm[im]) + it*2]   = zlm[im][(l-m)*NLAT_2 + (it+tm[im])*2];
					zlm[im][(l-m)*(NLAT_2-tm[im]) + it*2+1] = zlm[im][(l-m)*NLAT_2 + (it+tm[im])*2+1];
					if (vector) {
						dzlm[im][(l-m)*(NLAT_2-tm[im]) + it*2].t = dzlm[im][(l-m)*NLAT_2 + (it+tm[im])*2].t;
						dzlm[im][(l-m)*(NLAT_2-tm[im]) + it*2].p = dzlm[im][(l-m)*NLAT_2 + (it+tm[im])*2].p;
						dzlm[im][(l-m)*(NLAT_2-tm[im]) + it*2+1].t = dzlm[im][(l-m)*NLAT_2 + (it+tm[im])*2+1].t;
						dzlm[im][(l-m)*(NLAT_2-tm[im]) + it*2+1].p = dzlm[im][(l-m)*NLAT_2 + (it+tm[im])*2+1].p;
					}
				}
			  }
			  if (l==LMAX) {
				for (it=0; it<NLAT_2-tm[im]; it++) {
					zlm[im][(l-m)*(NLAT_2-tm[im]) + it]   = zlm[im][(l-m)*NLAT_2 + (it+tm[im])];
					if (vector) {
						dzlm[im][(l-m)*(NLAT_2-tm[im]) + it].t = dzlm[im][(l-m)*NLAT_2 + (it+tm[im])].t;
						dzlm[im][(l-m)*(NLAT_2-tm[im]) + it].p = dzlm[im][(l-m)*NLAT_2 + (it+tm[im])].p;
					}
				}
			  }
			}
		  }
		}
	}

/// Compression of dzlm for m=0, as .p is 0
	if ((vector) && (dzlm[0] != NULL)) {		// for sht_reg_poles there is no dzlm defined.
		im=0;	m=0;
		double* yg = (double *) dzlm[im];
		for (l=1; l<LMAX; l+=2) {		// l=0 is zero, so we start at l=1.
			for (it=0; it<NLAT_2; it++) {
				yg[(l-1)*NLAT_2 + it*2] = dzlm[im][(l-1)*NLAT_2 + it*2].t;	// l
				yg[(l-1)*NLAT_2 + it*2+1] = dzlm[im][(l-1)*NLAT_2 + it*2+1].t;	// l+1
			}
		}
		if (l==LMAX) {		// last l is stored right away, without interleaving.
			for (it=0; it<NLAT_2; it++) {
				yg[(l-1)*NLAT_2 + it] = dzlm[im][(l-1)*NLAT_2 + it].t;		// l (odd)
			}
		}
	}
	if ((zlm[0] != NULL) && (NLAT_2 & 1)) {
		zlm[0][NLAT_2] = 0.0;		// take care to write 0.0 for this sse2 padding value.
		if ( (MMAX==0) && (!(LMAX & 1)) ) {		// NLAT_2 odd + LMAX even
			zlm[0][(LMAX+1)*NLAT_2 +1] = 0.0;	// overflow for im=0
			zlm[0][(LMAX+1)*NLAT_2 +2] = 0.0;
		}
	}
}

// how to access the ylm matrix for im=0
#define YL0(it, l)  ( shtns->ylm[0][ (it)*(((LMAX+2)>>1)*2) + (l) ] )
#define DYL0(it, l) ( ((double*)(shtns->dylm[0]))[ (it)*(((LMAX+1)>>1)*2) + (l-1) ] )

#define YLM(it, l, im) ( (im==0) ? YL0(it,l) : shtns->ylm[im][(it)*(LMAX-(im)*MRES+1) + ((l)-(im)*MRES)] )
#define DTYLM(it, l, im) ( (im==0) ? ( (l==0) ? 0.0 : DYL0(it,l) ) : shtns->dylm[im][(it)*(LMAX-(im)*MRES+1) + ((l)-(im)*MRES)].t )
#define DPYLM(it, l, im) ( (im==0) ? 0.0 : shtns->dylm[im][(it)*(LMAX-(im)*MRES+1) + ((l)-(im)*MRES)].p )

/// \internal Precompute the matrix for SH synthesis.
static void init_SH_synth(shtns_cfg shtns)
{
	double *yl, *dyl;		// temp storage for derivative for legendre function values.
	long int it,im,m,l;
	double *ct = shtns->ct;
	double *st = shtns->st;
	int vector = (shtns->dylm != NULL);

	yl = (double*) malloc( 2*sizeof(double) *(LMAX+1) );
	dyl = yl + (LMAX+1);

	{	im=0;	m=0;
		for (it=0; it<NLAT_2; it++) {
			legendre_sphPlm_deriv_array_hp(shtns, LMAX, im, ct[it], st[it], yl, dyl);	// fixed im legendre functions lookup table.
			for (l=0; l<=LMAX; l++)   YL0(it, l) = yl[l];
			for (l=LMAX+1; l<=LMAX+3; l++)	YL0(it, l) = 0.0;	// allow overflow.
			if (vector) {
				for (l=1; l<=LMAX; l++)	 DYL0(it, l) = dyl[l];
				for (l=LMAX+1; l<=LMAX+2; l++)  DYL0(it, l) = 0.0;	// allow overflow.
			}
		}
	}
	for (im=1; im<=MMAX; im++) {
		double *ylm = shtns->ylm[im];
		struct DtDp* dylm = vector ? shtns->dylm[im] : NULL;
		m = im*MRES;
		for (it=0; it<NLAT_2; it++) {
			legendre_sphPlm_deriv_array_hp(shtns, LMAX, im, ct[it], st[it], yl, dyl);	// fixed im legendre functions lookup table.
			for (l=m; l<=LMAX; l++) {
				ylm[it*(LMAX-m+1) + (l-m)] = yl[l-m] * st[it];
				if (vector) {
					dylm[it*(LMAX-m+1) + (l-m)].t = dyl[l-m];
					dylm[it*(LMAX-m+1) + (l-m)].p = yl[l-m] *m;	// 1/sint(t) dYlm/dphi
				}
			}
		}
	}
	free(yl);
}


/// \internal Precompute matrices for SH synthesis and analysis, on a Gauss-Legendre grid.
static void init_SH_gauss(shtns_cfg shtns)
{
	long int it,im,m,l;
	int vector = (shtns->dylm != NULL);

	init_SH_synth(shtns);

// for analysis (decomposition, direct transform) : transpose and multiply by gauss weight and other normalizations.
// interleave l and l+1 : this stores data in the way it will be read.
	double **zlm = shtns->zlm;
	struct DtDp** dzlm = shtns->dzlm;
	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
		long int talign = 0;

		for (it=0;it<NLAT_2;it++) {
			double norm = shtns->wg[it];
			if ( (m>0) && (shtns->norm & SHT_REAL_NORM) )	norm *= 2;		// "Real" norm : zlm must be doubled for m>0
			long int l0 = m;
			if (m==0) {
				zlm[im][it] = YL0(it, 0) * norm;
				// les derivees sont nulles pour l=0
				l0++;
				talign = (NLAT_2&1);
			}
			for (l=l0; l<LMAX; l+=2) {
				double nz0 = norm;		double nz1 = norm;
				if (SHT_NORM == sht_schmidt) {
					nz0 *= (2*l+1);		nz1 *= (2*l+3);
				}
				zlm[im][(l-m)*NLAT_2 + it*2 +talign]    =  YLM(it, l, im) * nz0;
				zlm[im][(l-m)*NLAT_2 + it*2 +1 +talign] =  YLM(it, l+1, im) * nz1;
				if (vector) {
					nz0 *= shtns->l_2[l];		nz1 *= shtns->l_2[l+1];
					dzlm[im][(l-l0)*NLAT_2 + it*2].t = DTYLM(it, l, im) * nz0;
					dzlm[im][(l-l0)*NLAT_2 + it*2].p = DPYLM(it, l, im) * nz0;
					dzlm[im][(l-l0)*NLAT_2 + it*2+1].t = DTYLM(it, l+1, im)  * nz1;
					dzlm[im][(l-l0)*NLAT_2 + it*2+1].p = DPYLM(it, l+1, im) * nz1;
				}
			}
			if (l==LMAX) {		// last l is stored right away, without interleaving.
				double nz0 = norm;
				if (SHT_NORM == sht_schmidt)	nz0 *= (2*l+1);
				zlm[im][(l-m)*NLAT_2 + it +talign]    =  YLM(it, l, im) * nz0;
				if (vector) {
					nz0 *= shtns->l_2[l];
					dzlm[im][(l-l0)*NLAT_2 + it].t = DTYLM(it, l, im) * nz0;
					dzlm[im][(l-l0)*NLAT_2 + it].p = DPYLM(it, l, im) * nz0;
				}
			}
		}
	}
}

static void init_SH_dct_m(shtns_cfg shtns, double* is1, fftw_plan dct, fftw_plan idct, int analysis, const int m0, const int mstep)
{
	double *yk, *yk0, *dyk0, *yg;		// temp storage
	struct DtDp *dyg, *dyk;
	int it,im,m,l;
	double Z[2*NLAT_2] SSE;
	double dZt[2*NLAT_2] SSE;
	double dZp[2*NLAT_2] SSE;		// equally spaced theta points.

	double *st = shtns->st;
	double *st_1 = shtns->st_1;
	const int vector = (shtns->dylm != NULL);
	const int KMAX = LMAX+1;

	real iylm_fft_norm = 1.0;	// FFT/SHT normalization for zlm (4pi normalized)
 	if ((SHT_NORM != sht_fourpi)&&(SHT_NORM != sht_schmidt))  iylm_fft_norm = 4*M_PIl;	// FFT/SHT normalization for zlm (orthonormalized)
	iylm_fft_norm /= (2*NPHI*NLAT_2);

// Even/Odd symmetry : ylm is even or odd across equator, as l-m is even or odd => only NLAT_2 points required.
	// temp memory for ykm_dct.
	yk = (double *) malloc( sizeof(double) * (KMAX+1)*(LMAX+1) );
	dyk = (struct DtDp *) malloc( sizeof(struct DtDp)* (KMAX+1)*(LMAX+1) );
	if ((m0==0)&&(analysis)) {
		yk0 = (double *) malloc( sizeof(double) * (LMAX/2+1)*(2*NLAT_2) * 2 );		// temp for zlm_dct0
		dyk0 = yk0 + (LMAX/2+1)*(2*NLAT_2);
	}

	for (im=m0; im<=MMAX; im+=mstep) {
		m = im*MRES;
	// go to DCT space
		for (it=0;it<=KMAX;it+=2) {
			for(l=m; l<=LMAX; l++) {
				yk[(it/2)*(LMAX+1-m) + (l-m)] = 0.0;
				dyk[(it/2)*(LMAX+1-m) + (l-m)].t = 0.0;
				dyk[(it/2)*(LMAX+1-m) + (l-m)].p = 0.0;
			}
		}
		for (l=m; l<=LMAX; l++) {
			if (m & 1) {	// m odd
				for (it=0; it<NLAT_2; it++) {
					Z[it] = YLM(it, l, im) * st[it];	// P[l+1](x)	*st
					if (vector) {
						dZt[it] = DTYLM(it, l, im);	// P[l](x)	*1
						dZp[it] = DPYLM(it, l, im);		// P[l-1](x)	*1
					}
				}
			} else {	// m even
				for (it=0; it<NLAT_2; it++) {
					Z[it] = YLM(it, l, im);				// P[l](x)	*1
					if (vector) {
						dZt[it] = DTYLM(it, l, im) *st[it];	// P[l+1](x)	*st
						dZp[it] = YLM(it, l, im) * m;	// P[l](x)	*st
					}
				}
			}
			if ((l-m)&1) {	// odd
				for (it=NLAT_2; it<2*NLAT_2; it++) {
					Z[it]   = - Z[2*NLAT_2-it-1];	// reconstruct even/odd part
					dZt[it] =   dZt[2*NLAT_2-it-1];
					dZp[it] = - dZp[2*NLAT_2-it-1];
				}
			} else {	// even
				for (it=NLAT_2; it<2*NLAT_2; it++) {
					Z[it] =     Z[2*NLAT_2-it-1];	// reconstruct even/odd part
					dZt[it] = - dZt[2*NLAT_2-it-1];
					dZp[it] =   dZp[2*NLAT_2-it-1];
				}
			}
			fftw_execute_r2r(dct, Z, Z);
			fftw_execute_r2r(dct, dZt, dZt);
			fftw_execute_r2r(dct, dZp, dZp);
#if SHT_VERBOSE > 1
			if (LMAX <= 12)
			#pragma omp critical
			{
				printf("\nl=%d, m=%d ::\t", l,m);
				for(it=0;it<2*NLAT_2;it++) printf("%e ",Z[it]/(2*NLAT));
				printf("\n     dYt ::\t");
				for(it=0;it<2*NLAT_2;it++) printf("%e ",dZt[it]/(2*NLAT));
				printf("\n     dYp ::\t");
				for(it=0;it<2*NLAT_2;it++) printf("%e ",dZp[it]/(2*NLAT));
			}
#endif
			for (it=(l-m)&1; it<=l+1; it+=2) {
				yk[(it/2)*(LMAX+1-m) + (l-m)] = Z[it]/(2*NLAT);	// and transpose
				dyk[(it/2)*(LMAX+1-m) + (l-m)].p = dZp[it]/(2*NLAT);
			}
			for (it=(l+1-m)&1; it<=l+1; it+=2) {
				dyk[(it/2)*(LMAX+1-m) + (l-m)].t = dZt[it]/(2*NLAT);
			}
		}

	/* compute analysis coefficients (fast way)
	 * Wklm = int(Tk*Ylm) = int(Tk.sum(i,a_ilm*Ti)) = sum(i, a_ilm* int(Tk*Ti)) = sum(i, a_ilm*Jik)
	 * with Jik = int(Tk*Ti) = 1/(1-(k-i)^2) + 1/(1-(k+i)^2)
	*/
		if (analysis) {
	#if SHT_VERBOSE > 0
		if ((LMAX>126)&&(m0==0)) {		// only one thread prints this message.
			printf("computing weights m=%d\r",m);	fflush(stdout);
		}
	#endif
		for (l=m; l<=LMAX; l++) {
			unsigned k0,k1, k,i,d;
			double Jik0, Jik1, yy, dyp, dyt;
			double lnorm = iylm_fft_norm;

			if (SHT_NORM == sht_schmidt)	lnorm *= (2*l+1);		// Schmidt semi-normalization
			if ( (m>0) && (shtns->norm & SHT_REAL_NORM) )	lnorm *= 2;		// "real" norm : zlm must be doubled for m>0

			k0 = (l-m)&1;	k1 = 1-k0;
			for(k=0; k<NLAT; k++) {	Z[k] = 0.0;		dZt[k] = 0.0;	dZp[k] = 0.0; }
			for (i=0; i<=(l+1)/2; i++) {
				yy = yk[i*(LMAX+1-m) + (l-m)] * lnorm;
				dyp = dyk[i*(LMAX+1-m) + (l-m)].p * lnorm/(l*(l+1));
				dyt = dyk[i*(LMAX+1-m) + (l-m)].t * lnorm/(l*(l+1));
				if (i+k0==0) {	yy*=0.5;	dyp*=0.5; }
				if (i+k1==0) dyt*=0.5;
				for (k=0; k<NLAT_2; k++) {
					d = (k<i) ? i-k : k-i;
					Jik0 = is1[(i+k0)+k] + is1[d];
					Jik1 = is1[(i+k1)+k] + is1[d];
					Z[2*k+k0] += yy * Jik0;
					dZt[2*k+k1] += dyt * Jik1;
					if (m&1) dZp[2*k+k0] += dyp * Jik0;
				}
			}
#if SHT_VERBOSE > 1
		if (LMAX <= 12)
		#pragma omp critical
		{
			printf("\nl=%d, m=%d ::\t",l,m);
			for (k=0; k<(2*NLAT_2); k++) printf("%f ",Z[k]);
			printf("\n       dZt ::\t");
			for (k=0; k<(2*NLAT_2); k++) printf("%f ",dZt[k]);
			if (m&1) {
				printf("\n       dZp ::\t");
				for (k=0; k<(2*NLAT_2); k++) printf("%f ",dZp[k]);
			}
		}
#endif

			if (m == 0) {		// we store zlm in dct space for m=0
				if (k0==0) 	{
					yk0[((l-m)>>1)*(2*NLAT_2)] = Z[0]*0.5;         // store zlm_dct (k=0)
					for (k=1; k<(2*NLAT_2); k++) yk0[((l-m)>>1)*(2*NLAT_2) +k] = 0.0;		// zero out.
					k0=2;
				}
				for (k=k0; k<(2*NLAT_2); k+=2)
					yk0[((l-m)>>1)*(2*NLAT_2) +k] = Z[k];             // store zlm_dct
					
				if (l>0) {
					if (k1==0) 	{
						dyk0[((l-1-m)>>1)*(2*NLAT_2)] = dZt[0]*0.5;         // store dzlm_dct (k=0)
						for (k=1; k<(2*NLAT_2); k++) dyk0[((l-1-m)>>1)*(2*NLAT_2) +k] = 0.0;		// zero out.
						k1=2;
					}
					for (k=k1; k<(2*NLAT_2); k+=2)
						dyk0[((l-1-m)>>1)*(2*NLAT_2) +k] = dZt[k];             // store dzlm_dct
				}
			}

			fftw_execute_r2r(idct, Z, Z);	fftw_execute_r2r(idct, dZt, dZt);
			if (m == 0) {
				for (it=0; it<NLAT; it++) { dZp[it] = 0.0; 	dZt[it] *= st_1[it]; }
			} else if (m & 1) {	// m odd
				fftw_execute_r2r(idct, dZp, dZp);
				for (it=0; it<NLAT; it++) {	Z[it] *= st_1[it]; }
			} else {	// m even
				for (it=0; it<NLAT; it++) { dZp[it] = Z[it]*m/(l*(l+1)*st[it]); 	dZt[it] *= st_1[it]; }
			}

			long int l0 = (m==0) ? 1 : m;
			long int talign = (m==0)*(NLAT_2 & 1);
			long int sk = (l-l0)&1;
			if (l==0) {
				for (it=0; it<NLAT_2; it++) {
					shtns->zlm[im][it] =  Z[it];
				}
			} else if ((sk == 0)&&(l == LMAX)) {
				for (it=0; it<NLAT_2; it++) {
					shtns->zlm[im][(l-m)*NLAT_2 + it + talign] =  Z[it];
					if (vector) {
						shtns->dzlm[im][(l-l0)*NLAT_2 + it].p = dZp[it];
						shtns->dzlm[im][(l-l0)*NLAT_2 + it].t = dZt[it];
					}
				}
			} else {
				for (it=0; it<NLAT_2; it++) {
					shtns->zlm[im][(l-m-sk)*NLAT_2 + it*2 +sk + talign] = Z[it];
					if (vector) {
						shtns->dzlm[im][(l-l0-sk)*NLAT_2 + it*2 +sk].p = dZp[it];
						shtns->dzlm[im][(l-l0-sk)*NLAT_2 + it*2 +sk].t = dZt[it];
					}
				}
			}
		}
		}

		// Compact the coefficients for improved cache efficiency.
		yg = shtns->ykm_dct[im];
		for (it=0; it<= KMAX; it+=2) {
			l = (it < m) ? m : it-(m&1);
			while (l<=LMAX) {
				yg[0] = yk[(it/2)*(LMAX+1-m) + (l-m)];
				l++;	yg++;
			}
			if ((m==0) && ((LMAX & 1) == 0)) {	yg[0] = 0;		yg++;	}		// SSE2 padding.
		}
		if (MMAX==0) for (int j=0; j<3; j++) yg[j] = 0;		// allow some overflow.
		if (vector) {
			dyg = shtns->dykm_dct[im];
			for (it=0; it<= KMAX; it+=2) {
				l = (it-2 < m) ? m : it-2+(m&1);
				while (l<=LMAX) {
					dyg[0].t = dyk[(it/2)*(LMAX+1-m) + (l-m)].t;
					dyg[0].p = dyk[(it/2)*(LMAX+1-m) + (l-m)].p;
					l++;	dyg++;
				}
			}
			if (im == 0) {		// compact and reorder m=0 dylm because .p = 0 :
				dyg = shtns->dykm_dct[im];
				yg = (double *) shtns->dykm_dct[im];
				for (it=0; it<= KMAX; it+=2) {
					dyg++;
					for (l=it-1; l<=LMAX; l++) {
						if (l>0) {
							yg[0] = dyg[0].t;
							yg++;	dyg++;
						}
					}
					if (LMAX&1) {		// padding for SSE2 alignement.
						yg[0] = 0.0;	yg++;
					}
				}
			}
		}
	}
	
	// compact yk to zlm_dct0
	if ((m0==0)&&(analysis)) {
		long int klim = (LMAX * SHT_NL_ORDER) + 2;		// max k needed for nl-terms...
		klim = (klim/2)*2;		// must be even...
		if (klim > 2*NLAT_2) klim = 2*NLAT_2;		// but no more than 2*NLAT_2.
		shtns->klim = klim;		// store for use in codelets.
		yg = shtns->zlm_dct0;
		for (l=0; l<=LMAX; l+=2) {
			for (it=l; it<klim; it++) {	// for m=0, zl coeff with i<l are zeros.
				*yg = yk0[it];
				yg++;
			}
			yk0 += 2*NLAT_2;
		}
		if (vector) {
			yg = shtns->dzlm_dct0;
			for (l=1; l<=LMAX; l+=2) {
				for (it=l-1; it<klim; it++) {	// for m=0, dzl coeff with i<l-1 are zeros.
					*yg = dyk0[it];
					yg++;
				}
				dyk0 += 2*NLAT_2;
			}
		}
		free(yk0 - (2*NLAT_2)*(LMAX/2+1));
	}

	free(dyk);	free(yk);
}

/// \internal Computes the matrices required for SH transform on a regular grid (with or without DCT).
/// \param analysis : 0 => synthesis only.
static void init_SH_dct(shtns_cfg shtns, int analysis)
{
	fftw_plan dct, idct;
	long int it,im,m,l;
	long int sk, dsk;
	double Z[2*NLAT_2] SSE;
	double is1[NLAT];		// tabulate values for integrals.

	const long int marray_size = sizeof(void*)*(MMAX+1) + (MIN_ALIGNMENT-1);
	const int vector = (shtns->dylm != NULL);
	const int KMAX = LMAX+1;

	for(im=0, sk=0, dsk=0; im<=MMAX; im++) {	// how much memory to allocate for ykm_dct ?
		m = im*MRES;
		for (it=0; it<= KMAX; it+=2) {
			l = (it < m) ? m : it-(m&1);
			sk += LMAX+1 - l;
			if ((m==0) && ((LMAX & 1) ==0)) sk++;		// SSE padding for m=0
		}
		for (it=0; it<= KMAX; it+=2) {
			l = (it-2 < m) ? m : it-2+(m&1);
			dsk += LMAX+1 - l;
		}
	}
	if (MMAX == 0) sk+=3;	// allow some overflow.
	for (l=0, it=0; l<=LMAX; l+=2)	// how much memory for zlm_dct0 ?
		it += (2*NLAT_2 -l);
	for (l=1, im=0; l<=LMAX; l+=2)	// how much memory for dzlm_dct0 ?
		im += (2*NLAT_2 -l+1);

#if SHT_VERBOSE > 1
	printf("          Memory used for Ykm_dct matrices = %.3f Mb\n",sizeof(double)*(sk + 2.*dsk + it)/(1024.*1024.));
#endif
	shtns->ykm_dct = (double **) malloc( marray_size + sizeof(double)*sk );
	shtns->ykm_dct[0] = (double *) PTR_ALIGN( shtns->ykm_dct + (MMAX+1) );
	if (vector) {
		shtns->dykm_dct = (struct DtDp **) malloc( marray_size + sizeof(struct DtDp)*dsk );
		shtns->dykm_dct[0] = (struct DtDp *) PTR_ALIGN( shtns->dykm_dct + (MMAX+1) );
	}
	shtns->zlm_dct0 = (double *) VMALLOC( sizeof(double)* it );
	if (vector) {
		shtns->dzlm_dct0 = (double *) VMALLOC( sizeof(double)* im );
	}
	for (im=0; im<MMAX; im++) {
		m = im*MRES;
		for (it=0, sk=0; it<= KMAX; it+=2) {
			l = (it < m) ? m : it-(m&1);
			sk += LMAX+1 - l;
			if ((m==0) && ((LMAX & 1) ==0)) sk++;		// SSE padding for m=0
		}
		for (it=0, dsk=0; it<= KMAX; it+=2) {
			l = (it-2 < m) ? m : it-2+(m&1);
			dsk += LMAX+1 - l;
		}
		shtns->ykm_dct[im+1] = shtns->ykm_dct[im] + sk;
		if (vector)  shtns->dykm_dct[im+1] = shtns->dykm_dct[im] + dsk;
	}

#if SHT_VERBOSE > 1
	ticks tik0, tik1;
	tik0 = getticks();
#endif

	#ifdef OMP_FFTW
		fftw_plan_with_nthreads(1);
	#endif
	dct = fftw_plan_r2r_1d( 2*NLAT_2, Z, Z, FFTW_REDFT10, FFTW_MEASURE );	// quick and dirty dct.
	idct = fftw_plan_r2r_1d( 2*NLAT_2, Z, Z, FFTW_REDFT01, FFTW_MEASURE );	// quick and dirty idct.

	init_SH_synth(shtns);

// precomputation for scalar product of Chebychev polynomials.
	for(it=0; it<NLAT; it++)
		is1[it] = 1./(1. - 4.*it*it);

  #ifdef _OPENMP
	#pragma omp parallel
	{
		int n=omp_get_num_threads();
		int k=omp_get_thread_num();
		init_SH_dct_m(shtns, is1, dct, idct, analysis, k, n);
	}
  #else
	init_SH_dct_m(shtns, is1, dct, idct, analysis, 0, 1);
  #endif

#if SHT_VERBOSE > 1
	tik1 = getticks();
	printf("\n    ticks : %.3f\n", elapsed(tik1,tik0)/(NLM*NLAT*(MMAX+1)));
#endif
	fftw_destroy_plan(idct);	fftw_destroy_plan(dct);
}


/// \internal Generate a gauss grid (including weights)
static void grid_gauss(shtns_cfg shtns, double latdir)
{
	long int it;
	real iylm_fft_norm;
	real xg[NLAT], wgl[NLAT];	// gauss points and weights.
	const int overflow = 8*VSIZE2-1;

	shtns->grid = GRID_GAUSS;
	shtns->wg = VMALLOC((NLAT_2 +overflow) * sizeof(double));	// gauss weights, double precision.

	iylm_fft_norm = 1.0;	// FFT/SHT normalization for zlm (4pi normalized)
	if ((SHT_NORM != sht_fourpi)&&(SHT_NORM != sht_schmidt))  iylm_fft_norm = 4*M_PIl;	// FFT/SHT normalization for zlm (orthonormalized)
	iylm_fft_norm /= (2*NPHI);
#if SHT_VERBOSE > 0
	printf("        => using Gauss nodes\n");
	if (2*NLAT <= (SHT_NL_ORDER +1)*LMAX) printf("     !! Warning : Gauss-Legendre anti-aliasing condition 2*Nlat > %d*Lmax is not met.\n",SHT_NL_ORDER+1);
#endif
	gauss_nodes(xg,wgl,NLAT);	// generate gauss nodes and weights : ct = ]1,-1[ = cos(theta)
	for (it=0; it<NLAT; it++) {
		shtns->ct[it] = latdir * xg[it];
		shtns->st[it] = SQRT((1.-xg[it])*(1.+xg[it]));
		shtns->st_1[it] = 1.0/SQRT((1.-xg[it])*(1.+xg[it]));
	}
	for (it=0; it<NLAT_2; it++)
		shtns->wg[it] = wgl[it]*iylm_fft_norm;		// faster double-precision computations.
	if (NLAT & 1) {		// odd NLAT : adjust weigth of middle point.
		shtns->wg[NLAT_2-1] *= 0.5;
	}
	for (it=NLAT_2; it < NLAT_2 +overflow; it++) shtns->wg[it] = 0.0;		// padding for multi-way algorithm.

#if SHT_VERBOSE > 1
	printf(" NLAT=%d, NLAT_2=%d\n",NLAT,NLAT_2);
// TEST if gauss points are ok.
	double tmax = 0.0;
	for (it = 0; it<NLAT_2; it++) {
		double t = legendre_Pl(NLAT, shtns->ct[it]);
		if (t>tmax) tmax = t;
//		printf("i=%d, x=%12.12g, p=%12.12g\n",it,ct[it],t);
	}
	printf("          max zero at Gauss nodes for Pl[l=NLAT] : %g\n",tmax);
	if (NLAT_2 < 100) {
		printf("          Gauss nodes :");
		for (it=0;it<NLAT_2; it++)
			printf(" %g",shtns->ct[it]);
		printf("\n");
	}
#endif
}

/// \internal Generate an equi-spaced theta grid (Chebychev points, excluding poles) for Féjer-DCT SHT.
static void grid_dct(shtns_cfg shtns, double latdir)
{
	long int it;

	shtns->grid = GRID_REGULAR;
#if SHT_VERBOSE > 0
	printf("        => using equaly spaced nodes with DCT acceleration\n");
	if (NLAT <= SHT_NL_ORDER *LMAX)	printf("     !! Warning : DCT anti-aliasing condition Nlat > %d*Lmax is not met.\n",SHT_NL_ORDER);
	if (NLAT != fft_int(NLAT,7))	printf("     !! Warning : Nlat is not optimal for FFTW !\n");
#endif
	if (NLAT & 1) shtns_runerr("NLAT must be even (DCT)");
	if (NLAT <= LMAX+1) shtns_runerr("NLAT should be at least LMAX+2 (DCT)");

	for (it=0; it<NLAT; it++) {		// Chebychev points : equaly spaced but skipping poles.
		real th = M_PIl;	th = (th*(2*it+1))/(2*NLAT);
		shtns->ct[it] = latdir * COS(th);
		shtns->st[it] = SIN(th);
		shtns->st_1[it] = 1.0/SIN(th);
	}

#if SHT_VERBOSE > 1
	printf(" NLAT=%d, NLAT_2=%d\n",NLAT,NLAT_2);
	double tmax = 0.0;
	for (it=0;it<NLAT_2; it++) {
		double ct = shtns->ct[it];		double st = shtns->st[it];
		double t = fabs((ct*ct + st*st) -1.0);
		if (t > tmax) tmax=t;
	}
	printf(" max st^2 + ct^2 -1 = %g\n",tmax);
	if (NLAT_2 < 100) {
		printf("          DCT nodes :");
		for (it=0; it<NLAT_2; it++)
			printf(" %g",shtns->ct[it]);
		printf("\n");
	}
#endif
}

/// \internal Generate an equi-spaced theta grid including the poles, for synthesis only.
static void grid_equal_polar(shtns_cfg shtns, double latdir)
{
	long int j;

	shtns->grid = GRID_POLES;
#if SHT_VERBOSE > 0
	printf("        => using Equaly Spaced Nodes including poles\n");
#endif
// cos theta of latidunal points (equaly spaced in theta)
	double f = M_PIl/(NLAT-1);
	for (j=0; j<NLAT; j++) {
		shtns->ct[j] = latdir * cos(f*j);
		shtns->st[j] = sin(f*j);
		shtns->st_1[j] = 1.0/sin(f*j);
	}
#if SHT_VERBOSE > 0
	printf("     !! Warning : only synthesis (inverse transform) supported for this grid !\n");
#endif
}


/* TEST AND TIMING FUNCTIONS */

/// \internal return the max error for a back-and-forth SHT transform.
/// this function is used to internally measure the accuracy.
double SHT_error(shtns_cfg shtns, int vector)
{
	complex double *Tlm0, *Slm0, *Tlm, *Slm;
	double *Sh, *Th;
	double t, tmax, n2,  err;
	long int i, jj, nlm_cplx;
	
	srand( time(NULL) );	// init random numbers.
	
	Tlm0 = (complex double *) VMALLOC(sizeof(complex double)* NLM);
	Slm0 = (complex double *) VMALLOC(sizeof(complex double)* NLM);
	Slm = (complex double *) VMALLOC(sizeof(complex double)* NLM);
	Tlm = (complex double *) VMALLOC(sizeof(complex double)* NLM);
	Sh = (double *) VMALLOC( NSPAT_ALLOC(shtns) * sizeof(double) );
	Th = (double *) VMALLOC( NSPAT_ALLOC(shtns) * sizeof(double) );

// m = nphi/2 is also real if nphi is even.
	nlm_cplx = ( MMAX*2 == NPHI ) ? LiM(shtns, MRES*MMAX,MMAX) : NLM;
	t = 1.0 / (RAND_MAX/2);
	for (i=0; i<NLM; i++) {
		if ((i<=LMAX)||(i>=nlm_cplx)) {		// m=0 or m*2=nphi : real random data
			Slm0[i] = t*((double) (rand() - RAND_MAX/2));
			Tlm0[i] = t*((double) (rand() - RAND_MAX/2));
		} else {							// m>0 : complex random data
			Slm0[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
			Tlm0[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
		}
	}

	SH_to_spat(shtns, Slm0,Sh);		// scalar SHT
	spat_to_SH(shtns, Sh, Slm);
	for (i=0, tmax=0., n2=0., jj=0; i<NLM; i++) {		// compute error
		t = cabs(Slm[i] - Slm0[i]);
		n2 += t*t;
		if (t>tmax) { tmax = t; jj = i; }
	}
	err = tmax;
#if SHT_VERBOSE > 1
	printf("        scalar SH - poloidal   rms error = %.3g  max error = %.3g for l=%hu,lm=%d\n",sqrt(n2/NLM),tmax,shtns->li[jj],jj);
#endif

	if (vector) {
		Slm0[0] = 0.0; 	Tlm0[0] = 0.0;		// l=0, m=0 n'a pas de signification sph/tor
		SHsphtor_to_spat(shtns, Slm0, Tlm0, Sh, Th);		// vector SHT
		spat_to_SHsphtor(shtns, Sh, Th, Slm, Tlm);
		for (i=0, tmax=0., n2=0., jj=0; i<NLM; i++) {		// compute error
			t = cabs(Slm[i] - Slm0[i]);
			n2 += t*t;
			if (t>tmax) { tmax = t; jj = i; }
		}
		if (tmax > err) err = tmax;
	#if SHT_VERBOSE > 1
		printf("        vector SH - spheroidal rms error = %.3g  max error = %.3g for l=%hu,lm=%d\n",sqrt(n2/NLM),tmax,shtns->li[jj],jj);
	#endif
		for (i=0, tmax=0., n2=0., jj=0; i<NLM; i++) {		// compute error
			t = cabs(Tlm[i] - Tlm0[i]);
			n2 += t*t;
			if (t>tmax) { tmax = t; jj = i; }
		}
		if (tmax > err) err = tmax;
	#if SHT_VERBOSE > 1
		printf("                  - toroidal   rms error = %.3g  max error = %.3g for l=%hu,lm=%d\n",sqrt(n2/NLM),tmax,shtns->li[jj],jj);
	#endif
	}

	VFREE(Th);  VFREE(Sh);  VFREE(Tlm);  VFREE(Slm);  VFREE(Slm0);  VFREE(Tlm0);
	return(err);		// return max error.
}


#if SHT_VERBOSE == 1
  #define PRINT_DOT 	{	printf(".");	fflush(stdout);	}
#else
  #define PRINT_DOT (0);
#endif

/// \internal measure time used for a transform function
static double get_time(shtns_cfg shtns, int nloop, int npar, char* name, void *fptr, void *i1, void *i2, void *i3, void *o1, void *o2, void *o3, int l)
{
	double t;
	int i;
	ticks tik0, tik1;

	if (fptr == NULL) return(0.0);

	tik1 = getticks();
	for (i=0; i<nloop; i++) {
		switch(npar) {
			case 2: (*(pf2l)fptr)(shtns, i1,o1, l); break;			// l may be discarded.
			case 3: (*(pf3l)fptr)(shtns, i1,o1,o2, l); break;
			case 4: (*(pf4l)fptr)(shtns, i1,i2,o1,o2, l); break;
			default: (*(pf6l)fptr)(shtns, i1,i2,i3, o1,o2,o3, l); break;
		}
		if (i==0) tik0 = getticks();
	}
	if (nloop == 1) {
		t = elapsed(tik0, tik1);
	} else {
		tik1 = getticks();
		t = elapsed(tik1, tik0)/(nloop-1);		// discard first iteration.
	}
	#if SHT_VERBOSE > 1
		printf("  t(%s) = %.3g",name,t);	fflush(stdout);
	#endif
	return t;
}


/// \internal choose fastest between on-the-fly and gauss algorithms.
/// *nlp is the number of loops. If zero, it is set to a good value.
/// on_the_fly : 1 = skip all memory algorithm. 0 = include memory and on-the-fly. -1 = test only DCT.
/// returns time without dct / best time with dct (or 0 if no dct available).
static double choose_best_sht(shtns_cfg shtns, int* nlp, int vector, int dct_mtr)
{
	complex double *Qlm, *Slm, *Tlm;
	double *Qh, *Sh, *Th;
	int m, i, i0, minc, nloop, alg_end;
	int typ_lim = SHT_NTYP;		// time every type.
	double t0, t, tt, r;
	double tdct, tnodct;
	clock_t tcpu;
	int on_the_fly_only = (shtns->ylm == NULL);		// only on-the-fly.
	int otf_analys = (shtns->wg != NULL);			// on-the-fly analysis supported.

	if (NLAT < VSIZE2*4) return(0.0);			// on-the-fly not possible for NLAT_2 < 2*NWAY (overflow) and DCT not efficient for low NLAT.
	if ((dct_mtr != 0) && (shtns->ykm_dct == NULL)) return(0.0);		// no dct available : do nothing.

	m = 2*(NPHI/2+1) * NLAT * sizeof(double);
	i = sizeof(complex double)* NLM;
	if (i>m) m=i;
	Qh = (double *) VMALLOC(m);		Sh = (double *) VMALLOC(m);		Th = (double *) VMALLOC(m);
	Qlm = (complex double *) VMALLOC(m);	Slm = (complex double *) VMALLOC(m);	Tlm = (complex double *) VMALLOC(m);

	for (i=0;i<NLM;i++) {
		int l = shtns->li[i];
		Slm[i] = shtns->l_2[l] + 0.5*I*shtns->l_2[l];
		Tlm[i] = 0.5*shtns->l_2[l] + I*shtns->l_2[l];
	}

	#if SHT_VERBOSE > 0
		if (dct_mtr != 0)	printf("        finding optimal m-truncation for DCT synthesis");
		else 	printf("        finding optimal algorithm");
		fflush(stdout);
	#endif

	if (*nlp <= 0) {
		// find good nloop by requiring less than 3% difference between 2 consecutive timings.
		m=0;	nloop = 1;                     // number of loops to get timings.
		r = 0.0;	tt = 1.0;
		do {
			if ((r > 0.03)||(tt<0.1)) {
				m = 0;		nloop *= 3;
			} else 	m++;
			tcpu = clock();
			t0 = get_time(shtns, nloop, 2, "", sht_func[SHT_STD][SHT_TYP_SSY][SHT_FLY2], Slm, Tlm, Qlm, Sh, Th, Qh, LMAX);
			tcpu = clock() - tcpu;		tt = 1.e-6 * tcpu;
			if (tt >= SHT_TIME_LIMIT) break;			// we should not exceed 1 second
			t = get_time(shtns, nloop, 2, "", sht_func[SHT_STD][SHT_TYP_SSY][SHT_FLY2], Slm, Tlm, Qlm, Sh, Th, Qh, LMAX);
			r = fabs(2.0*(t-t0)/(t+t0));
			#if SHT_VERBOSE > 1
				printf(", nloop=%d, r=%g, m=%d (real time = %g s)\n",nloop,r,m,tt);
				if (tt >= 0.01) break;		// faster timing in debug mode.
			#endif
			PRINT_DOT
		} while((nloop<10000)&&(m < 3));
		*nlp = nloop;
	} else {
		nloop = *nlp;
	}
	#if SHT_VERBOSE > 1
		printf(" => nloop=%d (takes %g s)\n",nloop, tt);
	#endif
	if (vector == 0)	typ_lim = SHT_TYP_VSY;		// time only scalar transforms.
//	if (tt > 3.0)		typ_lim = SHT_TYP_VSY;		// time only scalar transforms.
//	if (tt > 10.0)	goto done;		// timing this will be too slow...

	int ityp = 0;	do {
		if ((dct_mtr != 0) && (ityp >= 4)) break;		// dct !=0 : only scalar and vector.
		if (ityp == 2) nloop = (nloop+1)/2;		// scalar ar done.
		t0 = 1e100;
		i0 = 0;
		if (MTR_DCT < 0)     i0 = SHT_MEM;		// skip dct.
		if (on_the_fly_only) i0 = SHT_SV;		// only on-the-fly (SV is then also on-the-fly)
		alg_end = SHT_NALG;
		if (shtns->nthreads <= 1) alg_end = SHT_OMP1;		// no OpenMP with 1 thread.
		if ((ityp&1) && (otf_analys == 0)) alg_end = SHT_FLY1;		// no on-the-fly analysis for regular grid.
		for (i=i0, m=0;	i<alg_end; i++) {
			if (sht_func[0][ityp][i] != NULL) m++;		// count number of algos
		}
		if (m >= 2) {		// don't time if there is only 1 algo !
			#if SHT_VERBOSE > 1
				printf("finding best %s ...",sht_type[ityp]);	fflush(stdout);
			#endif
			i = i0-1;		i0 = -1;
			while (++i < alg_end) {
				void *pf = sht_func[0][ityp][i];
				if (pf != NULL) {
					if (ityp&1) {	// analysis
						t = get_time(shtns, nloop, sht_npar[ityp], sht_name[i], pf, Sh, Th, Qh, Slm, Tlm, Qlm, LMAX);
					} else {
						t = get_time(shtns, nloop, sht_npar[ityp], sht_name[i], pf, Slm, Tlm, Qlm, Sh, Th, Qh, LMAX);
					}
					if (i < SHT_FLY1) t *= 1.03;	// 3% penality for memory based transforms.
				#ifdef _OPENMP
					if ((i >= SHT_OMP1)||(i == SHT_SV)) t *= 1.3;	// 30% penality for openmp transforms.
				#endif
					if (t < t0) {	i0 = i;		t0 = t;		PRINT_VERB("*");	}
				}
			}
			if (i0 >= 0) {
				for (int j=0; j<SHT_NVAR; j++) {
					shtns->fptr[j][ityp] = sht_func[j][ityp][i0];
					if (ityp == 4) shtns->fptr[j][ityp+1] = sht_func[j][ityp+1][i0];		// only one timing for both gradients variants.
				}
				PRINT_DOT
				#if SHT_VERBOSE > 1
					printf(" => %s\n",sht_name[i0]);
				#endif
			}
		}
		if (ityp == 4) ityp++;		// skip second gradient
	} while(++ityp < typ_lim);

	if (dct_mtr != 0) {		// find the best DCT timings...
 		#if SHT_VERBOSE > 1
			printf("finding best mtr_dct ...");	fflush(stdout);
		#endif
		minc = MMAX/20 + 1;             // don't test every single m.
		m = -1;		i = -1;		t0 = 0.0;		// reference = no dct.
		if (sht_func[SHT_STD][SHT_TYP_SSY][SHT_DCT] != NULL)
			t0 += get_time(shtns, *nlp, 2, "s", shtns->fptr[SHT_STD][SHT_TYP_SSY], Qlm, Slm, Tlm, Qh, Sh, Th, LMAX);
		if ( (sht_func[SHT_STD][SHT_TYP_VSY][SHT_DCT] != NULL) && (vector) )
			t0 += get_time(shtns, nloop, 4, "v", shtns->fptr[SHT_STD][SHT_TYP_VSY], Slm, Tlm, Qlm, Sh, Th, Qh, LMAX);
		tnodct = t0;
		for (m=0; m<=MMAX; m+=minc) {
			#if SHT_VERBOSE > 1
				printf("\n\tm=%d :",m);
			#endif
			Set_MTR_DCT(shtns, m);
			t = get_time(shtns, *nlp, 2, "sdct", sht_func[SHT_STD][SHT_TYP_SSY][SHT_DCT], Qlm, Slm, Tlm, Qh, Sh, Th, LMAX);
			if (vector)
				t += get_time(shtns, nloop, 4, "vdct", sht_func[SHT_STD][SHT_TYP_VSY][SHT_DCT], Slm, Tlm, Qlm, Sh, Th, Qh, LMAX);
			if (t < t0) {	t0 = t;		i = m;	PRINT_VERB("*"); }
			PRINT_DOT
		}
		tdct = t0;
		Set_MTR_DCT(shtns, i);		// the best DCT is chosen.
		#if SHT_VERBOSE > 0
			printf(" mtr_dct=%d  (%.1f%% performance gain)", MTR_DCT*MRES, 100.*(tnodct/tdct-1.));
		#endif
	}

done:
	#if SHT_VERBOSE > 0
		printf("\n");
	#endif
	VFREE(Tlm);	VFREE(Slm);	VFREE(Qlm);	VFREE(Th);	VFREE(Sh);	VFREE(Qh);
	if (dct_mtr > 0) {
		return(tnodct/tdct);
	} else	return(0.0);
}



#ifndef _HGID_
  #define _HGID_ "unknown"
#endif

void shtns_print_version() {
	printf("[" PACKAGE_STRING "] build " __DATE__ ", " __TIME__ ", id: " _HGID_ "\n");
}

void shtns_print_cfg(shtns_cfg shtns)
{
	printf("Lmax=%d, Mmax*Mres=%d, Mres=%d, Nlm=%d  [%d threads, ",LMAX, MMAX*MRES, MRES, NLM, shtns->nthreads);
	if (shtns->norm & SHT_REAL_NORM) printf("'real' norm, ");
	if (shtns->norm & SHT_NO_CS_PHASE) printf("no Condon-Shortley phase, ");
	if (SHT_NORM == sht_fourpi) printf("4.pi normalized]\n");
	else if (SHT_NORM == sht_schmidt) printf("Schmidt semi-normalized]\n");
	else printf("orthonormalized]\n");
	if (shtns->ct == NULL)	return;		// no grid is set

	switch(shtns->grid) {
		case GRID_GAUSS : printf("Gauss grid");	 break;
		case GRID_REGULAR : printf("Regular grid");	 break;
		case GRID_POLES : printf("Regular grid including poles");  break;
		default : printf("Unknown grid");
	}
	printf(" : Nlat=%d, Nphi=%d\n", NLAT, NPHI);
	printf("      ");
	for (int it=0; it<SHT_NTYP; it++)
		printf("%5s ",sht_type[it]);
	for (int iv=0; iv<SHT_NVAR; iv++) {
		printf("\n  %4s:",sht_var[iv]);
		for (int it=0; it<SHT_NTYP; it++) {
			if (shtns->fptr[iv][it] != NULL) {
				for (int ia=0; ia<SHT_NALG; ia++)
					if (sht_func[iv][it][ia] == shtns->fptr[iv][it]) {
						printf("%5s ",sht_name[ia]);	break;
					}
			} else  printf(" none ");
		}
	}
	printf("\n");
}

/// \internal returns 1 if val cannot fit in dest (unsigned)
#define IS_TOO_LARGE(val, dest) (sizeof(dest) >= sizeof(val)) ? 0 : ( ( val >= (1<<(8*sizeof(dest))) ) ? 1 : 0 )

/// \internal returns the size that must be allocated for an shtns_info.
#define SIZEOF_SHTNS_INFO(mmax) ( sizeof(struct shtns_info) + (mmax+1)*( sizeof(int)+sizeof(unsigned short) ) )

/* PUBLIC INITIALIZATION & DESTRUCTION */

/** \addtogroup init Initialization functions.
*/
//@{

/*! This sets the description of spherical harmonic coefficients.
 * It tells SHTns how to interpret spherical harmonic coefficient arrays, and it sets usefull arrays.
 * Returns the configuration to be passed to subsequent transform functions, which is basicaly a pointer to a \ref shtns_info struct.
 * \param lmax : maximum SH degree that we want to describe.
 * \param mmax : number of azimutal wave numbers.
 * \param mres : \c 2.pi/mres is the azimutal periodicity. \c mmax*mres is the maximum SH order.
 * \param norm : define the normalization of the spherical harmonics (\ref shtns_norm)
 * + optionaly disable Condon-Shortley phase (ex: \ref sht_schmidt | \ref SHT_NO_CS_PHASE)
 * + optionaly use a 'real' normalization (ex: \ref sht_fourpi | \ref SHT_REAL_NORM)
*/
shtns_cfg shtns_create(int lmax, int mmax, int mres, enum shtns_norm norm)
{
	shtns_cfg shtns, s2;
	int im, m, l, lm;
	int with_cs_phase = 1;		/// Condon-Shortley phase (-1)^m is used by default.
	double mpos_renorm = 1.0;	/// renormalization of m>0.
	int larrays_ok = 0;
	int legendre_ok = 0;
	int l_2_ok = 0;

//	if (lmax < 1) shtns_runerr("lmax must be larger than 1");
	if (lmax < 2) shtns_runerr("lmax must be at least 2");
	if (IS_TOO_LARGE(lmax, shtns->lmax)) shtns_runerr("lmax too large");
	if (mmax*mres > lmax) shtns_runerr("MMAX*MRES should not exceed LMAX");
	if (mres <= 0) shtns_runerr("MRES must be > 0");

	// allocate new setup and initialize some variables (used as flags) :
	shtns = malloc( SIZEOF_SHTNS_INFO(mmax) );
	if (shtns == NULL) return shtns;	// FAIL
	{
		void **p0 = (void**) &shtns->tm;	// first pointer in struct.
		void **p1 = (void**) &shtns->Y00_1;	// first non-pointer.
		while(p0 < p1)	 *p0++ = NULL;		// write NULL to every pointer.
		shtns->lmidx = (int*) (shtns + 1);		// lmidx is stored at the end of the struct...
		shtns->tm = (unsigned short*) (shtns->lmidx + (mmax+1));		// and tm just after.
		shtns->ct = NULL;	shtns->st = NULL;
		shtns->nphi = 0;	shtns->nlat = 0;	shtns->nlat_2 = 0;		shtns->nspat = 0;	// public data
	}

	// copy sizes.
	shtns->norm = norm;
	if (norm & SHT_NO_CS_PHASE)
		with_cs_phase = 0;
	if (norm & SHT_REAL_NORM)
		mpos_renorm = 0.5;		// normalization for 'real' spherical harmonics.

	shtns->mmax = mmax;		shtns->mres = mres;		shtns->lmax = lmax;
	shtns->nlm = nlm_calc(lmax, mmax, mres);
	shtns->nthreads = omp_threads;
	if (omp_threads > mmax+1) shtns->nthreads = mmax+1;	// limit the number of threads to mmax+1
	#if SHT_VERBOSE > 0
		shtns_print_version();
		printf("        ");		shtns_print_cfg(shtns);
	#endif

	s2 = sht_data;		// check if some data can be shared ...
	while(s2 != NULL) {
		if ((s2->mmax >= mmax) && (s2->mres == mres)) {
			if (s2->lmax == lmax) {		// we can reuse the l-related arrays (li + copy lmidx)
				shtns->li = s2->li;		shtns->mi = s2->mi;
				for (im=0; im<=mmax; im++)	shtns->lmidx[im] = s2->lmidx[im];
				larrays_ok = 1;
			}
			if ( (s2->lmax >= lmax) && (s2->norm == norm) ) {		// we can reuse the legendre tables.
				shtns->alm = s2->alm;		shtns->blm = s2->blm;
				legendre_ok = 1;
			}
		}
		if (s2->lmax >= lmax) {		// we can reuse l_2
			shtns->l_2 = s2->l_2;
			l_2_ok = 1;
		}
		s2 = s2->next;
	}
	if (larrays_ok == 0) {
		// alloc spectral arrays
		shtns->li = (unsigned short *) malloc( 2*NLM*sizeof(unsigned short) );	// NLM defined at runtime.
		shtns->mi = shtns->li + NLM;
		for (im=0, lm=0; im<=MMAX; im++) {	// init l-related arrays.
			m = im*MRES;
			shtns->lmidx[im] = lm -m;		// virtual pointer for l=0
			for (l=im*MRES;l<=LMAX;l++) {
				shtns->li[lm] = l;		shtns->mi[lm] = m;
				lm++;
			}
		}
		if (lm != NLM) shtns_runerr("unexpected error");
	}
	if (legendre_ok == 0) {	// this quickly precomputes some values for the legendre recursion.
		legendre_precomp(shtns, SHT_NORM, with_cs_phase, mpos_renorm);
	}
	if (l_2_ok == 0) {
		shtns->l_2 = (double *) malloc( (LMAX+1)*sizeof(double) );
		shtns->l_2[0] = 0.0;	// undefined for l=0 => replace with 0.
		real one = 1.0;
		for (l=1; l<=LMAX; l++)		shtns->l_2[l] = one/(l*(l+1));
	}

	switch(SHT_NORM) {
		case sht_schmidt:
			shtns->Y00_1 = 1.0;		shtns->Y10_ct = 1.0;
			break;
		case sht_fourpi:
			shtns->Y00_1 = 1.0;		shtns->Y10_ct = sqrt(1./3.);
			break;
		case sht_orthonormal:
		default:
			shtns->Y00_1 = sqrt(4.*M_PI);		shtns->Y10_ct = sqrt(4.*M_PI/3.);
//			Y11_st = sqrt(2.*M_PI/3.);		// orthonormal :  \f$ \sin\theta\cos\phi/(Y_1^1 + Y_1^{-1}) = -\sqrt{2 \pi /3} \f$
	}
	shtns->Y11_st = shtns->Y10_ct * sqrt(0.5/mpos_renorm);
	if (with_cs_phase)	shtns->Y11_st *= -1.0;		// correct Condon-Shortley phase

// save a pointer to this setup and return.
	shtns->next = sht_data;		// reference of previous setup (may be NULL).
	sht_data = shtns;			// keep track of new setup.
	return(shtns);
}

/// Copy a given config but allow a different (smaller) mmax and the possibility to enable/disable fft (beta).
shtns_cfg shtns_create_with_grid(shtns_cfg base, int mmax, int nofft)
{
	shtns_cfg shtns;

	if (mmax > base->mmax) return (NULL);		// fail if mmax larger than source config.

	shtns = malloc( SIZEOF_SHTNS_INFO(mmax) );
	memcpy(shtns, base, SIZEOF_SHTNS_INFO(mmax) );		// copy all
	shtns->lmidx = (int*) shtns+1;		// lmidx is stored at the end of the struct...
	shtns->tm = (unsigned short*) (shtns->lmidx + (mmax+1));		// ...and tm just after.

	if (mmax != shtns->mmax) {
		shtns->mmax = mmax;
		for (int im=0; im<=mmax; im++) {
			shtns->lmidx[im] = base->lmidx[im];
			shtns->tm[im] = base->tm[im];
		}
		if (mmax < shtns->mtr_dct) {
			shtns->idct = NULL;		// do not destroy the plan of the source.
			Set_MTR_DCT(shtns, mmax);		// adjut mtr_dct if required.
		}
		if (mmax == 0) {
			// TODO we may disable fft and replace with a phi-averaging function ...
			// ... then switch to axisymmetric functions :
			// init_sht_array_func(shtns);
			// choose_best_sht(shtns, &nloop, 0);
		}
	}
	if (nofft != 0) {
		shtns->ncplx_fft = -1;		// fft disabled.
	}

// save a pointer to this setup and return.
	shtns->next = sht_data;		// reference of previous setup (may be NULL).
	sht_data = shtns;			// keep track of new setup.
	return(shtns);
}

/// release all resources allocated by a grid.
void shtns_unset_grid(shtns_cfg shtns)
{
	if (ref_count(shtns, &shtns->wg) == 1)	VFREE(shtns->wg);
	shtns->wg = NULL;
	free_SH_dct(shtns);
	free_SHTarrays(shtns);
	shtns->nlat = 0;	shtns->nlat_2 = 0;
	shtns->nphi = 0;	shtns->nspat = 0;
}

/// release all resources allocated by a given shtns_cfg.
void shtns_destroy(shtns_cfg shtns)
{
	free_unused(shtns, &shtns->l_2);
	if (shtns->blm != shtns->alm)
		free_unused(shtns, &shtns->blm);
	free_unused(shtns, &shtns->alm);
	free_unused(shtns, &shtns->li);

	shtns_unset_grid(shtns);

	if (sht_data == shtns) {
		sht_data = shtns->next;		// forget shtns
	} else {
		shtns_cfg s2 = sht_data;
		while (s2 != NULL) {
			if (s2->next == shtns) {
				s2->next = shtns->next;		// forget shtns
				break;
			}
			s2 = s2->next;
		}
	}
	free(shtns);
}

/// clear all allocated memory (hopefully) and go back to 0 state.
void shtns_reset()
{
	while (sht_data != NULL) {
		shtns_destroy(sht_data);
	}
}

/*! Initialization of Spherical Harmonic transforms (backward and forward, vector and scalar, ...) of given size.
 * <b>This function must be called after \ref shtns_create and before any SH transform.</b> and sets all global variables and internal data.
 * returns the required number of doubles to be allocated for a spatial field.
 * \param shtns is the config created by \ref shtns_create for which the grid will be set.
 * \param nlat,nphi pointers to the number of latitudinal and longitudinal grid points respectively. If 0, they are set to optimal values.
 * \param nl_order defines the maximum SH degree to be resolved by analysis : lmax_analysis = lmax*nl_order. It is used to set an optimal and anti-aliasing nlat. If 0, the default SHT_DEFAULT_NL_ORDER is used.
 * \param flags allows to choose the type of transform (see \ref shtns_type) and the spatial data layout (see \ref spat)
 * \param eps polar optimization threshold : polar values of Legendre Polynomials below that threshold are neglected (for high m), leading to increased performance (a few percents)
 *  0 = no polar optimization;  1.e-14 = VERY safe;  1.e-10 = safe;  1.e-6 = aggresive, but still good accuracy.
*/
int shtns_set_grid_auto(shtns_cfg shtns, enum shtns_type flags, double eps, int nl_order, int *nlat, int *nphi)
{
	double t, mem;
	int im,m;
	int layout;
	int nloop = 0;
	int n_gauss = 0;
	int on_the_fly = 0;
	int quick_init = 0;
	int vector = !(flags & SHT_SCALAR_ONLY);
	int latdir = (flags & SHT_SOUTH_POLE_FIRST) ? -1 : 1;		// choose latitudinal direction (change sign of ct)
	int analys = 1;

	#if _GCC_VEC_
		if (*nlat & 1) shtns_runerr("Nlat must be even\n");
	#endif
	shtns_unset_grid(shtns);		// release grid if previously allocated.
	if (nl_order <= 0) nl_order = SHT_DEFAULT_NL_ORDER;
/*	shtns.lshift = 0;
	if (nl_order == 0) nl_order = SHT_DEFAULT_NL_ORDER;
	if (nl_order < 0) {	shtns.lshift = -nl_order;	nl_order = 1; }		// linear with a shift in l.
*/
	shtns->nspat = 0;
	shtns->nlorder = nl_order;
	shtns->mtr_dct = -1;		// dct switched off completely.
	layout = flags & 0xFFFF00;
	flags = flags & 255;	// clear higher bits.

	switch (flags) {
		case sht_gauss_fly :  flags = sht_gauss;  on_the_fly = 1;  break;
		case sht_quick_init : flags = sht_gauss;  quick_init = 1;  break;
		case sht_reg_poles :  analys = 0;         quick_init = 1;  break;
		default : break;
	}

	if (*nphi == 0) {
		*nphi = fft_int((nl_order+1)*MMAX+1, 7);		// required fft nodes
	}
	if (*nlat == 0) {
		n_gauss = ((nl_order+1)*LMAX)/2 +1;		// required gauss nodes
		n_gauss += (n_gauss&1);		// even is better.
		if (flags != sht_gauss) {
			m = fft_int(nl_order*LMAX+2, 7);		// required dct nodes
			*nlat = m + (m&1);		// even is better.
		} else *nlat = n_gauss;
	}

	mem = sht_mem_size(shtns->lmax, shtns->mmax, shtns->mres, *nlat);
	t=mem;	if (analys) t*=2;		if (vector) t*=3;
	#if SHT_VERBOSE > 1
		printf("Memory required for precomputed matrices (estimate) : %.3f Mb\n",t);
	#endif
	if ( t > SHTNS_MAX_MEMORY ) {		// huge transform has been requested
		on_the_fly = 1;
		if ( (flags == sht_reg_dct) || (flags == sht_reg_fast) ) shtns_runerr("Memory limit exceeded, try using sht_gauss or increase SHTNS_MAX_MEMORY in sht_config.h");
		if (flags != sht_reg_poles) {
			flags = sht_gauss;
			if (n_gauss > 0) *nlat = n_gauss;
		}
//		if (t > 10*SHTNS_MAX_MEMORY) quick_init =1;			// do not time such large transforms.
	}

	if (quick_init == 0) {		// do not waste too much time finding optimal fftw.
		shtns->fftw_plan_mode = FFTW_EXHAUSTIVE;		// defines the default FFTW planner mode.
	// fftw_set_timelimit(60.0);		// do not search plans for more than 1 minute (does it work well ???)
		if (*nphi > 512) shtns->fftw_plan_mode = FFTW_PATIENT;
		if (*nphi > 1024) shtns->fftw_plan_mode = FFTW_MEASURE;
	} else {
		shtns->fftw_plan_mode = FFTW_ESTIMATE;
		if ((mem < 1.0) && (SHT_VERBOSE < 2)) shtns->nthreads = 1;		// disable threads for small transforms (in quickinit mode).
		if ((VSIZE2 >= 4) && (*nlat >= VSIZE2*4)) on_the_fly = 1;		// with AVX, on-the-fly should be the default (faster).
		if ((shtns->nthreads > 1) && (*nlat >= VSIZE2*16)) on_the_fly = 1;		// force multi-thread transforms
	}

	if (flags == sht_auto) {
		if ( ((nl_order>=2)&&(MMAX*MRES > LMAX/2)) || (*nlat < SHT_MIN_NLAT_DCT) || (*nlat & 1) || (*nlat <= LMAX+1) ) {
			flags = sht_gauss;		// avoid computing DCT stuff when it is not expected to be faster.
			if (n_gauss > 0) *nlat = n_gauss;
		}
	}

	if (*nlat <= shtns->lmax) shtns_runerr("Nlat must be larger than Lmax");
	if (IS_TOO_LARGE(*nlat, shtns->nlat)) shtns_runerr("Nlat too large");
	if (IS_TOO_LARGE(*nphi, shtns->nphi)) shtns_runerr("Nphi too large");

	// copy to global variables.
	shtns->nphi = *nphi;
	shtns->nlat_2 = (*nlat+1)/2;	shtns->nlat = *nlat;

	planFFT(shtns, layout, on_the_fly);		// initialize fftw
	shtns->zlm_dct0 = NULL;		// used as a flag.
	init_sht_array_func(shtns);		// array of SHT functions is now set.

	if (flags == sht_reg_dct) {		// pure dct.
		alloc_SHTarrays(shtns, on_the_fly, vector, analys);		// allocate dynamic arrays
		grid_dct(shtns, latdir);		init_SH_dct(shtns, 1);
		OptimizeMatrices(shtns, eps);
		Set_MTR_DCT(shtns, MMAX);
	}
	if ((flags == sht_auto)||(flags == sht_reg_fast))
	{
		alloc_SHTarrays(shtns, on_the_fly, vector, analys);		// allocate dynamic arrays
		grid_dct(shtns, latdir);		init_SH_dct(shtns, 1);
		OptimizeMatrices(shtns, eps);
		if (NLAT >= SHT_MIN_NLAT_DCT) {			// dct requires large NLAT to perform well.
			t = choose_best_sht(shtns, &nloop, vector, 1);		// find optimal MTR_DCT.
			if ((n_gauss > 0)&&(flags == sht_auto)) t *= ((double) n_gauss)/NLAT;	// we can revert to gauss with a smaller nlat.
			if (t < MIN_PERF_IMPROVE_DCT) {
				Set_MTR_DCT(shtns, -1);		// turn off DCT.
			} else {
				t = SHT_error(shtns, vector);
				if (t > MIN_ACCURACY_DCT) {
				#if SHT_VERBOSE > 0
					printf("     !! Not enough accuracy (%.3g) => DCT disabled.\n",t);
				#endif
				#if SHT_VERBOSE < 2
					Set_MTR_DCT(shtns, -1);		// turn off DCT.
				#endif
				}
			}
		}
		if (MTR_DCT < 0) {			// free memory used by DCT and disables DCT.
			free_SH_dct(shtns);			// free now useless arrays.
			if (flags == sht_auto) {
				flags = sht_gauss;		// switch to gauss grid, even better accuracy.
		#if SHT_VERBOSE > 0
				printf("        => switching back to Gauss Grid\n");
		#endif
				for (im=1; im<=MMAX; im++) {	//	im >= 1
					m = im*MRES;
					shtns->ylm[im]  -= shtns->tm[im]*(LMAX-m+1);		// restore pointers altered by OptimizeMatrices().
					if (vector)  shtns->dylm[im] -= shtns->tm[im]*(LMAX-m+1);
				}
				if (n_gauss > 0) {		// we should use the optimal size for gauss-legendre
					free_SHTarrays(shtns);
					*nlat = n_gauss;
					shtns->nlat_2 = (*nlat+1)/2;	shtns->nlat = *nlat;
					planFFT(shtns, layout, on_the_fly);		// fft must be replanned because NLAT has changed.
				}
			}
		}
	}
	if (flags == sht_gauss)
	{
		alloc_SHTarrays(shtns, on_the_fly, vector, analys);		// allocate dynamic arrays
		grid_gauss(shtns, latdir);
		if (on_the_fly == 0) {
			init_SH_gauss(shtns);			// precompute matrices
			OptimizeMatrices(shtns, eps);
		}
	}
	if (flags == sht_reg_poles)
	{
		alloc_SHTarrays(shtns, on_the_fly, vector, 0);		// allocate dynamic arrays (no analysis)
		grid_equal_polar(shtns, latdir);
		if (on_the_fly == 0) {
			init_SH_synth(shtns);
			for (im=0; im<=MMAX; im++) shtns->tm[im] = 0;	// avoid problems with tm[im] modified ????
		}
	}

	if (on_the_fly == 1) {
  #if SHT_VERBOSE > 0
		printf("        + using on-the-fly transforms.\n");
  #endif
		if (NLAT < VSIZE2*4) shtns_runerr("on-the-fly only available for nlat>=32");		// avoid overflow with NLAT_2 < VSIZE2*2
		PolarOptimize(shtns, eps);
		set_sht_fly(shtns, 0);		// switch function pointers to "on-the-fly" functions.
	}

	if (quick_init == 0) {
		choose_best_sht(shtns, &nloop, vector, 0);
		if (on_the_fly == 0) free_unused_matrices(shtns);
		t = SHT_error(shtns, vector);		// compute SHT accuracy.
  #if SHT_VERBOSE > 0
		printf("        + SHT accuracy = %.3g\n",t);
  #endif
  #if SHT_VERBOSE < 2
		if (t > 1.e-3) shtns_runerr("bad SHT accuracy");		// stop if something went wrong (but not in debug mode)
  #endif
	}
//	set_sht_fly(shtns, SHT_TYP_VAN);
  #if SHT_VERBOSE > 1
	if (omp_threads > 1) printf(" nthreads = %d\n",shtns->nthreads);
  #endif
  #if SHT_VERBOSE > 0
	printf("        => " PACKAGE_NAME " is ready.\n");
  #endif
	return(shtns->nspat);	// returns the number of doubles to be allocated for a spatial field.
}


/*! Initialization of Spherical Harmonic transforms (backward and forward, vector and scalar, ...) of given size.
 * <b>This function must be called after \ref shtns_create and before any SH transform.</b> and sets all global variables.
 * returns the required number of doubles to be allocated for a spatial field.
 * \param shtns is the config created by shtns_create for which the grid will be set.
 * \param flags allows to choose the type of transform (see \ref shtns_type) and the spatial data layout (see \ref spat)
 * \param eps polar optimization threshold : polar values of Legendre Polynomials below that threshold are neglected (for high m), leading to increased performance (a few percents)
 * \param nlat,nphi respectively the number of latitudinal and longitudinal grid points.
 *  0 = no polar optimization;  1.e-14 = VERY safe;  1.e-10 = safe;  1.e-6 = aggresive, but still good accuracy.
*/
int shtns_set_grid(shtns_cfg shtns, enum shtns_type flags, double eps, int nlat, int nphi)
{
	if ((nlat == 0)||(nphi == 0)) shtns_runerr("nlat or nphi is zero !");
	return( shtns_set_grid_auto(shtns, flags, eps, 0, &nlat, &nphi) );
}

/*! Simple initialization of Spherical Harmonic transforms (backward and forward, vector and scalar, ...) of given size.
 * This function sets all global variables by calling \ref shtns_create followed by \ref shtns_set_grid, with the
 * default normalization and the default polar optimization (see \ref sht_config.h).
 * Returns the configuration to be passed to subsequent transform functions, which is basicaly a pointer to a \ref shtns_info struct.
 * \param lmax : maximum SH degree that we want to describe.
 * \param mmax : number of azimutal wave numbers.
 * \param mres : \c 2.pi/mres is the azimutal periodicity. \c mmax*mres is the maximum SH order.
 * \param nlat,nphi : respectively the number of latitudinal and longitudinal grid points.
 * \param flags allows to choose the type of transform (see \ref shtns_type) and the spatial data layout (see \ref spat)
*/
shtns_cfg shtns_init(enum shtns_type flags, int lmax, int mmax, int mres, int nlat, int nphi)
{
	shtns_cfg shtns = shtns_create(lmax, mmax, mres, SHT_DEFAULT_NORM);
	if (shtns != NULL)
		shtns_set_grid(shtns, flags, SHT_DEFAULT_POLAR_OPT, nlat, nphi);
	return shtns;
}

/** Enables OpenMP parallel transforms, if available (see \ref compil).
 Call before any initialization of shtns to use mutliple threads. Returns the actual number of threads.
 \li If num_threads > 0, specifies the maximum number of threads that should be used.
 \li If num_threads <= 0, maximum number of threads is automatically set to the number of processors.
 \li If num_threads == 1, openmp will be disabled. */
int shtns_use_threads(int num_threads)
{
#ifdef _OPENMP
	int procs = omp_get_num_procs();
	if (num_threads <= 0)  num_threads = procs;
	else if (num_threads > 4*procs) num_threads = 4*procs;		// limit the number of threads
	omp_threads = num_threads;
#endif
#ifdef OMP_FFTW
	fftw_init_threads();		// enable threads for FFTW.
#endif
	return omp_threads;
}

/// fill the given array with Gauss weights. returns the number of weights written, which
/// may be zero if the grid is not a Gauss grid.
int shtns_gauss_wts(shtns_cfg shtns, double *wts)
{
	int i = 0;
	if (shtns->wg) {
		double rescale = 2*NPHI;		// weights are stored with a rescaling that depends on SHT_NORM.
		if ((SHT_NORM != sht_fourpi)&&(SHT_NORM != sht_schmidt))  rescale *= 0.25/M_PI;

		do {
			wts[i] = shtns->wg[i] * rescale;
		} while(++i < shtns->nlat_2);
	}
	return i;
}

//@}


#ifdef SHT_F77_API

/* FORTRAN API */

/** \addtogroup fortapi Fortran API.
* Call from fortran without the trailing '_'.
* see the \link SHT_example.f Fortran example \endlink for a simple usage of SHTns from Fortran language.
*/
//@{

/// Enable threads
void shtns_use_threads_(int *num_threads)
{
	shtns_use_threads(*num_threads);
}
	
/// Initializes spherical harmonic transforms of given size using Gauss algorithm with default polar optimization.
void shtns_init_sh_gauss_(int *layout, int *lmax, int *mmax, int *mres, int *nlat, int *nphi)
{
	shtns_reset();
	shtns_cfg shtns = shtns_create(*lmax, *mmax, *mres, SHT_DEFAULT_NORM);
	shtns_set_grid(shtns, sht_gauss | *layout, SHT_DEFAULT_POLAR_OPT, *nlat, *nphi);
}

/// Initializes spherical harmonic transforms of given size using Fastest available algorithm and polar optimization.
void shtns_init_sh_auto_(int *layout, int *lmax, int *mmax, int *mres, int *nlat, int *nphi)
{
	shtns_reset();
	shtns_cfg shtns = shtns_create(*lmax, *mmax, *mres, SHT_DEFAULT_NORM);
	shtns_set_grid(shtns, sht_auto | *layout, SHT_DEFAULT_POLAR_OPT, *nlat, *nphi);
}

/// Initializes spherical harmonic transforms of given size using a regular grid and agressive optimizations.
void shtns_init_sh_reg_fast_(int *layout, int *lmax, int *mmax, int *mres, int *nlat, int *nphi)
{
	shtns_reset();
	shtns_cfg shtns = shtns_create(*lmax, *mmax, *mres, SHT_DEFAULT_NORM);
	shtns_set_grid(shtns, sht_reg_fast | *layout, 1.e-6, *nlat, *nphi);
}

/// Initializes spherical harmonic transform SYNTHESIS ONLY of given size using a regular grid including poles.
void shtns_init_sh_poles_(int *layout, int *lmax, int *mmax, int *mres, int *nlat, int *nphi)
{
	shtns_reset();
	shtns_cfg shtns = shtns_create(*lmax, *mmax, *mres, SHT_DEFAULT_NORM);
	shtns_set_grid(shtns, sht_reg_poles | *layout, 0, *nlat, *nphi);
}

/// Defines the size and convention of the transform.
/// Allow to choose the normalization and whether or not to include the Condon-Shortley phase.
/// \see shtns_create
void shtns_set_size_(int *lmax, int *mmax, int *mres, int *norm)
{
	shtns_reset();
	shtns_create(*lmax, *mmax, *mres, *norm);
}

/// Precompute matrices for synthesis and analysis.
/// Allow to choose polar optimization threshold and algorithm type.
/// \see shtns_set_grid
void shtns_precompute_(int *type, int *layout, double *eps, int *nlat, int *nphi)
{
	shtns_set_grid(sht_data, *type | *layout, *eps, *nlat, *nphi);
}

/// Same as shtns_precompute_ but choose optimal nlat and/or nphi.
/// \see shtns_set_grid_auto
void shtns_precompute_auto_(int *type, int *layout, double *eps, int *nl_order, int *nlat, int *nphi)
{
	shtns_set_grid_auto(sht_data, *type | *layout, *eps, *nl_order, nlat, nphi);
}

/// Clear everything
void shtns_reset_() {
	shtns_reset();
}

/// returns nlm, the number of complex*16 elements in an SH array.
/// call from fortran using \code call shtns_calc_nlm(nlm, lmax, mmax, mres) \endcode
void shtns_calc_nlm_(int *nlm, const int *const lmax, const int *const mmax, const int *const mres)
{
    *nlm = nlm_calc(*lmax, *mmax, *mres);
}

/// returns lm, the index in an SH array of mode (l,m).
/// call from fortran using \code call shtns_lmidx(lm, l, m) \endcode
void shtns_lmidx_(int *lm, const int *const l, const int *const m)
{
	unsigned im = *m;
	unsigned mres = sht_data->mres;
	if (mres > 1) {
		unsigned k = im % mres;
		im = im / mres;
		if (k) printf("wrong m");
	}
    *lm = LiM(sht_data, *l, im) + 1;	// convert to fortran convention index.
}

/// returns l and m, degree and order of an index in SH array lm.
/// call from fortran using \code call shtns_l_m(l, m, lm) \endcode
void shtns_l_m_(int *l, int *m, const int *const lm)
{
	*l = sht_data->li[*lm -1];	// convert from fortran convention index.
	*m = sht_data->mi[*lm -1];
}

/// fills the given array with the cosine of the co-latitude angle (NLAT real*8)
/// if no grid has been set, the first element will be set to zero.
void shtns_cos_array_(double *costh)
{
	if (sht_data->ct) {
		for (int i=0; i<sht_data->nlat; i++)
			costh[i] = sht_data->ct[i];
	} else costh[0] = 0.0;	// mark as invalid.
}

/// fills the given array with the gaussian quadrature weights ((NLAT+1)/2 real*8).
/// when there is no gaussian grid, the first element is set to zero.
void shtns_gauss_wts_(double *wts)
{
	int i = shtns_gauss_wts(sht_data, wts);
	if (i==0) wts[0] = 0;	// mark as invalid.
}

/** \name Point evaluation of Spherical Harmonics
Evaluate at a given point (\f$cos(\theta)\f$ and \f$\phi\f$) a spherical harmonic representation.
*/
//@{
/// \see SH_to_point for argument description
void shtns_sh_to_point_(double *spat, complex double *Qlm, double *cost, double *phi)
{
	*spat = SH_to_point(sht_data, Qlm, *cost, *phi);
}

/// \see SHqst_to_point for argument description
void shtns_qst_to_point_(double *vr, double *vt, double *vp,
		complex double *Qlm, complex double *Slm, complex double *Tlm, double *cost, double *phi)
{
	SHqst_to_point(sht_data, Qlm, Slm, Tlm, *cost, *phi, vr, vt, vp);
}
//@}

//@}

#endif
