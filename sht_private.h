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

/********************************************************************
 * SHTns : Spherical Harmonic Transform for numerical simulations.  *
 *    written by Nathanael Schaeffer / CNRS                         *
 ********************************************************************/

/// \internal \file sht_private.h private data used by sht transform functions.

#include <complex.h>
#include <math.h>
// FFTW la derivee d/dx = ik	(pas de moins !)
#include <fftw3.h>

#define SHTNS_PRIVATE

#include "sht_config.h"
#include "shtns.h"

// number of sht variants (std, ltr)
#define SHT_NVAR 2
#define SHT_STD 0
#define SHT_LTR 1
// number of sht types (scal synth, scal analys, vect synth, ...)
#define SHT_NTYP 8
#define SHT_TYP_SSY 0
#define SHT_TYP_SAN 1
#define SHT_TYP_VSY 2
#define SHT_TYP_VAN 3
#define SHT_TYP_GS1 4
#define SHT_TYP_GS2 5
#define SHT_TYP_3SY 6
#define SHT_TYP_3AN 7
// number of sht algorithms (hyb, fly1, ...)
#define SHT_NALG 7
#define SHT_HYB 0
#define SHT_FLY1 1
#define SHT_FLY2 2
#define SHT_FLY3 3
#define SHT_FLY4 4
#define SHT_FLY6 5
#define SHT_FLY8 6


struct DtDp {		// theta and phi derivatives stored together.
	double t, p;
};

// pointer to various function types
typedef void (*pf2)(shtns_cfg, void*, void*);
typedef void (*pf3)(shtns_cfg, void*, void*, void*);
typedef void (*pf4)(shtns_cfg, void*, void*, void*, void*);
typedef void (*pf6)(shtns_cfg, void*, void*, void*, void*, void*, void*);
typedef void (*pf2l)(shtns_cfg, void*, void*, int);
typedef void (*pf3l)(shtns_cfg, void*, void*, void*, int);
typedef void (*pf4l)(shtns_cfg, void*, void*, void*, void*, int);
typedef void (*pf6l)(shtns_cfg, void*, void*, void*, void*, void*, void*, int);

/// structure containing useful information about the SHT.
struct shtns_info {		// MUST start with "int nlm;"
/* PUBLIC PART (if modified, shtns.h should be modified acordingly) */
	int nlm;					///< total number of (l,m) spherical harmonics components.
	unsigned short lmax;		///< maximum degree (lmax) of spherical harmonics.
	unsigned short mmax;		///< maximum order (mmax*mres) of spherical harmonics.
	unsigned short mres;		///< the periodicity along the phi axis.
	unsigned short nphi;		///< number of spatial points in Phi direction (longitude)
	unsigned short nlat;		///< number of spatial points in Theta direction (latitude) ...
	unsigned short nlat_2;		///< ...and half of it (using (shtns.nlat+1)/2 allows odd shtns.nlat.)
	int *lmidx;					///< (virtual) index in SH array of given im : LiM(l,im) = lmidx[im] + l
	unsigned short *li;			///< degree l for given mode number : li[lm]
	double *el, *l2, *l_2;		///< l, l(l+1) and 1/(l(l+1)) arrays.
	double *ct, *st;			///< cos(theta) and sin(theta) arrays.
/* END OF PUBLIC PART */

	short sht_fft;				///< How to perform fft : 0=no fft, 1=in-place, 2=out-of-place.
	short mtr_dct;				///< m truncation for dct. -1 means no dct at all.
	unsigned short klim;		///< Limit to k for non-linear terms.
	unsigned short nlorder;		///< order of non-linear terms to be resolved by SH transform.

	unsigned short *tm;			// start theta value for SH (polar optimization : near the poles the legendre polynomials go to zero for high m's)
	double *wg;					// Gauss weights for Gauss-Legendre quadrature.
	double *st_1;				// 1/sin(theta);

	fftw_plan ifft, fft;			// plans for FFTW.

	double *al0;	double **alm;	// coefficient list for Legendre function recurrence (size 2*NLM)
	double *bl0;	double **blm;	// coefficient list for modified Legendre function recurrence for analysis (size 2*NLM)

	void* fptr[SHT_NVAR][SHT_NTYP];		// pointers to transform functions.

	double **ylm;		// matrix for inverse transform (synthesis)
	struct DtDp** dylm;	// theta and phi derivative of Ylm matrix
	double **zlm;		// matrix for direct transform (analysis)
	struct DtDp** dzlm;

	fftw_plan idct, dct_m0;			// (I)DCT for NPHI>1
	fftw_plan idct_r1, dct_r1;		// (I)DCT for axisymmetric case, NPHI=1

	double **ykm_dct;	// matrix for inverse transform (synthesis) using dct.
	struct DtDp** dykm_dct;	// theta and phi derivative of Ykm matrix.
	double *zlm_dct0;	// matrix for direct transform (analysis), only m=0
	double *dzlm_dct0;

	double Y00_1, Y10_ct, Y11_st;
	void *next;		// pointer to next sht_setup or NULL (records a chained list of SHT setup).
	int norm;		///< store the normalization of the Spherical Harmonics (enum \ref shtns_norm + \ref SHT_NO_CS_PHASE flag)
	// the end should be aligned on the size of int, to allow the storage of small arrays.
};

// define shortcuts to sizes + allow compile-time optimizations when SHT_AXISYM is defined.
#define LMAX shtns->lmax
#define NLAT shtns->nlat
#define NLAT_2 shtns->nlat_2
#ifndef SHT_AXISYM
  #define NPHI shtns->nphi
  #define MMAX shtns->mmax
  #define MRES shtns->mres
  #define SHT_FFT shtns->sht_fft
#else
  #define NPHI 1
  #define MMAX 0
  #define MRES 1
  #define SHT_FFT 0
#endif
#define NLM shtns->nlm
#define MTR_DCT shtns->mtr_dct
#define SHT_NL_ORDER shtns->nlorder

// SHT_NORM without CS_PHASE
#define SHT_NORM (shtns->norm & 0x0FF)

#ifndef M_PI
# define M_PI 3.1415926535897932384626433832795
#endif
#ifndef M_PIl
# define M_PIl 3.1415926535897932384626433832795L
#endif

// scale factor applied for LMAX larger than SHT_L_RESCALE. Allows accurate transforms up to l=2700 with 64 bit double precision.
#define SHT_LEG_SCALEF 1.10180321e-280
#define SHT_L_RESCALE 1536
// value for on-the-fly transforms is near the limit of double precsision, as there is an associated performance drop.
#define SHT_L_RESCALE_FLY 1792

/*
struct SHTdef {
	long int nlm;
	long int lmax,nlat,nlat_2;
	long int mmax,mres,nphi;
	long int mtr_dct;

	double *ct, *st, *st_1;		// cos(theta), sin(theta), 1/sin(theta);
	double *el, *l2, *l_2;		// l, l(l+1) and 1/(l(l+1))
	int *li;

	int *lmidx;		// (virtual) index in SH array of given im.
	int *tm;		// start theta value for SH (polar optimization : near the poles the legendre polynomials go to zero for high m's)

	double** ylm;		// matrix for inverse transform (synthesis)
	double** zlm;		// matrix for direct transform (analysis)
	double** ykm_dct;	// matrix for inverse transform (synthesis) using dct.
	double* zlm_dct0;	// matrix for direct transform (analysis), only m=0

	fftw_plan ifft, fft;	// plans for FFTW.
	fftw_plan idct, dctm0;
}

struct VSHTdef {
	long int nlm;
	long int lmax,nlat,nlat_2;
	long int mmax,mres,nphi;
	long int mtr_dct;

	double *ct, *st, *st_1;		// cos(theta), sin(theta), 1/sin(theta);
	double *el, *l2, *l_2;		// l, l(l+1) and 1/(l(l+1))
	int *li;

	int *lmidx;		// (virtual) index in SH array of given im.
	int *tm;		// start theta value for SH (polar optimization : near the poles the legendre polynomials go to zero for high m's)

	struct DtDp** dylm;	// theta and phi derivative of Ylm matrix
	struct DtDp** dzlm;
	struct DtDp** dykm_dct;	// theta and phi derivative of Ylm matrix

	fftw_plan ifft, fft;	// plans for FFTW.
	fftw_plan idct;

	spat_to_SH
	SH_to_spat
}
*/

/* for vectorization (SSE2) */

#define SSE __attribute__((aligned (16)))

#ifdef __INTEL_COMPILER
	#undef _GCC_VEC_
#endif

#if _GCC_VEC_ && __SSE2__
	typedef double s2d __attribute__ ((vector_size (16)));		// vector that should behave like a real scalar for complex number multiplication.
	typedef double v2d __attribute__ ((vector_size (16)));		// vector that contains a complex number
	#ifdef __SSE3__
		#include <pmmintrin.h>
		#warning "using GCC vector instructions (sse3)"
		#define addi(a,b) _mm_addsub_pd(a, _mm_shuffle_pd(b,b,1))		// a + I*b
		#define subadd(a,b) _mm_addsub_pd(a, b)		// [al-bl, ah+bh]
	#else
		#include <emmintrin.h>
		#warning "using GCC vector instructions (sse2)"
		#define addi(a,b) ( (a) + (_mm_shuffle_pd(b,b,1) * _mm_set_pd(1.0, -1.0)) )		// a + I*b		[note: _mm_set_pd(imag, real)) ]
		#define subadd(a,b) ( (a) + (b) * _mm_set_pd(1.0, -1.0) )		// [al-bl, ah+bh]
	#endif
	// vset(hi,lo) takes two doubles and pack them in a vector
	#define vset(hi, lo) _mm_set_pd(hi, lo)
	// vdup(x) takes a double and duplicate it to a vector of 2 doubles.
	#define vdup(x) _mm_set1_pd(x)
	// vxchg(a) exchange hi and lo component of vector a
	#define vxchg(a) _mm_shuffle_pd(a,a,1)
	#define vlo_to_cplx(a) _mm_shuffle_pd(a, vdup(0.0), 0)
	#define vhi_to_cplx(a) _mm_shuffle_pd(a, vdup(0.0), 1)
	#define vlo_to_dbl(a) __builtin_ia32_vec_ext_v2df (a, 0)
	#define vhi_to_dbl(a) __builtin_ia32_vec_ext_v2df (a, 1)
#else
	#undef _GCC_VEC_
	typedef double s2d;
	typedef complex double v2d;
	// typedef union {v2d v; complex double c; double d[2]; double r; } vcplx;
	#define vdup(x) (x)
	#define addi(a,b) ((a) + I*(b))
#endif

#define GLUE2(a,b) a##b
#define GLUE3(a,b,c) a##b##c

