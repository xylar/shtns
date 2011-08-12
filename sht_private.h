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

#include <stdlib.h>
#include <complex.h>
#include <math.h>
// FFTW la derivee d/dx = ik	(pas de moins !)
#include <fftw3.h>

#define SHTNS_PRIVATE

#include "sht_config.h"
#include "shtns.h"

// sht variants (std, ltr)
enum sht_variants { SHT_STD, SHT_LTR, SHT_M0, SHT_M0LTR, SHT_NVAR };
// sht types (scal synth, scal analys, vect synth, ...)
enum sht_types { SHT_TYP_SSY, SHT_TYP_SAN, SHT_TYP_VSY, SHT_TYP_VAN,
				SHT_TYP_GSP, SHT_TYP_GTO, SHT_TYP_3SY, SHT_TYP_3AN, SHT_NTYP };
// sht algorithms (hyb, fly1, ...)
enum sht_algos { SHT_DCT, SHT_MEM, SHT_SV, SHT_FLY1, SHT_FLY2, SHT_FLY3, SHT_FLY4, SHT_FLY6, SHT_FLY8, SHT_NALG };

// sht grids
enum sht_grids { GRID_NONE, GRID_GAUSS, GRID_REGULAR, GRID_POLES };

// pointer to various function types
typedef void (*pf2l)(shtns_cfg, void*, void*, long int);
typedef void (*pf3l)(shtns_cfg, void*, void*, void*, long int);
typedef void (*pf4l)(shtns_cfg, void*, void*, void*, void*, long int);
typedef void (*pf6l)(shtns_cfg, void*, void*, void*, void*, void*, void*, long int);

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
	int *lmidx;					///< (virtual) index in SH array of given im (size mmax+1) : LiM(l,im) = lmidx[im] + l
	unsigned short *li;			///< degree l for given mode number (size nlm) : li[lm] 
	double *ct, *st;			///< cos(theta) and sin(theta) arrays (size nlat)
	unsigned nspat;				///< number of real numbers that must be allocated in a spatial field.
/* END OF PUBLIC PART */

	int ncplx_fft;			///< number of complex numbers to allocate for the fft : -1 = no fft; 0 = in-place fft (no allocation).
	unsigned short *tm;		///< start theta value for SH (polar optimization : near the poles the legendre polynomials go to zero for high m's)
	double *wg;				///< Gauss weights for Gauss-Legendre quadrature.
	double *st_1;			///< 1/sin(theta);

	fftw_plan ifft, fft;		// plans for FFTW.

	/* Legendre function generation arrays */
	double *al0;	double **alm;	// coefficient list for Legendre function recurrence (size 2*NLM)
	double *bl0;	double **blm;	// coefficient list for modified Legendre function recurrence for analysis (size 2*NLM)
	double *l_2;	// array of size (LMAX+1) containing 1./l(l+1) for increasing integer l.

	void* fptr[SHT_NVAR][SHT_NTYP];		// pointers to transform functions.

	/* MEM matrices */
	double **ylm;		// matrix for inverse transform (synthesis)
	struct DtDp** dylm;	// theta and phi derivative of Ylm matrix
	double **zlm;		// matrix for direct transform (analysis)
	struct DtDp** dzlm;

	/* DCT stuff */
	short mtr_dct;			///< m truncation for dct. -1 means no dct at all.
	short nlorder;			///< order of non-linear terms to be resolved by SH transform.
	unsigned short klim;	///< Limit to k for non-linear terms (dct)
	short dummy;			// alignement.
	fftw_plan idct, dct_m0;			// (I)DCT for NPHI>1
	fftw_plan idct_r1, dct_r1;		// (I)DCT for axisymmetric case, NPHI=1
	double **ykm_dct;	// matrix for inverse transform (synthesis) using dct.
	struct DtDp** dykm_dct;	// theta and phi derivative of Ykm matrix.
	double *zlm_dct0;	// matrix for direct transform (analysis), only m=0
	double *dzlm_dct0;

	/* other misc informations */
	double Y00_1, Y10_ct, Y11_st;
	shtns_cfg next;		// pointer to next sht_setup or NULL (records a chained list of SHT setup).
	short norm;			// store the normalization of the Spherical Harmonics (enum \ref shtns_norm + \ref SHT_NO_CS_PHASE flag)
	short grid;			// store grid type.
	unsigned fftw_plan_mode;
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

struct DtDp {		// theta and phi derivatives stored together.
	double t, p;
};

#define GLUE2(a,b) a##b
#define GLUE3(a,b,c) a##b##c

// how to allocate memory aligned on 16 bytes ? this works only if fftw was compiled with --enable-sse2...
// in 64 bit system, malloc should be 16bytes aligned.
#define VMALLOC(s)	( (sizeof(void*) >= 8) ? malloc(s) : fftw_malloc(s) )
#define VFREE(s)	( (sizeof(void*) >= 8) ? free(s) : fftw_free(s) )
