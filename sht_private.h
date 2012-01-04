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
enum sht_variants { SHT_STD, SHT_LTR, SHT_NVAR };
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
	unsigned int nlm;			///< total number of (l,m) spherical harmonics components.
	unsigned short lmax;		///< maximum degree (lmax) of spherical harmonics.
	unsigned short mmax;		///< maximum order (mmax*mres) of spherical harmonics.
	unsigned short mres;		///< the periodicity along the phi axis.
	unsigned short nphi;		///< number of spatial points in Phi direction (longitude)
	unsigned short nlat;		///< number of spatial points in Theta direction (latitude) ...
	unsigned short nlat_2;		///< ...and half of it (using (shtns.nlat+1)/2 allows odd shtns.nlat.)
	int *lmidx;					///< (virtual) index in SH array of given im (size mmax+1) : LiM(l,im) = lmidx[im] + l
	unsigned short *li;			///< degree l for given mode number (size nlm) : li[lm] 
	double *ct, *st;			///< cos(theta) and sin(theta) arrays (size nlat)
	unsigned int nspat;			///< number of real numbers that must be allocated in a spatial field.
/* END OF PUBLIC PART */

	int ncplx_fft;			///< number of complex numbers to allocate for the fft : -1 = no fft; 0 = in-place fft (no allocation).
	unsigned short *tm;		///< start theta value for SH (polar optimization : near the poles the legendre polynomials go to zero for high m's)
	double *wg;				///< Gauss weights for Gauss-Legendre quadrature.
	double *st_1;			///< 1/sin(theta);

	fftw_plan ifft, fft;		// plans for FFTW.
	fftw_plan ifftc, fftc;

	/* Legendre function generation arrays */
	double **alm;	// coefficient list for Legendre function recurrence (size 2*NLM)
	double **blm;	// coefficient list for modified Legendre function recurrence for analysis (size 2*NLM)
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
	short fftc_mode;		///< how to perform the complex fft : -1 = no fft; 0 = interleaved/native; 1 = split/transpose.
	fftw_plan idct, dct_m0;			// (I)DCT
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

// define shortcuts to sizes.
#define NLM shtns->nlm
#define LMAX shtns->lmax
#define NLAT shtns->nlat
#define NLAT_2 shtns->nlat_2
#define NPHI shtns->nphi
#define MMAX shtns->mmax
#define MRES shtns->mres
#define MTR_DCT shtns->mtr_dct
#define SHT_NL_ORDER shtns->nlorder

// define index in alm/blm matrices
#define ALM_IDX(shtns, im) ( (im)*(2*shtns->lmax - ((im)-1)*shtns->mres) )

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
// set to a value close to the machine accuracy, it allows to speed-up on-the-fly SHTs with very large l (lmax > SHT_L_RESCALE_FLY).
#define SHT_ACCURACY 1.0e-15

/* for vectorization (SSE2) */

#define SSE __attribute__((aligned (16)))

#ifdef __INTEL_COMPILER
	#undef _GCC_VEC_
#endif

#if _GCC_VEC_ && __SSE2__
	#define VSIZE 2
	typedef double s2d __attribute__ ((vector_size (8*VSIZE)));		// vector that should behave like a real scalar for complex number multiplication.
	typedef double v2d __attribute__ ((vector_size (8*VSIZE)));		// vector that contains a complex number
	#ifdef __AVX__
		#define VSIZE2 4
		#include <immintrin.h>
		#warning "using GCC vector extensions (avx)"
		typedef double rnd __attribute__ ((vector_size (VSIZE2*8)));		// vector of 4 doubles.
		#define vall(x) _mm256_set1_pd(x)
		#define vread(mem, idx) _mm256_loadu_pd( ((double*)mem) + VSIZE2*(idx) )
		#define vlo(a) __builtin_ia32_vec_ext_v2df (_mm256_extractf128_pd(a,0), 0)
		#define S2D_STORE(mem, idx, ev, od) \
			((s2d*)mem)[(idx)*2] = _mm256_extractf128_pd(ev+od,0); \
			((s2d*)mem)[(idx)*2+1] = _mm256_extractf128_pd(ev+od,1); \
			((s2d*)mem)[NLAT_2-1 - (idx)*2] = _mm256_extractf128_pd(_mm256_permute_pd(ev-od,5), 0); \
			((s2d*)mem)[NLAT_2-2 - (idx)*2] = _mm256_extractf128_pd(_mm256_permute_pd(ev-od,5), 1);

		#define S2D_CSTORE(mem, idx, er, or, ei, oi)	{	\
			rnd aa = _mm256_permute_pd(ei+oi,5) + (er + or);		rnd bb = (er + or) - _mm256_permute_pd(ei+oi,5);	\
			((s2d*)mem)[(idx)*2] = _mm256_extractf128_pd(_mm256_shuffle_pd(bb, aa, 10 ),0);	\
			((s2d*)mem)[(idx)*2+1] = _mm256_extractf128_pd(_mm256_shuffle_pd(bb, aa, 10 ),1);	\
			((s2d*)mem)[(NPHI-2*im)*NLAT_2 + (idx)*2] = _mm256_extractf128_pd(_mm256_shuffle_pd(aa, bb, 10 ),0);	\
			((s2d*)mem)[(NPHI-2*im)*NLAT_2 + (idx)*2+1] = _mm256_extractf128_pd(_mm256_shuffle_pd(aa, bb, 10 ),1);	\
			aa = _mm256_permute_pd(er-or,5) + (ei - oi);		bb = _mm256_permute_pd(er-or,5) - (ei - oi);	\
			((s2d*)mem)[NLAT_2-1 -(idx)*2] = _mm256_extractf128_pd(_mm256_shuffle_pd(bb, aa, 10 ), 0);	\
			((s2d*)mem)[NLAT_2-2 -(idx)*2] = _mm256_extractf128_pd(_mm256_shuffle_pd(bb, aa, 10 ), 1);	\
			((s2d*)mem)[(NPHI+1-2*im)*NLAT_2 -1 -(idx)*2] = _mm256_extractf128_pd(_mm256_shuffle_pd(aa, bb, 10 ), 0);	\
			((s2d*)mem)[(NPHI+1-2*im)*NLAT_2 -2 -(idx)*2] = _mm256_extractf128_pd(_mm256_shuffle_pd(aa, bb, 10 ), 1);	}
	#else
		#define VSIZE2 2
		#ifdef __SSE3__
			#include <pmmintrin.h>
			#warning "using GCC vector extensions (sse3)"
		#else
			#include <emmintrin.h>
			#warning "using GCC vector extensions (sse2)"
		#endif
		typedef double rnd __attribute__ ((vector_size (VSIZE2*8)));		// vector of 2 doubles.
		#define vall(x) _mm_set1_pd(x)
		#define vread(mem, idx) ((s2d*)mem)[idx]
		#define vlo(a) __builtin_ia32_vec_ext_v2df (a, 0)
		#define S2D_STORE(mem, idx, ev, od)		((s2d*)mem)[idx] = ev+od;		((s2d*)mem)[NLAT_2-1 - (idx)] = vxchg(ev-od);
		#define S2D_CSTORE(mem, idx, er, or, ei, oi)	{	\
			rnd aa = vxchg(ei + oi) + (er + or);		rnd bb = (er + or) - vxchg(ei + oi);	\
			((s2d*)mem)[idx] = _mm_shuffle_pd(bb, aa, 2 );	\
			((s2d*)mem)[(NPHI-2*im)*NLAT_2 + (idx)] = _mm_shuffle_pd(aa, bb, 2 );	\
			aa = vxchg(er - or) + (ei - oi);		bb = vxchg(er - or) - (ei - oi);	\
			((s2d*)mem)[NLAT_2-1 -(idx)] = _mm_shuffle_pd(bb, aa, 2 );	\
			((s2d*)mem)[(NPHI+1-2*im)*NLAT_2 -1 -(idx)] = _mm_shuffle_pd(aa, bb, 2 );	}
	#endif
	#ifdef __SSE3__
		#define addi(a,b) _mm_addsub_pd(a, _mm_shuffle_pd(b,b,1))		// a + I*b
		#define subadd(a,b) _mm_addsub_pd(a, b)		// [al-bl, ah+bh]
	#else
		#define addi(a,b) ( (a) + (_mm_shuffle_pd(b,b,1) * _mm_set_pd(1.0, -1.0)) )		// a + I*b		[note: _mm_set_pd(imag, real)) ]
		#define subadd(a,b) ( (a) + (b) * _mm_set_pd(1.0, -1.0) )		// [al-bl, ah+bh]
	#endif
	// vset(lo, hi) takes two doubles and pack them in a vector
	#define vset(lo, hi) _mm_set_pd(hi, lo)
	// vdup(x) takes a double and duplicate it to a vector of 2 doubles.
	#define vdup(x) _mm_set1_pd(x)
	// vxchg(a) exchange hi and lo component of vector a
	#define vxchg(a) _mm_shuffle_pd(a,a,1)
	#define vlo_to_cplx(a) _mm_unpacklo_pd(a, vdup(0.0))
	#define vhi_to_cplx(a) _mm_unpackhi_pd(a, vdup(0.0))
	#ifdef __clang__
		// allow to compile with clang (llvm)
		#define vlo_to_dbl(a) (a)[0]
		#define vhi_to_dbl(a) (a)[1]
	#else
		// gcc extensions
		#define vlo_to_dbl(a) __builtin_ia32_vec_ext_v2df (a, 0)
		#define vhi_to_dbl(a) __builtin_ia32_vec_ext_v2df (a, 1)
	#endif
#else
	#undef _GCC_VEC_
	#define VSIZE 1
	#define VSIZE2 1
	typedef double s2d;
	typedef complex double v2d;
	typedef double rnd;
	#define vall(x) (x)
	#define vdup(x) (x)
	#define vxchg(x) (x)
	#define addi(a,b) ((a) + I*(b))
	#define vlo_to_dbl(a) (a)
	#define vhi_to_dbl(a) (a)
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

#define MIN_ALIGNMENT 16
/// align pointer on MIN_ALIGNMENT (must be a power of 2)
#define PTR_ALIGN(p) ((((size_t)(p)) + (MIN_ALIGNMENT-1)) & (~((size_t)(MIN_ALIGNMENT-1))))

// verbose printing
#if SHT_VERBOSE > 1
  #define PRINT_VERB(msg) printf(msg)
#else
  #define PRINT_VERB(msg) (0)
#endif

