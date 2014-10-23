/*
 * Copyright (c) 2010-2014 Centre National de la Recherche Scientifique.
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

/// \internal \file sht_private.h private data and options.

#include <stdlib.h>
#include <complex.h>
#include <math.h>
// FFTW la derivee d/dx = ik	(pas de moins !)
#include "fftw3/fftw3.h"

// config file generated by ./configure
#include "sht_config.h"

#define SHTNS_PRIVATE
#include "shtns.h"

#ifdef _OPENMP
  #include <omp.h>
#endif


/* BEGIN COMPILE-TIME SETTINGS */

/// defines the maximum amount of memory in megabytes that SHTns should use.
#define SHTNS_MAX_MEMORY 2048

/// Minimum performance improve for DCT in \ref sht_auto mode. If not atained, we may switch back to gauss.
#define MIN_PERF_IMPROVE_DCT 1.05

/// Try to enforce at least this accuracy for DCT in sht_auto mode.
#define MIN_ACCURACY_DCT 1.e-8

/// The default \ref opt_polar threshold (0 disabled, 1.e-6 is aggressive, 1.e-10 is safe, 1.e-14 is VERY safe)
#define SHT_DEFAULT_POLAR_OPT 1.e-10

/// The default \ref norm used by shtns_init
#define SHT_DEFAULT_NORM ( sht_orthonormal )
//#define SHT_DEFAULT_NORM ( sht_schmidt | SHT_NO_CS_PHASE )

/// The maximum order of non-linear terms to be resolved by SH transform by default.
/// 1 : no non-linear terms. 2 : quadratic non-linear terms (default), 3 : triadic, ...
/// must be larger or equal to 1.
#define SHT_DEFAULT_NL_ORDER 1

/// minimum NLAT to consider the use of DCT acceleration.
#define SHT_MIN_NLAT_DCT 64

/// time-limit for timing individual transforms (in seconds)
#define SHT_TIME_LIMIT 0.2

/* END COMPILE-TIME SETTINGS */

// sht variants (std, ltr)
enum sht_variants { SHT_STD, SHT_LTR, SHT_M, SHT_NVAR };
// sht types (scal synth, scal analys, vect synth, ...)
enum sht_types { SHT_TYP_SSY, SHT_TYP_SAN, SHT_TYP_VSY, SHT_TYP_VAN,
	SHT_TYP_GSP, SHT_TYP_GTO, SHT_TYP_3SY, SHT_TYP_3AN, SHT_NTYP };

// sht grids
enum sht_grids { GRID_NONE, GRID_GAUSS, GRID_REGULAR, GRID_POLES };

// pointer to various function types
typedef void (*pf2l)(shtns_cfg, void*, void*, long int);
typedef void (*pf3l)(shtns_cfg, void*, void*, void*, long int);
typedef void (*pf4l)(shtns_cfg, void*, void*, void*, void*, long int);
typedef void (*pf6l)(shtns_cfg, void*, void*, void*, void*, void*, void*, long int);
typedef void (*pf2ml)(shtns_cfg, int, void*, void*, long int);
typedef void (*pf3ml)(shtns_cfg, int, void*, void*, void*, long int);
typedef void (*pf4ml)(shtns_cfg, int, void*, void*, void*, void*, long int);
typedef void (*pf6ml)(shtns_cfg, int, void*, void*, void*, void*, void*, void*, long int);

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
	unsigned short *li;			///< degree l for given mode index (size nlm) : li[lm]
	unsigned short *mi;			///< order m for given mode index (size nlm) : mi[lm]
	double *ct, *st;			///< cos(theta) and sin(theta) arrays (size nlat)
	unsigned int nspat;			///< number of real numbers that must be allocated in a spatial field.
/* END OF PUBLIC PART */

	short fftc_mode;			///< how to perform the complex fft : -1 = no fft; 0 = interleaved/native; 1 = split/transpose.
	unsigned short nthreads;	///< number of threads (openmp).
	unsigned short *tm;			///< start theta value for SH (polar optimization : near the poles the legendre polynomials go to zero for high m's)
	int k_stride_a;				///< stride in theta direction
	int m_stride_a;				///< stride in phi direction (m)
	double *wg;					///< Gauss weights for Gauss-Legendre quadrature.
	double *st_1;				///< 1/sin(theta);

	fftw_plan ifft, fft;		// plans for FFTW.
	fftw_plan ifftc, fftc;

	/* Legendre function generation arrays */
	double *alm;	// coefficient list for Legendre function recurrence (size 2*NLM)
	double *blm;	// coefficient list for modified Legendre function recurrence for analysis (size 2*NLM)
	double *l_2;	// array of size (LMAX+1) containing 1./l(l+1) for increasing integer l.

	void* ftable[SHT_NVAR][SHT_NTYP];		// pointers to transform functions.

	/* MEM matrices */
	double **ylm;		// matrix for inverse transform (synthesis)
	struct DtDp** dylm;	// theta and phi derivative of Ylm matrix
	double **zlm;		// matrix for direct transform (analysis)
	struct DtDp** dzlm;

	int ncplx_fft;			///< number of complex numbers to allocate for the fft : -1 = no fft; 0 = in-place fft (no allocation).

	/* DCT stuff */
	short mtr_dct;			///< m truncation for dct. -1 means no dct at all.
	unsigned short klim;	///< Limit to k for non-linear terms (dct)
	fftw_plan idct, dct_m0;			// (I)DCT
	double **ykm_dct;	// matrix for inverse transform (synthesis) using dct.
	struct DtDp** dykm_dct;	// theta and phi derivative of Ykm matrix.
	double *zlm_dct0;	// matrix for direct transform (analysis), only m=0
	double *dzlm_dct0;

	/* other misc informations */
	unsigned char nlorder;	// order of non-linear terms to be resolved by SH transform.
	unsigned char grid;		// store grid type.
	short norm;				// store the normalization of the Spherical Harmonics (enum \ref shtns_norm + \ref SHT_NO_CS_PHASE flag)
	unsigned fftw_plan_mode;
	double Y00_1, Y10_ct, Y11_st;
	shtns_cfg next;		// pointer to next sht_setup or NULL (records a chained list of SHT setup).
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

// value for on-the-fly transforms is lower because it allows to optimize some more (don't compute l which are not significant).
#define SHT_L_RESCALE_FLY 1000
// set to a value close to the machine accuracy, it allows to speed-up on-the-fly SHTs with very large l (lmax > SHT_L_RESCALE_FLY).
#define SHT_ACCURACY 1.0e-20
// scale factor for extended range numbers (used in on-the-fly transforms to compute recurrence)
#define SHT_SCALE_FACTOR 2.9073548971824275622e+135
//#define SHT_SCALE_FACTOR 2.0370359763344860863e+90


#if _GCC_VEC_ == 0
	#undef _GCC_VEC_
#endif

/* are there vector extensions available ? */
#if !(defined __SSE2__ || defined __MIC__ || defined __VECTOR4DOUBLE__)
	#undef _GCC_VEC_
#endif
#ifdef __INTEL_COMPILER
	#if __INTEL_COMPILER < 1400
		#undef _GCC_VEC_
		#warning "no vector extensions available ! use gcc 4+ or icc 14+ for best performance."
	#endif
#endif
#ifdef __GNUC__
	#if __GNUC__ < 4
		#undef _GCC_VEC_
		#warning "no vector extensions available ! use gcc 4+ or icc 14+ for best performance."
	#endif
#endif

#if _GCC_VEC_ && __VECTOR4DOUBLE__
	// support Blue Gene/Q QPX vectors
	#define MIN_ALIGNMENT 32
	#define VSIZE 2
	typedef complex double v2d __attribute__((aligned (16)));		// vector that contains a complex number
	typedef double s2d __attribute__((aligned (16)));		// scalar number
	#define VSIZE2 4
	#define _SIMD_NAME_ "qpx"
	typedef vector4double rnd;		// vector of 4 doubles.
	#define vall(x) vec_splats(x)
	#define vread(mem, idx) vec_lda((idx)*32, ((double*)mem))
	#define vstor(mem, idx, v) vec_sta(v, (idx)*32, ((double*)mem))
	inline static double reduce_add(rnd a) {
		a += vec_perm(a, a, vec_gpci(02301));
		a += vec_perm(a, a, vec_gpci(01032));
		return( a[0] );
	}
	inline static v2d v2d_reduce(rnd a, rnd b) {
		a = vec_perm(a, b, vec_gpci(00426)) + vec_perm(a, b, vec_gpci(01537));
		a += vec_perm(a, a, vec_gpci(02301));
		return a[0] + I*a[1];
	}
	//#define v2d_reduce(a, b) ( reduce_add(a) +I* reduce_add(b) )

	#define S2D_STORE(mem, idx, ev, od) \
		vstor(mem, idx, ev+od); \
		vstor((double*)mem + NLAT-VSIZE2 - (idx)*VSIZE2, 0, vec_perm(ev-od, ev-od, vec_gpci(03210)));
	#define S2D_CSTORE(mem, idx, er, or, ei, oi)	{	\
		rnd aa = vec_perm(ei+oi, ei+oi, vec_gpci(01032)) + (er+or); \
		rnd bb = (er + or) - vec_perm(ei+oi, ei+oi, vec_gpci(01032)); \
		vstor(mem, idx, vec_perm(bb, aa, vec_gpci(00527))); \
		vstor(((double*)mem) + (NPHI-2*im)*NLAT, idx, vec_perm(aa, bb, vec_gpci(00527))); \
		aa = vec_perm(er-or, er-or, vec_gpci(01032)) + (ei-oi); \
		bb = vec_perm(er-or, er-or, vec_gpci(01032)) - (ei-oi); \
		vstor(((double*)mem) + NLAT, -(idx+1), vec_perm(bb, aa, vec_gpci(02705))); \
		vstor(((double*)mem) + (NPHI+1-2*im)*NLAT, -(idx+1), vec_perm(aa, bb, vec_gpci(02705))); }
	// TODO: S2D_CSTORE2 has not been tested and is probably wrong...
	#define S2D_CSTORE2(mem, idx, er, or, ei, oi)	{	\
		rnd aa = vec_perm(er+or, ei+oi, vec_gpci(00415)); \
		rnd bb = vec_perm(er+or, ei+oi, vec_gpci(02637)); \
		vstor(mem, idx*2, aa); \
		vstor(mem, idx*2+1, bb); \
		aa = vec_perm(er-or, ei-oi, vec_gpci(00415)); \
		bb = vec_perm(er-or, ei-oi, vec_gpci(02637)); \
		vstor(mem, NLAT_2-1-idx*2, aa); \
		vstor(mem, NLAT_2-2-idx*2, bb); }

	#define vdup(x) (x)

	#define vlo(a) (a[0])

	#define vcplx_real(a) creal(a)
	#define vcplx_imag(a) cimag(a)

	/*inline static void* VMALLOC(size_t s) {
		void* ptr;
		posix_memalign(&ptr, MIN_ALIGNMENT, s);
		return ptr;
	}*/
	#define VMALLOC(s)	malloc(s)
	#define VFREE(s)	free(s)
#endif

#if _GCC_VEC_ && __MIC__
	// these values must be adjusted for the larger vectors of the MIC
	#undef SHT_L_RESCALE_FLY
	#undef SHT_ACCURACY
	#define SHT_L_RESCALE_FLY 1800
	#define SHT_ACCURACY 1.0e-40

	#define MIN_ALIGNMENT 64
	#define VSIZE 2
	typedef complex double v2d __attribute__((aligned (16)));		// vector that contains a complex number
	typedef double s2d __attribute__((aligned (16)));		// scalar number
	#define VSIZE2 8
	#include <immintrin.h>
	#define _SIMD_NAME_ "mic"
	typedef double rnd __attribute__ ((vector_size (VSIZE2*8)));		// vector of 8 doubles.

	typedef union { rnd i; double v[8]; } vec_rnd;
	#define vall(x) ((rnd) _mm512_set1_pd(x))
	inline static rnd vread(double *mem, int idx) {		// unaligned load.
		rnd t;
		t = (rnd)_mm512_loadunpacklo_pd( t, (mem) + (idx)*VSIZE2 );
		t = (rnd)_mm512_loadunpackhi_pd( t, (mem) + (idx)*VSIZE2 + 64 );
		return t;
	}
	#define reduce_add(a) _mm512_reduce_add_pd(a)
	#define v2d_reduce(a, b) ( _mm512_reduce_add_pd(a) +I* _mm512_reduce_add_pd(b) )
	inline static void vstor(double *mem, int idx, rnd v) {		// unaligned store.
		_mm512_packstorelo_pd((mem) + (idx)*VSIZE2, v);
		_mm512_packstorehi_pd((mem) + (idx)*VSIZE2 + 64, v);
	}

	// could be simplified with scatter
	#define S2D_STORE(mem, idx, ev, od) \
		_mm512_store_pd( (double*)mem + (idx)*VSIZE2,   ev+od); \
		_mm512_store_pd( (double*)mem + NLAT_2-VSIZE2 - (idx)*VSIZE2, \
				_mm512_castsi512_pd(_mm512_shuffle_epi32(_mm512_permute4f128_epi32(_mm512_castpd_si512(ev-od), _MM_PERM_ABCD),_MM_PERM_BADC)));

	// could be simplified with scatter
	#define S2D_CSTORE(mem, idx, er, or, ei, oi)    {       \
		rnd aa = (rnd)_mm512_castsi512_pd(_mm512_shuffle_epi32(_mm512_castpd_si512(ei+oi), _MM_PERM_BADC)); \
		rnd bb = (er + or) - aa; \
		aa += er + or; \
		_mm512_store_pd( (double*)mem + (idx)*VSIZE2, _mm512_mask_mov_pd(bb, 170, aa) ); \
		_mm512_store_pd( (double*)mem + (NPHI-VSIZE2*im)*NLAT_2 + (idx)*VSIZE2, _mm512_mask_mov_pd(aa, 170, bb) ); \
		aa = (rnd)_mm512_castsi512_pd(_mm512_shuffle_epi32( _mm512_castpd_si512(er-or), _MM_PERM_BADC )); \
		bb = aa - (ei - oi);    \
		aa += ei - oi; \
		_mm512_store_pd( (double*)mem + NLAT_2-VSIZE2 - (idx)*VSIZE2, \
				_mm512_castsi512_pd(_mm512_shuffle_epi32(_mm512_permute4f128_epi32(_mm512_castpd_si512(_mm512_mask_mov_pd(bb,170,aa)), _MM_PERM_ABCD),_MM_PERM_BADC))); \
		_mm512_store_pd( (double*)mem + (NPHI+1-2*im)*NLAT_2-VSIZE2 - (idx)*VSIZE2, \
				_mm512_castsi512_pd(_mm512_shuffle_epi32(_mm512_permute4f128_epi32(_mm512_castpd_si512(_mm512_mask_mov_pd(bb,170,aa)), _MM_PERM_ABCD),_MM_PERM_BADC))); }

	#define vdup(x) (x)

	#define vlo(a) ((vec_rnd)a).v[0]

	#define vcplx_real(a) creal(a)
	#define vcplx_imag(a) cimag(a)

	#define VMALLOC(s)	_mm_malloc(s, MIN_ALIGNMENT)
	#define VFREE(s)	_mm_free(s)
#endif


#if _GCC_VEC_ && __SSE2__
	#define MIN_ALIGNMENT 16
	#define VSIZE 2
	typedef double s2d __attribute__ ((vector_size (8*VSIZE)));		// vector that should behave like a real scalar for complex number multiplication.
	typedef double v2d __attribute__ ((vector_size (8*VSIZE)));		// vector that contains a complex number
	#ifdef __AVX__
		#define VSIZE2 4
		#include <immintrin.h>
		#define _SIMD_NAME_ "avx"
		typedef double rnd __attribute__ ((vector_size (VSIZE2*8)));		// vector of 4 doubles.
		#define vall(x) ((rnd) _mm256_set1_pd(x))
		#define vread(mem, idx) ((rnd)_mm256_loadu_pd( ((double*)mem) + (idx)*4 ))
		#define vstor(mem, idx, v) _mm256_storeu_pd( ((double*)mem) + (idx)*4 , v)
		inline static double reduce_add(rnd a) {
			v2d t = (v2d)_mm256_castpd256_pd128(a) + (v2d)_mm256_extractf128_pd(a,1);
			return _mm_cvtsd_f64(t) + _mm_cvtsd_f64(_mm_unpackhi_pd(t,t));
		}
		inline static v2d v2d_reduce(rnd a, rnd b) {
			a = _mm256_hadd_pd(a, b);
			return (v2d)_mm256_castpd256_pd128(a) + (v2d)_mm256_extractf128_pd(a,1);
		}
		#define S2D_STORE(mem, idx, ev, od) \
			_mm256_storeu_pd(((double*)mem) + (idx)*4,   ev+od); \
			((s2d*)mem)[NLAT_2-1 - (idx)*2] = _mm256_castpd256_pd128(_mm256_shuffle_pd(ev-od, ev-od, 5)); \
			((s2d*)mem)[NLAT_2-2 - (idx)*2] = _mm256_extractf128_pd(_mm256_shuffle_pd(ev-od, ev-od, 5), 1);
		#define S2D_CSTORE(mem, idx, er, or, ei, oi)	{	\
			rnd aa = (rnd)_mm256_shuffle_pd(ei+oi,ei+oi,5) + (er + or);		rnd bb = (er + or) - (rnd)_mm256_shuffle_pd(ei+oi,ei+oi,5);	\
			_mm256_storeu_pd(((double*)mem) + (idx)*4, _mm256_shuffle_pd(bb, aa, 10 )); \
			_mm256_storeu_pd(((double*)mem) + (NPHI-2*im)*NLAT + (idx)*4, _mm256_shuffle_pd(aa, bb, 10 )); \
			aa = (rnd)_mm256_shuffle_pd(er-or,er-or,5) + (ei - oi);		bb = (rnd)_mm256_shuffle_pd(er-or,er-or,5) - (ei - oi);	\
			((s2d*)mem)[NLAT_2-1 -(idx)*2] = _mm256_castpd256_pd128(_mm256_shuffle_pd(bb, aa, 10 ));	\
			((s2d*)mem)[NLAT_2-2 -(idx)*2] = _mm256_extractf128_pd(_mm256_shuffle_pd(bb, aa, 10 ), 1);	\
			((s2d*)mem)[(NPHI+1-2*im)*NLAT_2 -1 -(idx)*2] = _mm256_castpd256_pd128(_mm256_shuffle_pd(aa, bb, 10 ));	\
			((s2d*)mem)[(NPHI+1-2*im)*NLAT_2 -2 -(idx)*2] = _mm256_extractf128_pd(_mm256_shuffle_pd(aa, bb, 10 ), 1);	}
		#define S2D_CSTORE2(mem, idx, er, or, ei, oi)	{	\
			rnd aa = (rnd)_mm256_unpacklo_pd(er+or, ei+oi);	rnd bb = (rnd)_mm256_unpackhi_pd(er+or, ei+oi);	\
			((s2d*)mem)[(idx)*4]   = _mm256_castpd256_pd128(aa);	\
			((s2d*)mem)[(idx)*4+1] = _mm256_castpd256_pd128(bb);	\
			((s2d*)mem)[(idx)*4+2] = _mm256_extractf128_pd(aa, 1);	\
			((s2d*)mem)[(idx)*4+3] = _mm256_extractf128_pd(bb, 1);	\
			aa = (rnd)_mm256_unpacklo_pd(er-or, ei-oi);	bb = (rnd)_mm256_unpackhi_pd(er-or, ei-oi);	\
			((s2d*)mem)[NLAT-1-(idx)*4] = _mm256_castpd256_pd128(aa);	\
			((s2d*)mem)[NLAT-2-(idx)*4] = _mm256_castpd256_pd128(bb);	\
			((s2d*)mem)[NLAT-3-(idx)*4] = _mm256_extractf128_pd(aa, 1);	\
			((s2d*)mem)[NLAT-4-(idx)*4] = _mm256_extractf128_pd(bb, 1);	}
	#else
		#define VSIZE2 2
		typedef double rnd __attribute__ ((vector_size (VSIZE2*8)));		// vector of 2 doubles.
		#ifdef __SSE3__
			#include <pmmintrin.h>
			#define _SIMD_NAME_ "sse3"
			inline static v2d v2d_reduce(v2d a, v2d b) {
				return _mm_hadd_pd(a,b);
			}
		#else
			#include <emmintrin.h>
			#define _SIMD_NAME_ "sse2"
			inline static v2d v2d_reduce(v2d a, v2d b) {
				v2d c = _mm_unpacklo_pd(a, b);		b = _mm_unpackhi_pd(a, b);
				return b + c;
			}
		#endif
		#define reduce_add(a) ( _mm_cvtsd_f64(a) + _mm_cvtsd_f64(_mm_unpackhi_pd(a,a)) )
		#define vall(x) ((rnd) _mm_set1_pd(x))
		#define vread(mem, idx) ((s2d*)mem)[idx]
		#define vstor(mem, idx, v) ((s2d*)mem)[idx] = v
		#define S2D_STORE(mem, idx, ev, od)		((s2d*)mem)[idx] = ev+od;		((s2d*)mem)[NLAT_2-1 - (idx)] = vxchg(ev-od);
		#define S2D_CSTORE(mem, idx, er, or, ei, oi)	{	\
			rnd aa = vxchg(ei + oi) + (er + or);		rnd bb = (er + or) - vxchg(ei + oi);	\
			((s2d*)mem)[idx] = _mm_shuffle_pd(bb, aa, 2 );	\
			((s2d*)mem)[(NPHI-2*im)*NLAT_2 + (idx)] = _mm_shuffle_pd(aa, bb, 2 );	\
			aa = vxchg(er - or) + (ei - oi);		bb = vxchg(er - or) - (ei - oi);	\
			((s2d*)mem)[NLAT_2-1 -(idx)] = _mm_shuffle_pd(bb, aa, 2 );	\
			((s2d*)mem)[(NPHI+1-2*im)*NLAT_2 -1 -(idx)] = _mm_shuffle_pd(aa, bb, 2 );	}
		#define S2D_CSTORE2(mem, idx, er, or, ei, oi)	{	\
			((s2d*)mem)[(idx)*2]   = _mm_unpacklo_pd(er+or, ei+oi);	\
			((s2d*)mem)[(idx)*2+1] = _mm_unpackhi_pd(er+or, ei+oi);	\
			((s2d*)mem)[NLAT-1-(idx)*2] = _mm_unpacklo_pd(er-or, ei-oi);	\
			((s2d*)mem)[NLAT-2-(idx)*2] = _mm_unpackhi_pd(er-or, ei-oi);	}
	#endif
	#ifdef __SSE3__
		#define addi(a,b) _mm_addsub_pd(a, _mm_shuffle_pd(b,b,1))		// a + I*b
		#define subadd(a,b) _mm_addsub_pd(a, b)		// [al-bl, ah+bh]
		//#define CMUL(a,b) _mm_addsub_pd(_mm_shuffle_pd(a,a,0)*b, _mm_shuffle_pd(a,a,3)*_mm_shuffle_pd(b,b,1))
	#else
		#define addi(a,b) ( (a) + (_mm_shuffle_pd(b,b,1) * _mm_set_pd(1.0, -1.0)) )		// a + I*b		[note: _mm_set_pd(imag, real)) ]
		#define subadd(a,b) ( (a) + (b) * _mm_set_pd(1.0, -1.0) )		// [al-bl, ah+bh]
	#endif

	// build mask (-0, -0) to change sign of both hi and lo values using xorpd
	#define SIGN_MASK_2  _mm_castsi128_pd(_mm_slli_epi64(_mm_cmpeq_epi16(_mm_set1_epi64x(0), _mm_set1_epi64x(0)), 63))
	// build mask (0, -0) to change sign of hi value using xorpd (used in CFFT_TO_2REAL)
	#define SIGN_MASK_HI  _mm_unpackhi_pd(vdup(0.0), SIGN_MASK_2 )
	// build mask (-0, 0) to change sign of lo value using xorpd
	#define SIGN_MASK_LO  _mm_unpackhi_pd(SIGN_MASK_2, vdup(0.0) )

	// vset(lo, hi) takes two doubles and pack them in a vector
	#define vset(lo, hi) _mm_set_pd(hi, lo)
	// vdup(x) takes a double and duplicate it to a vector of 2 doubles.
	#define vdup(x) ((s2d)_mm_set1_pd(x))
	// vxchg(a) exchange hi and lo component of vector a
	#define vxchg(a) ((v2d)_mm_shuffle_pd(a,a,1))
	#define vlo_to_cplx(a) _mm_unpacklo_pd(a, vdup(0.0))
	#define vhi_to_cplx(a) _mm_unpackhi_pd(a, vdup(0.0))
	#define vcplx_real(a) vlo_to_dbl(a)
	#define vcplx_imag(a) vhi_to_dbl(a)
	#ifdef __clang__
		// allow to compile with clang (llvm)
		#define vlo(a) (a)[0]
		#define vlo_to_dbl(a) (a)[0]
		#define vhi_to_dbl(a) (a)[1]
	#else
		// gcc extensions
		#ifdef __AVX__
			#define vlo(a) _mm_cvtsd_f64(_mm256_castpd256_pd128(a))
		#else
			#define vlo(a) _mm_cvtsd_f64(a)
		#endif
		#define vlo_to_dbl(a) _mm_cvtsd_f64(a)
		#define vhi_to_dbl(a) _mm_cvtsd_f64(_mm_unpackhi_pd(a,a))
	#endif

	// Allocate memory aligned on 16 bytes for SSE2 (fftw_malloc works only if fftw was compiled with --enable-sse2)
	// in 64 bit systems, malloc should be 16 bytes aligned anyway.
	#define VMALLOC(s)	( (sizeof(void*) >= 8) ? malloc(s) : _mm_malloc(s, MIN_ALIGNMENT) )
	#define VFREE(s)	( (sizeof(void*) >= 8) ? free(s) : _mm_free(s) )
#endif



#ifndef _GCC_VEC_
	#define MIN_ALIGNMENT 16
	#define VSIZE 1
	#define VSIZE2 1
	#define _SIMD_NAME_ "scalar"
	typedef double s2d;
	typedef complex double v2d;
	typedef double rnd;
	#define vread(mem, idx) ((double*)mem)[idx]
	#define vstor(mem, idx, v) ((double*)mem)[idx] = v;
	#define reduce_add(a) (a)
	#define v2d_reduce(a,b) ((a) +I*(b))	
	#define vlo(a) (a)
	#define vall(x) (x)
	#define vdup(x) (x)
	#define vxchg(x) (x)
	#define addi(a,b) ((a) + I*(b))
	#define vlo_to_dbl(a) (a)
	#define vhi_to_dbl(a) (a)
	#define vcplx_real(a) creal(a)
	#define vcplx_imag(a) cimag(a)

	// allocate memory aligned for FFTW. In 64 bit systems, malloc should be 16 bytes aligned.
	#define VMALLOC(s)	( (sizeof(void*) >= 8) ? malloc(s) : fftw_malloc(s) )
	#define VFREE(s)	( (sizeof(void*) >= 8) ? free(s) : fftw_free(s) )
#endif


#define SSE __attribute__((aligned (MIN_ALIGNMENT)))

/// align pointer on MIN_ALIGNMENT (must be a power of 2)
#define PTR_ALIGN(p) ((((size_t)(p)) + (MIN_ALIGNMENT-1)) & (~((size_t)(MIN_ALIGNMENT-1))))


struct DtDp {		// theta and phi derivatives stored together.
	double t, p;
};

#define GLUE2(a,b) a##b
#define GLUE3(a,b,c) a##b##c

// verbose printing
#if SHT_VERBOSE > 1
  #define PRINT_VERB(msg) printf(msg)
#else
  #define PRINT_VERB(msg) (0)
#endif

