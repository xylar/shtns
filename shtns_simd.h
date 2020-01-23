/*
 * Copyright (c) 2010-2020 Centre National de la Recherche Scientifique.
 * written by Nathanael Schaeffer (CNRS, ISTerre, Grenoble, France).
 * 
 * nathanael.schaeffer@univ-grenoble-alpes.fr
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
 
/****************************************************************************
 * SIMD macros and functions for processor agnostic vectorization of SHTns. *
 *    A subset also adapts to various vector-length (2, 4 or 8 doubles).    *
 *    Written by Nathanael Schaeffer / CNRS                                 *
 ****************************************************************************/

/// define _GCC_VEC_ to activate SIMD using gcc extensions.
#if _GCC_VEC_ == 0
	#undef _GCC_VEC_
#endif

/* are there supported vector extensions available ? */
#if !(defined __SSE2__ || defined __VSX__ )
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

#if _GCC_VEC_ && __VSX__
	// support VSX (IBM Power)
	#include <altivec.h>
	#define MIN_ALIGNMENT 16
	#define VSIZE 2
	//typedef double s2d __attribute__ ((vector_size (8*VSIZE)));		// vector that should behave like a real scalar for complex number multiplication.
	//typedef double v2d __attribute__ ((vector_size (8*VSIZE)));		// vector that contains a complex number
	//typedef double rnd __attribute__ ((vector_size (VSIZE2*8)));		// vector of 2 doubles.
	typedef __vector double s2d;
	typedef __vector double v2d;
	typedef __vector double rnd;
	#define VSIZE2 2
	#define _SIMD_NAME_ "vsx"
	#define vall(x) vec_splats((double)x)
	inline static s2d vreverse(s2d a) {
		const vector unsigned char perm = { 8,9,10,11,12,13,14,15, 0,1,2,3,4,5,6,7 };
		return vec_perm(a,a,perm);
	}
	#define vxchg(a) vreverse(a)
	#define vread(mem, idx) vec_ld((int)(idx)*16, ((const vector double*)(mem)))
	#define vstor(mem, idx, v) vec_st((v2d)v, (int)(idx)*16, ((vector double*)(mem)))
	inline static double reduce_add(rnd a) {
		rnd b = a + vec_mergel(a, a);
		return( vec_extract(b,0) );
	}
	inline static v2d v2d_reduce(rnd a, rnd b) {
		v2d c = vec_mergel(a, b);		b = vec_mergeh(a, b);
		return b + c;
	}
	inline static v2d addi(v2d a, v2d b) {
		const s2d mp = {-1.0, 1.0};
		return a + vxchg(b)*mp;
	}
	inline static s2d subadd(s2d a, s2d b) {
		const s2d mp = {-1.0, 1.0};
		return a + b*mp;
	}

	#define vdup(x) vec_splats((double)x)
	#define vlo(a) vec_extract(a, 0)
	#define vhi(a) vec_extract(a, 1)
	#define vcplx_real(a) vec_extract(a, 0)
	#define vcplx_imag(a) vec_extract(a, 1)
	inline static s2d vec_mix_lohi(s2d a, s2d b) {	// same as _mm_shuffle_pd(a,b,2)
		const vector unsigned char perm = {0,1,2,3,4,5,6,7, 24,25,26,27,28,29,30,31};
		return vec_perm(a,b,perm);
	}

		#define S2D_STORE(mem, idx, ev, od)		((s2d*)(mem))[idx] = ev+od;		((s2d*)(mem))[NLAT_2-1 - (idx)] = vxchg(ev-od);
		#define S2D_CSTORE(mem, idx, er, od, ei, oi)	{	\
			rnd aa = vxchg(ei + oi) + (er + od);		rnd bb = (er + od) - vxchg(ei + oi);	\
			((s2d*)mem)[idx] = vec_mix_lohi(bb, aa);	\
			((s2d*)mem)[(NPHI-2*im)*NLAT_2 + (idx)] = vec_mix_lohi(aa, bb);	\
			aa = vxchg(er - od) + (ei - oi);		bb = vxchg(er - od) - (ei - oi);	\
			((s2d*)mem)[NLAT_2-1 -(idx)] = vec_mix_lohi(bb, aa);	\
			((s2d*)mem)[(NPHI+1-2*im)*NLAT_2 -1 -(idx)] = vec_mix_lohi(aa, bb);	}
		#define S2D_CSTORE2(mem, idx, er, od, ei, oi)	{	\
			((s2d*)mem)[(idx)*2]   = vec_mergeh(er+od, ei+oi);	\
			((s2d*)mem)[(idx)*2+1] = vec_mergel(er+od, ei+oi);	\
			((s2d*)mem)[NLAT-1-(idx)*2] = vec_mergeh(er-od, ei-oi);	\
			((s2d*)mem)[NLAT-2-(idx)*2] = vec_mergel(er-od, ei-oi);	}
		inline static void S2D_STORE_4MAGIC(double* mem, long idx, rnd ev, rnd od) {
			s2d n = ev+od;		s2d s = ev-od;
			((s2d*)mem)[(idx)*2] = vec_mergeh(n, s);
			((s2d*)mem)[(idx)*2+1] = vec_mergel(n, s);
		}
		inline static void S2D_CSTORE_4MAGIC(double* mem, double* mem_m, long idx, rnd er, rnd od, rnd ei, rnd oi) {
			rnd nr = er + od;		rnd ni = ei + oi;
			rnd sr = er - od;		rnd si = ei - oi;
			((s2d*)mem)[idx*2]     = vec_mergeh(nr-si,  sr+ni);
			((s2d*)mem_m)[idx*2]   = vec_mergeh(nr+si,  sr-ni);
			((s2d*)mem)[idx*2+1]   = vec_mergel(nr-si,  sr+ni);
			((s2d*)mem_m)[idx*2+1] = vec_mergel(nr+si,  sr-ni);
		}
		inline static void S2D_CSTORE2_4MAGIC(double* mem, long idx, rnd er, rnd od, rnd ei, rnd oi) {
			((s2d*)mem)[idx*4]   = vec_mergeh(er+od, ei+oi);
			((s2d*)mem)[idx*4+1] = vec_mergeh(er-od, ei-oi);
			((s2d*)mem)[idx*4+2] = vec_mergel(er+od, ei+oi);
			((s2d*)mem)[idx*4+3] = vec_mergel(er-od, ei-oi);
		}
#endif



#if _GCC_VEC_ && __SSE2__
	#define VSIZE 2
	typedef double s2d __attribute__ ((vector_size (8*VSIZE)));		// vector that should behave like a real scalar for complex number multiplication.
	typedef double v2d __attribute__ ((vector_size (8*VSIZE)));		// vector that contains a complex number
	#ifdef __AVX__
		#include <immintrin.h>
		typedef double v4d __attribute__ ((vector_size (8*4)));		// vector that contains 2 complex numbers
		#define vall4(x) ((v4d) _mm256_set1_pd(x))
		#define vread4(mem, idx) ((v4d)_mm256_loadu_pd( ((double*)(mem)) + (idx)*4 ))
		#define vstor4(mem, idx, v) _mm256_storeu_pd( ((double*)(mem)) + (idx)*4 , v)
		#define vxchg(a) ((v2d)_mm_permute_pd(a,1))
		inline static v4d v2d_x2_to_v4d(v2d a, v2d b) {
			return (v4d) _mm256_insertf128_pd( _mm256_castpd128_pd256( a ), b, 1);
		}
		#ifdef __AVX2__
			inline static v4d vreverse4(v4d a) {		// reverse vector: [0,1,2,3] => [3,2,1,0]
				return (v4d) _mm256_permute4x64_pd(a, 0x1B);
			}
		#else
			inline static v4d vreverse4(v4d a) {		// reverse vector: [0,1,2,3] => [3,2,1,0]
				a = (v4d)_mm256_shuffle_pd(a, a, 5);		// [0,1,2,3] => [1,0,3,2]
				return (v4d)_mm256_permute2f128_pd(a,a, 1);	// => [3,2,1,0]
			}
		#endif
	#else
		#define vxchg(a) ((v2d)_mm_shuffle_pd(a,a,1))
	#endif
	#ifdef __AVX512F__
		#define MIN_ALIGNMENT 64
		#define VSIZE2 8
		// Allocate memory aligned on 64 bytes for AVX-512
		#define _SIMD_NAME_ "avx512"
		typedef double rnd __attribute__ ((vector_size (VSIZE2*8)));		// vector of 8 doubles.
		typedef double s4d __attribute__ ((vector_size (4*8)));				// vector of 4 doubles.
		#define vall(x) ((rnd) _mm512_set1_pd(x))
		#define vread(mem, idx) ((rnd)_mm512_loadu_pd( ((double*)(mem)) + (idx)*8 ))
		#define vstor(mem, idx, v) _mm512_storeu_pd( ((double*)(mem)) + (idx)*8 , v)
		inline static rnd vreverse(rnd a) {		// reverse vector: [0,1,2,3,4,5,6,7] => [7,6,5,4,3,2,1,0]
			a = _mm512_permute_pd(a,0x55);	// [1,0,3,2,5,4,7,6]
			return _mm512_shuffle_f64x2(a,a,0x1B);	// [7,6,5,4,3,2,1,0]
			//return (rnd) _mm512_permutexvar_pd(_mm512_set_epi64(0,1,2,3,4,5,6,7), a);		// same speed on KNL, but requires a constant to be loaded...
		}
		#define vdup_even(v) ((rnd)_mm512_movedup_pd(v))
		#define vdup_odd(v)	 ((rnd)_mm512_permute_pd(v,0xFF))
		#define vxchg_even_odd(v) ((rnd)_mm512_permute_pd(v,0x55))
		inline static rnd vneg_even_precalc(rnd v) {		// don't use in an intesive loop.
			return _mm512_fmaddsub_pd(vall(0.0), vall(0.0), v);
		}
		inline static rnd vneg_even_inloop(rnd v) {		// this needs to load a constant and may cause a cache miss, except in a loop where the constant is kept in register.
			static const unsigned long long neg0[2] = {0x8000000000000000ULL, 0};
			//return _mm512_xor_pd(v, _mm512_broadcast_f32x4(*(_m128*)neg0));	// requires __AVX512DQ__
			return (rnd) _mm512_castsi512_pd( _mm512_xor_epi64(_mm512_castpd_si512(v), _mm512_broadcast_i32x4(*(__m128i*)neg0)) );
		}
		inline static double reduce_add(rnd a) {
			return _mm512_reduce_add_pd(a);
		}
		inline static v2d v2d_reduce(rnd a, rnd b) {
			rnd x = (rnd)_mm512_unpackhi_pd(a, b) + (rnd)_mm512_unpacklo_pd(a, b);
			s4d y = (s4d)_mm512_castpd512_pd256(x) + (s4d)_mm512_extractf64x4_pd(x,1);
			return (v2d)_mm256_castpd256_pd128(y) + (v2d)_mm256_extractf128_pd(y,1);
		}
		/* Some AVX512 vector tricks:
		 * pair-wise exchange: _mm512_permute_pd(a, 0x55)	=> [1,0,3,2,5,4,7,6]
		 * exchange middle elements of quartets: _mm512_permutex_pd(a, 0xD8)  => [0,2,1,3,4,6,5,7]
		 * exchange middle elements of quartets + swap pairs: _mm512_permutex_pd(a, 0x8D)	=> [1,3,0,2,5,7,4,6]
		 * reverse vector of pairs: _mm512_shuffle_f64x2(a,a, 0x1B)	=> [6,7,4,5,2,3,0,1]
		 */
		#define S2D_STORE(mem, idx, ev, od) \
			_mm512_storeu_pd(((double*)mem) + (idx)*8,   ev+od); \
			_mm512_storeu_pd(((double*)mem) + NLAT-8 -(idx)*8,  vreverse(ev-od));

		void inline static cstore_north_south(double* mem, double* mem_m, long idx, long nlat, rnd er, rnd od, rnd ei, rnd oi) {
			rnd ni = ei+oi;		rnd sr = er-od;		rnd nr = er+od;		rnd si = ei-oi;
			ni = (rnd)_mm512_permute_pd(ni, 0x55);		sr = (rnd)_mm512_permute_pd(sr, 0x55);		// 
			const long ridx = nlat - (idx+1)*8;
			rnd aa = ni + nr;	rnd bb = nr - ni;		rnd cc = sr + si;	rnd dd = sr - si;
			nr = _mm512_mask_blend_pd((__mmask8) 0xAA, bb, aa );		ni = _mm512_mask_blend_pd((__mmask8) 0xAA, aa, bb );
			sr = _mm512_mask_blend_pd((__mmask8) 0xAA, dd, cc );		si = _mm512_mask_blend_pd((__mmask8) 0xAA, cc, dd );
			_mm512_storeu_pd(mem + idx*8, nr);
			_mm512_storeu_pd(mem_m + idx*8, ni);
			sr = _mm512_shuffle_f64x2(sr,sr, 0x1B);			si = _mm512_shuffle_f64x2(si,si, 0x1B);
			_mm512_storeu_pd(mem + ridx, sr);
			_mm512_storeu_pd(mem_m + ridx, si);
		}

		#define S2D_CSTORE(mem, idx, er, od, ei, oi)	{	\
			rnd ni = ei+oi;		rnd sr = er-od;		rnd nr = er+od;		rnd si = ei-oi;		\
			ni = (rnd)_mm512_permute_pd(ni, 0x55);		sr = (rnd)_mm512_permute_pd(sr, 0x55);	/* pair-wise exchange: [1,0,3,2,5,4,7,6] */ \
			rnd aa = ni + nr;	rnd bb = nr - ni;		rnd cc = sr + si;	rnd dd = sr - si;	\
			nr = _mm512_mask_blend_pd((__mmask8) 0xAA, bb, aa );		ni = _mm512_mask_blend_pd((__mmask8) 0xAA, aa, bb ); \
			sr = _mm512_mask_blend_pd((__mmask8) 0xAA, dd, cc );		si = _mm512_mask_blend_pd((__mmask8) 0xAA, cc, dd ); \
			_mm512_storeu_pd(((double*)mem) + (idx)*8, nr); \
			_mm512_storeu_pd(((double*)mem) + (NPHI-2*im)*NLAT + (idx)*8, ni); \
			sr = _mm512_shuffle_f64x2(sr,sr, 0x1B);			si = _mm512_shuffle_f64x2(si,si, 0x1B);	\
			_mm512_storeu_pd(((double*)mem) + NLAT-8 -(idx)*8, sr); \
			_mm512_storeu_pd(((double*)mem) + (NPHI+1-2*im)*NLAT-8 -(idx)*8, si); }

// This may replace the first half of the following macro:
//    rnd c = _mm512_permutex2var_pd(a, _mm512_set_epi64(11,3,10,2,9,1,8,0), b);
//    rnd d = _mm512_permutex2var_pd(a, _mm512_set_epi64(15,7,14,6,13,5,12,4), b);
		#define S2D_CSTORE2(mem, idx, er, od, ei, oi)	{	\
			rnd nr = er+od;		rnd ni = ei+oi; \
			rnd sr = er-od;		rnd si = ei-oi; \
			nr = (rnd)_mm512_permutex_pd(nr, 0xD8);		ni = (rnd)_mm512_permutex_pd(ni, 0xD8);	\
			sr = (rnd)_mm512_permutex_pd(sr, 0x8D);		si = (rnd)_mm512_permutex_pd(si, 0x8D);	\
			rnd aa = (rnd)_mm512_unpacklo_pd(nr, ni);	rnd bb = (rnd)_mm512_unpackhi_pd(nr, ni);	\
			nr = _mm512_shuffle_f64x2(aa,bb, 0x44);		ni = _mm512_shuffle_f64x2(aa,bb, 0xEE); \
			rnd cc = (rnd)_mm512_unpacklo_pd(sr, si);	rnd dd = (rnd)_mm512_unpackhi_pd(sr, si);	\
			sr = _mm512_shuffle_f64x2(dd,cc, 0xEE);		si = _mm512_shuffle_f64x2(dd,cc, 0x44); \
			_mm512_storeu_pd(((double*)mem) + (idx)*16,	nr);	\
			_mm512_storeu_pd(((double*)mem) + (idx)*16+8, ni);	\
			_mm512_storeu_pd(((double*)mem) + NLAT*2-16 - (idx)*16,	sr);	\
			_mm512_storeu_pd(((double*)mem) + NLAT*2-8  - (idx)*16, si); }

		inline static void S2D_STORE_4MAGIC(double* mem, long idx, rnd ev, rnd od) {
			rnd n = ev+od;		rnd s = ev-od;
			n =  _mm512_permutex_pd(n,0xD8);	s = _mm512_permutex_pd(s,0xD8);		// [0,2,1,3, 4,6,5,7]
			rnd a = _mm512_unpacklo_pd(n, s);	rnd b = _mm512_unpackhi_pd(n, s);	// [00,11, 44,55] and [22,33, 66,77]
			_mm512_storeu_pd(mem + idx*16,    _mm512_shuffle_f64x2(a,b, 0x44) );	// [00,11, 22,33]
			_mm512_storeu_pd(mem + idx*16 +8, _mm512_shuffle_f64x2(a,b, 0xEE) );	// [44,55, 66,77]
		}

		inline static void S2D_CSTORE_4MAGIC(double* mem, double* mem_m, long idx, rnd er, rnd od, rnd ei, rnd oi) {
			rnd nr = er + od;		rnd ni = ei + oi;
			rnd sr = er - od;		rnd si = ei - oi;
			er = nr - si;			ei = sr + ni;
			od = nr + si;			oi = sr - ni;
			er = _mm512_permutex_pd(er,0xD8);	ei = _mm512_permutex_pd(ei,0xD8);	// [0,2,1,3, 4,6,5,7]
			od = _mm512_permutex_pd(od,0xD8);	oi = _mm512_permutex_pd(oi,0xD8);
			nr = _mm512_unpacklo_pd(er,  ei);		// [00,11, 44,55]
			ni = _mm512_unpackhi_pd(er,  ei);		// [22,33, 66,77]
			sr = _mm512_unpacklo_pd(od,  oi);		// [00,11, 44,55]
			si = _mm512_unpackhi_pd(od,  oi);		// [22,33, 66,77]
			_mm512_storeu_pd(mem +   idx*16,    _mm512_shuffle_f64x2(nr,ni, 0x44) );		// [00,11,22,33]
			_mm512_storeu_pd(mem +   idx*16 +8, _mm512_shuffle_f64x2(nr,ni, 0xEE) );		// [44,55,66,77]
			_mm512_storeu_pd(mem_m + idx*16,    _mm512_shuffle_f64x2(sr,si, 0x44) );
			_mm512_storeu_pd(mem_m + idx*16 +8, _mm512_shuffle_f64x2(sr,si, 0xEE) );
		}

		inline static void S2D_CSTORE2_4MAGIC(double* mem, long idx, rnd er, rnd od, rnd ei, rnd oi) {
			rnd nr = er + od;		rnd ni = ei + oi;
			rnd sr = er - od;		rnd si = ei - oi;
			rnd aa = (rnd)_mm512_unpacklo_pd(nr, sr);	rnd bb = (rnd)_mm512_unpackhi_pd(nr, sr);	// [n0,s0,n2,s2,n4,s4,n6,s6] and [n1,s1,n3,s3,n5,s5,n7,s7]
			rnd cc = (rnd)_mm512_unpacklo_pd(ni, si);	rnd dd = (rnd)_mm512_unpackhi_pd(ni, si);	// same with imaginary part
			aa = _mm512_permutex_pd(aa, 0xD8);		bb = _mm512_permutex_pd(bb, 0xD8);	// [n0,n2,s0,s2,n4,n6,s4,s6] and [n1,n3,s1,s3,n5,n7,s5,s7]
			cc = _mm512_permutex_pd(cc, 0xD8);		dd = _mm512_permutex_pd(dd, 0xD8);	// same with imaginary part
			nr = (rnd)_mm512_unpacklo_pd(aa, cc);	// [n0,s0,n4,s4] packed real/imag
			sr = (rnd)_mm512_unpackhi_pd(aa, cc);	// [n2,s2,n6,s6] packed real/imag
			ni = (rnd)_mm512_unpacklo_pd(bb, dd);	// [n1,s1,n5,s5] packed real/imag
			si = (rnd)_mm512_unpackhi_pd(bb, dd);	// [n3,s3,n7,s7] packed real/imag
			_mm512_storeu_pd(mem + idx*32,     _mm512_shuffle_f64x2(nr,ni, 0x44) ); //[n0,s0,n1,s1]
			_mm512_storeu_pd(mem + idx*32 +8,  _mm512_shuffle_f64x2(sr,si, 0x44) ); //[n2,s2,n3,s3]
			_mm512_storeu_pd(mem + idx*32 +16, _mm512_shuffle_f64x2(nr,ni, 0xEE) ); //[n4,s4,n5,s5]
			_mm512_storeu_pd(mem + idx*32 +24, _mm512_shuffle_f64x2(sr,si, 0xEE) ); //[n6,s6,n7,s7]
		}
			
	#elif defined __AVX__
		#define MIN_ALIGNMENT 32
		#define VSIZE2 4
		#ifdef __AVX2__
			#define _SIMD_NAME_ "avx2"
		#else
			#define _SIMD_NAME_ "avx"
		#endif
		typedef double rnd __attribute__ ((vector_size (VSIZE2*8)));		// vector of 4 doubles.
		#define vall(x) ((rnd) _mm256_set1_pd(x))
		#define vread(mem, idx) ((rnd)_mm256_loadu_pd( ((double*)(mem)) + (idx)*4 ))
		#define vstor(mem, idx, v) _mm256_storeu_pd( ((double*)(mem)) + (idx)*4 , v)
		#define vreverse vreverse4
		#define vdup_even(v) ((rnd)_mm256_movedup_pd(v))
		#define vdup_odd(v)  ((rnd)_mm256_unpackhi_pd(v,v))
		#define vxchg_even_odd(v) ((rnd)_mm256_shuffle_pd(v,v,0x5))
		inline static rnd vneg_even_precalc(rnd v) {		// don't use in an intesive loop.
			return _mm256_addsub_pd(vall(0.0), v);
		}
		inline static rnd vneg_even_inloop(rnd v) {		// this needs to load a constant and may cause a cache miss, except in a loop where the constant is kept in register.
			static const unsigned long long neg0[2] = {0x8000000000000000ULL, 0};
			return _mm256_xor_pd(v, _mm256_broadcast_pd((v2d*)neg0));
		}
		//#define vneg_even(v) ((rnd)_mm256_xor_pd(v, _mm256_castsi256_pd( _mm256_setr_epi32(0,0x80000000, 0,0, 0,0x80000000, 0,0))))	// BUGGY ON GCC! DON'T USE!
		inline static double reduce_add(rnd a) {
			v2d t = (v2d)_mm256_castpd256_pd128(a) + (v2d)_mm256_extractf128_pd(a,1);
			return _mm_cvtsd_f64(t) + _mm_cvtsd_f64(_mm_unpackhi_pd(t,t));
		}
		inline static v2d v2d_reduce(rnd a, rnd b) {
			a = _mm256_hadd_pd(a, b);
			return (v2d)_mm256_castpd256_pd128(a) + (v2d)_mm256_extractf128_pd(a,1);
		}
		/*inline static v4d v4d_reduce(rnd a, rnd b, rnd c, rnd d) {
			a = _mm256_hadd_pd(a, b);
			c = _mm256_hadd_pd(c, d);
			v4d x = _mm256_blend_pd(a,c,0xc);			// _mm256_setr_pd(a[0],a[1],c[2],c[3]);
			v4d y = _mm256_permute2f128_pd(a,c,0x21);	// _mm256_setr_pd(a[2],a[3],c[0],c[1]);
			return x+y;
		}*/
		#define S2D_STORE(mem, idx, ev, od) \
			_mm256_storeu_pd(((double*)mem) + (idx)*4,   ev+od); \
			_mm256_storeu_pd(((double*)mem) + NLAT-4-(idx)*4,   vreverse(ev-od));

		void static cstore_north_south(double* mem, double* mem_m, long idx, long nlat, rnd er, rnd od, rnd ei, rnd oi) {
			rnd ni = ei+oi;		rnd sr = er-od;		rnd nr = er+od;		rnd si = ei-oi;
			const long ridx = nlat - (idx+1)*4;
			ni = (rnd)_mm256_shuffle_pd(ni,ni,5);		sr = (rnd)_mm256_shuffle_pd(sr,sr,5); 
			rnd aa = ni + nr;		rnd bb = nr - ni;
			rnd cc = sr + si;		rnd dd = sr - si;
			nr = _mm256_blend_pd (bb, aa, 10 );		ni = _mm256_blend_pd (aa, bb, 10 );
			sr = _mm256_blend_pd (dd, cc, 10 );		si = _mm256_blend_pd (cc, dd, 10 );
			_mm256_storeu_pd(mem + idx*4, nr);
			_mm256_storeu_pd(mem_m + idx*4, ni);
			sr = _mm256_permute2f128_pd(sr,sr,1);	si = _mm256_permute2f128_pd(si,si,1);
			_mm256_storeu_pd(mem + ridx,  sr );
			_mm256_storeu_pd(mem_m + ridx,  si );
		}

		#define S2D_CSTORE(mem, idx, er, od, ei, oi)	{	\
			rnd ni = ei+oi;		rnd sr = er-od;		rnd nr = er+od;		rnd si = ei-oi;  \
			ni = (rnd)_mm256_shuffle_pd(ni,ni,5);		sr = (rnd)_mm256_shuffle_pd(sr,sr,5);  \
			rnd aa = ni + nr;		rnd bb = nr - ni;	\
			rnd cc = sr + si;		rnd dd = sr - si;	\
			nr = _mm256_blend_pd (bb, aa, 10 );		ni = _mm256_blend_pd (aa, bb, 10 ); \
			sr = _mm256_blend_pd (dd, cc, 10 );		si = _mm256_blend_pd (cc, dd, 10 ); \
			_mm256_storeu_pd(((double*)mem) + (idx)*4, nr); \
			_mm256_storeu_pd(((double*)mem) + (NPHI-2*im)*NLAT + (idx)*4, ni); \
			sr = _mm256_permute2f128_pd(sr,sr,1);	si = _mm256_permute2f128_pd(si,si,1); \
			_mm256_storeu_pd(((double*)mem) + NLAT-4 - (idx)*4,  sr ); \
			_mm256_storeu_pd(((double*)mem) + (NPHI+1-2*im)*NLAT-4 - (idx)*4,  si ); }
			
		#define S2D_CSTORE2(mem, idx, er, od, ei, oi)	{	\
			rnd nr = er+od;		rnd ni = ei+oi;		rnd sr = er-od;		rnd si = ei-oi;  \
			rnd aa = (rnd)_mm256_unpacklo_pd(nr, ni);	rnd bb = (rnd)_mm256_unpackhi_pd(nr, ni);	\
			nr = _mm256_permute2f128_pd(aa, bb, 0x20);		ni = _mm256_permute2f128_pd(aa, bb, 0x31);  \
			rnd cc = (rnd)_mm256_unpacklo_pd(sr, si);	rnd dd = (rnd)_mm256_unpackhi_pd(sr, si);	\
			_mm256_storeu_pd(((double*)mem) + (idx)*8,    nr); \
			_mm256_storeu_pd(((double*)mem) + (idx)*8 +4, ni); \
			sr = _mm256_permute2f128_pd(dd, cc, 0x20);		si = _mm256_permute2f128_pd(dd, cc, 0x31);  \
			_mm256_storeu_pd(((double*)mem) + NLAT*2 -4 -(idx)*8, sr); \
			_mm256_storeu_pd(((double*)mem) + NLAT*2 -8 -(idx)*8, si); }

		inline static void S2D_STORE_4MAGIC(double* mem, long idx, rnd ev, rnd od) {
			rnd n = ev+od;		rnd s = ev-od;
			rnd a = _mm256_unpacklo_pd(n, s);	rnd b = _mm256_unpackhi_pd(n, s);
			_mm256_storeu_pd(((double*)mem) + (idx)*8,    _mm256_permute2f128_pd(a,b,0x20) );
			_mm256_storeu_pd(((double*)mem) + (idx)*8 +4, _mm256_permute2f128_pd(a,b,0x31) );
		}

		inline static void S2D_CSTORE_4MAGIC(double* mem, double* mem_m, long idx, rnd er, rnd od, rnd ei, rnd oi) {
			rnd nr = er + od;		rnd ni = ei + oi;
			rnd sr = er - od;		rnd si = ei - oi;
			rnd a0 = _mm256_unpacklo_pd(nr-si,  sr+ni);
			rnd b0 = _mm256_unpacklo_pd(nr+si,  sr-ni);
			rnd a1 = _mm256_unpackhi_pd(nr-si,  sr+ni);
			rnd b1 = _mm256_unpackhi_pd(nr+si,  sr-ni);
			_mm256_storeu_pd(mem   + idx*8,    _mm256_permute2f128_pd(a0,a1,0x20) );
			_mm256_storeu_pd(mem   + idx*8 +4, _mm256_permute2f128_pd(a0,a1,0x31) );
			_mm256_storeu_pd(mem_m + idx*8,    _mm256_permute2f128_pd(b0,b1,0x20) );
			_mm256_storeu_pd(mem_m + idx*8 +4, _mm256_permute2f128_pd(b0,b1,0x31) );
		}

		inline static void S2D_CSTORE2_4MAGIC(double* mem, long idx, rnd er, rnd od, rnd ei, rnd oi) {
			rnd aa = (rnd)_mm256_unpacklo_pd(er+od, ei+oi);	rnd bb = (rnd)_mm256_unpackhi_pd(er+od, ei+oi);
			rnd cc = (rnd)_mm256_unpacklo_pd(er-od, ei-oi);	rnd dd = (rnd)_mm256_unpackhi_pd(er-od, ei-oi);
			_mm256_storeu_pd(mem + idx*16,     _mm256_permute2f128_pd(aa,cc,0x20) );
			_mm256_storeu_pd(mem + idx*16 +4,  _mm256_permute2f128_pd(bb,dd,0x20) );
			_mm256_storeu_pd(mem + idx*16 +8,  _mm256_permute2f128_pd(aa,cc,0x31) );
			_mm256_storeu_pd(mem + idx*16 +12, _mm256_permute2f128_pd(bb,dd,0x31) );
		}
	#else
		#define MIN_ALIGNMENT 16
		#define VSIZE2 2
		typedef double rnd __attribute__ ((vector_size (VSIZE2*8)));		// vector of 2 doubles.
		// Allocate memory aligned on 16 bytes for SSE2 (fftw_malloc works only if fftw was compiled with --enable-sse2)
		// in 64 bit systems, malloc should be 16 bytes aligned anyway.
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
		inline static rnd vreverse(rnd a) {	return (rnd)_mm_shuffle_pd(a,a,1);	}
		#define vdup_even(v) ((rnd)_mm_unpacklo_pd(v,v))
		#define vdup_odd(v)  ((rnd)_mm_unpackhi_pd(v,v))
		#define vxchg_even_odd(v) ((rnd)_mm_shuffle_pd(v,v,0x1))
		inline static rnd vneg_even_precalc(rnd v) {		// don't use in an intesive loop.
			return _mm_addsub_pd(vall(0.0), v);
		}
		inline static rnd vneg_even_inloop(rnd v) {		// this needs to load a constant and may cause a cache miss, except in a loop where the constant is kept in register.
			static const unsigned long long neg0[2] = {0x8000000000000000ULL, 0};
			return _mm_xor_pd(v, (v2d*)neg0);
		}
		#define reduce_add(a) ( _mm_cvtsd_f64(a) + _mm_cvtsd_f64(_mm_unpackhi_pd(a,a)) )
		#define vall(x) ((rnd) _mm_set1_pd(x))
		#define vread(mem, idx) ((s2d*)(mem))[idx]
		#define vstor(mem, idx, v) ((s2d*)(mem))[idx] = v
		#define S2D_STORE(mem, idx, ev, od)		((s2d*)mem)[idx] = ev+od;		((s2d*)mem)[NLAT_2-1 - (idx)] = vxchg(ev-od);

		void inline static cstore_north_south(double* mem, double* mem_m, long idx, long nlat, rnd er, rnd od, rnd ei, rnd oi) {
			rnd aa = vxchg(ei + oi) + (er + od);		rnd bb = (er + od) - vxchg(ei + oi);
			((s2d*)mem)[idx]   = _mm_shuffle_pd(bb, aa, 2 );
			((s2d*)mem_m)[idx] = _mm_shuffle_pd(aa, bb, 2 );
			aa = vxchg(er - od) + (ei - oi);		bb = vxchg(er - od) - (ei - oi);
			((s2d*)mem)[(nlat>>1) -1 -idx]   = _mm_shuffle_pd(bb, aa, 2 );
			((s2d*)mem_m)[(nlat>>1) -1 -idx] = _mm_shuffle_pd(aa, bb, 2 );
		}

		#define S2D_CSTORE(mem, idx, er, od, ei, oi)	{	\
			rnd aa = vxchg(ei + oi) + (er + od);		rnd bb = (er + od) - vxchg(ei + oi);	\
			((s2d*)mem)[idx] = _mm_shuffle_pd(bb, aa, 2 );	\
			((s2d*)mem)[(NPHI-2*im)*NLAT_2 + (idx)] = _mm_shuffle_pd(aa, bb, 2 );	\
			aa = vxchg(er - od) + (ei - oi);		bb = vxchg(er - od) - (ei - oi);	\
			((s2d*)mem)[NLAT_2-1 -(idx)] = _mm_shuffle_pd(bb, aa, 2 );	\
			((s2d*)mem)[(NPHI+1-2*im)*NLAT_2 -1 -(idx)] = _mm_shuffle_pd(aa, bb, 2 );	}
		#define S2D_CSTORE2(mem, idx, er, od, ei, oi)	{	\
			((s2d*)mem)[(idx)*2]   = _mm_unpacklo_pd(er+od, ei+oi);	\
			((s2d*)mem)[(idx)*2+1] = _mm_unpackhi_pd(er+od, ei+oi);	\
			((s2d*)mem)[NLAT-1-(idx)*2] = _mm_unpacklo_pd(er-od, ei-oi);	\
			((s2d*)mem)[NLAT-2-(idx)*2] = _mm_unpackhi_pd(er-od, ei-oi);	}

		inline static void S2D_STORE_4MAGIC(double* mem, long idx, rnd ev, rnd od) {
			rnd n = ev+od;		rnd s = ev-od;
			((s2d*)mem)[idx*2]   = _mm_unpacklo_pd(n, s);	// [n0,s0]
			((s2d*)mem)[idx*2+1] = _mm_unpackhi_pd(n, s);	// [n1,s1]
		}

		inline static void S2D_CSTORE_4MAGIC(double* mem, double* mem_m, long idx, rnd er, rnd od, rnd ei, rnd oi) {
			rnd nr = er + od;		rnd ni = ei + oi;
			rnd sr = er - od;		rnd si = ei - oi;
			er = nr - si;		od = sr+ni;
			ei = nr + si;		oi = sr-ni;
			nr = _mm_unpacklo_pd(er,  od);
			sr = _mm_unpackhi_pd(er,  od);
			ni = _mm_unpacklo_pd(ei,  oi);
			si = _mm_unpackhi_pd(ei,  oi);
			((s2d*)mem)[idx*2]     = nr;
			((s2d*)mem)[idx*2+1]   = sr;
			((s2d*)mem_m)[idx*2]   = ni;
			((s2d*)mem_m)[idx*2+1] = si;
		}

		inline static void S2D_CSTORE2_4MAGIC(double* mem, long idx, rnd er, rnd od, rnd ei, rnd oi) {
			((s2d*)mem)[idx*4]   = _mm_unpacklo_pd(er+od, ei+oi);	// aa = north_ri[0]
			((s2d*)mem)[idx*4+1] = _mm_unpacklo_pd(er-od, ei-oi);	// cc = south_ri[0]
			((s2d*)mem)[idx*4+2] = _mm_unpackhi_pd(er+od, ei+oi);	// bb = north_ri[1]
			((s2d*)mem)[idx*4+3] = _mm_unpackhi_pd(er-od, ei-oi);	// dd = south_ri[1]
		}
	#endif
	#ifdef __SSE3__
		#define addi(a,b) _mm_addsub_pd(a, _mm_shuffle_pd((b),(b),1))		// a + I*b
		#define subadd(a,b) _mm_addsub_pd(a, b)		// [al-bl, ah+bh]
		//#define CMUL(a,b) _mm_addsub_pd(_mm_shuffle_pd(a,a,0)*b, _mm_shuffle_pd(a,a,3)*_mm_shuffle_pd(b,b,1))
	#else
		#define addi(a,b) ( (a) + (_mm_shuffle_pd((b),(b),1) * _mm_set_pd(1.0, -1.0)) )		// a + I*b		[note: _mm_set_pd(imag, real)) ]
		#define subadd(a,b) ( (a) + (b) * _mm_set_pd(1.0, -1.0) )		// [al-bl, ah+bh]
	#endif

	// build mask (-0, -0) to change sign of both hi and lo values using xorpd
	#define SIGN_MASK_2  _mm_castsi128_pd(_mm_slli_epi64(_mm_set1_epi64x(0), _mm_set1_epi64x(0)), 63))
	// build mask (0, -0) to change sign of hi value using xorpd (used in CFFT_TO_2REAL)
	#define SIGN_MASK_HI  _mm_unpackhi_pd(vdup(0.0), SIGN_MASK_2 )
	// build mask (-0, 0) to change sign of lo value using xorpd
	#define SIGN_MASK_LO  _mm_unpackhi_pd(SIGN_MASK_2, vdup(0.0) )

	// vset(lo, hi) takes two doubles and pack them in a vector
	#define vset(lo, hi) _mm_set_pd(hi, lo)
	// vdup(x) takes a double and duplicate it to a vector of 2 doubles.
	#define vdup(x) ((s2d)_mm_set1_pd(x))
	// vxchg(a) exchange hi and lo component of vector a
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
		#ifdef __AVX512F__
			#define vlo(a) _mm_cvtsd_f64(_mm512_castpd512_pd128(a))
			#define v2d_lo(a) _mm512_castpd512_pd128(a)
		#elif defined __AVX__
			#define vlo(a) _mm_cvtsd_f64(_mm256_castpd256_pd128(a))
			#define v2d_lo(a) _mm256_castpd256_pd128(a)
		#else
			#define vlo(a) _mm_cvtsd_f64(a)
			#define v2d_lo(a) (a)
		#endif
		#define vlo_to_dbl(a) _mm_cvtsd_f64(a)
		#define vhi_to_dbl(a) _mm_cvtsd_f64(_mm_unpackhi_pd(a,a))
	#endif
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

	#define S2D_STORE(mem, idx, ev, od)		((double*)mem)[idx] = ev+od;		((double*)mem)[NLAT-1 - (idx)] = ev-od;
	#define S2D_CSTORE(mem, idx, er, od, ei, oi)	mem[idx] = (er+od) + I*(ei+oi); 	mem[NLAT-1-(idx)] = (er-od) + I*(ei-oi);
	#define S2D_CSTORE2(mem, idx, er, od, ei, oi)	mem[idx] = (er+od) + I*(ei+oi); 	mem[NLAT-1-(idx)] = (er-od) + I*(ei-oi);
	void inline static cstore_north_south(double* mem, double* mem_m, long idx, long nlat, double er, double od, double ei, double oi) {
		((cplx*)mem)[idx] = (er+od) + I*(ei+oi);
		((cplx*)mem)[nlat-1-idx] = (er-od) + I*(ei-oi);
	}

	#define S2D_CSTOREX(mem, idx, v, er, od, ei, oi) { \
		double a0 = (er[(v)]  +od[(v)])   + (ei[(v)+1]+oi[(v)+1]); \
		double b0 = (er[(v)]  +od[(v)])   - (ei[(v)+1]+oi[(v)+1]); \
		double a1 = (er[(v)+1]+od[(v)+1]) + (ei[(v)]  +oi[(v)]); \
		double b1 = (er[(v)+1]+od[(v)+1]) - (ei[(v)]  +oi[(v)]); \
		((cplx*)mem)[(idx)] = b0 + I*a1; \
		((cplx*)mem)[(NPHI-2*im)*NLAT_2 + (idx)] = a0 + I*b1; \
		a1 = (er[(v)]  -od[(v)])   + (ei[(v)+1]-oi[(v)+1]); \
		b1 = (er[(v)]  -od[(v)])   - (ei[(v)+1]-oi[(v)+1]); \
		a0 = (er[(v)+1]-od[(v)+1]) + (ei[(v)]  -oi[(v)]); \
		b0 = (er[(v)+1]-od[(v)+1]) - (ei[(v)]  -oi[(v)]); \
		((cplx*)mem)[NLAT_2-1 -(idx)] = b0 + I*a1; \
		((cplx*)mem)[(NPHI+1-2*im)*NLAT_2 -1 -(idx)] = a0 + I*b1; }

	inline static void S2D_STORE_4MAGIC(double* mem, long idx, double ev, double od) {
		((cplx*)mem)[idx] = (ev+od) + I*(ev-od);
	}

	inline static void S2D_CSTORE_4MAGIC(double* mem, double* mem_m, long idx, double er, double od, double ei, double oi) {
		double nr = er + od;		double ni = ei + oi;
		double sr = er - od;		double si = ei - oi;
		er = nr - si;		od = sr+ni;
		ei = nr + si;		oi = sr-ni;		
		((cplx*)mem)[idx]   = er + I*od;
		((cplx*)mem_m)[idx] = ei + I*oi;
	}

	inline static void S2D_CSTORE2_4MAGIC(double* mem, long idx, double er, double od, double ei, double oi) {
		((cplx*)mem)[2*idx]   = (er+od) + I*(ei+oi);
		((cplx*)mem)[2*idx+1] = (er-od) + I*(ei-oi);
	}
#endif

/// Aligned malloc on 64 bytes (cache-line) that fits any vector size up to 512 bits.
#if _GCC_VEC_ && __SSE2__
	#define VMALLOC(s)	_mm_malloc(s, 64)
	#define VFREE(s)	_mm_free(s)
#else
	inline static void* VMALLOC(size_t s) {
		void* ptr = 0;		// return value will be zero on failure.
		posix_memalign(&ptr, 64, s);
		return ptr;
	}
	#define VFREE(s)	free(s)
#endif

#define SSE __attribute__((aligned (MIN_ALIGNMENT)))

/// align pointer on MIN_ALIGNMENT (must be a power of 2)
#define PTR_ALIGN(p) ((((size_t)(p)) + (MIN_ALIGNMENT-1)) & (~((size_t)(MIN_ALIGNMENT-1))))

