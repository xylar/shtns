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
 *    written by Nathanael Schaeffer / LGIT,CNRS                    *
 ********************************************************************/

// build interfaces to spherical harmonic transforms

#define GLUE2(a,b) a##b
#define GLUE3(a,b,c) a##b##c

#ifdef SUFFIX
  #define GEN(name,sfx) GLUE2(name,sfx)
  #define GENFLY(name,nw,sfx) GLUE3(name,nw,sfx)
#else
  #define GEN(name,sfx) name
  #define GENFLY(name,nw,sfx) GLUE2(name,nw)
#endif
#ifndef SUPARG
  #define SUPARG
#endif
#ifndef SUPARG2
  #define SUPARG2
#endif

typedef void (*pf2)(void*, void* SUPARG);
#ifndef SHT_AXISYM
typedef void (*pfg)(void*, void*, void* SUPARG);
#else
typedef void (*pfg)(void*, void* SUPARG);
#endif
typedef void (*pf4)(void*, void*, void*, void* SUPARG);
typedef void (*pf6)(void*, void*, void*, void*, void*, void* SUPARG);

/// \name scalar transforms
//@{

/// \b Scalar Spherical Harmonics Transform (analysis) : convert a spatial scalar field to its spherical harmonic representation.
void GEN(spat_to_SH_hyb,SUFFIX)(double *Vr, complex double *Qlm SUPARG)
{
	#include "spat_to_SH.c"
}
void GEN(SH_to_spat_hyb,SUFFIX)(complex double *Qlm, double *Vr SUPARG)
{
	#include "SH_to_spat.c"
}

#define SHT_NO_DCT
void GEN(spat_to_SH_nodct,SUFFIX)(double *Vr, complex double *Qlm SUPARG)
{
	#include "spat_to_SH.c"
}
void GEN(SH_to_spat_nodct,SUFFIX)(complex double *Qlm, double *Vr SUPARG)
{
	#include "SH_to_spat.c"
}
#undef SHT_NO_DCT

#define NWAY 1
#include "spat_to_SH_fly.c"
#include "SH_to_spat_fly.c"
#undef NWAY
#define NWAY 2
#include "spat_to_SH_fly.c"
#include "SH_to_spat_fly.c"
#undef NWAY
#define NWAY 3
#include "spat_to_SH_fly.c"
#include "SH_to_spat_fly.c"
#undef NWAY
#define NWAY 4
#include "spat_to_SH_fly.c"
#include "SH_to_spat_fly.c"
#undef NWAY
#define NWAY 6
#include "spat_to_SH_fly.c"
#include "SH_to_spat_fly.c"
#undef NWAY
#define NWAY 8
#include "spat_to_SH_fly.c"
#include "SH_to_spat_fly.c"
#undef NWAY

// function pointers.
pf2 GEN(SH_to_spat_ptr, SUFFIX) = (pf2) &GEN(SH_to_spat_hyb, SUFFIX);
pf2 GEN(spat_to_SH_ptr,SUFFIX) = (pf2) &GEN(spat_to_SH_hyb, SUFFIX);

/// Backward \b Scalar Spherical Harmonic Transform (synthesis).
void GEN(SH_to_spat,SUFFIX)(complex double *Qlm, double *Vr SUPARG)
{
	(*GEN(SH_to_spat_ptr,SUFFIX))(Qlm, Vr SUPARG2);
	return;
}

void GEN(spat_to_SH,SUFFIX)(double *Vr, complex double *Qlm SUPARG)
{
	(*GEN(spat_to_SH_ptr,SUFFIX))(Vr, Qlm SUPARG2);
	return;
}

//@}


/** \name 2D vector transforms
 * These functions use the toroidal/spheroidal representation for vectors on a sphere \f$ (v_\theta,v_\phi) \f$:
 * \f[ v_\theta = \frac{1}{\sin\theta} \frac{\partial T}{\partial \phi} + \frac{\partial S}{\partial \theta} \f]
 * \f[ v_\phi = \frac{1}{\sin\theta} \frac{\partial S}{\partial \phi} - \frac{\partial T}{\partial \theta} \f]
 * or \f$ \mathbf{v} = \nabla \times (T \mathbf{r}) + \nabla S \f$
 * where T and S are respectively the toroidal and spheroidal scalar.
 */

//@{

/// Backward \b Vector Spherical Harmonic Transform (synthesis).
void GEN(SHsphtor_to_spat_hyb,SUFFIX)(complex double *Slm, complex double *Tlm, double *Vt, double *Vp SUPARG)
{
#ifndef SHT_SCALAR_ONLY
	#include "SHst_to_spat.c"
#endif
}

/// \b Vector Spherical Harmonics Transform (analysis) : convert a spatial vector field (theta,phi components) to its spheroidal/toroidal spherical harmonic representation.
void GEN(spat_to_SHsphtor_hyb,SUFFIX)(double *Vt, double *Vp, complex double *Slm, complex double *Tlm SUPARG)
{
#ifndef SHT_SCALAR_ONLY
	#include "spat_to_SHst.c"
#endif
}

#define NWAY 1
#include "spat_to_SHst_fly.c"
#include "SHst_to_spat_fly.c"
#undef NWAY
#define NWAY 2
#include "spat_to_SHst_fly.c"
#include "SHst_to_spat_fly.c"
#undef NWAY
#define NWAY 3
#include "spat_to_SHst_fly.c"
#include "SHst_to_spat_fly.c"
#undef NWAY


/* GRADIENTS */
#define SHT_GRAD

#ifndef SHT_AXISYM
/// Spheroidal only synthesis.
void GEN(SHsph_to_spat_hyb,SUFFIX)(complex double *Slm, double *Vt, double *Vp SUPARG)
#else
/// Spheroidal m=0 only synthesis (results in theta component only).
void GEN(SHsph_to_spat_hyb,SUFFIX)(complex double *Slm, double *Vt SUPARG)
#endif
{
#ifndef SHT_SCALAR_ONLY
	#include "SHs_to_spat.c"
#endif
}

#ifndef SHT_AXISYM
/// Toroidal only synthesis.
void GEN(SHtor_to_spat_hyb,SUFFIX)(complex double *Tlm, double *Vt, double *Vp SUPARG)
#else
/// Toroidal m=0 only synthesis (results in phi component only).
void GEN(SHtor_to_spat_hyb,SUFFIX)(complex double *Tlm, double *Vp SUPARG)
#endif
{
#ifndef SHT_SCALAR_ONLY
	#include "SHt_to_spat.c"
#endif
}

#define NWAY 1
#include "SHs_to_spat_fly.c"
#include "SHt_to_spat_fly.c"
#undef NWAY
#define NWAY 2
#include "SHs_to_spat_fly.c"
#include "SHt_to_spat_fly.c"
#undef NWAY
#define NWAY 3
#include "SHs_to_spat_fly.c"
#include "SHt_to_spat_fly.c"
#undef NWAY
#define NWAY 4
#include "SHs_to_spat_fly.c"
#include "SHt_to_spat_fly.c"
#undef NWAY
#undef SHT_GRAD


// function pointers.
pf4 GEN(spat_to_SHsphtor_ptr,SUFFIX) = (pf4) &GEN(spat_to_SHsphtor_hyb,SUFFIX);
pf4 GEN(SHsphtor_to_spat_ptr,SUFFIX) = (pf4) &GEN(SHsphtor_to_spat_hyb,SUFFIX);
pfg GEN(SHsph_to_spat_ptr,SUFFIX) = (pfg) &GEN(SHsph_to_spat_hyb,SUFFIX);
pfg GEN(SHtor_to_spat_ptr,SUFFIX) = (pfg) &GEN(SHtor_to_spat_hyb,SUFFIX);

/// Backward \b Vector Spherical Harmonic Transform (synthesis).
void GEN(SHsphtor_to_spat,SUFFIX)(complex double *Slm, complex double *Tlm, double *Vt, double *Vp SUPARG)
{
	(*GEN(SHsphtor_to_spat_ptr,SUFFIX))(Slm, Tlm, Vt, Vp SUPARG2);
	return;
}

#ifndef SHT_AXISYM
/// Spheroidal only synthesis.
void GEN(SHsph_to_spat,SUFFIX)(complex double *Slm, double *Vt, double *Vp SUPARG)
{
	(*GEN(SHsph_to_spat_ptr,SUFFIX))(Slm, Vt, Vp SUPARG2);
	return;
}

/// Toroidal only synthesis.
void GEN(SHtor_to_spat,SUFFIX)(complex double *Tlm, double *Vt, double *Vp SUPARG)
{
	(*GEN(SHtor_to_spat_ptr,SUFFIX))(Tlm, Vt, Vp SUPARG2);
	return;
}
#else
/// Spheroidal m=0 only synthesis (results in theta component only).
void GEN(SHsph_to_spat,SUFFIX)(complex double *Slm, double *Vt SUPARG)
{
	(*GEN(SHsph_to_spat_ptr,SUFFIX))(Slm, Vt SUPARG2);
	return;
}

/// Toroidal m=0 only synthesis (results in phi component only).
void GEN(SHtor_to_spat,SUFFIX)(complex double *Tlm, double *Vp SUPARG)
{
	(*GEN(SHtor_to_spat_ptr,SUFFIX))(Tlm, Vp SUPARG2);
	return;
}
#endif


/// \b Vector Spherical Harmonics Transform (analysis) : convert a spatial vector field (theta,phi components) to its spheroidal/toroidal spherical harmonic representation.
void GEN(spat_to_SHsphtor,SUFFIX)(double *Vt, double *Vp, complex double *Slm, complex double *Tlm SUPARG)
{
	(*GEN(spat_to_SHsphtor_ptr,SUFFIX))(Vt, Vp, Slm, Tlm SUPARG2);
	return;
}

//@}

/** \name 3D vector transforms
 * 3D vectors can be handled as 2D vectors + a scalar radial component.
  * For a divergenceless 3D vector \f$ (v_r,v_\theta,v_\phi) \f$, the radial scalar \f$ Q = v_r \f$ and the spheroidal scalar S can be derived from the same poloidal scalar P :
 * \f[ Q = \frac{l(l+1)}{r} P \f]
 * \f[ S = \frac{1}{r} \frac{\partial \, rP}{\partial r} \f]
 * which corresponds to the poloidal/toroidal decomposition : \f$ \mathbf{v} = \nabla \times (T \mathbf{r}) + \nabla \times \nabla \times (P \mathbf{r}) \f$
*/
//@{


/// \b 3D Vector Spherical Harmonics Transform (analysis) : convert a 3D vector field (r,theta,phi components) to its radial/spheroidal/toroidal spherical harmonic representation.
/// This is basically a shortcut to call both spat_to_SH* and spat_to_SHsphtor* but may be significantly faster.
#define SHT_3COMP
void GEN(spat_to_SHqst_hyb,SUFFIX)(double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm SUPARG)
{
#ifndef SHT_SCALAR_ONLY
	#include "spat_to_SHqst.c"
#endif
}

#define NWAY 1
#include "spat_to_SHqst_fly.c"
#include "SHqst_to_spat_fly.c"
#undef NWAY
#define NWAY 2
#include "spat_to_SHqst_fly.c"
#include "SHqst_to_spat_fly.c"
#undef NWAY
#define NWAY 3
#include "spat_to_SHqst_fly.c"
#include "SHqst_to_spat_fly.c"
#undef NWAY
#undef SHT_3COMP

/// Backward \b 3D Vector Spherical Harmonic Transform (synthesis).
/// This is basically a shortcut to call both SH_to_spat* and SHsphtor_to spat* but may be significantly faster.

void GEN(SHqst_to_spat_hyb,SUFFIX)(complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp SUPARG)
{
#ifndef SHT_SCALAR_ONLY
	#define SHT_3COMP
	#include "SHqst_to_spat.c"
	#undef SHT_3COMP
#endif
}
// combining vector and scalar.
void GEN(SHqst_to_spat_2,SUFFIX)(complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp SUPARG)
{
	GEN(SH_to_spat,SUFFIX)(Qlm, Vr SUPARG2);
	GEN(SHsphtor_to_spat,SUFFIX)(Slm, Tlm, Vt, Vp SUPARG2);
}
void GEN(spat_to_SHqst_2,SUFFIX)(double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm SUPARG)
{
	GEN(spat_to_SH,SUFFIX)(Vr, Qlm SUPARG2);
	GEN(spat_to_SHsphtor,SUFFIX)(Vt, Vp, Slm, Tlm SUPARG2);
}

// function pointers.
pf6 GEN(spat_to_SHqst_ptr,SUFFIX) = (pf6) &GEN(spat_to_SHqst_hyb,SUFFIX);
pf6 GEN(SHqst_to_spat_ptr,SUFFIX) = (pf6) &GEN(SHqst_to_spat_hyb,SUFFIX);

void GEN(spat_to_SHqst,SUFFIX)(double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm SUPARG)
{
	(*GEN(spat_to_SHqst_ptr,SUFFIX))(Vr, Vt, Vp, Qlm, Slm, Tlm SUPARG2);
	return;
}

void GEN(SHqst_to_spat,SUFFIX)(complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp SUPARG)
{
	(*GEN(SHqst_to_spat_ptr,SUFFIX))(Qlm, Slm, Tlm, Vr, Vt, Vp SUPARG2);
	return;
}

//@}

typedef struct {
  pf2 fp;  char *name;	int flags;
} frec2;
typedef struct {
  pfg fp1;	pfg fp2;	char *name;		int flags;
} frecg;
typedef struct {
  pf4 fp;	char * name;	int flags;
} frec4;
typedef struct {
  pf6 fp;	char * name;	int flags;
} frec6;

/// use on-the-fly alogorithm (good guess without measuring)
void GEN(set_fly,SUFFIX)()
{
	GEN(SH_to_spat_ptr,SUFFIX) = (pf2) &GEN(SH_to_spat_fly2, SUFFIX);
	GEN(SHsph_to_spat_ptr,SUFFIX) = (pfg) &GEN(SHsph_to_spat_fly2,SUFFIX);
	GEN(SHtor_to_spat_ptr,SUFFIX) = (pfg) &GEN(SHtor_to_spat_fly2,SUFFIX);
	GEN(SHsphtor_to_spat_ptr,SUFFIX) = (pf4) &GEN(SHsphtor_to_spat_fly1,SUFFIX);
	GEN(SHqst_to_spat_ptr,SUFFIX) = (pf6) &GEN(SHqst_to_spat_fly1, SUFFIX);
	if ( (wg != NULL) && ((shtns.norm & SHT_REAL_NORM) == 0) && (SHT_NORM != sht_schmidt) ) {
		GEN(spat_to_SH_ptr,SUFFIX) = (pf2) &GEN(spat_to_SH_fly4, SUFFIX);
		GEN(spat_to_SHsphtor_ptr,SUFFIX) = (pf4) &GEN(spat_to_SHsphtor_fly2,SUFFIX);
		GEN(spat_to_SHqst_ptr,SUFFIX) = (pf6) &GEN(spat_to_SHqst_fly2, SUFFIX);
	}
}

#include "../cycle.h"
#include <time.h> 

#if SHT_VERBOSE > 1
  #define PRINT_VERB(msg) printf(msg)
#else
  #define PRINT_VERB(msg) (0)
#endif

double GEN(get_time_2,SUFFIX)(int nloop, char* name, void (*fptr)(void*, void* SUPARG), void *a, void *b SUPARG)
{
	int i;
	ticks tik0, tik1;
	
		(*fptr)(a,b SUPARG2);		// caching...
	tik0 = getticks();
	for (i=0; i<nloop; i++) {
		(*fptr)(a,b SUPARG2);
	}
	tik1 = getticks();
	double t = elapsed(tik1,tik0)/nloop;
	#if SHT_VERBOSE > 1
		printf("  t(%s) = %.3g",name,t);
	#endif
	return t;
}

double GEN(get_time_3,SUFFIX)(int nloop, char* name, void (*fptr)(void*, void*, void* SUPARG), void *a, void *b, void *c SUPARG)
{
	int i;
	ticks tik0, tik1;
	
		(*fptr)(a,b,c SUPARG2);		// caching...
	tik0 = getticks();
	for (i=0; i<nloop; i++) {
		(*fptr)(a,b,c SUPARG2);
	}
	tik1 = getticks();
	double t = elapsed(tik1,tik0)/nloop;
	#if SHT_VERBOSE > 1
		printf("  t(%s) = %.3g",name,t);
	#endif
	return t;
}

double GEN(get_time_4,SUFFIX)(int nloop, char* name, void (*fptr)(void*, void*, void*, void* SUPARG), void *a, void *b, void *c, void *d SUPARG)
{
	int i;
	ticks tik0, tik1;
	
		(*fptr)(a,b,c,d SUPARG2);		// caching...
	tik0 = getticks();
	for (i=0; i<nloop; i++) {
		(*fptr)(a,b,c,d SUPARG2);
	}
	tik1 = getticks();
	double t = elapsed(tik1,tik0)/nloop;
	#if SHT_VERBOSE > 1
		printf("  t(%s) = %.3g",name,t);
	#endif
	return t;
}

double GEN(get_time_6,SUFFIX)(int nloop, char* name, void (*fptr)(void*, void*, void*, void*, void*, void* SUPARG), void *a, void *b, void *c, void *d, void *e, void *f SUPARG)
{
	int i;
	ticks tik0, tik1;
	
		(*fptr)(a,b,c,d,e,f SUPARG2);		// caching...
	tik0 = getticks();
	for (i=0; i<nloop; i++) {
		(*fptr)(a,b,c,d,e,f SUPARG2);
	}
	tik1 = getticks();
	double t = elapsed(tik1,tik0)/nloop;
	#if SHT_VERBOSE > 1
		printf("  t(%s) = %.3g",name,t);
	#endif
	return t;
}

/// choose fastest between on-the-fly and gauss algorithms.
/// *nlp is the number of loops. If zero, it is set to a good value.
/// on_the_fly : 1 = skip all memory algorithm. 0 = include memory and on-the-fly. -1 = test only DCT.
double GEN(choose_best_sht,SUFFIX)(int* nlp, int on_the_fly SUPARG)
{
	complex double *Qlm, *Slm, *Tlm;
	double *Qh, *Sh, *Th;
	char *nb;
	int m, i, minc, nloop;
	int dct = 0;
	int analys = 1;		// check also analysis.
	double t0, t, tt;
	double tdct, tnodct;
	ticks tik0, tik1;
	clock_t tcpu;
	char ndct[20];

	frec2 scal_synth[] = {
		{ .fp = (pf2) GEN(SH_to_spat_hyb, SUFFIX),  .name = "hyb" },
		{ .fp = (pf2) GEN(SH_to_spat_fly1, SUFFIX), .name = "fly1" },
		{ .fp = (pf2) GEN(SH_to_spat_fly2, SUFFIX), .name = "fly2" },
		{ .fp = (pf2) GEN(SH_to_spat_fly3, SUFFIX), .name = "fly3" },
		{ .fp = (pf2) GEN(SH_to_spat_fly4, SUFFIX), .name = "fly4" },
		{ .fp = (pf2) GEN(SH_to_spat_fly6, SUFFIX), .name = "fly6" },
		{ .fp = (pf2) GEN(SH_to_spat_fly8, SUFFIX), .name = "fly8" },
		{ .fp = NULL, .name = NULL }
	};
	frec2 scal_analys[] = {
		{ .fp = (pf2) GEN(spat_to_SH_hyb, SUFFIX),  .name = "hyb" },
		{ .fp = (pf2) GEN(spat_to_SH_fly1, SUFFIX), .name = "fly1" },
		{ .fp = (pf2) GEN(spat_to_SH_fly2, SUFFIX), .name = "fly2" },
		{ .fp = (pf2) GEN(spat_to_SH_fly3, SUFFIX), .name = "fly3" },
		{ .fp = (pf2) GEN(spat_to_SH_fly4, SUFFIX), .name = "fly4" },
		{ .fp = (pf2) GEN(spat_to_SH_fly6, SUFFIX), .name = "fly6" },
		{ .fp = (pf2) GEN(spat_to_SH_fly8, SUFFIX), .name = "fly8" },
		{ .fp = NULL, .name = NULL }
	};

	frecg grad_synth[] = {
		{ .fp1 = (pfg) GEN(SHsph_to_spat_hyb, SUFFIX),  .fp2 = (pfg) GEN(SHtor_to_spat_hyb, SUFFIX),  .name = "hyb" },
		{ .fp1 = (pfg) GEN(SHsph_to_spat_fly1, SUFFIX), .fp2 = (pfg) GEN(SHtor_to_spat_fly1, SUFFIX), .name = "fly1" },
		{ .fp1 = (pfg) GEN(SHsph_to_spat_fly2, SUFFIX), .fp2 = (pfg) GEN(SHtor_to_spat_fly2, SUFFIX), .name = "fly2" },
		{ .fp1 = (pfg) GEN(SHsph_to_spat_fly3, SUFFIX), .fp2 = (pfg) GEN(SHtor_to_spat_fly3, SUFFIX), .name = "fly3" },
		{ .fp1 = (pfg) GEN(SHsph_to_spat_fly4, SUFFIX), .fp2 = (pfg) GEN(SHtor_to_spat_fly4, SUFFIX), .name = "fly4" },
		{ .fp1 = NULL, .fp2 =NULL, .name = NULL }
	};
	frec4 vect_synth[] = {
		{ .fp = (pf4) GEN(SHsphtor_to_spat_hyb, SUFFIX),   .name = "hyb" },
		{ .fp = (pf4) GEN(SHsphtor_to_spat_fly1, SUFFIX),  .name = "fly1" },
		{ .fp = (pf4) GEN(SHsphtor_to_spat_fly2, SUFFIX),  .name = "fly2" },
		{ .fp = (pf4) GEN(SHsphtor_to_spat_fly3, SUFFIX),  .name = "fly3" },
		{ .fp = NULL, .name = NULL }
	};
	frec4 vect_analys[] = {
		{ .fp = (pf4) GEN(spat_to_SHsphtor_hyb, SUFFIX),   .name = "hyb" },
		{ .fp = (pf4) GEN(spat_to_SHsphtor_fly1, SUFFIX),  .name = "fly1" },
		{ .fp = (pf4) GEN(spat_to_SHsphtor_fly2, SUFFIX),  .name = "fly2" },
		{ .fp = (pf4) GEN(spat_to_SHsphtor_fly3, SUFFIX),  .name = "fly3" },
		{ .fp = NULL, .name = NULL }
	};

	frec6 v3d_synth[] = {
		{ .fp = (pf6) GEN(SHqst_to_spat_hyb, SUFFIX),   .name = "hyb" },
		{ .fp = (pf6) GEN(SHqst_to_spat_2, SUFFIX),     .name = "s+v" },
		{ .fp = (pf6) GEN(SHqst_to_spat_fly1, SUFFIX),  .name = "fly1" },
		{ .fp = (pf6) GEN(SHqst_to_spat_fly2, SUFFIX),  .name = "fly2" },
		{ .fp = (pf6) GEN(SHqst_to_spat_fly3, SUFFIX),  .name = "fly3" },
		{ .fp = NULL, .name = NULL }
	};
	frec6 v3d_analys[] = {
		{ .fp = (pf6) GEN(spat_to_SHqst_hyb, SUFFIX),   .name = "hyb" },
		{ .fp = (pf6) GEN(spat_to_SHqst_2, SUFFIX),     .name = "s+v" },
		{ .fp = (pf6) GEN(spat_to_SHqst_fly1, SUFFIX),  .name = "fly1" },
		{ .fp = (pf6) GEN(spat_to_SHqst_fly2, SUFFIX),  .name = "fly2" },
		{ .fp = (pf6) GEN(spat_to_SHqst_fly3, SUFFIX),  .name = "fly3" },
		{ .fp = NULL, .name = NULL }
	};


	if (NLAT < 32) return(0.0);		// on-the-fly not possible for NLAT_2 < 2*NWAY (overflow) and DCT not efficient for low NLAT.
	if (on_the_fly == -1) {
		on_the_fly = 0;		dct = 1;		// choose mtr_dct.
	}
	if ( (wg == NULL) || (shtns.norm & SHT_REAL_NORM) || (SHT_NORM == sht_schmidt) )	analys = 0;		// on-the-fly analysis not supported.

	Qh = (double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Sh = (double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Th = (double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Qlm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Slm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Tlm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);

	for (i=0;i<NLM;i++) {
		Slm[i] = l_2[i] + 0.5*I*l_2[i];
		Tlm[i] = 0.5*l_2[i] + I*l_2[i];
	}

	if (*nlp <= 0) {
		// find good nloop by requiring less than 3% difference between 2 consecutive timings.
		m=0;	nloop = 10;                     // number of loops to get timings.
		do {
			tcpu = clock();
			t0 = GEN(get_time_2, SUFFIX)(nloop, "", GEN(SH_to_spat_ptr,SUFFIX), Slm, Sh SUPARG2);
			t = GEN(get_time_2, SUFFIX)(nloop, "", GEN(SH_to_spat_ptr,SUFFIX), Slm, Sh SUPARG2);
			tcpu = clock() - tcpu;
			double r = fabs(2.0*(t-t0)/(t+t0));
			if (r > 0.03) {
				m = 0;		nloop *= 3;
			} else 	m++;
			tt = 1.e-6 * tcpu;		// real time should not exceed 1 sec.
			#if SHT_VERBOSE > 1
				printf(", nloop=%d, t0=%g, t=%g, r=%g, m=%d (real time = %g s)\n",nloop,t0,t,r,m,tt);
			#endif
		} while((nloop<10000)&&(m < 3)&&(tt<0.35));
		*nlp = nloop;
	} else {
		nloop = *nlp;
	}
	#if SHT_VERBOSE > 1
		printf("nloop=%d\n",nloop);
	#endif

	// scalar
	#if SHT_VERBOSE > 1
		printf("finding best scalar synthesis ...");	fflush(stdout);
	#endif
	t0 = 1e100;		i = on_the_fly - 1;		// skip hybrid if on_the_fly == 1
	while ( scal_synth[++i].fp != NULL) {
		t = GEN(get_time_2, SUFFIX)(nloop, scal_synth[i].name, scal_synth[i].fp, Slm, Sh SUPARG2);
		if (i==0) t *= 1.0/MIN_PERF_IMPROVE_DCT;
		if (t < t0) {	GEN(SH_to_spat_ptr,SUFFIX) = scal_synth[i].fp;	t0 = t;	nb=scal_synth[i].name;	PRINT_VERB("*"); }
	}
	#if SHT_VERBOSE > 1
		printf(" => %s",nb);
	#endif
	if (analys) {
		#if SHT_VERBOSE > 1
			printf("\nfinding best scalar analysis ...");	fflush(stdout);
		#endif
		t0 = 1e100;		i = on_the_fly - 1;		// skip hybrid if on_the_fly == 1
		while ( scal_analys[++i].fp != NULL) {
			t = GEN(get_time_2, SUFFIX)(nloop, scal_analys[i].name, scal_analys[i].fp, Sh, Slm SUPARG2);
			if (i==0) t *= 1.0/MIN_PERF_IMPROVE_DCT;
			if (t < t0) {	GEN(spat_to_SH_ptr,SUFFIX) = scal_analys[i].fp;	t0 = t;	nb=scal_analys[i].name;	PRINT_VERB("*"); }
		}
		#if SHT_VERBOSE > 1
			printf(" => %s",nb);
		#endif
	}

	nloop = nloop/2;
	// vector
	#if SHT_VERBOSE > 1
		printf("\nfinding best vector synthesis ...");	fflush(stdout);
	#endif
	t0 = 1e100;		i = on_the_fly - 1;		// skip hybrid if on_the_fly == 1
	while ( vect_synth[++i].fp != NULL) {
		t = GEN(get_time_4, SUFFIX)(nloop, vect_synth[i].name, vect_synth[i].fp, Slm, Tlm, Sh, Th SUPARG2);
		if (i==0) t *= 1.0/MIN_PERF_IMPROVE_DCT;
		if (t < t0) {	GEN(SHsphtor_to_spat_ptr,SUFFIX) = vect_synth[i].fp;	t0 = t;	nb=vect_synth[i].name;	PRINT_VERB("*"); }
	}
	#if SHT_VERBOSE > 1
		printf(" => %s",nb);
	#endif
	if (analys) {
		#if SHT_VERBOSE > 1
			printf("\nfinding best vector analysis ...");	fflush(stdout);
		#endif
		t0 = 1e100;		i = on_the_fly - 1;		// skip hybrid if on_the_fly == 1
		while ( vect_analys[++i].fp != NULL) {
			t = GEN(get_time_4, SUFFIX)(nloop, vect_analys[i].name, vect_analys[i].fp, Sh, Th, Slm, Tlm SUPARG2);
			if (i==0) t *= 1.0/MIN_PERF_IMPROVE_DCT;
			if (t < t0) {	GEN(spat_to_SHsphtor_ptr,SUFFIX) = vect_analys[i].fp;	t0 = t;	nb=vect_analys[i].name;	PRINT_VERB("*"); }
		}
		#if SHT_VERBOSE > 1
			printf(" => %s",nb);
		#endif
	}

  if (dct == 0) {
	// 3d
	#if SHT_VERBOSE > 1
		printf("\nfinding best 3d synthesis ...");	fflush(stdout);
	#endif
	t0 = 1e100;		i = on_the_fly - 1;		// skip hybrid if on_the_fly == 1
	while ( v3d_synth[++i].fp != NULL) {
		t = GEN(get_time_6, SUFFIX)(nloop, v3d_synth[i].name, v3d_synth[i].fp, Qlm, Slm, Tlm, Qh, Sh, Th SUPARG2);
//		if (i<2) t *= 1.0/MIN_PERF_IMPROVE_DCT;
		if (t < t0) {	GEN(SHqst_to_spat_ptr,SUFFIX) = v3d_synth[i].fp;	t0 = t;	nb=v3d_synth[i].name;	PRINT_VERB("*"); }
	}
	#if SHT_VERBOSE > 1
		printf(" => %s",nb);
	#endif
	if (analys) {
		#if SHT_VERBOSE > 1
			printf("\nfinding best 3d analysis ...");	fflush(stdout);
		#endif
		t0 = 1e100;		i = on_the_fly - 1;		// skip hybrid if on_the_fly == 1
		while ( v3d_analys[++i].fp != NULL) {
			t = GEN(get_time_6, SUFFIX)(nloop, v3d_analys[i].name, v3d_analys[i].fp, Qh, Sh, Th, Qlm, Slm, Tlm SUPARG2);
//			if (i<2) t *= 1.0/MIN_PERF_IMPROVE_DCT;
			if (t < t0) {	GEN(spat_to_SHqst_ptr,SUFFIX) = v3d_analys[i].fp;	t0 = t; nb=v3d_analys[i].name;	PRINT_VERB("*"); }
		}
		#if SHT_VERBOSE > 1
			printf(" => %s",nb);
		#endif
	}
	
	// gradients
	#if SHT_VERBOSE > 1
		printf("\nfinding best gradient synthesis ...");	fflush(stdout);
	#endif
	t0 = 1e100;		i = on_the_fly - 1;		// skip hybrid if on_the_fly == 1
	while ( grad_synth[++i].fp1 != NULL) {
	  #ifndef SHT_AXISYM
		t = GEN(get_time_3, SUFFIX)(nloop, grad_synth[i].name, grad_synth[i].fp1, Slm, Sh, Th SUPARG2);
	  #else
		t = GEN(get_time_2, SUFFIX)(nloop, grad_synth[i].name, grad_synth[i].fp1, Slm, Sh SUPARG2);
	  #endif
		if (i==0) t *= 1.0/MIN_PERF_IMPROVE_DCT;
		if (t < t0) {	GEN(SHsph_to_spat_ptr,SUFFIX) = grad_synth[i].fp1;	GEN(SHtor_to_spat_ptr,SUFFIX) = grad_synth[i].fp2;	t0 = t;	nb=grad_synth[i].name;	PRINT_VERB("*"); }
	}
	#if SHT_VERBOSE > 1
		printf(" => %s",nb);
	#endif
  }

  #ifndef SUFFIX
	// DCT
	if (dct > 0) {		// find the best DCT timings...
        minc = MMAX/20 + 1;             // don't test every single m.
	#if SHT_VERBOSE > 1
		printf("\nfinding best dct synthesis ...");
	#endif
		m = -1;		i = -1;		// reference = no dct.
			t0 = GEN(get_time_2, SUFFIX)(nloop*2, "s", SH_to_spat_ptr, Qlm, Qh SUPARG2);
			t0 += GEN(get_time_4, SUFFIX)(nloop, "v", SHsphtor_to_spat_ptr, Slm, Tlm, Sh, Th SUPARG2);
			tnodct = t0;
		for (m=0; m<=MMAX; m+=minc) {
			#if SHT_VERBOSE > 1
				printf("\nm=%d  ",m);
			#endif
			Set_MTR_DCT(m);
			t = GEN(get_time_2, SUFFIX)(nloop*2, "sdct", (pf2) &SH_to_spat_hyb, Qlm, Qh SUPARG2);
			t += GEN(get_time_4, SUFFIX)(nloop, "vdct", (pf4) &SHsphtor_to_spat_hyb, Slm, Tlm, Sh, Th SUPARG2);
			if (t < t0) {	t0 = t;		i = m;	PRINT_VERB("*"); }
		}
		tdct = t0;
		Set_MTR_DCT(i);		// the best DCT is chosen.
	}
  #endif

	#if SHT_VERBOSE > 1
		printf("\n");
	#endif
	fftw_free(Tlm);	fftw_free(Slm);	fftw_free(Qlm);	fftw_free(Th);	fftw_free(Sh);	fftw_free(Qh);

	if (dct > 0) {
		return(tdct/tnodct);
	} else	return(0.0);
}

//@}

// Fortran 77 api
#ifdef SHT_F77_API

// Fortran API : Call from fortran without the trailing '_'
//@{

#ifdef SUFFIX
  #define GENF(name,sfx) _GENF(name,sfx)
  #define _GENF(a,b) shtns_##a##b##_
#else
  #define GENF(name,sfx) _GENF(name,sfx)
  #define _GENF(a,b) shtns_##a##_
#endif

#ifndef SUPARGF
  #define SUPARGF
#endif
#ifndef SUPARGF2
  #define SUPARGF2
#endif

/// \ingroup fortapi
void GENF(spat_to_sh,SUFFIX)(double *Vr, complex double *Qlm SUPARGF) {
	GEN(spat_to_SH,SUFFIX)(Vr, Qlm SUPARGF2);
}

/// \ingroup fortapi
void GENF(sh_to_spat,SUFFIX)(complex double *Qlm, double *Vr SUPARGF) {
	GEN(SH_to_spat,SUFFIX)(Qlm, Vr SUPARGF2);
}

/// \ingroup fortapi
void GENF(sphtor_to_spat,SUFFIX)(complex double *Slm, complex double *Tlm, double *Vt, double *Vp SUPARGF) {
	GEN(SHsphtor_to_spat,SUFFIX)(Slm, Tlm, Vt, Vp SUPARGF2);
}

#ifndef SHT_AXISYM
/// \ingroup fortapi
void GENF(sph_to_spat,SUFFIX)(complex double *Slm, double *Vt, double *Vp SUPARGF) {
	GEN(SHsph_to_spat,SUFFIX)(Slm, Vt, Vp SUPARGF2);
}

/// \ingroup fortapi
void GENF(tor_to_spat,SUFFIX)(complex double *Tlm, double *Vt, double *Vp SUPARGF) {
	GEN(SHtor_to_spat,SUFFIX)(Tlm, Vt, Vp SUPARGF2);
}
#else
/// \ingroup fortapi
void GENF(sph_to_spat,SUFFIX)(complex double *Slm, double *Vt SUPARGF) {
	GEN(SHsph_to_spat,SUFFIX)(Slm, Vt SUPARGF2);
}

/// \ingroup fortapi
void GENF(tor_to_spat,SUFFIX)(complex double *Tlm, double *Vp SUPARGF) {
	GEN(SHtor_to_spat,SUFFIX)(Tlm, Vp SUPARGF2);
}
#endif

/// \ingroup fortapi
void GENF(qst_to_spat,SUFFIX)(complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp SUPARGF)
{
	GEN(SHqst_to_spat,SUFFIX)(Qlm, Slm, Tlm, Vr, Vt, Vp SUPARGF2);
}

/// \ingroup fortapi
void GENF(spat_to_sphtor,SUFFIX)(double *Vt, double *Vp, complex double *Slm, complex double *Tlm SUPARGF) {
	GEN(spat_to_SHsphtor,SUFFIX)(Vt, Vp, Slm, Tlm SUPARGF2);
}

/// \ingroup fortapi
void GENF(spat_to_qst,SUFFIX)(double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm SUPARGF)
{
	GEN(spat_to_SHqst,SUFFIX)(Vr, Vt, Vp, Qlm, Slm, Tlm SUPARGF2);
}

//@}

#endif

#undef GEN
#undef GENFLY
#undef GLUE2
#undef GLUE3
#undef GENF
#undef _GENF
#undef SUFFIX
#undef SUPARG
#undef SUPARG2
#undef SUPARGF
#undef SUPARGF2
#undef SHT_AXISYM
#undef SHT_VAR_LTR
