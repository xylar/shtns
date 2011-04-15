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

#ifdef SUFFIX
  #define GEN(name,sfx) _GEN(name,sfx)
  #define _GEN(a,b) a##b
#else
  #define GEN(name,sfx) name
#endif
#ifndef SUPARG
  #define SUPARG
#endif
#ifndef SUPARG2
  #define SUPARG2
#endif

/// \name scalar transforms
//@{

/// \b Scalar Spherical Harmonics Transform (analysis) : convert a spatial scalar field to its spherical harmonic representation.
void GEN(spat_to_SH_hyb,SUFFIX)(double *Vr, complex double *Qlm SUPARG)
{
	#include "spat_to_SH.c"
}

void GEN(spat_to_SH_nodct,SUFFIX)(double *Vr, complex double *Qlm SUPARG)
{
	#define SHT_NO_DCT
	#include "spat_to_SH.c"
	#undef SHT_NO_DCT
}

void GEN(spat_to_SH_fly1,SUFFIX)(double *Vr, complex double *Qlm SUPARG)
{
	#define NWAY 1
	#include "spat_to_SH_fly.c"
	#undef NWAY
}
void GEN(spat_to_SH_fly2,SUFFIX)(double *Vr, complex double *Qlm SUPARG)
{
	#define NWAY 2
	#include "spat_to_SH_fly.c"
	#undef NWAY
}
void GEN(spat_to_SH_fly3,SUFFIX)(double *Vr, complex double *Qlm SUPARG)
{
	#define NWAY 3
	#include "spat_to_SH_fly.c"
	#undef NWAY
}
void GEN(spat_to_SH_fly4,SUFFIX)(double *Vr, complex double *Qlm SUPARG)
{
	#define NWAY 4
	#include "spat_to_SH_fly.c"
	#undef NWAY
}
void GEN(spat_to_SH_fly6,SUFFIX)(double *Vr, complex double *Qlm SUPARG)
{
	#define NWAY 6
	#include "spat_to_SH_fly.c"
	#undef NWAY
}
void GEN(spat_to_SH_fly8,SUFFIX)(double *Vr, complex double *Qlm SUPARG)
{
	#define NWAY 8
	#include "spat_to_SH_fly.c"
	#undef NWAY
}


void GEN(SH_to_spat_hyb,SUFFIX)(complex double *Qlm, double *Vr SUPARG)
{
	#include "SH_to_spat.c"
}

void GEN(SH_to_spat_nodct,SUFFIX)(complex double *Qlm, double *Vr SUPARG)
{
	#define SHT_NO_DCT
	#include "SH_to_spat.c"
	#undef SHT_NO_DCT
}

void GEN(SH_to_spat_fly1,SUFFIX)(complex double *Qlm, double *Vr SUPARG)
{
	#define NWAY 1
	#include "SH_to_spat_fly.c"
	#undef NWAY
}
void GEN(SH_to_spat_fly2,SUFFIX)(complex double *Qlm, double *Vr SUPARG)
{
	#define NWAY 2
	#include "SH_to_spat_fly.c"
	#undef NWAY
}
void GEN(SH_to_spat_fly3,SUFFIX)(complex double *Qlm, double *Vr SUPARG)
{
	#define NWAY 3
	#include "SH_to_spat_fly.c"
	#undef NWAY
}
void GEN(SH_to_spat_fly4,SUFFIX)(complex double *Qlm, double *Vr SUPARG)
{
	#define NWAY 4
	#include "SH_to_spat_fly.c"
	#undef NWAY
}
void GEN(SH_to_spat_fly6,SUFFIX)(complex double *Qlm, double *Vr SUPARG)
{
	#define NWAY 6
	#include "SH_to_spat_fly.c"
	#undef NWAY
}
void GEN(SH_to_spat_fly8,SUFFIX)(complex double *Qlm, double *Vr SUPARG)
{
	#define NWAY 8
	#include "SH_to_spat_fly.c"
	#undef NWAY
}

// function pointers.
void (*GEN(SH_to_spat_ptr,SUFFIX))(complex double*, double* SUPARG) = &GEN(SH_to_spat_hyb, SUFFIX);
void (*GEN(spat_to_SH_ptr,SUFFIX))(double*, complex double* SUPARG) = &GEN(spat_to_SH_hyb, SUFFIX);

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

void GEN(SHsphtor_to_spat_fly1,SUFFIX)(complex double *Slm, complex double *Tlm, double *Vt, double *Vp SUPARG)
{
	#define NWAY 1
	#include "SHst_to_spat_fly.c"
	#undef NWAY
}
void GEN(SHsphtor_to_spat_fly2,SUFFIX)(complex double *Slm, complex double *Tlm, double *Vt, double *Vp SUPARG)
{
	#define NWAY 2
	#include "SHst_to_spat_fly.c"
	#undef NWAY
}
void GEN(SHsphtor_to_spat_fly3,SUFFIX)(complex double *Slm, complex double *Tlm, double *Vt, double *Vp SUPARG)
{
	#define NWAY 3
	#include "SHst_to_spat_fly.c"
	#undef NWAY
}

#ifndef SHT_AXISYM
/// Spheroidal only synthesis.
void GEN(SHsph_to_spat_hyb,SUFFIX)(complex double *Slm, double *Vt, double *Vp SUPARG)
{
#ifndef SHT_SCALAR_ONLY
	#include "SHs_to_spat.c"
#endif
}

void GEN(SHsph_to_spat_fly1,SUFFIX)(complex double *Slm, double *Vt, double *Vp SUPARG)
{
	#define NWAY 1
	#include "SHs_to_spat_fly.c"
	#undef NWAY
}
void GEN(SHsph_to_spat_fly2,SUFFIX)(complex double *Slm, double *Vt, double *Vp SUPARG)
{
	#define NWAY 2
	#include "SHs_to_spat_fly.c"
	#undef NWAY
}
void GEN(SHsph_to_spat_fly4,SUFFIX)(complex double *Slm, double *Vt, double *Vp SUPARG)
{
	#define NWAY 4
	#include "SHs_to_spat_fly.c"
	#undef NWAY
}

/// Toroidal only synthesis.
void GEN(SHtor_to_spat_hyb,SUFFIX)(complex double *Tlm, double *Vt, double *Vp SUPARG)
{
#ifndef SHT_SCALAR_ONLY
	#include "SHt_to_spat.c"
#endif
}

void GEN(SHtor_to_spat_fly1,SUFFIX)(complex double *Tlm, double *Vt, double *Vp SUPARG)
{
	#define NWAY 1
	#include "SHt_to_spat_fly.c"
	#undef NWAY
}
void GEN(SHtor_to_spat_fly2,SUFFIX)(complex double *Tlm, double *Vt, double *Vp SUPARG)
{
	#define NWAY 2
	#include "SHt_to_spat_fly.c"
	#undef NWAY
}
void GEN(SHtor_to_spat_fly4,SUFFIX)(complex double *Tlm, double *Vt, double *Vp SUPARG)
{
	#define NWAY 4
	#include "SHt_to_spat_fly.c"
	#undef NWAY
}

#else
/// Spheroidal m=0 only synthesis (results in theta component only).
void GEN(SHsph_to_spat_hyb,SUFFIX)(complex double *Slm, double *Vt SUPARG)
{
#ifndef SHT_SCALAR_ONLY
	#include "SHs_to_spat.c"
#endif
}

void GEN(SHsph_to_spat_fly1,SUFFIX)(complex double *Slm, double *Vt SUPARG)
{
	#define NWAY 1
	#include "SHs_to_spat_fly.c"
	#undef NWAY
}
void GEN(SHsph_to_spat_fly2,SUFFIX)(complex double *Slm, double *Vt SUPARG)
{
	#define NWAY 2
	#include "SHs_to_spat_fly.c"
	#undef NWAY
}
void GEN(SHsph_to_spat_fly4,SUFFIX)(complex double *Slm, double *Vt SUPARG)
{
	#define NWAY 4
	#include "SHs_to_spat_fly.c"
	#undef NWAY
}

/// Toroidal m=0 only synthesis (results in phi component only).
void GEN(SHtor_to_spat_hyb,SUFFIX)(complex double *Tlm, double *Vp SUPARG)
{
#ifndef SHT_SCALAR_ONLY
	#include "SHt_to_spat.c"
#endif
}

void GEN(SHtor_to_spat_fly1,SUFFIX)(complex double *Tlm, double *Vp SUPARG)
{
	#define NWAY 1
	#include "SHt_to_spat_fly.c"
	#undef NWAY
}
void GEN(SHtor_to_spat_fly2,SUFFIX)(complex double *Tlm, double *Vp SUPARG)
{
	#define NWAY 2
	#include "SHt_to_spat_fly.c"
	#undef NWAY
}
void GEN(SHtor_to_spat_fly4,SUFFIX)(complex double *Tlm, double *Vp SUPARG)
{
	#define NWAY 4
	#include "SHt_to_spat_fly.c"
	#undef NWAY
}
#endif


/// \b Vector Spherical Harmonics Transform (analysis) : convert a spatial vector field (theta,phi components) to its spheroidal/toroidal spherical harmonic representation.
void GEN(spat_to_SHsphtor_hyb,SUFFIX)(double *Vt, double *Vp, complex double *Slm, complex double *Tlm SUPARG)
{
#ifndef SHT_SCALAR_ONLY
	#include "spat_to_SHst.c"
#endif
}

void GEN(spat_to_SHsphtor_fly1,SUFFIX)(double *Vt, double *Vp, complex double *Slm, complex double *Tlm SUPARG)
{
	#define NWAY 1
	#include "spat_to_SHst_fly.c"
	#undef NWAY
}
void GEN(spat_to_SHsphtor_fly2,SUFFIX)(double *Vt, double *Vp, complex double *Slm, complex double *Tlm SUPARG)
{
	#define NWAY 2
	#include "spat_to_SHst_fly.c"
	#undef NWAY
}
void GEN(spat_to_SHsphtor_fly3,SUFFIX)(double *Vt, double *Vp, complex double *Slm, complex double *Tlm SUPARG)
{
	#define NWAY 3
	#include "spat_to_SHst_fly.c"
	#undef NWAY
}

// function pointers.
void (*GEN(spat_to_SHsphtor_ptr,SUFFIX))(double*, double*, complex double*, complex double* SUPARG) = &GEN(spat_to_SHsphtor_hyb,SUFFIX);
void (*GEN(SHsphtor_to_spat_ptr,SUFFIX))(complex double*, complex double*, double*, double* SUPARG) = &GEN(SHsphtor_to_spat_hyb,SUFFIX);
#ifndef SHT_AXISYM
void (*GEN(SHsph_to_spat_ptr,SUFFIX))(complex double*, double*, double* SUPARG) = &GEN(SHsph_to_spat_hyb,SUFFIX);
void (*GEN(SHtor_to_spat_ptr,SUFFIX))(complex double*, double*, double* SUPARG) = &GEN(SHtor_to_spat_hyb,SUFFIX);
#else
void (*GEN(SHsph_to_spat_ptr,SUFFIX))(complex double*, double* SUPARG) = &GEN(SHsph_to_spat_hyb,SUFFIX);
void (*GEN(SHtor_to_spat_ptr,SUFFIX))(complex double*, double* SUPARG) = &GEN(SHtor_to_spat_hyb,SUFFIX);
#endif

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

void GEN(spat_to_SHqst_fly1,SUFFIX)(double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm SUPARG)
{
	#define NWAY 1
	#include "spat_to_SHqst_fly.c"
	#undef NWAY
}
void GEN(spat_to_SHqst_fly2,SUFFIX)(double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm SUPARG)
{
	#define NWAY 2
	#include "spat_to_SHqst_fly.c"
	#undef NWAY
}
void GEN(spat_to_SHqst_fly3,SUFFIX)(double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm SUPARG)
{
	#define NWAY 3
	#include "spat_to_SHqst_fly.c"
	#undef NWAY
}

/// Backward \b 3D Vector Spherical Harmonic Transform (synthesis).
/// This is basically a shortcut to call both SH_to_spat* and SHsphtor_to spat* but may be significantly faster.
void GEN(SHqst_to_spat_fly1,SUFFIX)(complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp SUPARG)
{
	#define NWAY 1
	#include "SHqst_to_spat_fly.c"
	#undef NWAY
}
void GEN(SHqst_to_spat_fly2,SUFFIX)(complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp SUPARG)
{
	#define NWAY 2
	#include "SHqst_to_spat_fly.c"
	#undef NWAY
}
void GEN(SHqst_to_spat_fly3,SUFFIX)(complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp SUPARG)
{
	#define NWAY 3
	#include "SHqst_to_spat_fly.c"
	#undef NWAY
}
#undef SHT_3COMP

void GEN(SHqst_to_spat_hyb,SUFFIX)(complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp SUPARG)
{
#ifndef SHT_SCALAR_ONLY
	if (MTR_DCT >= 0) {
		GEN(SH_to_spat,SUFFIX)(Qlm, Vr SUPARG2);
		GEN(SHsphtor_to_spat,SUFFIX)(Slm, Tlm, Vt, Vp SUPARG2);
	} else {
		#define SHT_3COMP
		#include "SHqst_to_spat.c"
		#undef SHT_3COMP
	}
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
void (*GEN(spat_to_SHqst_ptr,SUFFIX))(double*, double*, double*, complex double*, complex double*, complex double* SUPARG) = &GEN(spat_to_SHqst_hyb,SUFFIX);
void (*GEN(SHqst_to_spat_ptr,SUFFIX))(complex double*, complex double*, complex double*, double*, double*, double* SUPARG) = &GEN(SHqst_to_spat_hyb,SUFFIX);

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

/// use on-the-fly alogorithm (good guess without measuring)
void GEN(set_fly,SUFFIX)()
{
	GEN(SH_to_spat_ptr,SUFFIX) = &GEN(SH_to_spat_fly2, SUFFIX);
	GEN(SHsph_to_spat_ptr,SUFFIX) = &GEN(SHsph_to_spat_fly2,SUFFIX);
	GEN(SHtor_to_spat_ptr,SUFFIX) = &GEN(SHtor_to_spat_fly2,SUFFIX);
	GEN(SHsphtor_to_spat_ptr,SUFFIX) = &GEN(SHsphtor_to_spat_fly1,SUFFIX);
	GEN(SHqst_to_spat_ptr,SUFFIX) = &GEN(SHqst_to_spat_fly1, SUFFIX);
	if (wg != NULL) {
		GEN(spat_to_SH_ptr,SUFFIX) = &GEN(spat_to_SH_fly4, SUFFIX);
		GEN(spat_to_SHsphtor_ptr,SUFFIX) = &GEN(spat_to_SHsphtor_fly2,SUFFIX);
		GEN(spat_to_SHqst_ptr,SUFFIX) = &GEN(spat_to_SHqst_fly2, SUFFIX);
	}
}

#ifndef SUFFIX
#include "../cycle.h"

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
		printf("t(%s) = %.3g  ",name,t);
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
		printf("t(%s) = %.3g  ",name,t);
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
		printf("t(%s) = %.3g  ",name,t);
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
		printf("t(%s) = %.3g  ",name,t);
	#endif
	return t;
}

/// choose fastest between on-the-fly and gauss algorithms.
void GEN(choose_best_sht,SUFFIX)(int on_the_fly)
{
	complex double *Qlm, *Slm, *Tlm;
	double *Qh, *Sh, *Th;
	int m, i, minc, nloop;
	double t0, t0c, t, tt;
	ticks tik0, tik1;

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

	// find good nloop by requiring less than 3% difference between 2 consecutive timings.
	nloop = 10;                     // number of loops to get timings.
	m=0;
	do {
		t0 = get_time_2(nloop, "", SH_to_spat_ptr, Slm, Sh);
		t = get_time_2(nloop, "", SH_to_spat_ptr, Slm, Sh);
		double r = fabs(2.0*(t-t0)/(t+t0));
	#if SHT_VERBOSE > 1
		printf("nloop=%d, t0=%g, t=%g, r=%g, m=%d\n",nloop,t0,t,r,m);
	#endif
		if (r < 0.03) {
			m++;
		} else {
			nloop *= 3;
			m = 0;
		}
	} while((nloop<1000)&&(m < 3));
	#if SHT_VERBOSE > 1
		printf("nloop=%d\n",nloop);
	#endif

	// scalar
	#if SHT_VERBOSE > 1
		printf("finding best scalar synthesis ...\t");	fflush(stdout);
	#endif
	if (on_the_fly == 1) {
		t0c = 1e100;	t=t0c;
	} else {
		t0 = get_time_2(nloop, "hyb", &SH_to_spat_hyb, Slm, Sh);	t0c = t0/MIN_PERF_IMPROVE_DCT;		t=t0c;	// ref time.
	}
	tt = get_time_2(nloop, "fly1", &SH_to_spat_fly1, Slm, Sh);	if ((tt < t0c)&&(tt<t))	{	SH_to_spat_ptr = &SH_to_spat_fly1;	t = tt;	}
	tt = get_time_2(nloop, "fly2", &SH_to_spat_fly2, Slm, Sh);	if ((tt < t0c)&&(tt<t))	{	SH_to_spat_ptr = &SH_to_spat_fly2;	t = tt;	}
	tt = get_time_2(nloop, "fly3", &SH_to_spat_fly3, Slm, Sh);	if ((tt < t0c)&&(tt<t))	{	SH_to_spat_ptr = &SH_to_spat_fly3;	t = tt;	}
	tt = get_time_2(nloop, "fly4", &SH_to_spat_fly4, Slm, Sh);	if ((tt < t0c)&&(tt<t))	{	SH_to_spat_ptr = &SH_to_spat_fly4;	t = tt;	}
	tt = get_time_2(nloop, "fly6", &SH_to_spat_fly6, Slm, Sh);	if ((tt < t0c)&&(tt<t))	{	SH_to_spat_ptr = &SH_to_spat_fly6;	t = tt;	}
	tt = get_time_2(nloop, "fly8", &SH_to_spat_fly8, Slm, Sh);	if ((tt < t0c)&&(tt<t))	{	SH_to_spat_ptr = &SH_to_spat_fly8;	t = tt;	}
	if (wg != NULL) {
		#if SHT_VERBOSE > 1
			printf("\nfinding best scalar analysis ...\t");	fflush(stdout);
		#endif
		if (on_the_fly == 1) {
			t0c = 1e100;	t=t0c;
		} else {
			t0 = get_time_2(nloop, "hyb", &spat_to_SH_hyb, Sh, Slm);	t0c = t0/MIN_PERF_IMPROVE_DCT;		t=t0c;	//ref time.
		}
		tt = get_time_2(nloop, "fly1", &spat_to_SH_fly1, Sh, Slm);	if ((tt < t0c)&&(tt<t))	{	spat_to_SH_ptr = &spat_to_SH_fly1;	t = tt;	}
		tt = get_time_2(nloop, "fly2", &spat_to_SH_fly2, Sh, Slm);	if ((tt < t0c)&&(tt<t))	{	spat_to_SH_ptr = &spat_to_SH_fly2;	t = tt;	}
		tt = get_time_2(nloop, "fly3", &spat_to_SH_fly3, Sh, Slm);	if ((tt < t0c)&&(tt<t))	{	spat_to_SH_ptr = &spat_to_SH_fly3;	t = tt;	}
		tt = get_time_2(nloop, "fly4", &spat_to_SH_fly4, Sh, Slm);	if ((tt < t0c)&&(tt<t))	{	spat_to_SH_ptr = &spat_to_SH_fly4;	t = tt;	}
		tt = get_time_2(nloop, "fly6", &spat_to_SH_fly6, Sh, Slm);	if ((tt < t0c)&&(tt<t))	{	spat_to_SH_ptr = &spat_to_SH_fly6;	t = tt;	}
		tt = get_time_2(nloop, "fly8", &spat_to_SH_fly8, Sh, Slm);	if ((tt < t0c)&&(tt<t))	{	spat_to_SH_ptr = &spat_to_SH_fly8;	t = tt;	}
	}

	nloop = nloop/2;
	// gradients
	#if SHT_VERBOSE > 1
		printf("\nfinding best gradient synthesis ...\t");	fflush(stdout);
	#endif
	if (on_the_fly == 1) {
			t0c = 1e100;	t=t0c;
	} else {
  #ifndef SHT_AXISYM
		t0 = get_time_3(nloop, "hyb", &SHsph_to_spat_hyb, Slm, Sh, Th);	t0c = t0/MIN_PERF_IMPROVE_DCT;		t=t0c;	//ref time.
	}
	tt = get_time_3(nloop, "fly1", &SHsph_to_spat_fly1, Slm, Sh, Th);	if ((tt < t0c)&&(tt<t))	{	SHsph_to_spat_ptr = &SHsph_to_spat_fly1;	SHtor_to_spat_ptr = &SHtor_to_spat_fly1;	t = tt;	}
	tt = get_time_3(nloop, "fly2", &SHsph_to_spat_fly2, Slm, Sh, Th);	if ((tt < t0c)&&(tt<t))	{	SHsph_to_spat_ptr = &SHsph_to_spat_fly2;	SHtor_to_spat_ptr = &SHtor_to_spat_fly2;	t = tt;	}
	tt = get_time_3(nloop, "fly4", &SHsph_to_spat_fly4, Slm, Sh, Th);	if ((tt < t0c)&&(tt<t))	{	SHsph_to_spat_ptr = &SHsph_to_spat_fly4;	SHtor_to_spat_ptr = &SHtor_to_spat_fly4;	t = tt;	}
  #else
		t0 = get_time_2(nloop, "hyb", &SHsph_to_spat_hyb, Slm, Sh);	t0c = t0/MIN_PERF_IMPROVE_DCT;		t=t0c;	//ref time.
	}
	tt = get_time_2(nloop, "fly1", &SHsph_to_spat_fly1, Slm, Sh);	if ((tt < t0c)&&(tt<t))	{	SHsph_to_spat_ptr = &SHsph_to_spat_fly1;	SHtor_to_spat_ptr = &SHtor_to_spat_fly1;	t = tt;	}
	tt = get_time_2(nloop, "fly2", &SHsph_to_spat_fly2, Slm, Sh);	if ((tt < t0c)&&(tt<t))	{	SHsph_to_spat_ptr = &SHsph_to_spat_fly2;	SHtor_to_spat_ptr = &SHtor_to_spat_fly2;	t = tt;	}
	tt = get_time_2(nloop, "fly4", &SHsph_to_spat_fly4, Slm, Sh);	if ((tt < t0c)&&(tt<t))	{	SHsph_to_spat_ptr = &SHsph_to_spat_fly4;	SHtor_to_spat_ptr = &SHtor_to_spat_fly4;	t = tt;	}
  #endif

	// vector
	#if SHT_VERBOSE > 1
		printf("\nfinding best vector synthesis ...\t");	fflush(stdout);
	#endif
	if (on_the_fly == 1) {
			t0c = 1e100;	t=t0c;
	} else {
		t0 = get_time_4(nloop, "hyb", &SHsphtor_to_spat_hyb, Slm, Tlm, Sh, Th);	t0c = t0/MIN_PERF_IMPROVE_DCT;		t=t0c;	//ref time.
	}
	tt = get_time_4(nloop, "fly1", &SHsphtor_to_spat_fly1, Slm, Tlm, Sh, Th);	if ((tt < t0c)&&(tt<t))	{	SHsphtor_to_spat_ptr = &SHsphtor_to_spat_fly1;	t = tt; }
	tt = get_time_4(nloop, "fly2", &SHsphtor_to_spat_fly2, Slm, Tlm, Sh, Th);	if ((tt < t0c)&&(tt<t))	{	SHsphtor_to_spat_ptr = &SHsphtor_to_spat_fly2;	t = tt; }
	tt = get_time_4(nloop, "fly3", &SHsphtor_to_spat_fly3, Slm, Tlm, Sh, Th);	if ((tt < t0c)&&(tt<t))	{	SHsphtor_to_spat_ptr = &SHsphtor_to_spat_fly3;	t = tt; }
	if (wg != NULL) {
		#if SHT_VERBOSE > 1
			printf("\nfinding best vector analysis ...\t");	fflush(stdout);
		#endif
		if (on_the_fly == 1) {
			t0c = 1e100;	t=t0c;
		} else {
			t0 = get_time_4(nloop, "hyb", &spat_to_SHsphtor_hyb, Sh, Th, Slm, Tlm);	t0c = t0/MIN_PERF_IMPROVE_DCT;		t=t0c;	//ref time.
		}
		tt = get_time_4(nloop, "fly1", &spat_to_SHsphtor_fly1, Sh, Th, Slm, Tlm);	if ((tt < t0c)&&(tt<t))	{	spat_to_SHsphtor_ptr = &spat_to_SHsphtor_fly1;	t = tt; }
		tt = get_time_4(nloop, "fly2", &spat_to_SHsphtor_fly2, Sh, Th, Slm, Tlm);	if ((tt < t0c)&&(tt<t))	{	spat_to_SHsphtor_ptr = &spat_to_SHsphtor_fly2;	t = tt; }
		tt = get_time_4(nloop, "fly3", &spat_to_SHsphtor_fly3, Sh, Th, Slm, Tlm);	if ((tt < t0c)&&(tt<t))	{	spat_to_SHsphtor_ptr = &spat_to_SHsphtor_fly3;	t = tt; }
	}

	// 3d
	#if SHT_VERBOSE > 1
		printf("\nfinding best 3d synthesis ...\t");	fflush(stdout);
	#endif
	if (on_the_fly == 1) {
		t0c = 1e100;	t=t0c;
	} else {
		t0 = get_time_6(nloop, "hyb", &SHqst_to_spat_hyb, Qlm, Slm, Tlm, Qh, Sh, Th);	t0c = t0/MIN_PERF_IMPROVE_DCT;		t=t0c;	//ref time.
	}
	tt = get_time_6(nloop, "s+v", &SHqst_to_spat_2, Qlm, Slm, Tlm, Qh, Sh, Th);	if ((tt < t0c)&&(tt<t))	{	SHqst_to_spat_ptr = &SHqst_to_spat_2;	t = tt; }
	tt = get_time_6(nloop, "fly1", &SHqst_to_spat_fly1, Qlm, Slm, Tlm, Qh, Sh, Th);	if ((tt < t0c)&&(tt<t))	{	SHqst_to_spat_ptr = &SHqst_to_spat_fly1;	t = tt; }
	tt = get_time_6(nloop, "fly2", &SHqst_to_spat_fly2, Qlm, Slm, Tlm, Qh, Sh, Th);	if ((tt < t0c)&&(tt<t))	{	SHqst_to_spat_ptr = &SHqst_to_spat_fly2;	t = tt; }
	tt = get_time_6(nloop, "fly3", &SHqst_to_spat_fly3, Qlm, Slm, Tlm, Qh, Sh, Th);	if ((tt < t0c)&&(tt<t))	{	SHqst_to_spat_ptr = &SHqst_to_spat_fly3;	t = tt; }

	if (wg != NULL) {
		#if SHT_VERBOSE > 1
			printf("\nfinding best 3d analysis ...\t");	fflush(stdout);
		#endif
		if (on_the_fly == 1) {
			t0c = 1e100;	t=t0c;
		} else {
			t0 = get_time_6(nloop, "hyb", &spat_to_SHqst_hyb, Qh, Sh, Th, Qlm, Slm, Tlm);	t0c = t0/MIN_PERF_IMPROVE_DCT;		t=t0c;	//ref time.
		}
		tt = get_time_6(nloop, "s+v", &spat_to_SHqst_2, Qh, Sh, Th, Qlm, Slm, Tlm);	if ((tt < t0c)&&(tt<t))	{	spat_to_SHqst_ptr = &spat_to_SHqst_2;	t = tt; }
		tt = get_time_6(nloop, "fly1", &spat_to_SHqst_fly1, Qh, Sh, Th, Qlm, Slm, Tlm);	if ((tt < t0c)&&(tt<t))	{	spat_to_SHqst_ptr = &spat_to_SHqst_fly1;	t = tt; }
		tt = get_time_6(nloop, "fly2", &spat_to_SHqst_fly2, Qh, Sh, Th, Qlm, Slm, Tlm);	if ((tt < t0c)&&(tt<t))	{	spat_to_SHqst_ptr = &spat_to_SHqst_fly2;	t = tt; }
		tt = get_time_6(nloop, "fly3", &spat_to_SHqst_fly3, Qh, Sh, Th, Qlm, Slm, Tlm);	if ((tt < t0c)&&(tt<t))	{	spat_to_SHqst_ptr = &spat_to_SHqst_fly3;	t = tt; }
	}

	fftw_free(Tlm);	fftw_free(Slm);	fftw_free(Qlm);	fftw_free(Th);	fftw_free(Sh);	fftw_free(Qh);
}

#endif

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
#undef _GEN
#undef GENF
#undef _GENF
#undef SUFFIX
#undef SUPARG
#undef SUPARG2
#undef SUPARGF
#undef SUPARGF2
#undef SHT_AXISYM
#undef SHT_VAR_LTR
