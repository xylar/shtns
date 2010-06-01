/*
 * Copyright (c) 2010 Centre National de la Recherche Scientifique.
 * written by Nathanael Schaeffer (CNRS, LGIT, Grenoble, France).
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
void GEN(spat_to_SH,SUFFIX)(double *Vr, complex double *Qlm SUPARG)
{
	#include "spat_to_SH.c"
}

/// Backward \b Scalar Spherical Harmonic Transform (synthesis).
void GEN(SH_to_spat,SUFFIX)(complex double *Qlm, double *Vr SUPARG)
{
	#include "SH_to_spat.c"
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
void GEN(SHsphtor_to_spat,SUFFIX)(complex double *Slm, complex double *Tlm, double *Vt, double *Vp SUPARG)
{
	#include "SHst_to_spat.c"
}

#ifndef SHT_AXISYM
/// Spheroidal only synthesis.
void GEN(SHsph_to_spat,SUFFIX)(complex double *Slm, double *Vt, double *Vp SUPARG)
{
	#include "SHs_to_spat.c"
}

/// Toroidal only synthesis.
void GEN(SHtor_to_spat,SUFFIX)(complex double *Tlm, double *Vt, double *Vp SUPARG)
{
	#include "SHt_to_spat.c"
}
#else
/// Spheroidal m=0 only synthesis (results in theta component only).
void GEN(SHsph_to_spat,SUFFIX)(complex double *Slm, double *Vt SUPARG)
{
	#include "SHs_to_spat.c"
}

/// Toroidal m=0 only synthesis (results in phi component only).
void GEN(SHtor_to_spat,SUFFIX)(complex double *Tlm, double *Vp SUPARG)
{
	#include "SHt_to_spat.c"
}
#endif


/// \b Vector Spherical Harmonics Transform (analysis) : convert a spatial vector field (theta,phi components) to its spheroidal/toroidal spherical harmonic representation.
void GEN(spat_to_SHsphtor,SUFFIX)(double *Vt, double *Vp, complex double *Slm, complex double *Tlm SUPARG)
{
	#include "spat_to_SHst.c"
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
void GEN(spat_to_SHqst,SUFFIX)(double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm SUPARG)
{
	#define SHT_3COMP
	#include "spat_to_SHqst.c"
	#undef SHT_3COMP
}

/// Backward \b 3D Vector Spherical Harmonic Transform (synthesis).
/// This is basically a shortcut to call both SH_to_spat* and SHsphtor_to spat* but may be significantly faster.
void GEN(SHqst_to_spat,SUFFIX)(complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp SUPARG)
{
	if (MTR_DCT >= 0) {
		GEN(SH_to_spat,SUFFIX)(Qlm, Vr SUPARG2);
		GEN(SHsphtor_to_spat,SUFFIX)(Slm, Tlm, Vt, Vp SUPARG2);
	} else {
		#define SHT_3COMP
		#include "SHqst_to_spat.c"
		#undef SHT_3COMP
	}
}
//@}

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
//    GEN(SH_to_spat,SUFFIX)(Qlm, Vr SUPARGF2);
//    GEN(SHsphtor_to_spat,SUFFIX)(Slm, Tlm, Vt, Vp SUPARGF2);
}

/// \ingroup fortapi
void GENF(spat_to_sphtor,SUFFIX)(double *Vt, double *Vp, complex double *Slm, complex double *Tlm SUPARGF) {
	GEN(spat_to_SHsphtor,SUFFIX)(Vt, Vp, Slm, Tlm SUPARGF2);
}

/// \ingroup fortapi
void GENF(spat_to_qst,SUFFIX)(double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm SUPARGF)
{
    GEN(spat_to_SHsphtor,SUFFIX)(Vt, Vp, Slm, Tlm SUPARGF2);
    GEN(spat_to_SH,SUFFIX)(Vr, Qlm SUPARGF2);
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
