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
 *    written by Nathanael Schaeffer / LGIT,CNRS                    *
 ********************************************************************/

// build interfaces to spherical harmonic transforms

#ifdef SUFFIX
  #define GEN(name,sfx) GLUE2(name,sfx)
  #define GEN3(name,nw,sfx) GLUE3(name,nw,sfx)
#else
  #define GEN(name,sfx) name
  #define GEN3(name,nw,sfx) GLUE2(name,nw)
#endif
#ifndef SUPARG
  #define SUPARG
#endif
#ifndef SUPARG2
  #define SUPARG2 , shtns->lmax
#endif

#ifdef _OPENMP
  #ifndef SHT_AXISYM
    #define ADD_OPENMP
  #endif
#endif


/// \name scalar transforms
//@{

/* Scalar Spherical Harmonics Transform : from spherical harmonic representation to a spatial grid (and reverse) */
#define ID_NME hyb
#include "spat_to_SH.c"
#include "SH_to_spat.c"
#include "SHst_to_spat.c"
#include "spat_to_SHst.c"
#undef ID_NME

// fly are compiled only once, with SHT_VAR_LTR
#ifdef SHT_VAR_LTR
	#define NWAY 1
	#include "spat_to_SHst_fly.c"
	#include "SHst_to_spat_fly.c"
	#undef NWAY
	#define NWAY 2
	#include "spat_to_SH_fly.c"
	#include "SH_to_spat_fly.c"
	#include "spat_to_SHst_fly.c"
	#include "SHst_to_spat_fly.c"
	#undef NWAY
	#define NWAY 3
	#include "spat_to_SH_fly.c"
	#include "SH_to_spat_fly.c"
	#include "spat_to_SHst_fly.c"
	#include "SHst_to_spat_fly.c"
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
  #ifdef ADD_OPENMP
	#define NWAY 1
	#include "spat_to_SHst_omp.c"
	#include "SHst_to_spat_omp.c"
	#undef NWAY
	#define NWAY 2
	#include "spat_to_SH_omp.c"
	#include "SH_to_spat_omp.c"
	#include "spat_to_SHst_omp.c"
	#include "SHst_to_spat_omp.c"
	#undef NWAY
	#define NWAY 3
	#include "spat_to_SH_omp.c"
	#include "SH_to_spat_omp.c"
	#include "spat_to_SHst_omp.c"
	#include "SHst_to_spat_omp.c"
	#undef NWAY
	#define NWAY 4
	#include "spat_to_SH_omp.c"
	#include "SH_to_spat_omp.c"
	#undef NWAY
	#define NWAY 6
	#include "spat_to_SH_omp.c"
	#include "SH_to_spat_omp.c"
	#undef NWAY
	#define NWAY 8
	#include "spat_to_SH_omp.c"
	#include "SH_to_spat_omp.c"
	#undef NWAY
  #endif
#endif

#ifdef IVAR
/// Backward \b Scalar Spherical Harmonic Transform (synthesis).
void GEN(SH_to_spat,SUFFIX)(shtns_cfg shtns, complex double *Qlm, double *Vr SUPARG)
{
	((pf2l) shtns->fptr[IVAR][SHT_TYP_SSY])(shtns, Qlm, Vr SUPARG2);
	return;
}

void GEN(spat_to_SH,SUFFIX)(shtns_cfg shtns, double *Vr, complex double *Qlm SUPARG)
{
	((pf2l) shtns->fptr[IVAR][SHT_TYP_SAN])(shtns, Vr, Qlm SUPARG2);
	return;
}

/// Backward \b Vector Spherical Harmonic Transform (synthesis).
void GEN(SHsphtor_to_spat,SUFFIX)(shtns_cfg shtns, complex double *Slm, complex double *Tlm, double *Vt, double *Vp SUPARG)
{
	((pf4l) shtns->fptr[IVAR][SHT_TYP_VSY])(shtns, Slm, Tlm, Vt, Vp SUPARG2);
	return;
}

/// \b Vector Spherical Harmonics Transform (analysis) : convert a spatial vector field (theta,phi components) to its spheroidal/toroidal spherical harmonic representation.
void GEN(spat_to_SHsphtor,SUFFIX)(shtns_cfg shtns, double *Vt, double *Vp, complex double *Slm, complex double *Tlm SUPARG)
{
	((pf4l) shtns->fptr[IVAR][SHT_TYP_VAN])(shtns, Vt, Vp, Slm, Tlm SUPARG2);
	return;
}
#endif

//@}


/* GRADIENTS */
#define SHT_GRAD

#define ID_NME hyb
#include "SHs_to_spat.c"
#include "SHt_to_spat.c"
#undef ID_NME

// fly are compiled only once, with SHT_VAR_LTR
#ifdef SHT_VAR_LTR
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
  #ifdef ADD_OPENMP
	#define NWAY 1
	#include "SHs_to_spat_omp.c"
	#include "SHt_to_spat_omp.c"
	#undef NWAY
	#define NWAY 2
	#include "SHs_to_spat_omp.c"
	#include "SHt_to_spat_omp.c"
	#undef NWAY
	#define NWAY 3
	#include "SHs_to_spat_omp.c"
	#include "SHt_to_spat_omp.c"
	#undef NWAY
	#define NWAY 4
	#include "SHs_to_spat_omp.c"
	#include "SHt_to_spat_omp.c"
	#undef NWAY
  #endif
#endif

#undef SHT_GRAD

#ifdef IVAR
#ifndef SHT_AXISYM
/// Spheroidal only synthesis.
void GEN(SHsph_to_spat,SUFFIX)(shtns_cfg shtns, complex double *Slm, double *Vt, double *Vp SUPARG)
{
	((pf3l) shtns->fptr[IVAR][SHT_TYP_GSP])(shtns, Slm, Vt, Vp SUPARG2);
	return;
}

/// Toroidal only synthesis.
void GEN(SHtor_to_spat,SUFFIX)(shtns_cfg shtns, complex double *Tlm, double *Vt, double *Vp SUPARG)
{
	((pf3l) shtns->fptr[IVAR][SHT_TYP_GTO])(shtns, Tlm, Vt, Vp SUPARG2);
	return;
}
#else
/// Spheroidal m=0 only synthesis (results in theta component only).
void GEN(SHsph_to_spat,SUFFIX)(shtns_cfg shtns, complex double *Slm, double *Vt SUPARG)
{
	((pf3l) shtns->fptr[IVAR][SHT_TYP_GSP])(shtns, Slm, Vt, NULL SUPARG2);
	return;
}

/// Toroidal m=0 only synthesis (results in phi component only).
void GEN(SHtor_to_spat,SUFFIX)(shtns_cfg shtns, complex double *Tlm, double *Vp SUPARG)
{
	((pf3l) shtns->fptr[IVAR][SHT_TYP_GTO])(shtns, Tlm, NULL, Vp SUPARG2);
	return;
}
#endif
#endif


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

#define ID_NME hyb
#include "spat_to_SHqst.c"
// hybrid SHqst_to_spat possible only of axisymmetric transform
  #ifdef SHT_AXISYM
  #include "SHqst_to_spat.c"
  #endif
#undef ID_NME

// fly are compiled only once, with SHT_VAR_LTR
#ifdef SHT_VAR_LTR
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
  #ifdef ADD_OPENMP
	#define NWAY 1
	#include "spat_to_SHqst_omp.c"
	#include "SHqst_to_spat_omp.c"
	#undef NWAY
	#define NWAY 2
	#include "spat_to_SHqst_omp.c"
	#include "SHqst_to_spat_omp.c"
	#undef NWAY
	#define NWAY 3
	#include "spat_to_SHqst_omp.c"
	#include "SHqst_to_spat_omp.c"
	#undef NWAY
  #endif
#endif
#undef SHT_3COMP

#ifdef IVAR
// combining vector and scalar.
void GEN(SHqst_to_spat_2,SUFFIX)(shtns_cfg shtns, complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp SUPARG)
{
	((pf2l) shtns->fptr[IVAR][SHT_TYP_SSY])(shtns, Qlm, Vr SUPARG2);
	((pf4l) shtns->fptr[IVAR][SHT_TYP_VSY])(shtns, Slm, Tlm, Vt, Vp SUPARG2);
}
void GEN(spat_to_SHqst_2,SUFFIX)(shtns_cfg shtns, double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm SUPARG)
{
	((pf2l) shtns->fptr[IVAR][SHT_TYP_SAN])(shtns, Vr, Qlm SUPARG2);
	((pf4l) shtns->fptr[IVAR][SHT_TYP_VAN])(shtns, Vt, Vp, Slm, Tlm SUPARG2);
}

void GEN(spat_to_SHqst,SUFFIX)(shtns_cfg shtns, double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm SUPARG)
{
	((pf6l) shtns->fptr[IVAR][SHT_TYP_3AN])(shtns, Vr, Vt, Vp, Qlm, Slm, Tlm SUPARG2);
	return;
}

void GEN(SHqst_to_spat,SUFFIX)(shtns_cfg shtns, complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp SUPARG)
{
	((pf6l) shtns->fptr[IVAR][SHT_TYP_3SY])(shtns, Qlm, Slm, Tlm, Vr, Vt, Vp SUPARG2);
	return;
}
#endif


/* functions without dct, at the end for SHT_NO_DCT must no interfere with others */
#define SHT_NO_DCT
#define ID_NME mem
#include "spat_to_SH.c"
#include "SH_to_spat.c"
#include "SHst_to_spat.c"
#include "spat_to_SHst.c"
/* nodct and gradients */
#define SHT_GRAD
#include "SHs_to_spat.c"
#include "SHt_to_spat.c"
#undef SHT_GRAD
/* nodct 3 components */
#define SHT_3COMP
#include "spat_to_SHqst.c"
#include "SHqst_to_spat.c"
#undef SHT_3COMP
#undef ID_NME

/* FUNCTION POINTER ARRAY */
void* GEN(sht_array, SUFFIX)[SHT_NALG][SHT_NTYP] = {
/* hyb */	{ GEN(SH_to_spat_hyb, SUFFIX), GEN(spat_to_SH_hyb, SUFFIX), GEN(SHsphtor_to_spat_hyb, SUFFIX), GEN(spat_to_SHsphtor_hyb, SUFFIX), 
#ifdef SHT_AXISYM
				GEN(SHsph_to_spat_hyb, SUFFIX), GEN(SHtor_to_spat_hyb, SUFFIX), GEN(SHqst_to_spat_hyb, SUFFIX), GEN(spat_to_SHqst_hyb, SUFFIX) },
#else
				GEN(SHsph_to_spat_hyb, SUFFIX), GEN(SHtor_to_spat_hyb, SUFFIX), NULL, GEN(spat_to_SHqst_hyb, SUFFIX) },
#endif
/* mem */	{ GEN(SH_to_spat_mem, SUFFIX), GEN(spat_to_SH_mem, SUFFIX), GEN(SHsphtor_to_spat_mem, SUFFIX), GEN(spat_to_SHsphtor_mem, SUFFIX), 
				GEN(SHsph_to_spat_mem, SUFFIX), GEN(SHtor_to_spat_mem, SUFFIX), GEN(SHqst_to_spat_mem, SUFFIX), GEN(spat_to_SHqst_mem, SUFFIX) },
#ifdef IVAR
/* s+v */	{ NULL, NULL, NULL, NULL, NULL, NULL, GEN(SHqst_to_spat_2, SUFFIX), GEN(spat_to_SHqst_2, SUFFIX) },
#else
/* s+v */	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
#endif
#ifdef SHT_VAR_LTR
/* fly1 */	{ NULL, NULL, GEN(SHsphtor_to_spat_fly1, SUFFIX), GEN(spat_to_SHsphtor_fly1, SUFFIX), 
				GEN(SHsph_to_spat_fly1, SUFFIX), GEN(SHtor_to_spat_fly1, SUFFIX), GEN(SHqst_to_spat_fly1, SUFFIX), GEN(spat_to_SHqst_fly1, SUFFIX) },
/* fly2 */	{ GEN(SH_to_spat_fly2, SUFFIX), GEN(spat_to_SH_fly2, SUFFIX), GEN(SHsphtor_to_spat_fly2, SUFFIX), GEN(spat_to_SHsphtor_fly2, SUFFIX), 
				GEN(SHsph_to_spat_fly2, SUFFIX), GEN(SHtor_to_spat_fly2, SUFFIX), GEN(SHqst_to_spat_fly2, SUFFIX), GEN(spat_to_SHqst_fly2, SUFFIX) },
/* fly3 */	{ GEN(SH_to_spat_fly3, SUFFIX), GEN(spat_to_SH_fly3, SUFFIX), GEN(SHsphtor_to_spat_fly3, SUFFIX), GEN(spat_to_SHsphtor_fly3, SUFFIX), 
				GEN(SHsph_to_spat_fly3, SUFFIX), GEN(SHtor_to_spat_fly3, SUFFIX), GEN(SHqst_to_spat_fly3, SUFFIX), GEN(spat_to_SHqst_fly3, SUFFIX) },
/* fly4 */	{ GEN(SH_to_spat_fly4, SUFFIX), GEN(spat_to_SH_fly4, SUFFIX), NULL, NULL, 
				GEN(SHsph_to_spat_fly4, SUFFIX), GEN(SHtor_to_spat_fly4, SUFFIX), NULL, NULL },
/* fly6 */	{ GEN(SH_to_spat_fly6, SUFFIX), GEN(spat_to_SH_fly6, SUFFIX), NULL, NULL, 
				NULL, NULL, NULL, NULL },
/* fly8 */	{ GEN(SH_to_spat_fly8, SUFFIX), GEN(spat_to_SH_fly8, SUFFIX), NULL, NULL, 
				NULL, NULL, NULL, NULL },
  #ifdef ADD_OPENMP
/* omp1 */	{ NULL, NULL, GEN(SHsphtor_to_spat_omp1, SUFFIX), GEN(spat_to_SHsphtor_omp1, SUFFIX), 
				GEN(SHsph_to_spat_omp1, SUFFIX), GEN(SHtor_to_spat_omp1, SUFFIX), GEN(SHqst_to_spat_omp1, SUFFIX), GEN(spat_to_SHqst_omp1, SUFFIX) },
/* omp2 */	{ GEN(SH_to_spat_omp2, SUFFIX), GEN(spat_to_SH_omp2, SUFFIX), GEN(SHsphtor_to_spat_omp2, SUFFIX), GEN(spat_to_SHsphtor_omp2, SUFFIX), 
				GEN(SHsph_to_spat_omp2, SUFFIX), GEN(SHtor_to_spat_omp2, SUFFIX), GEN(SHqst_to_spat_omp2, SUFFIX), GEN(spat_to_SHqst_omp2, SUFFIX) },
/* omp3 */	{ GEN(SH_to_spat_omp3, SUFFIX), GEN(spat_to_SH_omp3, SUFFIX), GEN(SHsphtor_to_spat_omp3, SUFFIX), GEN(spat_to_SHsphtor_omp3, SUFFIX), 
				GEN(SHsph_to_spat_omp3, SUFFIX), GEN(SHtor_to_spat_omp3, SUFFIX), GEN(SHqst_to_spat_omp3, SUFFIX), GEN(spat_to_SHqst_omp3, SUFFIX) },
/* omp4 */	{ GEN(SH_to_spat_omp4, SUFFIX), GEN(spat_to_SH_omp4, SUFFIX), NULL, NULL, 
				GEN(SHsph_to_spat_omp4, SUFFIX), GEN(SHtor_to_spat_omp4, SUFFIX), NULL, NULL },
/* omp6 */	{ GEN(SH_to_spat_omp6, SUFFIX), GEN(spat_to_SH_omp6, SUFFIX), NULL, NULL, 
				NULL, NULL, NULL, NULL },
/* omp8 */	{ GEN(SH_to_spat_omp8, SUFFIX), GEN(spat_to_SH_omp8, SUFFIX), NULL, NULL, 
				NULL, NULL, NULL, NULL }
  #else
/* omp1 */	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
/* omp2 */	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
/* omp3 */	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
/* omp4 */	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
/* omp6 */	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
/* omp8 */	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL }
  #endif
#endif
};


//@}

// Fortran 77 api
#if defined(SHT_F77_API) && defined(IVAR)

extern shtns_cfg sht_data;

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
	GEN(spat_to_SH,SUFFIX)(sht_data, Vr, Qlm SUPARGF2);
}

/// \ingroup fortapi
void GENF(sh_to_spat,SUFFIX)(complex double *Qlm, double *Vr SUPARGF) {
	GEN(SH_to_spat,SUFFIX)(sht_data, Qlm, Vr SUPARGF2);
}

/// \ingroup fortapi
void GENF(sphtor_to_spat,SUFFIX)(complex double *Slm, complex double *Tlm, double *Vt, double *Vp SUPARGF) {
	GEN(SHsphtor_to_spat,SUFFIX)(sht_data, Slm, Tlm, Vt, Vp SUPARGF2);
}

#ifndef SHT_AXISYM
/// \ingroup fortapi
void GENF(sph_to_spat,SUFFIX)(complex double *Slm, double *Vt, double *Vp SUPARGF) {
	GEN(SHsph_to_spat,SUFFIX)(sht_data, Slm, Vt, Vp SUPARGF2);
}

/// \ingroup fortapi
void GENF(tor_to_spat,SUFFIX)(complex double *Tlm, double *Vt, double *Vp SUPARGF) {
	GEN(SHtor_to_spat,SUFFIX)(sht_data, Tlm, Vt, Vp SUPARGF2);
}
#else
/// \ingroup fortapi
void GENF(sph_to_spat,SUFFIX)(complex double *Slm, double *Vt SUPARGF) {
	GEN(SHsph_to_spat,SUFFIX)(sht_data, Slm, Vt SUPARGF2);
}

/// \ingroup fortapi
void GENF(tor_to_spat,SUFFIX)(complex double *Tlm, double *Vp SUPARGF) {
	GEN(SHtor_to_spat,SUFFIX)(sht_data, Tlm, Vp SUPARGF2);
}
#endif

/// \ingroup fortapi
void GENF(qst_to_spat,SUFFIX)(complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp SUPARGF)
{
	GEN(SHqst_to_spat,SUFFIX)(sht_data, Qlm, Slm, Tlm, Vr, Vt, Vp SUPARGF2);
}

/// \ingroup fortapi
void GENF(spat_to_sphtor,SUFFIX)(double *Vt, double *Vp, complex double *Slm, complex double *Tlm SUPARGF) {
	GEN(spat_to_SHsphtor,SUFFIX)(sht_data, Vt, Vp, Slm, Tlm SUPARGF2);
}

/// \ingroup fortapi
void GENF(spat_to_qst,SUFFIX)(double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm SUPARGF)
{
	GEN(spat_to_SHqst,SUFFIX)(sht_data, Vr, Vt, Vp, Qlm, Slm, Tlm SUPARGF2);
}

//@}

#endif

#undef GEN
#undef GEN3
#undef GENF
#undef _GENF
#undef SUFFIX
#undef SUPARG
#undef SUPARG2
#undef SUPARGF
#undef SUPARGF2
#undef SHT_AXISYM
#undef SHT_VAR_LTR
#undef IVAR
#undef ADD_OPENMP
