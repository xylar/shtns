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
  #define GEN(name,sfx) GLUE2(name,sfx)
  #define GEN3(name,nw,sfx) GLUE3(name,nw,sfx)
#else
  #define GEN(name,sfx) name
  #define GEN3(name,nw,sfx) GLUE2(name,nw)
#endif
#ifndef SUPARG
  #define SUPARG
  #define PF2 pf2
  #define PF3 pf3
  #define PF4 pf4
  #define PF6 pf6
#else
  #define PF2 pf2l
  #define PF3 pf3l
  #define PF4 pf4l
  #define PF6 pf6l
#endif
#ifndef SUPARG2
  #define SUPARG2
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

#define SHT_NO_DCT
#define ID_NME nodct
#include "spat_to_SH.c"
#include "SH_to_spat.c"
#undef ID_NME
#undef SHT_NO_DCT

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

/// Backward \b Scalar Spherical Harmonic Transform (synthesis).
void GEN(SH_to_spat,SUFFIX)(shtns_cfg shtns, complex double *Qlm, double *Vr SUPARG)
{
	((PF2) shtns->fptr[IVAR][SHT_TYP_SSY])(shtns, Qlm, Vr SUPARG2);
	return;
}

void GEN(spat_to_SH,SUFFIX)(shtns_cfg shtns, double *Vr, complex double *Qlm SUPARG)
{
	((PF2) shtns->fptr[IVAR][SHT_TYP_SAN])(shtns, Vr, Qlm SUPARG2);
	return;
}

/// Backward \b Vector Spherical Harmonic Transform (synthesis).
void GEN(SHsphtor_to_spat,SUFFIX)(shtns_cfg shtns, complex double *Slm, complex double *Tlm, double *Vt, double *Vp SUPARG)
{
	((PF4) shtns->fptr[IVAR][SHT_TYP_VSY])(shtns, Slm, Tlm, Vt, Vp SUPARG2);
	return;
}

/// \b Vector Spherical Harmonics Transform (analysis) : convert a spatial vector field (theta,phi components) to its spheroidal/toroidal spherical harmonic representation.
void GEN(spat_to_SHsphtor,SUFFIX)(shtns_cfg shtns, double *Vt, double *Vp, complex double *Slm, complex double *Tlm SUPARG)
{
	((PF4) shtns->fptr[IVAR][SHT_TYP_VAN])(shtns, Vt, Vp, Slm, Tlm SUPARG2);
	return;
}

//@}


/* GRADIENTS */
#define SHT_GRAD

#define ID_NME hyb
#include "SHs_to_spat.c"
#include "SHt_to_spat.c"
#undef ID_NME

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

#ifndef SHT_AXISYM
/// Spheroidal only synthesis.
void GEN(SHsph_to_spat,SUFFIX)(shtns_cfg shtns, complex double *Slm, double *Vt, double *Vp SUPARG)
{
	((PF3) shtns->fptr[IVAR][SHT_TYP_GSP])(shtns, Slm, Vt, Vp SUPARG2);
	return;
}

/// Toroidal only synthesis.
void GEN(SHtor_to_spat,SUFFIX)(shtns_cfg shtns, complex double *Tlm, double *Vt, double *Vp SUPARG)
{
	((PF3) shtns->fptr[IVAR][SHT_TYP_GTO])(shtns, Tlm, Vt, Vp SUPARG2);
	return;
}
#else
/// Spheroidal m=0 only synthesis (results in theta component only).
void GEN(SHsph_to_spat,SUFFIX)(shtns_cfg shtns, complex double *Slm, double *Vt SUPARG)
{
	((PF2) shtns->fptr[IVAR][SHT_TYP_GSP])(shtns, Slm, Vt SUPARG2);
	return;
}

/// Toroidal m=0 only synthesis (results in phi component only).
void GEN(SHtor_to_spat,SUFFIX)(shtns_cfg shtns, complex double *Tlm, double *Vp SUPARG)
{
	((PF2) shtns->fptr[IVAR][SHT_TYP_GTO])(shtns, Tlm, Vp SUPARG2);
	return;
}
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
#include "SHqst_to_spat.c"
#undef ID_NME

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

// combining vector and scalar.
void GEN(SHqst_to_spat_2,SUFFIX)(shtns_cfg shtns, complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp SUPARG)
{
	GEN(SH_to_spat,SUFFIX)(shtns, Qlm, Vr SUPARG2);
	GEN(SHsphtor_to_spat,SUFFIX)(shtns, Slm, Tlm, Vt, Vp SUPARG2);
}
void GEN(spat_to_SHqst_2,SUFFIX)(shtns_cfg shtns, double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm SUPARG)
{
	GEN(spat_to_SH,SUFFIX)(shtns, Vr, Qlm SUPARG2);
	GEN(spat_to_SHsphtor,SUFFIX)(shtns, Vt, Vp, Slm, Tlm SUPARG2);
}

void GEN(spat_to_SHqst,SUFFIX)(shtns_cfg shtns, double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm SUPARG)
{
	((PF6) shtns->fptr[IVAR][SHT_TYP_3AN])(shtns, Vr, Vt, Vp, Qlm, Slm, Tlm SUPARG2);
	return;
}

void GEN(SHqst_to_spat,SUFFIX)(shtns_cfg shtns, complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp SUPARG)
{
	((PF6) shtns->fptr[IVAR][SHT_TYP_3SY])(shtns, Qlm, Slm, Tlm, Vr, Vt, Vp SUPARG2);
	return;
}

/* FUNCTION POINTER ARRAY */
void* GEN(sht_array, SUFFIX)[SHT_NTYP][SHT_NALG] = {
		{ GEN(SH_to_spat_hyb, SUFFIX), NULL, NULL, GEN(SH_to_spat_fly2, SUFFIX),
		  GEN(SH_to_spat_fly3, SUFFIX), GEN(SH_to_spat_fly4, SUFFIX), GEN(SH_to_spat_fly6, SUFFIX), GEN(SH_to_spat_fly8, SUFFIX) },
		{ GEN(spat_to_SH_hyb, SUFFIX), NULL, NULL, GEN(spat_to_SH_fly2, SUFFIX),
		  GEN(spat_to_SH_fly3, SUFFIX), GEN(spat_to_SH_fly4, SUFFIX), GEN(spat_to_SH_fly6, SUFFIX), GEN(spat_to_SH_fly8, SUFFIX) },
		{ GEN(SHsphtor_to_spat_hyb, SUFFIX), NULL, GEN(SHsphtor_to_spat_fly1, SUFFIX), GEN(SHsphtor_to_spat_fly2, SUFFIX),
		  GEN(SHsphtor_to_spat_fly3, SUFFIX), NULL, NULL, NULL },
		{ GEN(spat_to_SHsphtor_hyb, SUFFIX), NULL, GEN(spat_to_SHsphtor_fly1, SUFFIX), GEN(spat_to_SHsphtor_fly2, SUFFIX),
		  GEN(spat_to_SHsphtor_fly3, SUFFIX), NULL, NULL, NULL },
		{ GEN(SHsph_to_spat_hyb, SUFFIX), NULL, GEN(SHsph_to_spat_fly1, SUFFIX), GEN(SHsph_to_spat_fly2, SUFFIX),
		  GEN(SHsph_to_spat_fly3, SUFFIX), GEN(SHsph_to_spat_fly4, SUFFIX), NULL, NULL },
		{ GEN(SHtor_to_spat_hyb, SUFFIX), NULL, GEN(SHtor_to_spat_fly1, SUFFIX), GEN(SHtor_to_spat_fly2, SUFFIX),
		  GEN(SHtor_to_spat_fly3, SUFFIX), GEN(SHtor_to_spat_fly4, SUFFIX), NULL, NULL },
		{ GEN(SHqst_to_spat_hyb, SUFFIX), GEN(SHqst_to_spat_2, SUFFIX), GEN(SHqst_to_spat_fly1, SUFFIX),
		  GEN(SHqst_to_spat_fly2, SUFFIX), GEN(SHqst_to_spat_fly3, SUFFIX), NULL, NULL, NULL },
		{ GEN(spat_to_SHqst_hyb, SUFFIX), GEN(spat_to_SHqst_2, SUFFIX), GEN(spat_to_SHqst_fly1, SUFFIX),
		  GEN(spat_to_SHqst_fly2, SUFFIX), GEN(spat_to_SHqst_fly3, SUFFIX), NULL, NULL, NULL }
};

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
#undef PF2
#undef PF3
#undef PF4
#undef PF6
#undef IVAR
