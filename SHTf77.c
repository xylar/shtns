/********************************************************************
 * SHTns : Spherical Harmonic Transform for numerical simulations.  *
 *    written by Nathanael Schaeffer / LGIT,CNRS                    *
 ********************************************************************/

/// \file SHTf77.c interface to Fortran language (compatible with gfortran).
/// see the \link SHT_example.f Fortran example \endlink for a simple usage of SHTns from Fortran language.

/** \example SHT_example.f
  \brief A Fortran example program that performs backward and forward Spherical Harmonic Transforms using SHTns.
  Compile using : \code make SHT_fort_ex \endcode
*/


#include <complex.h>
#include <math.h>
#include <stdio.h>
#include "SHT.h"

/** \name Initialization
Call from fortran using : \code
call shtns_init_sh_*( lmax, mmax, mres, nlat, nphi) \endcode
(without the trailing '_')
\see init_SH for argument description
*/
//@{
/// Initializes spherical harmonic transforms of given size using Gauss algorithm and no approximation
void shtns_init_sh_gauss_(int *layout, int *lmax, int *mmax, int *mres, int *nlat, int *nphi)
{
    shtns_init(sht_gauss | *layout, 0, *lmax, *mmax, *mres, *nlat, *nphi);
}

/// Initializes spherical harmonic transforms of given size using Fastest available algorithm and polar optimization.
void shtns_init_sh_auto_(int *layout, int *lmax, int *mmax, int *mres, int *nlat, int *nphi)
{
    shtns_init(sht_auto | *layout, 1.e-10, *lmax, *mmax, *mres, *nlat, *nphi);
}

/// Initializes spherical harmonic transforms of given size using a regular grid and agressive optimizations.
void shtns_init_sh_reg_fast_(int *layout, int *lmax, int *mmax, int *mres, int *nlat, int *nphi)
{
    shtns_init(sht_reg_fast | *layout, 1.e-6, *lmax, *mmax, *mres, *nlat, *nphi);
}

/// Initializes spherical harmonic transform SYNTHESIS ONLY of given size using a regular grid including poles.
void shtns_init_sh_poles_(int *layout, int *lmax, int *mmax, int *mres, int *nlat, int *nphi)
{
    shtns_init(sht_reg_poles | *layout, 0, *lmax, *mmax, *mres, *nlat, *nphi);
}

/// Defines the size and convention of the transform.
/// Allow to choose the normalization and whether or not to include the Condon-Shortley phase.
void shtns_set_size_(int *lmax, int *mmax, int *mres, int *norm, int *cs_phase)
{
	if (*cs_phase)
		shtns_set_size(*lmax, *mmax, *mres, *norm);
	else
		shtns_set_size(*lmax, *mmax, *mres, *norm | SHT_NO_CS_PHASE);
}

/// Precompute matrices for synthesis and analysis.
/// Allow to choose polar optimization threshold and algorithm type.
void shtns_precompute_(int *flags, double *eps, int *nlat, int *nphi)
{
	shtns_precompute(*flags, *eps, *nlat, *nphi);
}
//@}

/// returns nlm, the number of complex*16 elements in an SH array.
/// call from fortran using \code call shtns_get_nlm(nlm, lmax, mmax, mres) \endcode
void shtns_calc_nlm_(int *nlm, int *lmax, int *mmax, int *mres)
{
    *nlm = nlm_calc(*lmax, *mmax, *mres);
}

/// returns lm, the index in an SH array of mode (l,m).
/// call from fortran using \code call shtns_get_lmidx(lm, l, m) \endcode
void shtns_lmidx_(int *lm, int *l, int *m)
{
    *lm = LM(*l, *m) + 1;
}

/// fills the given array with the cosine of the co-latitude angle (NLAT real*8)
void shtns_cos_array_(double *costh)
{
	int i;	
	for (i=0; i<shtns.nlat; i++)
		costh[i] = ct[i];
}


/** \name Scalar transforms
Call from fortran using : \code
call shtns_spat_to_sh( Brf, Qlm )
call shtns_sh_to_spat( Qlm, Brf ) \endcode
\see spat_to_SH for argument description
*/
//@{
void shtns_spat_to_sh_(double *Vr, complex double *Qlm)
{
    spat_to_SH(Vr, Qlm);
}

void shtns_sh_to_spat_(complex double *Qlm, double *Vr)
{
    SH_to_spat(Qlm, Vr);
}
//@}


/** \name Surface vector transforms
\see spat_to_SHsphtor for argument description
*/
//@{
void shtns_spat_to_sphtor_(double *Vt, double *Vp, complex double *Slm, complex double *Tlm)
{
    spat_to_SHsphtor(Vt, Vp, Slm, Tlm);
}

void shtns_sphtor_to_spat_(complex double *Slm, complex double *Tlm, double *Vt, double *Vp)
{
    SHsphtor_to_spat(Slm, Tlm, Vt, Vp);
}

void shtns_sph_to_spat_(complex double *Slm, double *Vt, double *Vp)
{
    SHsph_to_spat(Slm, Vt, Vp);
}

void shtns_tor_to_spat_(complex double *Tlm, double *Vt, double *Vp)
{
    SHtor_to_spat(Tlm, Vt, Vp);
}
//@}

/** \name 3D vector transforms
\see spat_to_SHsphtor for argument description
*/
//@{
void shtns_spat_to_qst_(double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm)
{
    spat_to_SHsphtor(Vt, Vp, Slm, Tlm);
    spat_to_SH(Vr, Qlm);
}

void shtns_qst_to_spat_(complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp)
{
    SH_to_spat(Qlm, Vr);
    SHsphtor_to_spat(Slm, Tlm, Vt, Vp);
}
//@}

/** \name Point evaluation of Spherical Harmonics
Evaluate at a given point (\f$cos(\theta)\f$ and \f$\phi\f$) a spherical harmonic representation.
\see spat_to_SH and for argument description
*/
//@{
void shtns_sh_to_point_(double *spat, complex double *Qlm, double *cost, double *phi)
{
	*spat = SH_to_point(Qlm, *cost, *phi);
}

void shtns_qst_to_point_(double *vr, double *vt, double *vp,
		complex double *Qlm, complex double *Slm, complex double *Tlm, double *cost, double *phi)
{
	SHqst_to_point(Qlm, Slm, Tlm, *cost, *phi, vr, vt, vp);
}
//@}

/** \name Axisymmetric Spherical Harmonic transform
*/
//@{
void shtns_spat_to_sh_m0_(double *Vr, complex double *Qlm) {
	spat_to_SH_m0(Vr, Qlm);
}
void shtns_sh_to_spat_m0_(complex double *Qlm, double *Vr) {
    SH_to_spat_m0(Qlm, Vr);
}

void shtns_spat_to_sphtor_m0_(double *Vt, double *Vp, complex double *Slm, complex double *Tlm) {
    spat_to_SHsphtor_m0(Vt, Vp, Slm, Tlm);
}
void shtns_sphtor_to_spat_m0_(complex double *Slm, complex double *Tlm, double *Vt, double *Vp) {
    SHsphtor_to_spat_m0(Slm, Tlm, Vt, Vp);
}

void shtns_sph_to_spat_m0_(complex double *Slm, double *Vt) {
    SHsph_to_spat_m0(Slm, Vt);
}
void shtns_tor_to_spat_m0_(complex double *Tlm, double *Vp) {
    SHtor_to_spat_m0(Tlm, Vp);
}
//@}
