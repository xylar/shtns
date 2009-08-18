/********************************************************************
 * SHTns : Spherical Harmonic Transform for numerical simulations.  *
 *    written by Nathanael Schaeffer / LGIT,CNRS                    *
 ********************************************************************/

/// \file SHTf77.c interface to Fortran language (compatible with gfortran)
/// \example SHT_example.f
/// \brief A Fortran example program that performs backward and forward Spherical Harmonic Transforms using SHTns
/// compile using : make SHT_fort_ex


#include <complex.h>
#include <math.h>
#include "SHT.h"

/** \name Initialization
Call from fortran using : \code
call shtns_init_sh_*( lmax, mmax, mres, nlat, nphi) \endcode
\see init_SH for argument description
*/
//@{
/// Initializes spherical harmonic transforms of given size using Gauss algorithm and no approximation
void shtns_init_sh_gauss_(int *lmax, int *mmax, int *mres, int *nlat, int *nphi)
{
    init_SH(sht_gauss, 0, *lmax, *mmax, *mres, *nlat, *nphi);
}

/// Initializes spherical harmonic transforms of given size using Fastest available algorithm and polar optimization.
void shtns_init_sh_auto_(int *lmax, int *mmax, int *mres, int *nlat, int *nphi)
{
    init_SH(sht_auto, 1.e-6, *lmax, *mmax, *mres, *nlat, *nphi);
}

/// Initializes spherical harmonic transforms of given size using a regular grid and agressive optimizations.
void shtns_init_sh_reg_fast_(int *lmax, int *mmax, int *mres, int *nlat, int *nphi)
{
    init_SH(sht_reg_fast, 1.e-6, *lmax, *mmax, *mres, *nlat, *nphi);
}

/// Initializes spherical harmonic transform SYNTHESIS ONLY of given size using a regular grid including poles.
void shtns_init_sh_poles_(int *lmax, int *mmax, int *mres, int *nlat, int *nphi)
{
    init_SH(sht_reg_poles, 0, *lmax, *mmax, *mres, *nlat, *nphi);
}

/// returns the number of complex*16 elements in an SH array
void shtns_get_nlm_(int *nlm, int *lmax, int *mmax, int *mres)
{
    *nlm = nlm_calc(*lmax, *mmax, *mres);
}

/// returns the number of real*8 that must be allocated for spatial data
void shtns_get_nspat_alloc_(int *nspat, int* nlat, int* nphi)
{
    *nspat = (*nlat)*(*nphi/2+1)*2;
}
//@}

/** \name Scalar transforms
Call from fortran using : \code
call shtns_spat_to_sh( Brf, Qlm )
call shtns_sh_to_spat( Qlm, Brf ) \endcode
\see spat_to_SH and for argument description
*/
//@{
void shtns_spat_to_sh_(complex double *BrF, complex double *Qlm)
{
    spat_to_SH(BrF, Qlm);
}

void shtns_sh_to_spat_(complex double *Qlm, complex double *BrF)
{
    SH_to_spat(Qlm, BrF);
}
//@}
