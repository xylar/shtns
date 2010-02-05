/********************************************************************
 * SHTns : Spherical Harmonic Transform for numerical simulations.  *
 *    written by Nathanael Schaeffer / LGIT,CNRS                    *
 ********************************************************************/

// global variables definitions
#include "sht_private.h"

#define SHT_NO_DCT
#define LTR LMAX
#define MTR MMAX
#undef SHT_VAR_LTR

/// \addtogroup shteo SHT functions with assumed equatorial symmetry (processes only one hemisphere).
//@{
void SHeo_to_spat(complex double *Qlm, double *Vr, int parity)
{
	#include "SHT/SHeo_to_spat.c"
}

void spat_to_SHeo(double *Vr, complex double *Qlm, int parity)
{
	#include "SHT/spat_to_SHeo.c"
}

void SHeo_sphtor_to_spat(complex double *Slm, complex double *Tlm, double *Vt, double *Vp, int parity)
{
	#include "SHT/SHeost_to_spat.c"
}

void spat_to_SHeo_sphtor(double *Vt, double *Vp, complex double *Slm, complex double *Tlm, int parity)
{
	#include "SHT/spat_to_SHeost.c"
}

//@}

#ifdef SHT_F77_API

// TODO

#endif
