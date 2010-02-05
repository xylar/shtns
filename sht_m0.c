/********************************************************************
 * SHTns : Spherical Harmonic Transform for numerical simulations.  *
 *    written by Nathanael Schaeffer / LGIT,CNRS                    *
 ********************************************************************/

// \file sht_m0.c SHT functions restricted to axisymmetric data (m=0)

// global variables definitions
#include "sht_private.h"

// axisymmetric transform
#define SHT_AXISYM
// truncation at LMAX and MMAX
#define LTR LMAX
#define MTR MMAX
#undef SHT_VAR_LTR

#define SUFFIX _m0

/** \addtogroup shtm0 SHT functions restricted to axisymmetric data (m=0)
 * these work for any MMAX, and will only transform m=0, to/from arrays with NLAT contiguous theta-values.
 */
//@{
#include "SHT/sht_generic.c"
