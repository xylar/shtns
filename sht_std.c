/********************************************************************
 * SHTns : Spherical Harmonic Transform for numerical simulations.  *
 *    written by Nathanael Schaeffer / LGIT,CNRS                    *
 ********************************************************************/

// global variables definitions
#include "sht_private.h"

// truncation at LMAX and MMAX
#define LTR LMAX
#define MTR MMAX
#undef SHT_VAR_LTR

/** \addtogroup shtbase Basic Spherical Harmonic Transform functions.
 * Spherical Harmonic transforms for sizes set by \ref shtns_init, \ref shtns_set_size and \ref shtns_precompute
*/
//@{
#include "SHT/sht_generic.c"
