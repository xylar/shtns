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
 * This is the definition
*/
//@{
#include "SHT/sht_generic.c"
