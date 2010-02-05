/********************************************************************
 * SHTns : Spherical Harmonic Transform for numerical simulations.  *
 *    written by Nathanael Schaeffer / LGIT,CNRS                    *
 ********************************************************************/

// \file sht_ltr.c SHT functions with variable harmonic degree truncation.

// global variables definitions
#include "sht_private.h"

// truncation at LMAX and MMAX
#define LTR ltr
#define MTR MMAX
#define SHT_VAR_LTR

#define SUFFIX _l
#define SUPARG , int ltr
#define SUPARG2 , ltr
#define SUPARGF , int *ltr
#define SUPARGF2 , *ltr

/** \addtogroup shtltr SHT functions with variable harmonic degree (l) truncation.
 * All functions in this group have a last argument ltr that is the maximum degree of spherical harmonic that is taken into account.
 * 
 * For synthesis, coefficients wiht l>ltr are ignored, but for analysis coefficient with l>ltr are set to zero.
*/
//@{
#include "SHT/sht_generic.c"

