/********************************************************************
 * SHTns : Spherical Harmonic Transform for numerical simulations.  *
 *    written by Nathanael Schaeffer / LGIT,CNRS                    *
 ********************************************************************/

// global variables definitions
#include "sht_private.h"

// axisymmetric transform
#define SHT_AXISYM
// truncation at ltr and MMAX
#define LTR ltr
#define MTR MMAX
#define SHT_VAR_LTR

#define SUFFIX _m0l
#define SUPARG , int ltr
#define SUPARG2 , ltr
#define SUPARGF , int *ltr
#define SUPARGF2 , *ltr

/** \addtogroup shtm0ltr Axisymmetric SHT functions with variable harmonic degree (l) truncation.
 * 
 * these work for any MMAX, and will only transform m=0, to/from arrays with NLAT contiguous theta-values.
 * 
 * All functions in this group have a last argument ltr that is the maximum degree of spherical harmonic that is taken into account.
 * 
 * For synthesis, coefficients wiht l>ltr are ignored, but for analysis coefficient with l>ltr are set to zero.
*/
//@{
#include "SHT/sht_generic.c"

