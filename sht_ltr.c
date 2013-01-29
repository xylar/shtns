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
 *    written by Nathanael Schaeffer / CNRS                         *
 ********************************************************************/

// global variables definitions
#include "sht_private.h"

// truncation at MMAX
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
#define IVAR SHT_LTR
#include "SHT/sht_generic.c"

