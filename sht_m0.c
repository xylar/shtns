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

// \file sht_m0.c SHT functions restricted to axisymmetric data (m=0)

// global variables definitions
#include "sht_private.h"

// axisymmetric transform
#define SHT_AXISYM
// truncation at LMAX and MMAX
#define MTR MMAX
#undef SHT_VAR_LTR

#define SUFFIX _m0

/** \addtogroup shtm0 SHT functions restricted to axisymmetric data (m=0)
 * these work for any MMAX, and will only transform m=0, to/from arrays with NLAT contiguous theta-values.
 */
//@{
#include "SHT/sht_generic.c"
