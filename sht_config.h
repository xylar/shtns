/*
 * Copyright (c) 2010-2011 Centre National de la Recherche Scientifique.
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

/// \file sht_config.h compile-time configuration.

/// 0:no output, 1:output info to stdout, 2:more output (debug info), 3:also print fftw plans.
#define SHT_VERBOSE 2

/// defines the maximum amount of memory in megabytes that SHTns should use.
#define SHTNS_MAX_MEMORY 2048

/// Compile the \ref fortapi
//#define SHT_F77_API


/// Minimum performance improve for DCT in \ref sht_auto mode. If not atained, we switch back to gauss.
#define MIN_PERF_IMPROVE_DCT 0.95
/// Try to enforce at least this accuracy for DCT in sht_auto mode.
#define MIN_ACCURACY_DCT 1.e-8

/// The default \ref opt_polar threshold (0 disabled, 1.e-6 is aggressive, 1.e-10 is safe, 1.e-14 is VERY safe)
#define SHT_DEFAULT_POLAR_OPT 1.e-10

/// The default \ref norm used by shtns_init
#define SHT_DEFAULT_NORM ( sht_orthonormal )
//#define SHT_DEFAULT_NORM ( sht_schmidt | SHT_NO_CS_PHASE )

/// The maximum order of non-linear terms to be resolved by SH transform by default.
/// 1 : no non-linear terms. 2 : quadratic non-linear terms (default), 3 : triadic, ...
/// must be larger or equal to 1.
#define SHT_DEFAULT_NL_ORDER 1

/// I compile with GCC 4 or later, and I would like fast vectorized code (if SSE2 is supported) !
/// (set to zero to disable, may be useful for calling from Fortran)
#define _GCC_VEC_ 1

// if SHT_SCALAR_ONLY is defined, it will disable the vector transform (which saves some memory)
//#define SHT_SCALAR_ONLY

// if SHT_AXISYM is defined, an axisymmetric-only transform will be compiled.
//#define SHT_AXISYM

// if SHT_NO_DCT is defined, no DCT support will be compiled. (Gauss-legendre only).
//#define SHT_NO_DCT
