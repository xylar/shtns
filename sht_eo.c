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
#ifndef SHT_SCALAR_ONLY
	#include "SHT/SHeost_to_spat.c"
#endif
}

void spat_to_SHeo_sphtor(double *Vt, double *Vp, complex double *Slm, complex double *Tlm, int parity)
{
#ifndef SHT_SCALAR_ONLY
	#include "SHT/spat_to_SHeost.c"
#endif
}

//@}

#ifdef SHT_F77_API

// TODO

#endif
