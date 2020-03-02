/*
 * Copyright (c) 2010-2020 Centre National de la Recherche Scientifique.
 * written by Nathanael Schaeffer (CNRS, ISTerre, Grenoble, France).
 * 
 * nathanael.schaeffer@univ-grenoble-alpes.fr
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

#include "sht_private.h"

#define MTR MMAX
#define SHT_VAR_LTR

#define GEN(name,sfx) GLUE2(name,sfx)
#define GEN3(name,nw,sfx) GLUE3(name,nw,sfx)

// genaral case, hi lmax
#undef SUFFIX
#define SUFFIX _l
#define HI_LLIM

  #ifdef _GCC_VEC_
	#define NWAY 1
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
  #endif
	#define NWAY 2
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 6
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 8
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY

#define SHT_GRAD
	#define NWAY 2
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
	#ifdef _GCC_VEC_
	#define NWAY 1
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
	#endif
	#define NWAY 4
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
#undef SHT_GRAD

#define SHT_3COMP
  #ifdef _GCC_VEC_
	#define NWAY 1
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
  #endif
	#define NWAY 2
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 6
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 8
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
#undef SHT_3COMP

// genaral case, low lmax
#undef SUFFIX
#define SUFFIX _l
#undef HI_LLIM

  #ifdef _GCC_VEC_
	#define NWAY 1
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
  #endif
	#define NWAY 2
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 6
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 8
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY

#define SHT_GRAD
	#define NWAY 2
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
	#ifdef _GCC_VEC_
	#define NWAY 1
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
	#endif
	#define NWAY 4
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
#undef SHT_GRAD

#define SHT_3COMP
  #ifdef _GCC_VEC_
	#define NWAY 1
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
  #endif
	#define NWAY 2
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 6
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 8
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
#undef SHT_3COMP


// axisymmetric
#define SHT_AXISYM
#undef SUFFIX
#define SUFFIX _m0l

  #ifdef _GCC_VEC_
	#define NWAY 1
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
  #endif
	#define NWAY 2
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/SH_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 6
	#include "SHT/SH_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 8
	#include "SHT/SH_to_spat_kernel.c"
	#undef NWAY

#define SHT_GRAD
	#define NWAY 2
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
	#ifdef _GCC_VEC_
	#define NWAY 1
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
	#endif
	#define NWAY 4
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
#undef SHT_GRAD

#define SHT_3COMP
  #ifdef _GCC_VEC_
	#define NWAY 1
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
  #endif
	#define NWAY 2
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
#undef SHT_3COMP
