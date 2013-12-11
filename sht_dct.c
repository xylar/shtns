#include "sht_private.h"

#define MTR MMAX

/* Hybrid (DCT enabled) functions */
#define ID_NME hyb

// standard
#undef SHT_AXISYM
#undef SHT_VAR_LTR
#define SUFFIX
  #define GEN(name,sfx) name
  #define GEN3(name,nw,sfx) GLUE2(name,nw)

#include "SHT/spat_to_SH.c"
#include "SHT/SH_to_spat.c"
#include "SHT/SHst_to_spat.c"
#include "SHT/spat_to_SHst.c"
#define SHT_GRAD
#include "SHT/SHs_to_spat.c"
#include "SHT/SHt_to_spat.c"
#undef SHT_GRAD
#define SHT_3COMP
#include "SHT/spat_to_SHqst.c"
#undef SHT_3COMP

// axisymmetric
#define SHT_AXISYM
#undef SHT_VAR_LTR
#undef SUFFIX
#define SUFFIX _m0
  #undef GEN
  #undef GEN3
  #define GEN(name,sfx) GLUE2(name,sfx)
  #define GEN3(name,nw,sfx) GLUE3(name,nw,sfx)

#include "SHT/spat_to_SH.c"
#include "SHT/SH_to_spat.c"
#include "SHT/SHst_to_spat.c"
#include "SHT/spat_to_SHst.c"
#define SHT_GRAD
#include "SHT/SHs_to_spat.c"
#include "SHT/SHt_to_spat.c"
#undef SHT_GRAD
#define SHT_3COMP
#include "SHT/spat_to_SHqst.c"
#include "SHT/SHqst_to_spat.c"
#undef SHT_3COMP


// standard, l-truncation
#undef SHT_AXISYM
#define SHT_VAR_LTR
#undef SUFFIX
#define SUFFIX _l

// scalar
#include "SHT/spat_to_SH.c"
#include "SHT/SH_to_spat.c"
// vector
#include "SHT/SHst_to_spat.c"
#include "SHT/spat_to_SHst.c"
// gradients
#define SHT_GRAD
#include "SHT/SHs_to_spat.c"
#include "SHT/SHt_to_spat.c"
#undef SHT_GRAD
// 3 components
#define SHT_3COMP
#include "SHT/spat_to_SHqst.c"
#undef SHT_3COMP


// axisymm, l-truncation
#define SHT_AXISYM
#define SHT_VAR_LTR
#undef SUFFIX
#define SUFFIX _m0l

// scalar
#include "SHT/spat_to_SH.c"
#include "SHT/SH_to_spat.c"
// vector
#include "SHT/SHst_to_spat.c"
#include "SHT/spat_to_SHst.c"
// gradients
#define SHT_GRAD
#include "SHT/SHs_to_spat.c"
#include "SHT/SHt_to_spat.c"
#undef SHT_GRAD
// 3 components
#define SHT_3COMP
#include "SHT/spat_to_SHqst.c"
#include "SHT/SHqst_to_spat.c"
#undef SHT_3COMP

#undef ID_NME

void* fdct[SHT_NTYP] = { SH_to_spat_hyb, spat_to_SH_hyb, SHsphtor_to_spat_hyb, spat_to_SHsphtor_hyb,
		SHsph_to_spat_hyb, SHtor_to_spat_hyb, NULL, spat_to_SHqst_hyb };
void* fdct_l[SHT_NTYP] = { SH_to_spat_hyb_l, spat_to_SH_hyb_l, SHsphtor_to_spat_hyb_l, spat_to_SHsphtor_hyb_l,
		SHsph_to_spat_hyb_l, SHtor_to_spat_hyb_l, NULL, spat_to_SHqst_hyb_l };
void* fdct_m0[SHT_NTYP] = { SH_to_spat_hyb_m0, spat_to_SH_hyb_m0, SHsphtor_to_spat_hyb_m0, spat_to_SHsphtor_hyb_m0,
		SHsph_to_spat_hyb_m0, SHtor_to_spat_hyb_m0, SHqst_to_spat_hyb_m0, spat_to_SHqst_hyb_m0 };
void* fdct_m0l[SHT_NTYP] = { SH_to_spat_hyb_m0l, spat_to_SH_hyb_m0l, SHsphtor_to_spat_hyb_m0l, spat_to_SHsphtor_hyb_m0l,
		SHsph_to_spat_hyb_m0l, SHtor_to_spat_hyb_m0l, SHqst_to_spat_hyb_m0l, spat_to_SHqst_hyb_m0l };

