/**
Definition file for SHTns.
**/

/// GLOBAL VARIABLES : these are set by the call to init_SH. ///

extern long int LMAX;	// maximum degree (LMAX) of spherical harmonics.
extern long int NLAT;	// number of spatial points in Theta direction (latitude).
#ifndef SHT_AXISYM
  extern long int MMAX,MRES;	// maximum order (MMAX*MRES) of spherical harmonics. MRES is the periodicity along the phi coord.
  extern long int NPHI;	// number of spatial points in Phi direction (longitude)
#else
  #define MMAX 0
  #define NPHI 1
  #define MRES 1
#endif
extern long int NLM;	// total number of (l,m) spherical harmonics components.

extern double *ct, *st, *st_1;	// cos(theta), sin(theta), 1/sin(theta);
extern double *el, *l2, *l_2;	// l, l(l+1) and 1/(l(l+1))
extern int *li;

extern long int *lmidx;		// (virtual) index in SH array of given im.

/// MACROS ///

/** ACCESS TO SPHERICAL HARMONICS COMPONENTS **
  LM(l,m) : macro returning array index for given l and m.
  LiM(l,im) : macro returning array index for given l and im.
  LM_LOOP( action ) : macro that performs "action" for every (l,m), with l and lm set, but neither m nor im.
  el[NLM], l2[NLM], l_2[NLM] : floating point arrays containing l, l*(l+1) and 1/(l*(l+1))
**/
// LM(l,m) : index in the Spherical Harmonic coefficient array [ (l,m) space ]
#define LiM(l,im) ( lmidx[im] + l )
#define LM(l,m) ( lmidx[(m)/MRES] + l )

//LM_LOOP : loop over all lm's and perform "action". only lm is defined, neither l nor m.
#define LM_LOOP( action ) for (lm=0; lm<NLM; lm++) { action }

// LM_L_LOOP : loop over all (l,im) and perform "action"  : l and lm are defined. (but NOT m and im)
#define LM_L_LOOP( action ) for (lm=0; lm<NLM; lm++) { l=li[lm]; { action } }

#ifndef M_PI
  #define M_PI 3.1415926535897932384626433832795
#endif

// useful values for some basic spherical harmonic representations
// Y00_1 = 1/Y00 = spherical harmonic representation of 1 (l=0,m=0)
#define Y00_1 sqrt(4.*M_PI)
// Y10_ct = spherical harmonic representation of cos(theta) (l=1,m=0)
#define Y10_ct sqrt(4.*M_PI/3.)

enum shtns_type {
	sht_gauss,	// use gaussian grid and quadrature. highest accuracy.
	sht_auto,	// use a regular grid if dct is faster with goog accuracy, otherwise defaults to gauss.
	sht_reg_fast,	// use fastest algorithm, on a regular grid, mixing dct and regular quadrature.
	sht_reg_dct,	// use pure dct algorithm, on a regular grid.
	sht_reg_poles	// use a synthesis only algo including poles, not suitable for computations.
};

/// FUNCTIONS ///

// initialize everything.
void init_SH(enum shtns_type flags, double eps, int lmax, int mmax, int mres, int nlat, int nphi);
long int nlm_calc(long int lmax, long int mmax, long int mres);
void Set_MTR_DCT(int m);
int Get_MTR_DCT();

void spat_to_SH(complex double *BrF, complex double *Qlm);
void SH_to_spat(complex double *Qlm, complex double *BrF);

#define SH_to_grad_spat(S,Gt,Gp) SHsph_to_spat(S, Gt, Gp)
void SHsphtor_to_spat(complex double *Slm, complex double *Tlm, complex double *BtF, complex double *BpF);
void SHsph_to_spat(complex double *Slm, complex double *BtF, complex double *BpF);
void SHtor_to_spat(complex double *Tlm, complex double *BtF, complex double *BpF);

void spat_to_SHsphtor(complex double *BtF, complex double *BpF, complex double *Slm, complex double *Tlm);

double SH_to_point(complex double *Qlm, double cost, double phi);
void SHqst_to_point(complex double *Qlm, complex double *Slm, complex double *Tlm, double cost, double phi,
					   double *vr, double *vt, double *vp);
					   
void spat_to_SH_l(complex double *BrF, complex double *Qlm, int LTR);
void SH_to_spat_l(complex double *Qlm, complex double *BrF, int LTR);

#define SH_to_grad_spat_l(S,Gt,Gp,ltr) SHsph_to_spat(S, Gt, Gp, ltr)
void SHsphtor_to_spat_l(complex double *Slm, complex double *Tlm, complex double *BtF, complex double *BpF, int LTR);
void SHsph_to_spat_l(complex double *Slm, complex double *BtF, complex double *BpF, int LTR);
void SHtor_to_spat_l(complex double *Tlm, complex double *BtF, complex double *BpF, int LTR);
void spat_to_SHsphtor_l(complex double *BtF, complex double *BpF, complex double *Slm, complex double *Tlm, int LTR);
