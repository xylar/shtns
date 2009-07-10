/** \file SHT.h
 \brief SHT.h is the definition file for SHTns : include this file in your source code to use SHTns.
**/

// GLOBAL VARIABLES : do not modify, they are set by the call to init_SH. //

/*! \name size variables
 * The spherical harmonic coefficients are stored in a one-dimensional array of size NLM (the SH array).
 * The following variables are used to record the size of the physical space fields and its spherical harmonic description.
 * 
 * <b>DO NOT MODIFY</b> these variables ! They are set once for all by \ref init_SH.
 */
//@{
extern long int LMAX;	///< maximum degree (LMAX) of spherical harmonics.
extern long int NLAT;	///< number of spatial points in Theta direction (latitude).
#ifndef SHT_AXISYM
  extern long int MRES;	///< \c 2.pi/MRES is the periodicity along the phi coord. (MRES=1 means no assumed periodicity).
  extern long int MMAX;	///< maximum number of azimutal modes. \c (MMAX*MRES) is the maximum order of spherical harmonics. 
  extern long int NPHI;	///< number of spatial points in Phi direction (longitude).
#else
  #define MMAX 0
  #define NPHI 1
  #define MRES 1
#endif
extern long int NLM;	///< total number of (l,m) spherical harmonics components (complex double).
//@}

/// total number of 'doubles' required for a spatial field (includes FFTW padding).
/// only the first NLAT*NPHI are real spatial data, the remaining is used by the Fourier Transform. more info : \ref spat_data
#define NSPAT_ALLOC (NLAT*(NPHI/2+1)*2)

/*! \name physical space coordinates arrays
 * functions of the co-latitude theta for latitudinal grid points [0..NLAT-1]
 */
//@{
extern double *ct;	///< ct[i] = cos(theta)
extern double *st;	///< st[i] = sin(theta)
extern double *st_1;	///< st_1[i] = 1/sin(theta)
//@}
/*! \name spherical harmonic space arrays
 * useful combinations of the \b degree \b l of the spherical harmonic stored at offset \c lm in the SH array.
 */
//@{
extern int *li;		///< li[lm] is the \b degree \b l of the spherical harmonic coefficient stored in lm position (integer array)
extern double *el;	///< el[lm] is the \b degree \b l of the spherical harmonic coefficient stored in lm position (double array)
extern double *l2;	///< l2[lm] = l(l+1) 
extern double *l_2;	///< l_2[lm] = 1/(l(l+1))
//@}

/// required by some macros. Do not use directly.
extern long int *lmidx;

// MACROS //

/*! \name Access to spherical harmonic components
 * The following macros give access to single spherical harmonic coefficient or perform loops spanning all of them.
**/
//@{
///LiM(l,im) : macro returning array index for given l and im.
#define LiM(l,im) ( lmidx[im] + l )
/// LM(l,m) : macro returning array index for given l and m.
#define LM(l,m) ( lmidx[(m)/MRES] + l )
/// LM_LOOP( action ) : macro that performs "action" for every (l,m), with lm set, but neither l, m nor im.
/// \c lm must be a declared int and is the loop counter and the SH array index. more info : \ref spec_data
#define LM_LOOP( action ) for (lm=0; lm<NLM; lm++) { action }
/// LM_L_LOOP : loop over all (l,im) and perform "action"  : l and lm are defined (but NOT m and im).
/// \c lm and \c m must be declared int's. \c lm is the loop counter and SH array index, while \c l is the SH degree. more info : \ref spec_data
#define LM_L_LOOP( action ) for (lm=0; lm<NLM; lm++) { l=li[lm]; { action } }
//@}

#ifndef M_PI
  #define M_PI 3.1415926535897932384626433832795
#endif

// useful values for some basic spherical harmonic representations
/// Y00_1 = \f$ 1/Y00 = \sqrt{4 \pi} \f$ spherical harmonic representation of 1 (l=0,m=0)
#define Y00_1 sqrt(4.*M_PI)
/// Y10_ct = \f$ 1/Y10 = \sqrt{4 \pi /3} \f$ spherical harmonic representation of cos(theta) (l=1,m=0)
#define Y10_ct sqrt(4.*M_PI/3.)
/// Y11_st = \f$ 1/Y11 = \sqrt{2 \pi /3} \f$ spherical harmonic representation of sin(theta)*cos(phi) (l=1,m=1)
#define Y11_st sqrt(2.*M_PI/3.)

/// different SHT types and algorithms
enum shtns_type {
	sht_gauss,	///< use gaussian grid and quadrature. highest accuracy.
	sht_auto,	///< use a regular grid if dct is faster with goog accuracy, otherwise defaults to gauss.
	sht_reg_fast,	///< use fastest algorithm, on a regular grid, mixing dct and regular quadrature.
	sht_reg_dct,	///< use pure dct algorithm, on a regular grid.
	sht_reg_poles	///< use a synthesis only algo including poles, not suitable for computations.
};

// FUNCTIONS //

/// Initializes spherical harmonic transforms of given size, and sets all global variables.
void init_SH(enum shtns_type flags, double eps, int lmax, int mmax, int mres, int nlat, int nphi);

/// compute number of spherical harmonics modes (l,m) for given size parameters. Does not require a previous call to init_SH
long int nlm_calc(long int lmax, long int mmax, long int mres);
void Set_MTR_DCT(int m);
int Get_MTR_DCT();

/// \name Scalar transforms
//@{
void spat_to_SH(complex double *BrF, complex double *Qlm);
void SH_to_spat(complex double *Qlm, complex double *BrF);
//@}

/// \name Vector transforms
//@{
void SHsphtor_to_spat(complex double *Slm, complex double *Tlm, complex double *BtF, complex double *BpF);
void SHsph_to_spat(complex double *Slm, complex double *BtF, complex double *BpF);
void SHtor_to_spat(complex double *Tlm, complex double *BtF, complex double *BpF);

void spat_to_SHsphtor(complex double *BtF, complex double *BpF, complex double *Slm, complex double *Tlm);
//@}
#define SH_to_grad_spat(S,Gt,Gp) SHsph_to_spat(S, Gt, Gp)

/// \name Evalution of a SH representation at a given point in physical space.
//@{
double SH_to_point(complex double *Qlm, double cost, double phi);
void SHqst_to_point(complex double *Qlm, complex double *Slm, complex double *Tlm, double cost, double phi,
					   double *vr, double *vt, double *vp);
//@}

/// \name Truncated transforms at given degree l
//@{
void spat_to_SH_l(complex double *BrF, complex double *Qlm, int LTR);
void SH_to_spat_l(complex double *Qlm, complex double *BrF, int LTR);

void SHsphtor_to_spat_l(complex double *Slm, complex double *Tlm, complex double *BtF, complex double *BpF, int LTR);
void SHsph_to_spat_l(complex double *Slm, complex double *BtF, complex double *BpF, int LTR);
void SHtor_to_spat_l(complex double *Tlm, complex double *BtF, complex double *BpF, int LTR);
void spat_to_SHsphtor_l(complex double *BtF, complex double *BpF, complex double *Slm, complex double *Tlm, int LTR);
//@}
#define SH_to_grad_spat_l(S,Gt,Gp,ltr) SHsph_to_spat(S, Gt, Gp, ltr)
