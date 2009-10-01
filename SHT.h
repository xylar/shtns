/** \file SHT.h
 \brief SHT.h is the definition file for SHTns : include this file in your source code to use SHTns.
**/

// GLOBAL VARIABLES : do not modify, they are set by the call to init_SH. //

struct sht_sze {
	long int lmax;			///< maximum degree (lmax) of spherical harmonics.
	long int mmax;			///< maximum order (mmax*mres) of spherical harmonics.
	long int mres;			///< the periodicity along the phi axis.
	long int nlat;			///< number of spatial points in Theta direction (latitude) ...
	long int nlat_2;		///< ...and half of it (using (shtns.nlat+1)/2 allows odd shtns.nlat.)
	long int nphi;			///< number of spatial points in Phi direction (longitude)
	long int nlm;			///< total number of (l,m) spherical harmonics components.
	long int mtr_dct;		///< m truncation for dct. -1 means no dct at all.
};
extern struct sht_sze shtns;

/// total number of 'doubles' required for a spatial field (includes FFTW padding).
/// only the first shtns.nlat*shtns.nphi are real spatial data, the remaining is used by the Fourier Transform. more info : \ref spat_data
#define NSPAT_ALLOC (shtns.nlat*(shtns.nphi/2+1)*2)

/*! \name physical space coordinates arrays
 * functions of the co-latitude theta for latitudinal grid points [0..shtns.nlat-1]
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
#define LM(l,m) ( lmidx[(m)/shtns.mres] + l )
/// LM_LOOP( action ) : macro that performs "action" for every (l,m), with lm set, but neither l, m nor im.
/// \c lm must be a declared int and is the loop counter and the SH array index. more info : \ref spec_data
#define LM_LOOP( action ) for (lm=0; lm<shtns.nlm; lm++) { action }
/// LM_L_LOOP : loop over all (l,im) and perform "action"  : l and lm are defined (but NOT m and im).
/// \c lm and \c m must be declared int's. \c lm is the loop counter and SH array index, while \c l is the SH degree. more info : \ref spec_data
#define LM_L_LOOP( action ) for (lm=0; lm<shtns.nlm; lm++) { l=li[lm]; { action } }
//@}

#ifndef M_PI
  #define M_PI 3.1415926535897932384626433832795
#endif

// useful values for some basic spherical harmonic representations
/// Y00_1 = \f$ 1/Y_0^0 = \sqrt{4 \pi} \f$ spherical harmonic representation of 1 (l=0,m=0)
#define Y00_1 sqrt(4.*M_PI)
/// Y10_ct = \f$ \cos\theta/Y_1^0 = \sqrt{4 \pi /3} \f$ spherical harmonic representation of cos(theta) (l=1,m=0)
#define Y10_ct sqrt(4.*M_PI/3.)
/// Y11_st = \f$ \sin\theta\cos\phi/(Y_1^1 + Y_1^{-1}) = \sqrt{2 \pi /3} \f$ spherical harmonic representation of sin(theta)*cos(phi) (l=1,m=1)
#define Y11_st -sqrt(2.*M_PI/3.)

/// different SHT types and algorithms
enum shtns_type {
	sht_gauss,	///< use <b>gaussian grid</b> and quadrature. highest accuracy.
	sht_auto,	///< use a regular grid if dct is faster with goog accuracy, otherwise defaults to gauss.
	sht_reg_fast,	///< use fastest algorithm, on a <b>regular grid</b>, mixing dct and regular quadrature.
	sht_reg_dct,	///< use pure dct algorithm, on a <b>regular grid</b>.
	sht_quick_init, ///< gauss grid, with minimum initialization time (useful for pre/post-processing)
	sht_reg_poles	///< use a <b>synthesis only</b> algo <b>including poles</b>, not suitable for computations. \ref shtns.nlat odd is supported even if \link compil SHT_shtns.nlat_EVEN \endlink is defined, useful for vizualisation.
};

#define SHT_NATIVE_LAYOUT 0		///< Tells shtns_init to use \ref native
#define SHT_THETA_CONTIGUOUS 256	///< use \ref theta_fast
#define SHT_PHI_CONTIGUOUS 256*2	///< use \ref phi_fast

// FUNCTIONS //

/// Initializes spherical harmonic transforms of given size, and sets all global variables.
int shtns_init(enum shtns_type flags, double eps, int lmax, int mmax, int mres, int nlat, int nphi);

/// compute number of spherical harmonics modes (l,m) for given size parameters. Does not require a previous call to init_SH
long int nlm_calc(long int lmax, long int mmax, long int mres);
void Set_MTR_DCT(int m);
int Get_MTR_DCT();

/// \name Scalar transforms
//@{
void spat_to_SH(double *Vr, complex double *Qlm);
void SH_to_spat(complex double *Qlm, double *Vr);
//@}

/// \name Vector transforms
//@{
void SHsphtor_to_spat(complex double *Slm, complex double *Tlm, double *Vt, double *Vp);
void SHsph_to_spat(complex double *Slm, double *Vt, double *Vp);
void SHtor_to_spat(complex double *Tlm, double *Vt, double *Vp);

void spat_to_SHsphtor(double *Vt, double *Vp, complex double *Slm, complex double *Tlm);
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
void spat_to_SH_l(double *Vr, complex double *Qlm, int LTR);
void SH_to_spat_l(complex double *Qlm, double *Vr, int LTR);

void SHsphtor_to_spat_l(complex double *Slm, complex double *Tlm, double *Vt, double *Vp, int LTR);
void SHsph_to_spat_l(complex double *Slm, double *Vt, double *Vp, int LTR);
void SHtor_to_spat_l(complex double *Tlm, double *Vt, double *Vp, int LTR);
void spat_to_SHsphtor_l(double *Vt, double *Vp, complex double *Slm, complex double *Tlm, int LTR);
//@}
#define SH_to_grad_spat_l(S,Gt,Gp,ltr) SHsph_to_spat(S, Gt, Gp, ltr)
