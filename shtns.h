/** \file shtns.h
 \brief shtns.h is the definition file for SHTns : include this file in your source code to use SHTns.
**/

/// different Spherical Harmonic normalizations.
/// see also section \ref norm for details.
enum shtns_norm {
	sht_orthonormal,	///< orthonormalized spherical harmonics (default).
	sht_fourpi,			///< Geodesy and spectral analysis : 4.pi normalization.
	sht_schmidt			///< Schmidt semi-normalized : 4.pi/(2l+1)
};

/// different SHT types and algorithms
enum shtns_type {
	sht_gauss,	///< use <b>gaussian grid</b> and quadrature. highest accuracy.
	sht_auto,	///< use a regular grid if dct is faster with goog accuracy, otherwise defaults to gauss.
	sht_reg_fast,	///< use fastest algorithm, on a <b>regular grid</b>, mixing dct and regular quadrature.
	sht_reg_dct,	///< use pure dct algorithm, on a <b>regular grid</b>.
	sht_quick_init, ///< gauss grid, with minimum initialization time (useful for pre/post-processing)
	sht_reg_poles	///< use a <b>synthesis only</b> algo <b>including poles</b>, not suitable for computations. Useful for vizualisation.
};

/// structure containing useful information about the SHT.
struct sht_sze {
	int lmax;		///< maximum degree (lmax) of spherical harmonics.
	int mmax;		///< maximum order (mmax*mres) of spherical harmonics.
	int mres;		///< the periodicity along the phi axis.
	int nlm;		///< total number of (l,m) spherical harmonics components.

	int nlat;		///< number of spatial points in Theta direction (latitude) ...
	int nphi;		///< number of spatial points in Phi direction (longitude)

	int *lmidx;		///< (virtual) index in SH array of given im.

	int mtr_dct;	///< m truncation for dct. -1 means no dct at all.
	int nlat_2;		///< ...and half of it (using (shtns.nlat+1)/2 allows odd shtns.nlat.)
	int sht_fft;	///< How to perform fft : 0=no fft, 1=in-place, 2=out-of-place.
	int klim;		///< Limit to k for non-linear terms.

	int norm;		///< store the normalization of the Spherical Harmonics (enum \ref shtns_norm + \ref SHT_NO_CS_PHASE flag)
};

// GLOBAL VARIABLES : do not modify, they are set by the call to init_SH. //
extern struct sht_sze shtns;

/*! \name physical space coordinates arrays
 * functions of the co-latitude theta for latitudinal grid points [0..shtns.nlat-1]
 *///@{
extern double *ct;	///< ct[i] = cos(theta)
extern double *st;	///< st[i] = sin(theta)
extern double *st_1;	///< st_1[i] = 1/sin(theta)
//@}

/*! \name spherical harmonic space arrays
 * useful combinations of the \b degree \b l of the spherical harmonic stored at offset \c lm in the SH array.
 *///@{
extern int *li;		///< li[lm] is the \b degree \b l of the spherical harmonic coefficient stored in lm position (integer array)
extern double *el;	///< el[lm] is the \b degree \b l of the spherical harmonic coefficient stored in lm position (double array)
extern double *l2;	///< l2[lm] = l(l+1) 
extern double *l_2;	///< l_2[lm] = 1/(l(l+1))
//@}

// MACROS //

/*! \name Access to spherical harmonic components
 * The following macros give access to single spherical harmonic coefficient or perform loops spanning all of them.
**///@{
///LiM(l,im) : macro returning array index for given l and im.
#define LiM(l,im) ( shtns.lmidx[im] + l )
/// LM(l,m) : macro returning array index for given l and m.
#define LM(l,m) ( shtns.lmidx[(m)/shtns.mres] + l )
/// LM_LOOP( action ) : macro that performs "action" for every (l,m), with lm set, but neither l, m nor im.
/// \c lm must be a declared int and is the loop counter and the SH array index. more info : \ref spec_data
#define LM_LOOP( action ) for (lm=0; lm<shtns.nlm; lm++) { action }
/// LM_L_LOOP : loop over all (l,im) and perform "action"  : l and lm are defined (but NOT m and im).
/// \c lm and \c m must be declared int's. \c lm is the loop counter and SH array index, while \c l is the SH degree. more info : \ref spec_data
#define LM_L_LOOP( action ) for (lm=0; lm<shtns.nlm; lm++) { l=li[lm]; { action } }
//@}

#define SHT_NATIVE_LAYOUT 0			///< Tells shtns_init to use \ref native
#define SHT_THETA_CONTIGUOUS 256	///< use \ref theta_fast
#define SHT_PHI_CONTIGUOUS (256*2)	///< use \ref phi_fast

#define SHT_NO_CS_PHASE (256*4)	///< don't include Condon-Shortley phase (add to last argument of \ref shtns_set_size)
#define SHT_REAL_NORM (256*8)	///< use a "real" normalization. (add to last argument of \ref shtns_set_size)

/// total number of 'doubles' required for a spatial field (includes FFTW reserved space).
/// only the first shtns.nlat*shtns.nphi are real spatial data, the remaining is used by the Fourier Transform. more info : \ref spat
#define NSPAT_ALLOC (shtns.nlat*(shtns.nphi/2+1)*2)

// FUNCTIONS //

/// \name initialization
//@{
/// Simple initialization of the spherical harmonic transforms of given size. Calls \ref shtns_set_size and \ref shtns_precompute.
int shtns_init(enum shtns_type flags, int lmax, int mmax, int mres, int nlat, int nphi);
/// defines the sizes of the spectral description. Use for advanced initialization.
int shtns_set_size(int lmax, int mmax, int mres, enum shtns_norm norm);
/// precompute everything for a given spatial grid. Use for advanced initialization, after \ref shtns_set_size.
int shtns_precompute(enum shtns_type flags, double eps, int nlat, int nphi);
/// precompute everything and choose the optimal nlat and nphi for a given non-linear order.
int shtns_precompute_auto(enum shtns_type flags, double eps, int nl_order, int *nlat, int *nphi);

/// compute number of spherical harmonics modes (l,m) for given size parameters. Does not require a previous call to init_SH
int nlm_calc(int lmax, int mmax, int mres);
void Set_MTR_DCT(int m);
int Get_MTR_DCT();
//@}

/// \name special values
//@{
double sh00_1();	///< returns the spherical harmonic representation of 1 (l=0,m=0)
double sh10_ct();	///< returns the spherical harmonic representation of cos(theta) (l=1,m=0)
double sh11_st();	///< returns the spherical harmonic representation of sin(theta)*cos(phi) (l=1,m=1)
//@}

/// \name Scalar transforms
//@{
void spat_to_SH(double *Vr, complex double *Qlm);
void SH_to_spat(complex double *Qlm, double *Vr);
//@}

/// \name Vector transforms
//@{
void spat_to_SHsphtor(double *Vt, double *Vp, complex double *Slm, complex double *Tlm);
void SHsphtor_to_spat(complex double *Slm, complex double *Tlm, double *Vt, double *Vp);
void SHsph_to_spat(complex double *Slm, double *Vt, double *Vp);
void SHtor_to_spat(complex double *Tlm, double *Vt, double *Vp);
//@}
/// Compute the spatial representation of the gradient of a scalar SH field. Alias for \ref SHsph_to_spat
#define SH_to_grad_spat(S,Gt,Gp) SHsph_to_spat(S, Gt, Gp)

/// \name Local and partial evalutions of a SH representation :
/// Does not require a call to \ref shtns_precompute
//@{
double SH_to_point(complex double *Qlm, double cost, double phi);
void SHqst_to_point(complex double *Qlm, complex double *Slm, complex double *Tlm,
					double cost, double phi, double *vr, double *vt, double *vp);

void SHqst_to_lat(complex double *Qlm, complex double *Slm, complex double *Tlm, double cost,
					double *vr, double *vt, double *vp, int nphi, int ltr, int mtr);
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
/// Compute the spatial representation of the gradient of a scalar SH field. Alias for \ref SHsph_to_spat_l
#define SH_to_grad_spat_l(S,Gt,Gp,ltr) SHsph_to_spat_l(S, Gt, Gp, ltr)

/*! \name Axisymmetric transforms m=0 only
 * these work for any MMAX, and will only transform m=0, to/from arrays with NLAT contiguous theta-values.
 *///@{
void spat_to_SH_m0(double *Vr, complex double *Qlm);
void SH_to_spat_m0(complex double *Qlm, double *Vr);
void SHsphtor_to_spat_m0(complex double *Slm, complex double *Tlm, double *Vt, double *Vp);
void SHsph_to_spat_m0(complex double *Slm, double *Vt);
void SHtor_to_spat_m0(complex double *Tlm, double *Vp);
void spat_to_SHsphtor_m0(double *Vt, double *Vp, complex double *Slm, complex double *Tlm);
//@}
/// Compute the spatial representation of the gradient of a scalar SH field. Alias for \ref SHsph_to_spat_m0
#define SH_to_grad_spat_m0(S,Gt) SHsph_to_spat_m0(S, Gt)

/*! \name SHT transforms with assumed equatorial symmetry (parity = 0 or 1)
 * these work with (NLAT+1)/2 latitudinal points, and do not overwrite SH coefficients of other parity
 *///@{
void SHeo_to_spat(complex double *Qlm, double *Vr, int parity);
void spat_to_SHeo(double *Vr, complex double *Qlm, int parity);
void SHeo_sphtor_to_spat(complex double *Slm, complex double *Tlm, double *Vt, double *Vp, int parity);
void spat_to_SHeo_sphtor(double *Vt, double *Vp, complex double *Slm, complex double *Tlm, int parity);
//@}
