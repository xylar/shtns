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

/** \file shtns.h
 \brief shtns.h is the definition file for SHTns : include this file in your source code to use SHTns.
**/

typedef struct shtns_info* shtns_cfg;

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
	sht_reg_poles,	///< use a <b>synthesis only</b> algo <b>including poles</b>, not suitable for computations. Useful for vizualisation.
	sht_gauss_fly	///< legendre polynomials are recomputed on-the-fly for each transform (may be faster on some machines, saves memory and bandwidth).
};

#ifndef SHTNS_PRIVATE
struct shtns_info {		// allow read-only access to some data (useful for optimization and helper macros)
	const int nlm;					///< total number of (l,m) spherical harmonics components.
	const unsigned short lmax;		///< maximum degree (lmax) of spherical harmonics.
	const unsigned short mmax;		///< maximum order (mmax*mres) of spherical harmonics.
	const unsigned short mres;		///< the periodicity along the phi axis.
	const unsigned short nphi;		///< number of spatial points in Phi direction (longitude)
	const unsigned short nlat;		///< number of spatial points in Theta direction (latitude) ...
	const unsigned short nlat_2;		///< ...and half of it (using (shtns.nlat+1)/2 allows odd shtns.nlat.)
	const int *const lmidx;					///< (virtual) index in SH array of given im : LiM(l,im) = lmidx[im] + l
	const unsigned short *const li;			///< degree l for given mode number : li[lm]
	const double *const el;			///< l, l(l+1) and 1/(l(l+1)) arrays.
	const double *const l2;			
	const double *const l_2;		
	const double *const ct;			///< cos(theta) and sin(theta) arrays.
	const double *const st;			
};
#endif

// MACROS //

/*! \name Access to spherical harmonic components
 * The following macros give access to single spherical harmonic coefficient or perform loops spanning all of them.
**///@{
///LiM(l,im) : macro returning array index for given l and im.
#define LiM(shtns, l,im) ( shtns->lmidx[im] + l )
/// LM(l,m) : macro returning array index for given l and m.
#define LM(shtns, l,m) ( shtns->lmidx[(m)/shtns->mres] + l )
/// LM_LOOP( action ) : macro that performs "action" for every (l,m), with lm set, but neither l, m nor im.
/// \c lm must be a declared int and is the loop counter and the SH array index. more info : \ref spec_data
#define LM_LOOP( shtns, action ) { int lm=0; do { action } while(++lm < shtns->nlm); }
/// LM_L_LOOP : loop over all (l,im) and perform "action"  : l and lm are defined (but NOT m and im).
/// \c lm and \c m must be declared int's. \c lm is the loop counter and SH array index, while \c l is the SH degree. more info : \ref spec_data
#define LM_L_LOOP( shtns, action ) { int lm=0; do { int l=shtns->li[lm]; action } while(++lm < shtns->nlm); }
//@}

#define SHT_NATIVE_LAYOUT 0			///< Tells shtns_init to use \ref native
#define SHT_THETA_CONTIGUOUS 256	///< use \ref theta_fast
#define SHT_PHI_CONTIGUOUS (256*2)	///< use \ref phi_fast

#define SHT_NO_CS_PHASE (256*4)	///< don't include Condon-Shortley phase (add to last argument of \ref shtns_set_size)
#define SHT_REAL_NORM (256*8)	///< use a "real" normalization. (add to last argument of \ref shtns_set_size)

/// total number of 'doubles' required for a spatial field (includes FFTW reserved space).
/// only the first shtns.nlat*shtns.nphi are real spatial data, the remaining is used by the Fourier Transform. more info : \ref spat
#define NSPAT_ALLOC(shtns) (shtns->nlat*(shtns->nphi/2+1)*2)

// HELPER MACROS //

/// phi angle value in degrees for given index ip.
#define PHI_DEG(shtns, ip) (360./(shtns->nphi*shtns->mres))*(ip)
/// phi angle value in radians for given index ip.
#define PHI_RAD(shtns, ip) (2.*M_PI/(shtns->nphi*shtns->mres))*(ip)


// FUNCTIONS //

/// \name initialization
//@{
/// Simple initialization of the spherical harmonic transforms of given size. Calls \ref shtns_set_size and \ref shtns_precompute.
shtns_cfg shtns_init(enum shtns_type flags, int lmax, int mmax, int mres, int nlat, int nphi);
/// defines the sizes of the spectral description. Use for advanced initialization.
shtns_cfg shtns_create(int lmax, int mmax, int mres, enum shtns_norm norm);
/// precompute everything for a given spatial grid. Use for advanced initialization, after \ref shtns_set_size.
int shtns_set_grid(shtns_cfg, enum shtns_type flags, double eps, int nlat, int nphi);
/// precompute everything and choose the optimal nlat and nphi for a given non-linear order.
int shtns_set_grid_auto(shtns_cfg, enum shtns_type flags, double eps, int nl_order, int *nlat, int *nphi);

/// compute number of spherical harmonics modes (l,m) for given size parameters. Does not require a previous call to init_SH
int nlm_calc(int lmax, int mmax, int mres);
void Set_MTR_DCT(shtns_cfg, int m);
int Get_MTR_DCT(shtns_cfg);
//@}

/// \name special values
//@{
double sh00_1(shtns_cfg);	///< returns the spherical harmonic representation of 1 (l=0,m=0)
double sh10_ct(shtns_cfg);	///< returns the spherical harmonic representation of cos(theta) (l=1,m=0)
double sh11_st(shtns_cfg);	///< returns the spherical harmonic representation of sin(theta)*cos(phi) (l=1,m=1)
double shlm_e1(shtns_cfg, int l, int m);		///< returns the l,m SH coefficient corresponding to unit energy.
//@}

/// \name Scalar transforms
//@{
void spat_to_SH(shtns_cfg, double *Vr, complex double *Qlm);
void SH_to_spat(shtns_cfg, complex double *Qlm, double *Vr);
//@}

/// \name Vector transforms
//@{
void spat_to_SHsphtor(shtns_cfg, double *Vt, double *Vp, complex double *Slm, complex double *Tlm);
void SHsphtor_to_spat(shtns_cfg, complex double *Slm, complex double *Tlm, double *Vt, double *Vp);
void SHsph_to_spat(shtns_cfg, complex double *Slm, double *Vt, double *Vp);
void SHtor_to_spat(shtns_cfg, complex double *Tlm, double *Vt, double *Vp);
//@}
/// Compute the spatial representation of the gradient of a scalar SH field. Alias for \ref SHsph_to_spat
#define SH_to_grad_spat(shtns, S,Gt,Gp) SHsph_to_spat(shtns, S, Gt, Gp)

/// \name 3D transforms (combine Scalar and Vector)
//@{
void spat_to_SHqst(shtns_cfg, double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm);
void SHqst_to_spat(shtns_cfg, complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp);
//@}

/// \name Local and partial evalutions of a SH representation :
/// Does not require a call to \ref shtns_precompute
//@{
double SH_to_point(shtns_cfg, complex double *Qlm, double cost, double phi);
void SHqst_to_point(shtns_cfg, complex double *Qlm, complex double *Slm, complex double *Tlm,
					double cost, double phi, double *vr, double *vt, double *vp);

void SHqst_to_lat(shtns_cfg, complex double *Qlm, complex double *Slm, complex double *Tlm, double cost,
					double *vr, double *vt, double *vp, int nphi, int ltr, int mtr);
//@}

/// \name Truncated transforms at given degree l
//@{
void spat_to_SH_l(shtns_cfg, double *Vr, complex double *Qlm, int LTR);
void SH_to_spat_l(shtns_cfg, complex double *Qlm, double *Vr, int LTR);

void SHsphtor_to_spat_l(shtns_cfg, complex double *Slm, complex double *Tlm, double *Vt, double *Vp, int LTR);
void SHsph_to_spat_l(shtns_cfg, complex double *Slm, double *Vt, double *Vp, int LTR);
void SHtor_to_spat_l(shtns_cfg, complex double *Tlm, double *Vt, double *Vp, int LTR);
void spat_to_SHsphtor_l(shtns_cfg, double *Vt, double *Vp, complex double *Slm, complex double *Tlm, int LTR);

void spat_to_SHqst_l(shtns_cfg, double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm, int LTR);
void SHqst_to_spat_l(shtns_cfg, complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp, int LTR);
//@}
/// Compute the spatial representation of the gradient of a scalar SH field. Alias for \ref SHsph_to_spat_l
#define SH_to_grad_spat_l(shtns, S,Gt,Gp,ltr) SHsph_to_spat_l(shtns, S, Gt, Gp, ltr)

/*! \name Axisymmetric transforms m=0 only
 * these work for any MMAX, and will only transform m=0, to/from arrays with NLAT contiguous theta-values.
 *///@{
void spat_to_SH_m0(shtns_cfg, double *Vr, complex double *Qlm);
void SH_to_spat_m0(shtns_cfg, complex double *Qlm, double *Vr);
void SHsphtor_to_spat_m0(shtns_cfg, complex double *Slm, complex double *Tlm, double *Vt, double *Vp);
void SHsph_to_spat_m0(shtns_cfg, complex double *Slm, double *Vt);
void SHtor_to_spat_m0(shtns_cfg, complex double *Tlm, double *Vp);
void spat_to_SHsphtor_m0(shtns_cfg, double *Vt, double *Vp, complex double *Slm, complex double *Tlm);
//@}
/// Compute the spatial representation of the gradient of a scalar SH field. Alias for \ref SHsph_to_spat_m0
#define SH_to_grad_spat_m0(shtns, S,Gt) SHsph_to_spat_m0(shtns, S, Gt)
