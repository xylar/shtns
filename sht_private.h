/********************************************************************
 * SHTns : Spherical Harmonic Transform for numerical simulations.  *
 *    written by Nathanael Schaeffer / LGIT,CNRS                    *
 ********************************************************************/

/// \file sht_private.h private data used by sht transform functions.

#include <complex.h>
#include <math.h>
// FFTW la derivee d/dx = ik	(pas de moins !)
#include <fftw3.h>

#include "sht_config.h"

#include "SHT.h"

// define shortcuts to sizes + allow compile-time optimizations when SHT_AXISYM is defined.
#define LMAX shtns.lmax
#define NLAT shtns.nlat
#define NLAT_2 shtns.nlat_2
#ifndef SHT_AXISYM
  #define NPHI shtns.nphi
  #define MMAX shtns.mmax
  #define MRES shtns.mres
  #define SHT_FFT shtns.sht_fft
#else
  #define NPHI 1
  #define MMAX 0
  #define MRES 1
  #define SHT_FFT 0
#endif
#define NLM shtns.nlm
#define MTR_DCT shtns.mtr_dct

#ifndef M_PI
# define M_PI 3.1415926535897932384626433832795
#endif
#ifndef M_PIl
# define M_PIl 3.1415926535897932384626433832795L
#endif

/*
struct SHTdef {
	long int nlm;
	long int lmax,nlat,nlat_2;
	long int mmax,mres,nphi;
	long int mtr_dct;

	double *ct, *st, *st_1;		// cos(theta), sin(theta), 1/sin(theta);
	double *el, *l2, *l_2;		// l, l(l+1) and 1/(l(l+1))
	int *li;

	int *lmidx;		// (virtual) index in SH array of given im.
	int *tm;		// start theta value for SH (polar optimization : near the poles the legendre polynomials go to zero for high m's)

	double** ylm;		// matrix for inverse transform (synthesis)
	double** zlm;		// matrix for direct transform (analysis)
	double** ykm_dct;	// matrix for inverse transform (synthesis) using dct.
	double* zlm_dct0;	// matrix for direct transform (analysis), only m=0

	fftw_plan ifft, fft;	// plans for FFTW.
	fftw_plan idct, dctm0;
}

struct VSHTdef {
	long int nlm;
	long int lmax,nlat,nlat_2;
	long int mmax,mres,nphi;
	long int mtr_dct;

	double *ct, *st, *st_1;		// cos(theta), sin(theta), 1/sin(theta);
	double *el, *l2, *l_2;		// l, l(l+1) and 1/(l(l+1))
	int *li;

	int *lmidx;		// (virtual) index in SH array of given im.
	int *tm;		// start theta value for SH (polar optimization : near the poles the legendre polynomials go to zero for high m's)

	struct DtDp** dylm;	// theta and phi derivative of Ylm matrix
	struct DtDp** dzlm;
	struct DtDp** dykm_dct;	// theta and phi derivative of Ylm matrix

	fftw_plan ifft, fft;	// plans for FFTW.
	fftw_plan idct;

	spat_to_SH
	SH_to_spat
}
*/

struct DtDp {		// theta and phi derivatives stored together.
	double t, p;
};

extern int *tm;			// start theta value for SH (polar optimization : near the poles the legendre polynomials go to zero for high m's)
extern double** ylm;		// matrix for inverse transform (synthesis)
extern struct DtDp** dylm;	// theta and phi derivative of Ylm matrix
extern double** zlm;		// matrix for direct transform (analysis)
extern struct DtDp** dzlm;

extern double** ykm_dct;	// matrix for inverse transform (synthesis) using dct.
extern struct DtDp** dykm_dct;	// theta and phi derivative of Ylm matrix
extern double* zlm_dct0;	// matrix for direct transform (analysis), only m=0
extern double* dzlm_dct0;

extern fftw_plan ifft, fft;		// plans for FFTW.
extern fftw_plan idct, dct_m0;			// (I)DCT for NPHI>1
extern fftw_plan idct_r1, dct_r1;		// (I)DCT for axisymmetric case, NPHI=1
extern fftw_plan ifft_eo, fft_eo;		// for half the size (given parity)

#define SSE __attribute__((aligned (16)))

