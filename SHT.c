/*
 * Copyright (c) 2010 Centre National de la Recherche Scientifique.
 * written by Nathanael Schaeffer (CNRS, LGIT, Grenoble, France).
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

/** \internal \file SHT.c
 * \brief main source file for SHTns.
 * This files contains initialization code and also some partial transforms (point or latitudinal evaluations)
 */

#include <stdio.h>
#include <stdlib.h>

// global variables definitions
#include "sht_private.h"

// cycle counter from FFTW
#include "cycle.h"

/// sizes for SHT, in a structure to avoid problems with conflicting names.
struct sht_sze shtns;

double *al0;	double **alm;	// coefficient list for Legendre function recurrence (size 2*NLM)
double **dlm;	// coefficient list for Legendre function and theta derivative recurrence (size 4*NLM)

double *ct, *st, *st_1;		// cos(theta), sin(theta), 1/sin(theta);
double *wg;					// gauss weights (if current grid is a gauss grid).

int *tm;			// start theta value for SH (polar optimization : near the poles the legendre polynomials go to zero for high m's)
double** ylm;		// matrix for inverse transform (synthesis)
double** zlm;		// matrix for direct transform (analysis)
#ifndef SHT_SCALAR_ONLY
struct DtDp** dylm;	// theta and phi derivative of Ylm matrix
struct DtDp** dzlm;
#endif

double** ykm_dct;	// matrix for inverse transform (synthesis) using dct.
double* zlm_dct0;	// matrix for direct transform (analysis), only m=0
#ifndef SHT_SCALAR_ONLY
struct DtDp** dykm_dct;	// theta and phi derivative of Ylm matrix
double* dzlm_dct0;
#endif

double *el, *l2, *l_2;		// l, l(l+1) and 1/(l(l+1))
int *li = NULL;				// used as flag for shtns_set_size

fftw_plan ifft, fft;			// plans for FFTW.
fftw_plan idct, dct_m0;			// (I)DCT for NPHI>1
fftw_plan idct_r1, dct_r1;		// (I)DCT for axisymmetric case, NPHI=1
fftw_plan ifft_eo, fft_eo;		// for half the size (given parity)
unsigned fftw_plan_mode = FFTW_EXHAUSTIVE;		// defines the default FFTW planner mode.



/// Abort program with error message.
void shtns_runerr(const char * error_text)
{
	printf("*** [SHTns] Run-time error : %s\n",error_text);
	exit(1);
}

double Y00_1, Y10_ct, Y11_st;
/// returns the l=0, m=0 SH coefficient corresponding to a uniform value of 1.
double sh00_1() {
	return Y00_1;
}
/// returns the l=1, m=0 SH coefficient corresponding to cos(theta).
double sh10_ct() {
	return Y10_ct;
}
/// returns the l=1, m=1 SH coefficient corresponding to sin(theta).cos(phi).
double sh11_st() {
	return Y11_st;
}
/// returns the l,m SH coefficient corresponding to unit energy.
double shlm_e1(int l, int m) {
	double x = Y00_1/sqrt(4.*M_PI);
	if (SHT_NORM == sht_schmidt) x *= sqrt(2*l+1);
	if ((m!=0)&&((shtns.norm & SHT_REAL_NORM)==0)) x *= sqrt(0.5);
	return(x);
}


/*  LEGENDRE FUNCTIONS  */
#include "sht_legendre.c"


/// return the smallest power of 2 larger than n.
int next_power_of_2(int n)
{
	int f = 1;
	if ( (n<=0) || (n>(1<<(sizeof(int)*8-2))) ) return 0;
	while (f<n) f*=2;
	return f;
}

/// find the closest integer that is larger than n and that contains only factors up to fmax.
/// fmax is 7 for optimal FFTW fourier transforms.
/// return only even integers for n>fmax.
int fft_int(int n, int fmax)
{
	int k,f;

	if (n<=fmax) return n;
	if (fmax<2) return 0;
	if (fmax==2) return next_power_of_2(n);

	n -= 2-(n&1);		// only even n
	do {
		n+=2;	f=2;
		while ((2*f <= n) && ((n&f)==0)) f *= 2;		// no divisions for factor 2.
		k=3;
		while ((k<=fmax) &&  (k*f <= n)) {
			while ((k*f <= n) && (n%(k*f)==0)) f *= k;
			k+=2;
		}
	} while (f != n);

	k = next_power_of_2(n);			// what is the closest power of 2 ?
	if ((k-n)*33 < n) return k;		// rather choose power of 2 if not too far (3%)

	return n;
}

/*
	SHT FUNCTIONS
*/

// truncation at LMAX and MMAX
#define LTR LMAX
#define MTR MMAX

/** \addtogroup local Local and partial evaluation of SH fields.
 * These do only require a call to \ref shtns_set_size, but not to \ref shtns_precompute.
*/
//@{

/// Evaluate scalar SH representation \b Qlm at physical point defined by \b cost = cos(theta) and \b phi
double SH_to_point(complex double *Qlm, double cost, double phi)
{
	double yl[LMAX+1];
	double vr;
	complex double *Ql;
	long int l,m,im;

	vr = 0.0;
	m=0;	im=0;
		legendre_sphPlm_array(LTR, im, cost, &yl[m]);
		for (l=m; l<=LTR; l++)
			vr += yl[l] * creal( Qlm[l] );
	if (MTR>0) {
		complex double eip, eimp;
		eip = cos(phi*MRES) + I*sin(phi*MRES);	eimp = 2.0;
		for (im=1; im<=MTR; im++) {
			m = im*MRES;
			legendre_sphPlm_array(LTR, im, cost, &yl[m]);
//			eimp = 2.*(cos(m*phi) + I*sin(m*phi));
			eimp *= eip;			// not so accurate, but it should be enough for rendering uses.
			Ql = &Qlm[LiM(0,im)];	// virtual pointer for l=0 and im
			for (l=m; l<=LTR; l++)
				vr += yl[l] * creal( Ql[l]*eimp );
		}
	}
	return vr;
}

/// Evaluate vector SH representation \b Qlm at physical point defined by \b cost = cos(theta) and \b phi
void SHqst_to_point(complex double *Qlm, complex double *Slm, complex double *Tlm, double cost, double phi,
					   double *vr, double *vt, double *vp)
{
	double yl[LMAX+1];
	double dtyl[LMAX+1];
	double vst, vtt, vsp, vtp, vr0, vrm;
	complex double *Ql, *Sl, *Tl;
	long int l,m,im;

	const double sint = sqrt((1.-cost)*(1.+cost));
	vst = 0.; vtt = 0.; vsp = 0.; vtp =0.; vr0 = 0.; vrm = 0.;
	m=0;	im=0;
		legendre_sphPlm_deriv_array(LTR, im, cost, sint, &yl[m], &dtyl[m]);
		for (l=m; l<=LTR; l++) {
			vr0 += yl[l] * creal( Qlm[l] );
			vst += dtyl[l] * creal( Slm[l] );
			vtt += dtyl[l] * creal( Tlm[l] );
		}
	if (MTR>0) {
		complex double eip, eimp, imeimp;
		eip = cos(phi*MRES) + I*sin(phi*MRES);	eimp = 2.0;
		for (im=1; im<=MTR; im++) {
			m = im*MRES;
			legendre_sphPlm_deriv_array(LTR, im, cost, sint, &yl[m], &dtyl[m]);
//			eimp = 2.*(cos(m*phi) + I*sin(m*phi));
			eimp *= eip;		// not so accurate, but it should be enough for rendering uses.
			imeimp = I*m*eimp;
			Ql = &Qlm[LiM(0,im)];	Sl = &Slm[LiM(0,im)];	Tl = &Tlm[LiM(0,im)];
			for (l=m; l<=LTR; l++) {
				vrm += yl[l] * creal( Ql[l]*eimp );
				vst += dtyl[l] * creal(Sl[l]*eimp);
				vtt += dtyl[l] * creal(Tl[l]*eimp);
				vsp += yl[l] * creal(Sl[l]*imeimp);
				vtp += yl[l] * creal(Tl[l]*imeimp);
			}
		}
		vrm *= sint;
	}
	*vr = vr0 + vrm;
	*vt = vtp + vst;	// Bt = I.m/sint *T  + dS/dt
	*vp = vsp - vtt;	// Bp = I.m/sint *S  - dT/dt
}
//@}
	
#undef LTR
#undef MTR


/*
	SYNTHESIS AT A GIVEN LATITUDE
	(does not require a previous call to shtns_precompute)
*/

fftw_plan ifft_lat = NULL;		///< fftw plan for SHqst_to_lat
long int nphi_lat = 0;			///< nphi of previous SHqst_to_lat
double* ylm_lat = NULL;
double* dylm_lat;
double ct_lat = 2.0;
double st_lat;

/// synthesis at a given latitude, on nphi equispaced longitude points.
/// vr, vt, and vp arrays must have nphi+2 doubles allocated (fftw requirement).
/// It does not require a previous call to shtns_precompute
/// \ingroup local
void SHqst_to_lat(complex double *Qlm, complex double *Slm, complex double *Tlm, double cost,
					double *vr, double *vt, double *vp, int nphi, int ltr, int mtr)
{
	complex double vst, vtt, vsp, vtp, vrr;
	complex double *vrc, *vtc, *vpc;
	long int m, l, j;

	if (ltr > LMAX) ltr=LMAX;
	if (mtr > MMAX) mtr=MMAX;
	if (mtr*MRES > ltr) mtr=ltr/MRES;
	if (mtr*2*MRES >= nphi) mtr = (nphi-1)/(2*MRES);

	vrc = (complex double *) vr;
	vtc = (complex double *) vt;
	vpc = (complex double *) vp;

	if ((nphi != nphi_lat)||(ifft_lat == NULL)) {
		if (ifft_lat != NULL) fftw_destroy_plan(ifft_lat);
		ifft_lat = fftw_plan_dft_c2r_1d(nphi, vrc, vr, FFTW_ESTIMATE);
		nphi_lat = nphi;
	}
	if (ylm_lat == NULL) {
		ylm_lat = (double *) malloc(sizeof(double)* NLM*2);
		dylm_lat = ylm_lat + NLM;
	}
	if (cost != ct_lat) {		// don't recompute if same latitude (ie equatorial disc rendering)
		st_lat = sqrt((1.-cost)*(1.+cost));	// sin(theta)
		for (m=0,j=0; m<=mtr; m++) {
			legendre_sphPlm_deriv_array(ltr, m, cost, st_lat, &ylm_lat[j], &dylm_lat[j]);
			j += LMAX -m*MRES +1;
		}
	}

	for (m = 0; m<nphi/2+1; m++) {	// init with zeros
		vrc[m] = 0.0;	vtc[m] = 0.0;	vpc[m] = 0.0;
	}
	j=0;
	m=0;
		vrr=0;	vtt=0;	vst=0;
		for(l=m; l<=ltr; l++, j++) {
			vrr += ylm_lat[j] * creal(Qlm[j]);
			vst += dylm_lat[j] * creal(Slm[j]);
			vtt += dylm_lat[j] * creal(Tlm[j]);
		}
		j += (LMAX-ltr);
		vrc[m] = vrr;
		vtc[m] =  vst;	// Vt =   dS/dt
		vpc[m] = -vtt;	// Vp = - dT/dt
	for (m=MRES; m<=mtr*MRES; m+=MRES) {
		vrr=0;	vtt=0;	vst=0;	vsp=0;	vtp=0;
		for(l=m; l<=ltr; l++, j++) {
			vrr += ylm_lat[j] * Qlm[j];
			vst += dylm_lat[j] * Slm[j];
			vtt += dylm_lat[j] * Tlm[j];
			vsp += ylm_lat[j] * Slm[j];
			vtp += ylm_lat[j] * Tlm[j];
		}
		j+=(LMAX-ltr);
		vrc[m] = vrr*st_lat;
		vtc[m] = I*m*vtp + vst;	// Vt = I.m/sint *T  + dS/dt
		vpc[m] = I*m*vsp - vtt;	// Vp = I.m/sint *S  - dT/dt
	}
	fftw_execute_dft_c2r(ifft_lat,vrc,vr);
	fftw_execute_dft_c2r(ifft_lat,vtc,vt);
	fftw_execute_dft_c2r(ifft_lat,vpc,vp);
//	free(ylm_lat);
}


/*
	INITIALIZATION FUNCTIONS
*/

void alloc_SHTarrays()
{
	long int im,m, l0;
	long int lstride;

	l0 = ((NLAT+1)>>1)*2;		// round up to even
	ct = (double *) fftw_malloc(sizeof(double) * l0*3);
	st = ct + l0;		st_1 = ct + 2*l0;
	ylm = (double **) fftw_malloc( sizeof(double *) * (MMAX+1)*3 );
	zlm = ylm + (MMAX+1);		ykm_dct = ylm + (MMAX+1)*2;
  #ifndef SHT_SCALAR_ONLY
	dylm = (struct DtDp **) fftw_malloc( sizeof(struct DtDp *) * (MMAX+1)*3);
	dzlm = dylm + (MMAX+1);		dykm_dct = dylm + (MMAX+1)*2;
  #endif

// Allocate legendre functions lookup tables.
	lstride = (LMAX+1);		lstride += (lstride&1);		// even stride.
	ylm[0] = (double *) fftw_malloc(sizeof(double)* (NLM-(LMAX+1)+lstride)*NLAT_2);
  #ifndef SHT_SCALAR_ONLY
	dylm[0] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* NLM*NLAT_2);
  #endif
	zlm[0] = (double *) fftw_malloc(sizeof(double)* (NLM*NLAT_2 + (NLAT_2 & 1)));
  #ifndef SHT_SCALAR_ONLY
	dzlm[0] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* (NLM-1)*NLAT_2);		// remove l=0
  #endif
	for (im=0; im<MMAX; im++) {
		m = im*MRES;	l0 = (m==0) ? 1 : m;
		if (im>0) lstride = (LMAX+1 -m);
		ylm[im+1] = ylm[im] + NLAT_2*lstride;
		zlm[im+1] = zlm[im] + NLAT_2*(LMAX+1 -m) + ((m==0)*(NLAT_2&1));
  #ifndef SHT_SCALAR_ONLY
		dylm[im+1] = dylm[im] + NLAT_2*(LMAX+1 -m);
		dzlm[im+1] = dzlm[im] + NLAT_2*(LMAX+1 -l0);
  #endif
	}
#if SHT_VERBOSE > 1
	printf("          Memory used for Ylm and Zlm matrices = %.3f Mb x2\n",3.0*sizeof(double)*NLM*NLAT_2/(1024.*1024.));
#endif
}

void free_SHTarrays()
{
#ifndef SHT_SCALAR_ONLY
	fftw_free(dzlm[0]);		
#endif	
	fftw_free(zlm[0]);
#ifndef SHT_SCALAR_ONLY
	fftw_free(dylm[0]);
#endif
	fftw_free(ylm[0]);
#ifndef SHT_SCALAR_ONLY
	fftw_free(dylm);
#endif
	fftw_free(ylm);		fftw_free(ct);
}

/// \code return (mmax+1)*(lmax+1) - mres*(mmax*(mmax+1))/2; \endcode */
/// \ingroup init
int nlm_calc(int lmax, int mmax, int mres)
{
	if (mmax*mres > lmax) mmax = lmax/mres;
	return( (mmax+1)*(lmax+1) - mres*(mmax*(mmax+1))/2 );	// this is wrong if lmax < mmax*mres
}

/*
long int nlm_calc_eo(long int lmax, long int mmax, long int mres) {
	long int im,l,lm;
	for (im=0, lm=0; im<=mmax; im++) {
		if (im*mres <= lmax) lm += (lmax+2-im*mres)/2;
	}
	return lm;
}
*/


/// Generates an equi-spaced theta grid including the poles, for synthesis only.
void EqualPolarGrid()
{
	int j;
	double f;

#if SHT_VERBOSE > 0
	printf("        => using Equaly Spaced Nodes including poles\n");
#endif
// cos theta of latidunal points (equaly spaced in theta)
	f = M_PI/(NLAT-1.0);
	for (j=0; j<NLAT; j++) {
		ct[j] = cos(f*j);
		st[j] = sin(f*j);
		st_1[j] = 1.0/st[j];
	}
#if SHT_VERBOSE > 0
	printf("     !! Warning : only synthesis (inverse transform) supported for this grid !\n");
#endif
}


/// initialize FFTs using FFTW.
/// \param[in] layout defines the spatial layout (see \ref spat).
/// returns the number of double to be allocated for a spatial field.
/// \todo we should time out-of-place and in-place, and chose the fastest.
/// \todo even/odd transform do not work with out-of-place FFT (segfault).
int planFFT(int layout)
{
	complex double *ShF;
	double *Sh;
	int nfft, ncplx, nreal;
	int theta_inc, phi_inc, phi_embed;

	if (NPHI <= 2*MMAX) shtns_runerr("the sampling condition Nphi > 2*Mmax is not met.");

	switch (layout) {
		case SHT_NATIVE_LAYOUT : 	theta_inc=1;  phi_inc=NLAT;  phi_embed=2*(NPHI/2+1);  break;
		case SHT_THETA_CONTIGUOUS :	theta_inc=1;  phi_inc=NLAT;  phi_embed=NPHI;  break;
		default :
		case SHT_PHI_CONTIGUOUS :	phi_inc=1;  theta_inc=NPHI;  phi_embed=NPHI;  break;
	}

  if (NPHI>1) {
	SHT_FFT = 1;		// yes, do some fft
	nfft = NPHI;
	ncplx = NPHI/2 +1;
	nreal = phi_embed;
	if ((theta_inc != 1)||(phi_inc != NLAT)||(nreal < 2*ncplx)) {
		SHT_FFT = 2;		// we need to do the fft out-of-place.
	}

#if SHT_VERBOSE > 0
	printf("        => using FFTW : Mmax=%d, Nphi=%d, Nlat=%d  (data layout : phi_inc=%d, theta_inc=%d, phi_embed=%d)\n",MMAX,NPHI,NLAT,phi_inc,theta_inc,phi_embed);
	if (NPHI <= (SHT_NL_ORDER+1)*MMAX)	printf("     !! Warning : anti-aliasing condition Nphi > %d*Mmax is not met !\n", SHT_NL_ORDER+1);
	if (NPHI != fft_int(NPHI,7))		printf("     !! Warning : Nphi is not optimal for FFTW !\n");
	if (SHT_FFT > 1) printf("        ** out-of-place fft **\n");
#endif

// Allocate dummy Spatial Fields.
	if (SHT_FFT > 1) {
		ShF = (complex double *) fftw_malloc(ncplx * NLAT * sizeof(complex double));
		Sh = (double *) fftw_malloc(ncplx * NLAT * sizeof(complex double));
	} else {
		ShF = (complex double *) fftw_malloc(ncplx * NLAT * sizeof(complex double));
		Sh = (double *) ShF;
	}

// IFFT : unnormalized.  FFT : must be normalized.
		ifft = fftw_plan_many_dft_c2r(1, &nfft, NLAT, ShF, &ncplx, NLAT, 1, Sh, &nreal, phi_inc, theta_inc, fftw_plan_mode);
		if (ifft == NULL) shtns_runerr("[FFTW] ifft planning failed !");
		fft = fftw_plan_many_dft_r2c(1, &nfft, NLAT, Sh, &nreal, phi_inc, theta_inc, ShF, &ncplx, NLAT, 1, fftw_plan_mode);
		if (fft == NULL) shtns_runerr("[FFTW] fft planning failed !");

		ifft_eo = fftw_plan_many_dft_c2r(1, &nfft, NLAT_2, ShF, &ncplx, NLAT_2, 1, Sh, &nreal, (phi_inc+1)/2, theta_inc, fftw_plan_mode);
		if (ifft_eo == NULL) shtns_runerr("[FFTW] ifft planning failed !");
		fft_eo = fftw_plan_many_dft_r2c(1, &nfft, NLAT_2, Sh, &nreal, (phi_inc+1)/2, theta_inc, ShF, &ncplx, NLAT_2, 1, fftw_plan_mode);
		if (fft_eo == NULL) shtns_runerr("[FFTW] fft planning failed !");

#if SHT_VERBOSE > 2
	printf(" *** fft plan :\n");
	fftw_print_plan(fft);
	printf("\n *** ifft plan :\n");
	fftw_print_plan(ifft);
	printf("\n");
#endif

//	fft_norm = 1.0/nfft;
	if (SHT_FFT > 1) fftw_free(Sh);
	fftw_free(ShF);
  } else {
	if (theta_inc != 1) shtns_runerr("only contiguous spatial data is supported for Nphi=1");
#if SHT_VERBOSE > 0
	printf("        => no fft : Mmax=0, Nphi=1, Nlat=%d\n",NLAT);
#endif
	SHT_FFT = 0;	// no fft.
  }
	dct_m0 = NULL;	idct = NULL;		// set dct plans to uninitialized.
	dct_r1 = NULL;	idct_r1 = NULL;
	return(phi_embed * NLAT);
}

#ifndef SHT_NO_DCT
/// initialize DCTs using FFTW. Must be called if MTR_DCT is changed.
void planDCT()
{
	double *Sh;
	int ndct = NLAT;
	fftw_r2r_kind r2r_kind;
	fftw_iodim dims, hdims[2];
	double Sh0[NLAT] SSE;				// temp storage on the stack, aligned.
	
// real NPHI=1, allocate only once since it does not change.
	if ((dct_r1 == NULL)||(idct_r1 == NULL)) {
		Sh = (double *) fftw_malloc( NLAT * sizeof(double) );
		if (dct_r1 == NULL) {
			r2r_kind = FFTW_REDFT10;
			dct_r1 = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, 1, NLAT, Sh, &ndct, 1, NLAT, &r2r_kind, fftw_plan_mode);
		}
		if (idct_r1 == NULL) {
			r2r_kind = FFTW_REDFT01;
			idct_r1 = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, 1, NLAT, Sh, &ndct, 1, NLAT, &r2r_kind, fftw_plan_mode);
		}
		fftw_free(Sh);
		if ((dct_r1 == NULL)||(idct_r1 == NULL))
			shtns_runerr("[FFTW] (i)dct_r1 planning failed !");
#if SHT_VERBOSE > 2
			printf(" *** idct_r1 plan :\n");	fftw_print_plan(idct_r1);
			printf("\n *** dct_r1 plan :\n");	fftw_print_plan(dct_r1);	printf("\n");
#endif
	}

#ifndef SHT_AXISYM
	if (idct != NULL) fftw_destroy_plan(idct);
	// Allocate dummy Spatial Fields.
	Sh = (double *) fftw_malloc((NPHI/2 +1) * NLAT*2 * sizeof(double));

	dims.n = NLAT;	dims.is = 2;	dims.os = 2;		// real and imaginary part.
	hdims[0].n = MTR_DCT+1;	hdims[0].is = 2*NLAT; 	hdims[0].os = 2*NLAT;
	hdims[1].n = 2;			hdims[1].is = 1; 	hdims[1].os = 1;

	if (NPHI>1) {		// complex data for NPHI>1, recompute as it does depend on MTR_DCT
		r2r_kind = FFTW_REDFT01;
		idct = fftw_plan_guru_r2r(1, &dims, 2, hdims, Sh, Sh, &r2r_kind, fftw_plan_mode);
		if (idct == NULL)
			shtns_runerr("[FFTW] idct planning failed !");
#if SHT_VERBOSE > 2
			printf(" *** idct plan :\n");	fftw_print_plan(idct);	printf("\n");
#endif
		if (dct_m0 == NULL) {
			r2r_kind = FFTW_REDFT10;
//			dct_m0 = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, 2, 2*NLAT, Sh, &ndct, 2, 2*NLAT, &r2r_kind, fftw_plan_mode);
			dct_m0 = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, 2, 2*NLAT, Sh0, &ndct, 1, NLAT, &r2r_kind, fftw_plan_mode);	// out-of-place.
			if (dct_m0 == NULL)
				shtns_runerr("[FFTW] dct_m0 planning failed !");
#if SHT_VERBOSE > 2
				printf(" *** dct_m0 plan :\n");		fftw_print_plan(dct_m0);	printf("\n");
#endif
		}
	} else {	// NPHI==1
		if (dct_m0 == NULL) {
			r2r_kind = FFTW_REDFT10;
			dct_m0 = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, 1, NLAT, Sh0, &ndct, 1, NLAT, &r2r_kind, fftw_plan_mode);	// out-of-place.
			if (dct_m0 == NULL)
				shtns_runerr("[FFTW] dct_m0 planning failed !");
#if SHT_VERBOSE > 2
				printf(" *** dct_m0 plan :\n");		fftw_print_plan(dct_m0);	printf("\n");
#endif
		}
	}
	fftw_free(Sh);
#endif
}
#endif

/// SET MTR_DCT and updates fftw_plan for DCT's
void Set_MTR_DCT(int m)
{
#ifndef SHT_NO_DCT
	if ((zlm_dct0 == NULL)||(m == MTR_DCT)) return;
	if ( m < 0 ) {	// don't use dct
		MTR_DCT = -1;
	} else {
		if (m>MMAX) m=MMAX;
		MTR_DCT = m;
		planDCT();
	}
#endif
}

int Get_MTR_DCT() {
	return MTR_DCT;
}

/// TIMINGS
	double get_time(int m, int nloop, double *Sh, double *Th, complex double *Slm, complex double *Tlm)
	{
		int i;
		ticks tik0, tik1;

		Set_MTR_DCT(m);
			SH_to_spat(Tlm,Sh);			// caching...
			SHsphtor_to_spat(Slm,Tlm,Sh,Th);
		tik0 = getticks();
		for (i=0; i<nloop; i++) {
			SH_to_spat(Tlm,Sh);
			SHsphtor_to_spat(Slm,Tlm,Sh,Th);
		}
		tik1 = getticks();
	#if SHT_VERBOSE > 1
		printf("m=%d - ticks : %.3f\n", m*MRES, elapsed(tik1,tik0)/(nloop*NLM*NLAT));
	#endif
		return elapsed(tik1,tik0);
	}

double Find_Optimal_SHT()
{
	complex double *Slm, *Tlm;
	double *Sh, *Th;
	int m, i, minc, nloop;
	double t, tsum, tsum0;

	Sh = (double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Th = (double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Slm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Tlm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);

	t = 1.0 / (RAND_MAX/2);		// some random data to play with.
	for (i=0;i<NLM;i++) {
		Slm[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
		Tlm[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
	}

        minc = MMAX/20 + 1;             // don't test every single m.
        nloop = 10;                     // number of loops to get time.
        if (NLM*NLAT > 1024*1024)
                nloop = 1 + (10*1024*1024)/(NLM*NLAT);          // don't use too many loops for large transforms.
#if SHT_VERBOSE > 1
	printf("\nminc = %d, nloop = %d\n",minc,nloop);
#endif

	m = -1;
	tsum0 = get_time(m, nloop, Sh, Th, Slm, Tlm);
	tsum=tsum0;	i = m;
	t = get_time(m, nloop, Sh, Th, Slm, Tlm);
	if (t<tsum0) tsum0 = t;		// recheck m=-1
	for (m=0; m<=MMAX; m+=minc) {
		t = get_time(m, nloop, Sh, Th, Slm, Tlm);
		if ((m==-1)&&(t<tsum0)) tsum0 = t;
		if (t < tsum) {	tsum = t;	i = m; }
	}
	m=-1;	// recheck m=-1;
		t = get_time(m, nloop, Sh, Th, Slm, Tlm);
		if (t<tsum0) tsum0 = t;
		if (t < tsum) { tsum = t;	i = m; }
		t = get_time(m, nloop, Sh, Th, Slm, Tlm);	// twice
		if (t<tsum0) tsum0 = t;
		if (t < tsum) { tsum = t;	i = m; }

	free(Tlm);	free(Slm);	fftw_free(Th);		fftw_free(Sh);

	Set_MTR_DCT(i);
	return(tsum/tsum0);	// returns ratio of "optimal" time over "no_dct" time
}


/// Perform some optimization on the SHT matrices.
void OptimizeMatrices(double eps)
{
	double *yg;
	int im,m,l,it;

/// POLAR OPTIMIZATION : analyzing coefficients, some can be safely neglected.
	if (eps > 0.0) {
		tm[0] = 0;
		for (im=1;im<=MMAX;im++) {
			m = im*MRES;
			tm[im] = NLAT_2;
			for (l=m;l<=LMAX;l++) {
				it=0;
				while( fabs(ylm[im][it*(LMAX-m+1) + (l-m)]) < eps ) { it++; }
				if (tm[im] > it) tm[im] = it;
			}
		}
#if SHT_VERBOSE > 0
		printf("        + polar optimization threshold = %.1e\n",eps);
#endif
#if SHT_VERBOSE > 1
		printf("          tm[im]=");
		for (im=0;im<=MMAX;im++)
			printf(" %d",tm[im]);
		printf("\n");
#endif
		for (im=1; im<=MMAX; im++) {	//	im >= 1
		  if (tm[im] > 0) {		// we can remove the data corresponding to polar values.
			m = im*MRES;
			ylm[im]  += tm[im]*(LMAX-m+1);		// shift pointers (still one block for each m)
#ifndef SHT_SCALAR_ONLY
			dylm[im] += tm[im]*(LMAX-m+1);
#endif
			if (zlm[0] != NULL) {
			  for (l=m; l<LMAX; l+=2) {
				for (it=0; it<NLAT_2-tm[im]; it++) {	// copy data to avoid cache misses.
					zlm[im][(l-m)*(NLAT_2-tm[im]) + it*2]   = zlm[im][(l-m)*NLAT_2 + (it+tm[im])*2];
					zlm[im][(l-m)*(NLAT_2-tm[im]) + it*2+1] = zlm[im][(l-m)*NLAT_2 + (it+tm[im])*2+1];
#ifndef SHT_SCALAR_ONLY
					dzlm[im][(l-m)*(NLAT_2-tm[im]) + it*2].t = dzlm[im][(l-m)*NLAT_2 + (it+tm[im])*2].t;
					dzlm[im][(l-m)*(NLAT_2-tm[im]) + it*2].p = dzlm[im][(l-m)*NLAT_2 + (it+tm[im])*2].p;
					dzlm[im][(l-m)*(NLAT_2-tm[im]) + it*2+1].t = dzlm[im][(l-m)*NLAT_2 + (it+tm[im])*2+1].t;
					dzlm[im][(l-m)*(NLAT_2-tm[im]) + it*2+1].p = dzlm[im][(l-m)*NLAT_2 + (it+tm[im])*2+1].p;
#endif
				}
			  }
			  if (l==LMAX) {
				for (it=0; it<NLAT_2-tm[im]; it++) {
					zlm[im][(l-m)*(NLAT_2-tm[im]) + it]   = zlm[im][(l-m)*NLAT_2 + (it+tm[im])];
#ifndef SHT_SCALAR_ONLY
					dzlm[im][(l-m)*(NLAT_2-tm[im]) + it].t = dzlm[im][(l-m)*NLAT_2 + (it+tm[im])].t;
					dzlm[im][(l-m)*(NLAT_2-tm[im]) + it].p = dzlm[im][(l-m)*NLAT_2 + (it+tm[im])].p;
#endif
				}
			  }
			}
		  }
		}
	} else {
		for (im=0;im<=MMAX;im++)	tm[im] = 0;
	}

#ifndef SHT_SCALAR_ONLY
/// Compression of dylm and dzlm for m=0, as .p is 0
	im=0;	m=0;
	yg = (double *) dylm[im];
	long int lstride = LMAX + (LMAX & 1);	// we need an even stride here, for SSE2 alignement.
	for (it=0; it<NLAT_2; it++) {
		for (l=1; l<=LMAX; l++)		// start at l=1, as l=0 is null.
			yg[it*lstride + (l-1)] = dylm[im][it*(LMAX+1) + l].t;
		if (LMAX & 1) yg[it*lstride + LMAX] = 0.0;	// add some zero here for vectorization.
	}
	if (dzlm[0] != NULL) {		// for sht_reg_poles there is no dzlm defined.
		yg = (double *) dzlm[im];
		for (l=1; l<LMAX; l+=2) {		// l=0 is zero, so we start at l=1.
			for (it=0; it<NLAT_2; it++) {
				yg[(l-1)*NLAT_2 + it*2] = dzlm[im][(l-1)*NLAT_2 + it*2].t;	// l
				yg[(l-1)*NLAT_2 + it*2+1] = dzlm[im][(l-1)*NLAT_2 + it*2+1].t;	// l+1
			}
		}
		if (l==LMAX) {		// last l is stored right away, without interleaving.
			for (it=0; it<NLAT_2; it++) {
				yg[(l-1)*NLAT_2 + it] = dzlm[im][(l-1)*NLAT_2 + it].t;		// l (odd)
			}
		}
	}
#endif
}


/// Precompute the matrix for SH synthesis.
void init_SH_synth()
{
	double dtylm[LMAX+1];		// temp storage for derivative : d(P_l^m(x))/dx
	long int it,im,m,l;
	long int lstride;

	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
		lstride = LMAX-m+1;
		if (im==0) lstride += (lstride&1);		// even stride for m=0.
		for (it=0; it<NLAT_2; it++) {
			legendre_sphPlm_deriv_array(LMAX, im, ct[it], st[it], ylm[im] + it*lstride, dtylm);	// fixed im legendre functions lookup table.
			for (l=m; l<=LMAX; l++) {
  #ifndef SHT_SCALAR_ONLY
				dylm[im][it*(LMAX-m+1) + (l-m)].t = dtylm[l-m];
				dylm[im][it*(LMAX-m+1) + (l-m)].p = ylm[im][it*lstride + (l-m)] *m;	// 1/sint(t) dYlm/dphi
  #endif
				if (m>0) ylm[im][it*lstride + (l-m)] *= st[it];
			}
			if ((im==0) && ((LMAX&1) == 0)) ylm[im][it*lstride + LMAX+1] = 0.0;	// add one zero for padding.
		}
	}
}


/// Precompute matrices for SH synthesis and analysis, on a Gauss-Legendre grid.
void init_SH_gauss()
{
	double t,tmax;
	long int it,im,m,l;
	long double iylm_fft_norm;
	long double xg[NLAT], wgl[NLAT];	// gauss points and weights.
	
	wg = malloc(NLAT_2 * sizeof(double));	// gauss weights, double precision.

 	if ((SHT_NORM == sht_fourpi)||(SHT_NORM == sht_schmidt)) {
 		iylm_fft_norm = 0.5/NPHI;		// FFT/SHT normalization for zlm (4pi normalized)
	} else {
 		iylm_fft_norm = 2.0*M_PIl/NPHI;		// FFT/SHT normalization for zlm (orthonormalized)
	}
#if SHT_VERBOSE > 0
	printf("        => using Gauss nodes\n");
	if (2*NLAT <= (SHT_NL_ORDER +1)*LMAX) printf("     !! Warning : Gauss-Legendre anti-aliasing condition 2*Nlat > %d*Lmax is not met.\n",SHT_NL_ORDER+1);
#endif
	gauss_nodes(xg,wgl,NLAT);	// generate gauss nodes and weights : ct = ]1,-1[ = cos(theta)
	for (it=0; it<NLAT; it++) {
		ct[it] = xg[it];
		st[it] = sqrtl((1.-xg[it])*(1.+xg[it]));
		st_1[it] = 1.0/sqrtl((1.-xg[it])*(1.+xg[it]));
	}
	for (it=0; it<NLAT_2; it++)
		wg[it] = wgl[it]*iylm_fft_norm;		// faster double-precision computations.

	if (NLAT & 1) {		// odd NLAT : adjust weigth of middle point.
		wg[NLAT_2-1] *= 0.5;
	}

#if SHT_VERBOSE > 1
	printf(" NLAT=%d, NLAT_2=%d\n",NLAT,NLAT_2);
// TEST if gauss points are ok.
	tmax = 0.0;
	for (it = 0; it<NLAT_2; it++) {
		t = legendre_Pl(NLAT, ct[it]);
		if (t>tmax) tmax = t;
//		printf("i=%d, x=%12.12g, p=%12.12g\n",it,ct[it],t);
	}
	printf("          max zero at Gauss nodes for Pl[l=NLAT] : %g\n",tmax);
	if (NLAT_2 < 100) {
		printf("          Gauss nodes :");
		for (it=0;it<NLAT_2; it++)
			printf(" %g",ct[it]);
		printf("\n");
	}
#endif

	init_SH_synth();

// for analysis (decomposition, direct transform) : transpose and multiply by gauss weight and other normalizations.
// interleave l and l+1 : this stores data in the way it will be read.
	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
		long int talign = 0;
		long int lstride = LMAX-m+1;
		if (im==0) lstride += (lstride&1);		// even stride for m=0.

		for (it=0;it<NLAT_2;it++) {
			double nz0, nz1, norm;
			norm = wg[it];
			if ( (m>0) && (shtns.norm & SHT_REAL_NORM) )	norm *= 2;		// "Real" norm : zlm must be doubled for m>0
			nz0 = norm;		nz1 = norm;
			long int l0 = m;
			if (m==0) {
				zlm[im][it] = ylm[im][it*lstride] * norm;
				// les derivees sont nulles pour l=0
				l0++;
				talign = (NLAT_2&1);
			}
			for (l=l0; l<LMAX; l+=2) {
				if (SHT_NORM == sht_schmidt) {
					nz0 = norm*(2*l+1);		nz1 = norm*(2*l+3);
				}
				zlm[im][(l-m)*NLAT_2 + it*2 +talign]    =  ylm[im][it*lstride + (l-m)]   * nz0;
				zlm[im][(l-m)*NLAT_2 + it*2 +1 +talign] =  ylm[im][it*lstride + (l+1-m)] * nz1;
  #ifndef SHT_SCALAR_ONLY
				dzlm[im][(l-l0)*NLAT_2 + it*2].t = dylm[im][it*(LMAX-m+1) + (l-m)].t * nz0 /(l*(l+1));
				dzlm[im][(l-l0)*NLAT_2 + it*2].p = dylm[im][it*(LMAX-m+1) + (l-m)].p * nz0 /(l*(l+1));
				dzlm[im][(l-l0)*NLAT_2 + it*2+1].t = dylm[im][it*(LMAX-m+1) + (l+1-m)].t * nz1 /((l+1)*(l+2));
				dzlm[im][(l-l0)*NLAT_2 + it*2+1].p = dylm[im][it*(LMAX-m+1) + (l+1-m)].p * nz1 /((l+1)*(l+2));
  #endif
			}
			if (l==LMAX) {		// last l is stored right away, without interleaving.
				if (SHT_NORM == sht_schmidt)
					nz0 = norm*(2*l+1);
				zlm[im][(l-m)*NLAT_2 + it +talign]    =  ylm[im][it*lstride + (l-m)]   * nz0;
  #ifndef SHT_SCALAR_ONLY
				dzlm[im][(l-l0)*NLAT_2 + it].t = dylm[im][it*(LMAX-m+1) + (l-m)].t * nz0 /(l*(l+1));
				dzlm[im][(l-l0)*NLAT_2 + it].p = dylm[im][it*(LMAX-m+1) + (l-m)].p * nz0 /(l*(l+1));
  #endif
			}
		}
	}
}

/*
eval_sht(double *val, complex double *ql, double *al, double *bl, long int lmax, double ct)
{
	complex double t1, t2;
	long int l;

// basic "naive"
	t2 = ql[lmax];								// al[l] = (2l-1)/(l-m);	bl[l] = -(l-1+m)/(l-m);
	t1 = ql[lmax-1] + (al[lmax]*ct)*t2;
	for (l=lmax-2; l>m; l-=2) {		// ops: 10*, 8+ (double)
		t2 = ql[l]   + (al[l+1]*ct)*t1 + bl[l+2]*t2;		// 5*, 4+
		t1 = ql[l-1] + (al[l]*ct)*t2   + bl[l+1]*t1;		// 5*, 4+
	}
	if (l==m)
		t1 = ql[l] + (al[l+1]*ct)*t1 + bl[l+2]*t2;		// 5*, 4+
	*val = t1 * y_mm(ct);		// y_mm   = sqrt((2.0+1.0/m)/(4.0*M_PI)) * sgn * exp(lnpre); (cf gsl)

// advanced "even/odd" (saves some op, and reduces loop dependancy)
	te2 = 0.0;
	to2 = ql[lmax];
	te1 = ql[lmax-1];
	to1 = al[lmax]*ct*to2;
	for (l=lmax-2; l>m; l-=2) {		// ops: 18*, 12+  (instead of 20*, 16+, + cache miss)
		te2 =           (al[l+1]*ct)*te1 + bl[l+2]*te2;			// 4*, 2+
		to2 = ql[l]   + (al[l+1]*ct)*to1 + bl[l+2]*to2;			// 5*, 4+
		te1 = ql[l-1] + (al[l]*ct)*te2   + bl[l+1]*te1;			// 4*, 4+
		to1 =           (al[l]*ct)*to2   + bl[l+1]*to1;			// 5*, 2+
	}		// compare to 4*, 4+ of the hyb_SH_to_spat...
	if (l==m) {		// only 1 left, + exchange even and odd.
		te2 =         al[l+1]*ct*te1 + bl[l+2]*te2;		// 5*, 2+
		to2 = ql[l] + al[l+1]*ct*to1 + bl[l+2]*to2;		// 5*, 4+
		te1 = to2;
		to1 = te2;
	}
	*vale = te1 * y_mm(ct);		// y_mm   = sqrt((2.0+1.0/m)/(4.0*M_PI)) * sgn * exp(lnpre); (cf gsl)
	*valo = to1 * y_mm(ct);		// y_mm   = sqrt((2.0+1.0/m)/(4.0*M_PI)) * sgn * exp(lnpre); (cf gsl)
	    // P_{\ell}^{\ell}(x) = (-1)^l (2\ell-1)!! (1- x^2)^{(l/2)}
}

/// Clenshaw algorithm to evaluate a partial Chebyshev series dct[] up to degree n at any ct = cos(theta)
/// works for n>=2. see http://en.wikipedia.org/wiki/Clenshaw_algorithm
inline eval_dct(double *val, double *dct, long int n, double ct)
{
	double t1, t2, tmp;
	double ct2 = 2.*ct;
	long int i;

	t2 = dct[n];
	t1 = dct[n-1] + ct2*t2;
	for (i=n-2; i>0; i--) {
		tmp = t1;
		t1 = dct[i] + ct2*t1 - t2;		// 1*, 2+
		t2 = tmp;
	}
	*val = dct[0] + ct*t1 - t2;
}

inline eval_dct_cplx(complex double *val, complex double *dct, long int n, double ct)
{
	complex double t1, t2, tmp;
	double ct2 = 2.*ct;
	long int i;

	t2 = dct[n];
	t1 = dct[n-1] + ct2*t2;
	for (i=n-2; i>0; i--) {
		tmp = t1;
		t1 = dct[i] + ct2*t1 - t2;		// 2*, 4+
		t2 = tmp;
	}
	*val = dct[0] + ct*t1 - t2;

// unrolled version
	t2 = dct[n];
	t1 = dct[n-1] + ct2*t2;
	for (i=n-2; i>1; i-=2) {
		t2 = dct[i]   + ct2*t1 - t2;		// 2*, 4+
		t1 = dct[i-1] + ct2*t2 - t1;		// 2*, 4+
	}
	if (i != 0) {
		tmp = t1;
		t1 = dct[i] + ct2*t1 - t2;		// 2*, 4+
		t2 = tmp;
	}
	*val = dct[0] + ct*t1 - t2;

// even/odd acceleration (saves 2+, reduces loop dependancy)
	t2e = 0.0;
	t2o = dct[n];
	t1e = dct[n-1];
	t1o = ct2*t2o;
	for (i=n-2; i>1; i-=2) {	// total : 8*, 12+
		t2e =            ct2*t1e - t2e;		// 2*, 2+
		t2o = dct[i]   + ct2*t1o - t2o;		// 2*, 4+
		t1e = dct[i-1] + ct2*t2e - t1e;		// 2*, 4+
		t1o =            ct2*t2o - t1o;		// 2*, 2+
	}
	if (i != 0) {	// i == 1
		t2e =          ct2*t1e - t2e;		// 2*, 2+
		t2o = dct[i] + ct2*t1o - t2o;		// 2*, 4+
		t1e = dct[0] + ct*t2e - t1e;		// 2*, 4+
		t1o =          ct*t2o - t1o;		// 2*, 2+
	} else {	// i == 0: exchange even and odd.
		t2e =          ct*t1e - t2e;
		t2o = dct[0] + ct*t1o - t2o;
		t1e = t2o;
		t1o = t2e;
	}
	*vale = t1e;
	*valo = t1o;
}
*/

void free_SH_dct()
{
#ifndef SHT_SCALAR_ONLY
	fftw_free(dzlm_dct0);
#endif
	fftw_free(zlm_dct0);
#ifndef SHT_SCALAR_ONLY
	fftw_free(dykm_dct[0]);
#endif
	fftw_free(ykm_dct[0]);
}

#ifndef SHT_NO_DCT
/// Computes the matrices required for SH transform on a regular grid (with or without DCT).
/// \param analysis : 0 => synthesis only.
void init_SH_dct(int analysis)
{
	fftw_plan dct, idct;
	double *yk, *yk0, *dyk0, *yg;		// temp storage
	struct DtDp *dyg, *dyk;
	double dtylm[LMAX+1];		// temp storage for derivative : d(P_l^m(x))/dx
	double iylm_fft_norm;
	long int it,im,m,l;
	long int sk, dsk;
	double Z[2*NLAT_2], dZt[2*NLAT_2], dZp[2*NLAT_2];		// equally spaced theta points.
	double is1[NLAT];		// tabulate values for integrals.

	if ((SHT_NORM == sht_fourpi)||(SHT_NORM == sht_schmidt)) {
 		iylm_fft_norm = 0.5/(NPHI*NLAT_2);	// FFT/DCT/SHT normalization for zlm (4pi)
	} else {
 		iylm_fft_norm = 2.0*M_PI/(NPHI*NLAT_2);	// FFT/DCT/SHT normalization for zlm (orthonormal)
	}

#if SHT_VERBOSE > 0
	printf("        => using equaly spaced nodes with DCT acceleration\n");
	if (NLAT <= SHT_NL_ORDER *LMAX)	printf("     !! Warning : DCT anti-aliasing condition Nlat > %d*Lmax is not met.\n",SHT_NL_ORDER);
	if (NLAT != fft_int(NLAT,7))	printf("     !! Warning : Nlat is not optimal for FFTW !\n");
#endif
	if ((NLAT_2)*2 <= LMAX+1) shtns_runerr("NLAT_2*2 should be at least LMAX+2 (DCT)");
	if (NLAT & 1) shtns_runerr("NLAT must be even (DCT)");
	for (it=0; it<NLAT; it++) {	// Chebychev points : equaly spaced but skipping poles.
		long double th = M_PIl*(it+0.5)/NLAT;
		ct[it] = cosl(th);
		st[it] = sinl(th);
		st_1[it] = 1.0/sinl(th);
	}
#if SHT_VERBOSE > 1
	{
	double tsum, t;
	printf(" NLAT=%d, NLAT_2=%d\n",NLAT,NLAT_2);
	if (NLAT_2 < 100) {
		printf("          DCT nodes :");
		tsum = 0.0;
		for (it=0;it<NLAT_2; it++) {
			printf(" %g",ct[it]);
			t = fabs(ct[it]*ct[it] + st[it]*st[it] -1.0);
			if (t > tsum) tsum=t;
		}
	}
	printf("\n max st^2 + ct^2 -1 = %lg\n",tsum);
	}
#endif

#define KMAX (LMAX+1)

	for(im=0, sk=0, dsk=0; im<=MMAX; im++) {	// how much memory to allocate for ykm_dct ?
		m = im*MRES;
		for (it=0; it<= KMAX; it+=2) {
			l = (it < m) ? m : it-(m&1);
			sk += LMAX+1 - l;
			if ((m==0) && ((LMAX & 1) ==0)) sk++;		// SSE padding for m=0
		}
		for (it=0; it<= KMAX; it+=2) {
			l = (it-2 < m) ? m : it-2+(m&1);
			dsk += LMAX+1 - l;
		}
	}
	for (l=0, it=0; l<=LMAX; l+=2)	// how much memory for zlm_dct0 ?
		it += (2*NLAT_2 -l);
	for (l=1, im=0; l<=LMAX; l+=2)	// how much memory for dzlm_dct0 ?
		im += (2*NLAT_2 -l+1);

#if SHT_VERBOSE > 1
	printf("          Memory used for Ykm_dct matrices = %.3f Mb\n",sizeof(double)*(sk + 2.*dsk + it)/(1024.*1024.));
#endif
	ykm_dct[0] = (double *) fftw_malloc(sizeof(double)* sk);
#ifndef SHT_SCALAR_ONLY
	dykm_dct[0] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* dsk);
#endif
	zlm_dct0 = (double *) fftw_malloc( sizeof(double)* it );
#ifndef SHT_SCALAR_ONLY
	dzlm_dct0 = (double *) fftw_malloc( sizeof(double)* im );
#endif
	for (im=0; im<MMAX; im++) {
		m = im*MRES;
		for (it=0, sk=0; it<= KMAX; it+=2) {
			l = (it < m) ? m : it-(m&1);
			sk += LMAX+1 - l;
			if ((m==0) && ((LMAX & 1) ==0)) sk++;		// SSE padding for m=0
		}
		for (it=0, dsk=0; it<= KMAX; it+=2) {
			l = (it-2 < m) ? m : it-2+(m&1);
			dsk += LMAX+1 - l;
		}
		ykm_dct[im+1] = ykm_dct[im] + sk;
#ifndef SHT_SCALAR_ONLY
		dykm_dct[im+1] = dykm_dct[im] + dsk;
#endif
	}

	dct = fftw_plan_r2r_1d( 2*NLAT_2, Z, Z, FFTW_REDFT10, FFTW_ESTIMATE );	// quick and dirty dct.
	idct = fftw_plan_r2r_1d( 2*NLAT_2, Z, Z, FFTW_REDFT01, FFTW_ESTIMATE );	// quick and dirty idct.

#if SHT_VERBOSE > 1
	ticks tik0, tik1;
	tik0 = getticks();
#endif

// precomputation for scalar product of Chebychev polynomials.
	for(it=0; it<NLAT; it++)
		is1[it] = 1./(1. - 4.*it*it);

// Even/Odd symmetry : ylm is even or odd across equator, as l-m is even or odd => only NLAT_2 points required.
	// temp memory for ykm_dct.
	yk = (double *) malloc( sizeof(double) * (KMAX+1)*(LMAX+1) );
	dyk = (struct DtDp *) malloc( sizeof(struct DtDp)* (KMAX+1)*(LMAX+1) );
	if (analysis) {
		yk0 = (double *) malloc( sizeof(double) * (LMAX/2+1)*(2*NLAT_2) * 2 );		// temp for zlm_dct0
		dyk0 = yk0 + (LMAX/2+1)*(2*NLAT_2);
	}

	init_SH_synth();

	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
		long int lstride = LMAX-m+1;
		if (im==0) lstride += (lstride&1);		// even stride for m=0.

	// go to DCT space
		for (it=0;it<=KMAX;it+=2) {
			for(l=m; l<=LMAX; l++) {
				yk[(it/2)*(LMAX+1-m) + (l-m)] = 0.0;
				dyk[(it/2)*(LMAX+1-m) + (l-m)].t = 0.0;
				dyk[(it/2)*(LMAX+1-m) + (l-m)].p = 0.0;
			}
		}
		for (l=m; l<=LMAX; l++) {
			if (m & 1) {	// m odd
				for (it=0; it<NLAT_2; it++) {
					Z[it] = ylm[im][it*lstride + (l-m)] * st[it];	// P[l+1](x)	*st
#ifndef SHT_SCALAR_ONLY
					dZt[it] = dylm[im][it*(LMAX-m+1) + (l-m)].t;	// P[l](x)	*1
					dZp[it] = dylm[im][it*(LMAX-m+1) + (l-m)].p;		// P[l-1](x)	*1
#endif
				}
			} else {	// m even
				for (it=0; it<NLAT_2; it++) {
					Z[it] = ylm[im][it*lstride + (l-m)];		// P[l](x)	*1
#ifndef SHT_SCALAR_ONLY
					dZt[it] = dylm[im][it*(LMAX-m+1) + (l-m)].t *st[it];	// P[l+1](x)	*st
					dZp[it] = ylm[im][it*lstride + (l-m)] * m;	// P[l](x)	*st
#endif
				}
			}
			if ((l-m)&1) {	// odd
				for (it=NLAT_2; it<2*NLAT_2; it++) {
					Z[it]   = - Z[2*NLAT_2-it-1];	// reconstruct even/odd part
					dZt[it] =   dZt[2*NLAT_2-it-1];
					dZp[it] = - dZp[2*NLAT_2-it-1];
				}
			} else {	// even
				for (it=NLAT_2; it<2*NLAT_2; it++) {
					Z[it] =     Z[2*NLAT_2-it-1];	// reconstruct even/odd part
					dZt[it] = - dZt[2*NLAT_2-it-1];
					dZp[it] =   dZp[2*NLAT_2-it-1];
				}
			}
			fftw_execute(dct);
			fftw_execute_r2r(dct, dZt, dZt);
			fftw_execute_r2r(dct, dZp, dZp);
#if SHT_VERBOSE > 1
			if (LMAX <= 12) {
				printf("\nl=%d, m=%d ::\t", l,m);
				for(it=0;it<2*NLAT_2;it++) printf("%e ",Z[it]/(2*NLAT));
				printf("\n     dYt ::\t");
				for(it=0;it<2*NLAT_2;it++) printf("%e ",dZt[it]/(2*NLAT));
				printf("\n     dYp ::\t");
				for(it=0;it<2*NLAT_2;it++) printf("%e ",dZp[it]/(2*NLAT));
			}
#endif
			for (it=(l-m)&1; it<=l+1; it+=2) {
				yk[(it/2)*(LMAX+1-m) + (l-m)] = Z[it]/(2*NLAT);	// and transpose
				dyk[(it/2)*(LMAX+1-m) + (l-m)].p = dZp[it]/(2*NLAT);
			}
			for (it=(l+1-m)&1; it<=l+1; it+=2) {
				dyk[(it/2)*(LMAX+1-m) + (l-m)].t = dZt[it]/(2*NLAT);
			}
		}

	/* compute analysis coefficients (fast way)
	 * Wklm = int(Tk*Ylm) = int(Tk.sum(i,a_ilm*Ti)) = sum(i, a_ilm* int(Tk*Ti)) = sum(i, a_ilm*Jik)
	 * with Jik = int(Tk*Ti) = 1/(1-(k-i)^2) + 1/(1-(k+i)^2)
	*/
		if (analysis) {
	#if SHT_VERBOSE > 0
		printf("computing weights m=%d\r",m);	fflush(stdout);
	#endif
		for (l=m; l<=LMAX; l++) {
			long int k0,k1, k,i,d;
			double Jik, yy, dyy;
			double lnorm = iylm_fft_norm;

			if (SHT_NORM == sht_schmidt)	lnorm *= (2*l+1);		// Schmidt semi-normalization
			if ( (m>0) && (shtns.norm & SHT_REAL_NORM) )	lnorm *= 2;		// "real" norm : zlm must be doubled for m>0

			k0 = (l-m)&1;	k1 = 1-k0;
			for(k=0; k<NLAT; k++) {	Z[k] = 0.0;		dZt[k] = 0.0;	dZp[k] = 0.0; }
			for (i=k0; i<=l+1; i+=2) {		// i+k even
				yy = yk[(i/2)*(LMAX+1-m) + (l-m)] * lnorm;
				dyy = dyk[(i/2)*(LMAX+1-m) + (l-m)].p * lnorm/(l*(l+1));
				if (i==0) {	yy*=0.5;	dyy*=0.5; }
				for (k=k0; k<NLAT; k+=2) {
					d = (k<i) ? i-k : k-i;		// d=|i-k|
					Jik = is1[(i+k)/2] + is1[d/2];
					Z[k] += yy * Jik;
					if (m&1) dZp[k] += dyy * Jik;
				}
			}
			if (l != 0) {
				for (i=k1; i<=l+1; i+=2) {		// i+k even
					yy = dyk[(i/2)*(LMAX+1-m) + (l-m)].t * lnorm/(l*(l+1));
					if (i==0) yy*=0.5;
					for (k=k1; k<NLAT; k+=2) {
						d = (k<i) ? i-k : k-i;		// d=|i-k|
						Jik = is1[(i+k)/2] + is1[d/2];
						dZt[k] += yy * Jik;
					}
				}
			}
#if SHT_VERBOSE > 1
		if (LMAX <= 12) {
			printf("\nl=%d, m=%d ::\t",l,m);
			for (k=0; k<(2*NLAT_2); k++) printf("%f ",Z[k]);
			printf("\n       dZt ::\t");
			for (k=0; k<(2*NLAT_2); k++) printf("%f ",dZt[k]);
			if (m&1) {
				printf("\n       dZp ::\t");
				for (k=0; k<(2*NLAT_2); k++) printf("%f ",dZp[k]);
			}
		}
#endif

			if (m == 0) {		// we store zlm in dct space for m=0
				if (k0==0) 	{
					yk0[((l-m)>>1)*(2*NLAT_2)] = Z[0]*0.5;         // store zlm_dct (k=0)
					for (k=1; k<(2*NLAT_2); k++) yk0[((l-m)>>1)*(2*NLAT_2) +k] = 0.0;		// zero out.
					k0=2;
				}
				for (k=k0; k<(2*NLAT_2); k+=2)
					yk0[((l-m)>>1)*(2*NLAT_2) +k] = Z[k];             // store zlm_dct
					
				if (l>0) {
					if (k1==0) 	{
						dyk0[((l-1-m)>>1)*(2*NLAT_2)] = dZt[0]*0.5;         // store dzlm_dct (k=0)
						for (k=1; k<(2*NLAT_2); k++) dyk0[((l-1-m)>>1)*(2*NLAT_2) +k] = 0.0;		// zero out.
						k1=2;
					}
					for (k=k1; k<(2*NLAT_2); k+=2)
						dyk0[((l-1-m)>>1)*(2*NLAT_2) +k] = dZt[k];             // store dzlm_dct
				}
			}

			fftw_execute_r2r(idct, Z, Z);	fftw_execute_r2r(idct, dZt, dZt);
			if (m == 0) {
				for (it=0; it<NLAT; it++) { dZp[it] = 0.0; 	dZt[it] *= st_1[it]; }
			} else if (m & 1) {	// m odd
				fftw_execute_r2r(idct, dZp, dZp);
				for (it=0; it<NLAT; it++) {	Z[it] *= st_1[it]; }
			} else {	// m even
				for (it=0; it<NLAT; it++) { dZp[it] = Z[it]*m/(l*(l+1)*st[it]); 	dZt[it] *= st_1[it]; }
			}

			long int l0 = (m==0) ? 1 : m;
			long int talign = (m==0)*(NLAT_2 & 1);
			sk = (l-l0)&1;
			if (l==0) {
				for (it=0; it<NLAT_2; it++) {
					zlm[im][it] =  Z[it];
				}
			} else if ((sk == 0)&&(l == LMAX)) {
				for (it=0; it<NLAT_2; it++) {
					zlm[im][(l-m)*NLAT_2 + it + talign] =  Z[it];
	#ifndef SHT_SCALAR_ONLY
					dzlm[im][(l-l0)*NLAT_2 + it].p = dZp[it];
					dzlm[im][(l-l0)*NLAT_2 + it].t = dZt[it];
	#endif
				}
			} else {
				for (it=0; it<NLAT_2; it++) {
					zlm[im][(l-m-sk)*NLAT_2 + it*2 +sk + talign] = Z[it];
	#ifndef SHT_SCALAR_ONLY
					dzlm[im][(l-l0-sk)*NLAT_2 + it*2 +sk].p = dZp[it];
					dzlm[im][(l-l0-sk)*NLAT_2 + it*2 +sk].t = dZt[it];
	#endif
				}
			}
		}
		}

		// Compact the coefficients for improved cache efficiency.
		yg = ykm_dct[im];
		for (it=0; it<= KMAX; it+=2) {
			l = (it < m) ? m : it-(m&1);
			while (l<=LMAX) {
				yg[0] = yk[(it/2)*(LMAX+1-m) + (l-m)];
				l++;	yg++;
			}
			if ((m==0) && ((LMAX & 1) == 0)) {	yg[0] = 0;		yg++;	}		// SSE2 padding.
		}
#ifndef SHT_SCALAR_ONLY		
		dyg = dykm_dct[im];
		for (it=0; it<= KMAX; it+=2) {
			l = (it-2 < m) ? m : it-2+(m&1);
			while (l<=LMAX) {
				dyg[0].t = dyk[(it/2)*(LMAX+1-m) + (l-m)].t;
				dyg[0].p = dyk[(it/2)*(LMAX+1-m) + (l-m)].p;
				l++;	dyg++;
			}
		}
		if (im == 0) {		// compact and reorder m=0 dylm because .p = 0 :
			dyg = dykm_dct[im];
			yg = (double *) dykm_dct[im];
			for (it=0; it<= KMAX; it+=2) {
				dyg++;
				for (l=it-1; l<=LMAX; l++) {
					if (l>0) {
						yg[0] = dyg[0].t;
						yg++;	dyg++;
					}
				}
				if (LMAX&1) {		// padding for SSE2 alignement.
					yg[0] = 0.0;	yg++;
				}
			}
		}
#endif
	}
	
	// compact yk to zlm_dct0
	if (analysis) {
		long int klim = (LMAX * SHT_NL_ORDER) + 2;		// max k needed for nl-terms...
		klim = (klim/2)*2;		// must be even...
		if (klim > 2*NLAT_2) klim = 2*NLAT_2;		// but no more than 2*NLAT_2.
		shtns.klim = klim;		// store for use in codelets.
		yg = zlm_dct0;
		for (l=0; l<=LMAX; l+=2) {
			for (it=l; it<klim; it++) {	// for m=0, zl coeff with i<l are zeros.
				*yg = yk0[it];
				yg++;
			}
			yk0 += 2*NLAT_2;
		}
	#ifndef SHT_SCALAR_ONLY
		yg = dzlm_dct0;
		for (l=1; l<=LMAX; l+=2) {
			for (it=l-1; it<klim; it++) {	// for m=0, dzl coeff with i<l-1 are zeros.
				*yg = dyk0[it];
				yg++;
			}
			dyk0 += 2*NLAT_2;
		}
	#endif
		free(yk0 - (2*NLAT_2)*(LMAX/2+1));
	}

#if SHT_VERBOSE > 1
	tik1 = getticks();
	printf("\n    ticks : %.3f\n", elapsed(tik1,tik0)/(NLM*NLAT*(MMAX+1)));
#endif
	free(dyk);	free(yk);
	fftw_destroy_plan(idct);	fftw_destroy_plan(dct);
}
#endif

/// return the max error for a back-and-forth SHT transform.
/// this function is used to internally measure the accuracy.
double SHT_error()
{
	complex double *Tlm0, *Slm0, *Tlm, *Slm;
	double *Sh, *Th;
	double t, tmax, n2,  err;
	long int i, jj, nlm_cplx;
	
	srand( time(NULL) );	// init random numbers.
	
	Tlm0 = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Slm0 = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Slm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Tlm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Sh = (double *) fftw_malloc( NSPAT_ALLOC * sizeof(double) );
	Th = (double *) fftw_malloc( NSPAT_ALLOC * sizeof(double) );

// m = nphi/2 is also real if nphi is even.
	nlm_cplx = ( MMAX*2 == NPHI ) ? LiM(MRES*MMAX,MMAX) : NLM;
	t = 1.0 / (RAND_MAX/2);
	for (i=0; i<NLM; i++) {
		if ((i<=LMAX)||(i>=nlm_cplx)) {		// m=0 or m*2=nphi : real random data
			Slm0[i] = t*((double) (rand() - RAND_MAX/2));
			Tlm0[i] = t*((double) (rand() - RAND_MAX/2));
		} else {							// m>0 : complex random data
			Slm0[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
			Tlm0[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
		}
	}

	SH_to_spat(Slm0,Sh);		// scalar SHT
	spat_to_SH(Sh, Slm);
	for (i=0, tmax=0., n2=0., jj=0; i<NLM; i++) {		// compute error
		t = cabs(Slm[i] - Slm0[i]);
		n2 += t*t;
		if (t>tmax) { tmax = t; jj = i; }
	}
	err = tmax;
#if SHT_VERBOSE > 1
	printf("        scalar SH - poloidal   rms error = %.3g  max error = %.3g for l=%d,lm=%d\n",sqrt(n2/NLM),tmax,li[jj],jj);
#endif

	Slm0[0] = 0.0; 	Tlm0[0] = 0.0;		// l=0, m=0 n'a pas de signification sph/tor
	SHsphtor_to_spat(Slm0, Tlm0, Sh, Th);		// vector SHT
	spat_to_SHsphtor(Sh, Th, Slm, Tlm);
	for (i=0, tmax=0., n2=0., jj=0; i<NLM; i++) {		// compute error
		t = cabs(Slm[i] - Slm0[i]);
		n2 += t*t;
		if (t>tmax) { tmax = t; jj = i; }
	}
	if (tmax > err) err = tmax;
#if SHT_VERBOSE > 1
	printf("        vector SH - spheroidal rms error = %.3g  max error = %.3g for l=%d,lm=%d\n",sqrt(n2/NLM),tmax,li[jj],jj);
#endif
	for (i=0, tmax=0., n2=0., jj=0; i<NLM; i++) {		// compute error
		t = cabs(Tlm[i] - Tlm0[i]);
		n2 += t*t;
		if (t>tmax) { tmax = t; jj = i; }
	}
	if (tmax > err) err = tmax;
#if SHT_VERBOSE > 1
	printf("                  - toroidal   rms error = %.3g  max error = %.3g for l=%d,lm=%d\n",sqrt(n2/NLM),tmax,li[jj],jj);
#endif
	return(err);		// return max error.
}

#ifndef _HGID_
  #define _HGID_ "unknown"
#endif

/** \addtogroup init Initialization functions.
*/
//@{

/*! This sets the description of spherical harmonic coefficients.
 * It tells SHTns how to interpret spherical harmonic coefficient arrays, and it sets usefull arrays.
 * returns the number of modes (complex double) to describe a scalar field.
 * \param lmax : maximum SH degree that we want to describe.
 * \param mmax : number of azimutal wave numbers.
 * \param mres : \c 2.pi/mres is the azimutal periodicity. \c mmax*mres is the maximum SH order.
 * \param norm : define the normalization of the spherical harmonics (\ref shtns_norm)
 * + optionaly disable Condon-Shortley phase (ex: sht_schmidt | SHT_NO_CS_PHASE)
 * + optionaly use a 'real' normalization (ex: sht_fourpi | SHT_REAL_NORM)
*/
int shtns_set_size(int lmax, int mmax, int mres, enum shtns_norm norm)
{
	int im, m, l, lm;
	int with_cs_phase = 1;		/// Condon-Shortley phase (-1)^m is used by default.
	double mpos_renorm = 1.0;	/// renormalization of m>0.

//	if (lmax < 1) shtns_runerr("lmax must be larger than 1");
	if (lmax < 2) shtns_runerr("lmax must be at least 2");
	if (li != NULL) {
		if ( (lmax != LMAX)||(mmax != MMAX)||(mres != MRES) )
			shtns_runerr("different size already set");
		return(NLM);
	}

	// copy to global variables.
	shtns.norm = norm;
	if (norm & SHT_NO_CS_PHASE)
		with_cs_phase = 0;
	if (norm & SHT_REAL_NORM)
		mpos_renorm = 0.5;		// normalization for 'real' spherical harmonics.

#ifdef SHT_AXISYM
	shtns.mmax = 0;		shtns.mres = 2;		shtns.nphi = 1;
	if (mmax != 0) shtns_runerr("axisymmetric version : only Mmax=0 allowed");
#else
	MMAX = mmax;	MRES = mres;
#endif
	LMAX = lmax;
	NLM = nlm_calc(LMAX,MMAX,MRES);
#if SHT_VERBOSE > 0
	printf("[SHTns] build " __DATE__ ", " __TIME__ ", id: " _HGID_ "\n");
	printf("        Lmax=%d, Mmax*Mres=%d, Mres=%d, Nlm=%d  [",LMAX,MMAX*MRES,MRES,NLM);
	if (norm & SHT_REAL_NORM) printf("'real' norm, ");
	if (!with_cs_phase) printf("no Condon-Shortley phase, ");
	if (SHT_NORM == sht_fourpi) printf("4.pi normalized]\n");
	else if (SHT_NORM == sht_schmidt) printf("Schmidt semi-normalized]\n");
	else printf("orthonormalized]\n");
#endif
	if (MMAX*MRES > LMAX) shtns_runerr("MMAX*MRES should not exceed LMAX");
	if (MRES <= 0) shtns_runerr("MRES must be > 0");

	// alloc spectral arrays
	shtns.lmidx = (int *) fftw_malloc(sizeof(int) * (MMAX+1)*2);
	tm =  shtns.lmidx + (MMAX+1);
	el = (double *) fftw_malloc( 3*NLM*sizeof(double) + NLM*sizeof(int) );	// NLM defined at runtime.
	l2 = el + NLM;	l_2 = el + 2*NLM;
	li = (int *) (el + 3*NLM);

	for (im=0, lm=0; im<=MMAX; im++) {		// init l-related arrays.
		m = im*MRES;
		shtns.lmidx[im] = lm -m;		// virtual pointer for l=0
		for (l=im*MRES;l<=LMAX;l++) {
			el[lm] = l;	l2[lm] = l*(l+1.0);	l_2[lm] = 1.0/(l*(l+1.0));
			li[lm] = l;
			lm++;
		}
	}
	l_2[0] = 0.0;	// undefined for l=0 => replace with 0.	
	if (lm != NLM) shtns_runerr("unexpected error");

	// this quickly precomputes some values for the legendre recursion.
	legendre_precomp(SHT_NORM, with_cs_phase, mpos_renorm);

	switch(SHT_NORM) {
		case sht_schmidt:
			Y00_1 = 1.0;		Y10_ct = 1.0;
			break;
		case sht_fourpi:
			Y00_1 = 1.0;		Y10_ct = sqrt(1./3.);
			break;
		case sht_orthonormal:
		default:
			Y00_1 = sqrt(4.*M_PI);		Y10_ct = sqrt(4.*M_PI/3.);
//			Y11_st = sqrt(2.*M_PI/3.);		// orthonormal :  \f$ \sin\theta\cos\phi/(Y_1^1 + Y_1^{-1}) = -\sqrt{2 \pi /3} \f$
	}
	Y11_st = Y10_ct * sqrt(0.5/mpos_renorm);
	if (with_cs_phase)	Y11_st *= -1.0;		// correct Condon-Shortley phase

	return(NLM);
}


/*! Initialization of Spherical Harmonic transforms (backward and forward, vector and scalar, ...) of given size.
 * <b>This function must be called after \ref shtns_set_size and before any SH transform.</b> and sets all global variables and internal data.
 * returns the required number of doubles to be allocated for a spatial field.
 * \param nlat,nphi pointers to the number of latitudinal and longitudinal grid points respectively. If 0, they are set to optimal values.
 * \param nl_order defines the maximum SH degree to be resolved by analysis : lmax_analysis = lmax*nl_order. It is used to set an optimal and anti-aliasing nlat. If 0, the default SHT_DEFAULT_NL_ORDER is used.
 * \param flags allows to choose the type of transform (see \ref shtns_type) and the spatial data layout (see \ref spat)
 * \param eps polar optimization threshold : polar values of Legendre Polynomials below that threshold are neglected (for high m), leading to increased performance (a few percents)
 *  0 = no polar optimization;  1.e-14 = VERY safe;  1.e-10 = safe;  1.e-6 = aggresive, but still good accuracy.
*/
int shtns_precompute_auto(enum shtns_type flags, double eps, int nl_order, int *nlat, int *nphi)
{
	double t;
	int im,m,l,lm;
	int layout, nspat;
	int n_gauss = 0;

	if (nl_order <= 0) nl_order = SHT_DEFAULT_NL_ORDER;
/*	shtns.lshift = 0;
	if (nl_order == 0) nl_order = SHT_DEFAULT_NL_ORDER;
	if (nl_order < 0) {	shtns.lshift = -nl_order;	nl_order = 1; }		// linear with a shift in l.
*/
	shtns.nlorder = nl_order;
	shtns.mtr_dct = -1;
	layout = flags & 0xFFFF00;
	flags = flags & 255;	// clear higher bits.

	if ((flags == sht_reg_poles)||(flags == sht_quick_init))
		fftw_plan_mode = FFTW_ESTIMATE;		// quick fftw init

	if (*nphi == 0) {
		*nphi = fft_int((nl_order+1)*MMAX+1, 7);		// required fft nodes
	}
	if (*nlat == 0) {
		n_gauss = ((nl_order+1)*LMAX)/2 +1;		// required gauss nodes
		n_gauss += (n_gauss&1);		// even is better.
		if ((flags == sht_auto)||(flags == sht_reg_fast)) {
			m = fft_int(nl_order*LMAX+1, 7);		// required dct nodes
			*nlat = m + (m&1);		// even is better.
		} else *nlat = n_gauss;
	}

	// copy to global variables.
#ifdef SHT_AXISYM
	shtns.nphi = 1;
	if (*nphi != 1) shtns_runerr("axisymmetric version : only Nphi=1 allowed");
#else
	NPHI = *nphi;
#endif
	NLAT_2 = (*nlat+1)/2;	NLAT = *nlat;
	if ((NLAT_2)*2 <= LMAX) shtns_runerr("NLAT_2*2 should be at least LMAX+1");

	alloc_SHTarrays();		// allocate dynamic arrays
	nspat = planFFT(layout);		// initialize fftw
	zlm_dct0 = NULL;		// used as a flag.

#ifndef SHT_NO_DCT
	if (flags == sht_reg_dct) {		// pure dct.
		init_SH_dct(1);
		OptimizeMatrices(0.0);
		Set_MTR_DCT(MMAX);
	}
	if ((flags == sht_auto)||(flags == sht_reg_fast))
	{
		init_SH_dct(1);
		OptimizeMatrices(eps);
  #if SHT_VERBOSE > 0
		printf("finding optimal MTR_DCT ...\r");	fflush(stdout);
  #endif
		t = Find_Optimal_SHT();
  #if SHT_VERBOSE > 0
		printf("        + optimal MTR_DCT = %d  (%.1f%% performance gain)\n", MTR_DCT*MRES, 100.*(1/t-1));
  #endif
  		if ((n_gauss > 0)&&(flags == sht_auto)) t *= ((double) NLAT)/n_gauss;	// we can revert to gauss with a smaller nlat.
		if (t > MIN_PERF_IMPROVE_DCT) {
			Set_MTR_DCT(-1);		// turn off DCT.
		} else {
			t = SHT_error();
			if (t > MIN_ACCURACY_DCT) {
  #if SHT_VERBOSE > 0
				printf("     !! Not enough accuracy (%.3g) => turning off DCT.\n",t);
  #endif
  #if SHT_VERBOSE < 2
				Set_MTR_DCT(-1);		// turn off DCT.
  #endif
			}
		}
  #if SHT_VERBOSE < 2
		if (MTR_DCT == -1) {			// free memory used by DCT and disables DCT.
			free_SH_dct();			// free now useless arrays.
			zlm_dct0 = NULL;		// this disables DCT completely.
			if (idct != NULL) fftw_destroy_plan(idct);	// free unused dct plans
			if (dct_m0 != NULL) fftw_destroy_plan(dct_m0);
			if (flags == sht_auto) {
				flags = sht_gauss;		// switch to gauss grid, even better accuracy.
	#if SHT_VERBOSE > 0
				printf("        => switching back to Gauss Grid\n");
	#endif
				for (im=1; im<=MMAX; im++) {	//	im >= 1
					m = im*MRES;
					ylm[im]  -= tm[im]*(LMAX-m+1);		// restore pointers altered by OptimizeMatrices().
#ifndef SHT_SCALAR_ONLY
					dylm[im] -= tm[im]*(LMAX-m+1);
#endif
				}
				if (n_gauss > 0) {		// we should use the optimal size for gauss-legendre
					if (NPHI>1) {
						fftw_destroy_plan(fft_eo);  fftw_destroy_plan(ifft_eo);		fftw_destroy_plan(fft);  fftw_destroy_plan(ifft);
					}
					free_SHTarrays();
					*nlat = n_gauss;
					NLAT_2 = (*nlat+1)/2;	NLAT = *nlat;
					alloc_SHTarrays();
					nspat = planFFT(layout);		// fft must be replanned because NLAT has changed.
				}
			}
		}
  #endif
	}
	if (flags == sht_gauss)
#else
	if ((flags == sht_gauss)||(flags == sht_auto)||(flags == sht_reg_fast)||(flags == sht_reg_dct))
#endif
	{
		MTR_DCT = -1;		// we do not use DCT !!!
		init_SH_gauss();
		OptimizeMatrices(eps);
	}
	if (flags == sht_reg_poles)
	{
		MTR_DCT = -1;		// we do not use DCT !!!
	  #ifndef SHT_SCALAR_ONLY
		fftw_free(dzlm[0]);		dzlm[0] = NULL;
	  #endif
		fftw_free(zlm[0]);		zlm[0] = NULL;		// no inverse transform, mark as unused.
		EqualPolarGrid();
		init_SH_synth();
		OptimizeMatrices(0.0);
	}
	if (flags == sht_quick_init)
	{
		MTR_DCT = -1;		// we do not use DCT !!!
		init_SH_gauss();
		OptimizeMatrices(eps);
	}

	if ((flags != sht_reg_poles)&&(flags != sht_quick_init)) {
		t = SHT_error();		// compute SHT accuracy.
  #if SHT_VERBOSE > 0
		printf("        + SHT accuracy = %.3g\n",t);
  #endif
  #if SHT_VERBOSE < 2
		if (t > 1.e-3) shtns_runerr("bad SHT accuracy");		// stop if something went wrong (but not in debug mode)
  #endif
	}
	return(nspat);	// returns the number of doubles to be allocated for a spatial field.
}


/*! Initialization of Spherical Harmonic transforms (backward and forward, vector and scalar, ...) of given size.
 * <b>This function must be called after \ref shtns_set_size and before any SH transform.</b> and sets all global variables.
 * returns the required number of doubles to be allocated for a spatial field.
 * \param nlat,nphi respectively the number of latitudinal and longitudinal grid points.
 * \param flags allows to choose the type of transform (see \ref shtns_type) and the spatial data layout (see \ref spat)
 * \param eps polar optimization threshold : polar values of Legendre Polynomials below that threshold are neglected (for high m), leading to increased performance (a few percents)
 *  0 = no polar optimization;  1.e-14 = VERY safe;  1.e-10 = safe;  1.e-6 = aggresive, but still good accuracy.
*/
int shtns_precompute(enum shtns_type flags, double eps, int nlat, int nphi)
{
	if ((nlat == 0)||(nphi == 0)) shtns_runerr("nlat or nphi is zero !");
	return( shtns_precompute_auto(flags, eps, 0, &nlat, &nphi) );
}


/*! Simple initialization of Spherical Harmonic transforms (backward and forward, vector and scalar, ...) of given size.
 * This function sets all global variables by calling \ref shtns_set_size followed by \ref shtns_precompute, with the
 * default normalization and the default polar optimization (see \ref sht_config.h).
 * returns the number of modes to describe a scalar field.
 * \param lmax : maximum SH degree that we want to describe.
 * \param mmax : number of azimutal wave numbers.
 * \param mres : \c 2.pi/mres is the azimutal periodicity. \c mmax*mres is the maximum SH order.
 * \param nlat,nphi : respectively the number of latitudinal and longitudinal grid points.
 * \param flags allows to choose the type of transform (see \ref shtns_type) and the spatial data layout (see \ref spat)
*/
int shtns_init(enum shtns_type flags, int lmax, int mmax, int mres, int nlat, int nphi)
{
	shtns_set_size(lmax, mmax, mres, SHT_DEFAULT_NORM);
	return shtns_precompute(flags, SHT_DEFAULT_POLAR_OPT, nlat, nphi);
}

//@}


#ifdef SHT_F77_API

/** \addtogroup fortapi Fortran API.
* Call from fortran without the trailing '_'.
* see the \link SHT_example.f Fortran example \endlink for a simple usage of SHTns from Fortran language.
*/
//@{
	
/// Initializes spherical harmonic transforms of given size using Gauss algorithm and no approximation
void shtns_init_sh_gauss_(int *layout, int *lmax, int *mmax, int *mres, int *nlat, int *nphi)
{
	shtns_set_size(*lmax, *mmax, *mres, SHT_DEFAULT_NORM);
	shtns_precompute(sht_gauss | *layout, 0, *nlat, *nphi);
}

/// Initializes spherical harmonic transforms of given size using Fastest available algorithm and polar optimization.
void shtns_init_sh_auto_(int *layout, int *lmax, int *mmax, int *mres, int *nlat, int *nphi)
{
	shtns_set_size(*lmax, *mmax, *mres, SHT_DEFAULT_NORM);
	shtns_precompute(sht_auto | *layout, 1.e-10, *nlat, *nphi);
}

/// Initializes spherical harmonic transforms of given size using a regular grid and agressive optimizations.
void shtns_init_sh_reg_fast_(int *layout, int *lmax, int *mmax, int *mres, int *nlat, int *nphi)
{
	shtns_set_size(*lmax, *mmax, *mres, SHT_DEFAULT_NORM);
	shtns_precompute(sht_reg_fast | *layout, 1.e-6, *nlat, *nphi);
}

/// Initializes spherical harmonic transform SYNTHESIS ONLY of given size using a regular grid including poles.
void shtns_init_sh_poles_(int *layout, int *lmax, int *mmax, int *mres, int *nlat, int *nphi)
{
	shtns_set_size(*lmax, *mmax, *mres, SHT_DEFAULT_NORM);
	shtns_precompute(sht_reg_poles | *layout, 0, *nlat, *nphi);
}

/// Defines the size and convention of the transform.
/// Allow to choose the normalization and whether or not to include the Condon-Shortley phase.
/// \see shtns_set_size
void shtns_set_size_(int *lmax, int *mmax, int *mres, int *norm)
{
	shtns_set_size(*lmax, *mmax, *mres, *norm);
}

/// Precompute matrices for synthesis and analysis.
/// Allow to choose polar optimization threshold and algorithm type.
/// \see shtns_precompute
void shtns_precompute_(int *type, int *layout, double *eps, int *nlat, int *nphi)
{
	shtns_precompute(*type | *layout, *eps, *nlat, *nphi);
}

/// Same as shtns_precompute_ but choose optimal nlat and/or nphi.
void shtns_precompute_auto_(int *type, int *layout, double *eps, int *nl_order, int *nlat, int *nphi)
{
	shtns_precompute_auto(*type | *layout, *eps, *nl_order, nlat, nphi);
}

/// returns nlm, the number of complex*16 elements in an SH array.
/// call from fortran using \code call shtns_calc_nlm(nlm, lmax, mmax, mres) \endcode
void shtns_calc_nlm_(int *nlm, int *lmax, int *mmax, int *mres)
{
    *nlm = nlm_calc(*lmax, *mmax, *mres);
}

/// returns lm, the index in an SH array of mode (l,m).
/// call from fortran using \code call shtns_lmidx(lm, l, m) \endcode
void shtns_lmidx_(int *lm, int *l, int *m)
{
    *lm = LM(*l, *m) + 1;
}

/// fills the given array with the cosine of the co-latitude angle (NLAT real*8)
void shtns_cos_array_(double *costh)
{
	int i;	
	for (i=0; i<shtns.nlat; i++)
		costh[i] = ct[i];
}

/** \name Point evaluation of Spherical Harmonics
Evaluate at a given point (\f$cos(\theta)\f$ and \f$\phi\f$) a spherical harmonic representation.
*/
//@{
/// \see SH_to_point for argument description
void shtns_sh_to_point_(double *spat, complex double *Qlm, double *cost, double *phi)
{
	*spat = SH_to_point(Qlm, *cost, *phi);
}

/// \see SHqst_to_point for argument description
void shtns_qst_to_point_(double *vr, double *vt, double *vp,
		complex double *Qlm, complex double *Slm, complex double *Tlm, double *cost, double *phi)
{
	SHqst_to_point(Qlm, Slm, Tlm, *cost, *phi, vr, vt, vp);
}
//@}

//@}

#endif
