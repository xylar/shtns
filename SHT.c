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
 *    written by Nathanael Schaeffer / CNRS                         *
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

// chained list of sht_setup : start with NULL
shtns_cfg sht_data = NULL;

/// Abort program with error message.
void shtns_runerr(const char * error_text)
{
	printf("*** [SHTns] Run-time error : %s\n",error_text);
	exit(1);
}

/// returns the l=0, m=0 SH coefficient corresponding to a uniform value of 1.
double sh00_1(shtns_cfg shtns) {
	return shtns->Y00_1;
}
/// returns the l=1, m=0 SH coefficient corresponding to cos(theta).
double sh10_ct(shtns_cfg shtns) {
	return shtns->Y10_ct;
}
/// returns the l=1, m=1 SH coefficient corresponding to sin(theta).cos(phi).
double sh11_st(shtns_cfg shtns) {
	return shtns->Y11_st;
}
/// returns the l,m SH coefficient corresponding to unit energy.
double shlm_e1(shtns_cfg shtns, int l, int m) {
	double x = shtns->Y00_1/sqrt(4.*M_PI);
	if (SHT_NORM == sht_schmidt) x *= sqrt(2*l+1);
	if ((m!=0)&&((shtns->norm & SHT_REAL_NORM)==0)) x *= sqrt(0.5);
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



/// \code return (mmax+1)*(lmax+1) - mres*(mmax*(mmax+1))/2; \endcode */
/// \ingroup init
int nlm_calc(int lmax, int mmax, int mres)
{
	if (mmax*mres > lmax) mmax = lmax/mres;
	return( (mmax+1)*(lmax+1) - mres*(mmax*(mmax+1))/2 );	// this is wrong if lmax < mmax*mres
}

/// returns an aproximation of the memory usage in mega bytes.
/// \ingroup init
double sht_mem_size(int lmax, int mmax, int mres, int nlat)
{
	double s = 1./(1024*1024);
	s *= (nlat+1) * sizeof(double) * nlm_calc(lmax, mmax, mres);
  #ifndef SHT_SCALAR_ONLY
	s *= 3.0;		// scalar + vector arrays.
  #endif
	return s;
}

/*
	SHT FUNCTIONS
*/

// truncation at LMAX and MMAX
#define LTR LMAX
#define MTR MMAX

/** \addtogroup local Local and partial evaluation of SH fields.
 * These do only require a call to \ref shtns_create, but not to \ref shtns_set_grid.
*/
//@{

/// Evaluate scalar SH representation \b Qlm at physical point defined by \b cost = cos(theta) and \b phi
double SH_to_point(shtns_cfg shtns, complex double *Qlm, double cost, double phi)
{
	double yl[LMAX+1];
	double vr;
	complex double *Ql;
	long int l,m,im;

	vr = 0.0;
	m=0;	im=0;
		legendre_sphPlm_array(shtns, LTR, im, cost, &yl[m]);
		for (l=m; l<=LTR; l++)
			vr += yl[l] * creal( Qlm[l] );
	if (MTR>0) {
		complex double eip, eimp;
		eip = cos(phi*MRES) + I*sin(phi*MRES);	eimp = 2.0;
		for (im=1; im<=MTR; im++) {
			m = im*MRES;
			legendre_sphPlm_array(shtns, LTR, im, cost, &yl[m]);
//			eimp = 2.*(cos(m*phi) + I*sin(m*phi));
			eimp *= eip;			// not so accurate, but it should be enough for rendering uses.
			Ql = &Qlm[LiM(shtns, 0,im)];	// virtual pointer for l=0 and im
			for (l=m; l<=LTR; l++)
				vr += yl[l] * creal( Ql[l]*eimp );
		}
	}
	return vr;
}

/// Evaluate vector SH representation \b Qlm at physical point defined by \b cost = cos(theta) and \b phi
void SHqst_to_point(shtns_cfg shtns, complex double *Qlm, complex double *Slm, complex double *Tlm, double cost, double phi,
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
		legendre_sphPlm_deriv_array(shtns, LTR, im, cost, sint, &yl[m], &dtyl[m]);
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
			legendre_sphPlm_deriv_array(shtns, LTR, im, cost, sint, &yl[m], &dtyl[m]);
//			eimp = 2.*(cos(m*phi) + I*sin(m*phi));
			eimp *= eip;		// not so accurate, but it should be enough for rendering uses.
			imeimp = I*m*eimp;
			l = LiM(shtns, 0,im);
			Ql = &Qlm[l];	Sl = &Slm[l];	Tl = &Tlm[l];
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
	(does not require a previous call to shtns_set_grid)
*/

fftw_plan ifft_lat = NULL;		///< fftw plan for SHqst_to_lat
long int nphi_lat = 0;			///< nphi of previous SHqst_to_lat
double* ylm_lat = NULL;
double* dylm_lat;
double ct_lat = 2.0;
double st_lat;

/// synthesis at a given latitude, on nphi equispaced longitude points.
/// vr, vt, and vp arrays must have nphi+2 doubles allocated (fftw requirement).
/// It does not require a previous call to shtns_set_grid
/// \ingroup local
void SHqst_to_lat(shtns_cfg shtns, complex double *Qlm, complex double *Slm, complex double *Tlm, double cost,
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
			legendre_sphPlm_deriv_array(shtns, ltr, m, cost, st_lat, &ylm_lat[j], &dylm_lat[j]);
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

/// allocate arrays for SHT related to a given grid.
void alloc_SHTarrays(shtns_cfg shtns, int on_the_fly)
{
	long int im,m, l0;
	long int lstride;

	l0 = ((NLAT+1)>>1)*2;		// round up to even
	shtns->ct = (double *) fftw_malloc( sizeof(double) * l0*3 );			/// ct[] (including st and st_1)
	shtns->st = shtns->ct + l0;		shtns->st_1 = shtns->ct + 2*l0;

  if (on_the_fly == 0) {
	shtns->ylm = (double **) fftw_malloc( sizeof(double *) * (MMAX+1)*3 );		/// ylm[] (including zlm and ykm_dct)
	shtns->zlm = shtns->ylm + (MMAX+1);		shtns->ykm_dct = shtns->ylm + (MMAX+1)*2;
  #ifndef SHT_SCALAR_ONLY
	shtns->dylm = (struct DtDp **) fftw_malloc( sizeof(struct DtDp *) * (MMAX+1)*3);		/// dylm[] (including dzlm and dykm_dct)
	shtns->dzlm = shtns->dylm + (MMAX+1);		shtns->dykm_dct = shtns->dylm + (MMAX+1)*2;
  #endif

// Allocate legendre functions lookup tables.
	lstride = (LMAX+1);		lstride += (lstride&1);		// even stride.
	shtns->ylm[0] = (double *) fftw_malloc(sizeof(double)* (NLM-(LMAX+1)+lstride)*NLAT_2);		/// ylm[][]
  #ifndef SHT_SCALAR_ONLY
	shtns->dylm[0] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* NLM*NLAT_2);				/// dylm[][]
  #endif
	shtns->zlm[0] = (double *) fftw_malloc(sizeof(double)* (NLM*NLAT_2 + (NLAT_2 & 1)));		/// zlm[][]
  #ifndef SHT_SCALAR_ONLY
	shtns->dzlm[0] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* (NLM-1)*NLAT_2);		// remove l=0
  #endif
	for (im=0; im<MMAX; im++) {
		m = im*MRES;	l0 = (m==0) ? 1 : m;
		if (im>0) lstride = (LMAX+1 -m);
		shtns->ylm[im+1] = shtns->ylm[im] + NLAT_2*lstride;
		shtns->zlm[im+1] = shtns->zlm[im] + NLAT_2*(LMAX+1 -m) + ((m==0)*(NLAT_2&1));
  #ifndef SHT_SCALAR_ONLY
		shtns->dylm[im+1] = shtns->dylm[im] + NLAT_2*(LMAX+1 -m);
		shtns->dzlm[im+1] = shtns->dzlm[im] + NLAT_2*(LMAX+1 -l0);
  #endif
	}
#if SHT_VERBOSE > 1
	printf("          Memory used for Ylm and Zlm matrices = %.3f Mb x2\n",3.0*sizeof(double)*NLM*NLAT_2/(1024.*1024.));
#endif
  }
}

void free_SHTarrays(shtns_cfg shtns)
{
	if (shtns->fft != NULL)		fftw_destroy_plan(shtns->fft);
	if (shtns->ifft != NULL)	fftw_destroy_plan(shtns->ifft);
	shtns->fft = NULL;		shtns->ifft = NULL;		shtns->sht_fft = 0;

	if (shtns->dylm != NULL) {
		if (shtns->dzlm[0] != NULL) fftw_free(shtns->dzlm[0]);
		shtns->dzlm[0] = NULL;
		if (shtns->dylm != NULL) fftw_free(shtns->dylm[0]);
		shtns->dylm[0] = NULL;
		fftw_free(shtns->dylm);		shtns->dylm = NULL;
	}
	if (shtns->ylm != NULL) {
		if (shtns->zlm[0] != NULL) fftw_free(shtns->zlm[0]);
		shtns->zlm[0] = NULL;
		if (shtns->ylm[0] != NULL) fftw_free(shtns->ylm[0]);
		shtns->ylm[0] = NULL;
		fftw_free(shtns->ylm);		shtns->ylm = NULL;
	}
	if (shtns->ct != NULL) {
		fftw_free(shtns->ct);		shtns->ct = NULL;
	}
}

/// Generates an equi-spaced theta grid including the poles, for synthesis only.
void EqualPolarGrid(shtns_cfg shtns)
{
	int j;
	double f;
	
#if SHT_VERBOSE > 0
	printf("        => using Equaly Spaced Nodes including poles\n");
#endif
// cos theta of latidunal points (equaly spaced in theta)
	f = M_PI/(NLAT-1.0);
	for (j=0; j<NLAT; j++) {
		shtns->ct[j] = cos(f*j);
		shtns->st[j] = sin(f*j);
		shtns->st_1[j] = 1.0/(shtns->st[j]);
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
int planFFT(shtns_cfg shtns, int layout)
{
	double cost_ip, cost_oop;
	complex double *ShF;
	double *Sh;
	fftw_plan fft2, ifft2, fft, ifft;
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
#endif

// Allocate dummy Spatial Fields.
	ShF = (complex double *) fftw_malloc(ncplx * NLAT * sizeof(complex double));
	Sh = (double *) fftw_malloc(ncplx * NLAT * sizeof(complex double));
	fft = NULL;		ifft = NULL;

// IFFT : unnormalized.  FFT : must be normalized.
	cost_ip = 0.0;		cost_oop = 0.0;
	if (SHT_FFT == 1) {		// in-place FFT allowed
		ifft2 = fftw_plan_many_dft_c2r(1, &nfft, NLAT, ShF, &ncplx, NLAT, 1, (double*) ShF, &nreal, phi_inc, theta_inc, shtns->fftw_plan_mode);
		if (ifft2 != NULL) {
			fft2 = fftw_plan_many_dft_r2c(1, &nfft, NLAT, (double*) ShF, &nreal, phi_inc, theta_inc, ShF, &ncplx, NLAT, 1, shtns->fftw_plan_mode);
			if (fft2 != NULL) {
//				cost_ip = fftw_cost(ifft2) + fftw_cost(fft2);
			} else {
				fftw_destroy_plan(ifft2);	ifft2 = NULL;	SHT_FFT = 2;
			}
		} else SHT_FFT = 2;
	}
	if ((SHT_FFT > 1) || (cost_ip > 0.0)) {		// out-of-place FFT
		ifft = fftw_plan_many_dft_c2r(1, &nfft, NLAT, ShF, &ncplx, NLAT, 1, Sh, &nreal, phi_inc, theta_inc, shtns->fftw_plan_mode);
		if (ifft == NULL) shtns_runerr("[FFTW] ifft planning failed !");
		fft = fftw_plan_many_dft_r2c(1, &nfft, NLAT, Sh, &nreal, phi_inc, theta_inc, ShF, &ncplx, NLAT, 1, shtns->fftw_plan_mode);
		if (fft == NULL) shtns_runerr("[FFTW] fft planning failed !");
//		cost_oop = fftw_cost(ifft) + fftw_cost(fft);
	}
	if (SHT_FFT == 1) {
		if ((cost_oop >= cost_ip) || (cost_oop == 0.0)) {		// switch to in-place transforms.
			if (fft != NULL) fftw_destroy_plan(fft);
			if (ifft != NULL) fftw_destroy_plan(ifft);
			fft = fft2;		ifft = ifft2;
		} else {		// out-of-place is faster
			fftw_destroy_plan(fft2);	fftw_destroy_plan(ifft2);
			SHT_FFT = 2;		// switch to out-of-place.
		}
	}
	shtns->fft = fft;		shtns->ifft = ifft;

#if SHT_VERBOSE > 0
	if (SHT_FFT > 1) printf("        ** out-of-place fft **\n");
#endif
#if SHT_VERBOSE > 2
	printf("out-of-place cost = %g    in-place cost = %g\n",cost_oop, cost_ip);
	printf(" *** fft plan :\n");
	fftw_print_plan(fft);
	printf("\n *** ifft plan :\n");
	fftw_print_plan(ifft);
	printf("\n");
#endif

	fftw_free(Sh);		fftw_free(ShF);
  } else {
	if (theta_inc != 1) shtns_runerr("only contiguous spatial data is supported for Nphi=1");
#if SHT_VERBOSE > 0
	printf("        => no fft : Mmax=0, Nphi=1, Nlat=%d\n",NLAT);
#endif
	SHT_FFT = 0;	// no fft.
  }
	shtns->dct_m0 = NULL;	shtns->idct = NULL;		// set dct plans to uninitialized.
	shtns->dct_r1 = NULL;	shtns->idct_r1 = NULL;
	return(phi_embed * NLAT);
}

#ifndef SHT_NO_DCT
/// initialize DCTs using FFTW. Must be called if MTR_DCT is changed.
void planDCT(shtns_cfg shtns)
{
	double *Sh;
	int ndct = NLAT;
	fftw_r2r_kind r2r_kind;
	fftw_iodim dims, hdims[2];
	double Sh0[NLAT] SSE;				// temp storage on the stack, aligned.
	
// real NPHI=1, allocate only once since it does not change.
	if ((shtns->dct_r1 == NULL)||(shtns->idct_r1 == NULL)) {
		Sh = (double *) fftw_malloc( NLAT * sizeof(double) );
		if (shtns->dct_r1 == NULL) {
			r2r_kind = FFTW_REDFT10;
			shtns->dct_r1 = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, 1, NLAT, Sh, &ndct, 1, NLAT, &r2r_kind, shtns->fftw_plan_mode);
		}
		if (shtns->idct_r1 == NULL) {
			r2r_kind = FFTW_REDFT01;
			shtns->idct_r1 = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, 1, NLAT, Sh, &ndct, 1, NLAT, &r2r_kind, shtns->fftw_plan_mode);
		}
		fftw_free(Sh);
		if ((shtns->dct_r1 == NULL)||(shtns->idct_r1 == NULL))
			shtns_runerr("[FFTW] (i)dct_r1 planning failed !");
#if SHT_VERBOSE > 2
			printf(" *** idct_r1 plan :\n");	fftw_print_plan(shtns->idct_r1);
			printf("\n *** dct_r1 plan :\n");	fftw_print_plan(shtns->dct_r1);	printf("\n");
#endif
	}

#ifndef SHT_AXISYM
	if (shtns->idct != NULL) fftw_destroy_plan(shtns->idct);
	// Allocate dummy Spatial Fields.
	Sh = (double *) fftw_malloc((NPHI/2 +1) * NLAT*2 * sizeof(double));

	dims.n = NLAT;	dims.is = 2;	dims.os = 2;		// real and imaginary part.
	hdims[0].n = MTR_DCT+1;	hdims[0].is = 2*NLAT; 	hdims[0].os = 2*NLAT;
	hdims[1].n = 2;			hdims[1].is = 1; 	hdims[1].os = 1;

	if (NPHI>1) {		// complex data for NPHI>1, recompute as it does depend on MTR_DCT
		r2r_kind = FFTW_REDFT01;
		shtns->idct = fftw_plan_guru_r2r(1, &dims, 2, hdims, Sh, Sh, &r2r_kind, shtns->fftw_plan_mode);
		if (shtns->idct == NULL)
			shtns_runerr("[FFTW] idct planning failed !");
#if SHT_VERBOSE > 2
			printf(" *** idct plan :\n");	fftw_print_plan(shtns->idct);	printf("\n");
#endif
		if (shtns->dct_m0 == NULL) {
			r2r_kind = FFTW_REDFT10;
//			shtns->dct_m0 = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, 2, 2*NLAT, Sh, &ndct, 2, 2*NLAT, &r2r_kind, shtns->fftw_plan_mode);
			shtns->dct_m0 = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, 2, 2*NLAT, Sh0, &ndct, 1, NLAT, &r2r_kind, shtns->fftw_plan_mode);	// out-of-place.
			if (shtns->dct_m0 == NULL)
				shtns_runerr("[FFTW] dct_m0 planning failed !");
#if SHT_VERBOSE > 2
				printf(" *** dct_m0 plan :\n");		fftw_print_plan(shtns->dct_m0);	printf("\n");
#endif
		}
	} else {	// NPHI==1
		if (shtns->dct_m0 == NULL) {
			r2r_kind = FFTW_REDFT10;
			shtns->dct_m0 = fftw_plan_many_r2r(1, &ndct, 1, Sh, &ndct, 1, NLAT, Sh0, &ndct, 1, NLAT, &r2r_kind, shtns->fftw_plan_mode);	// out-of-place.
			if (shtns->dct_m0 == NULL)
				shtns_runerr("[FFTW] dct_m0 planning failed !");
#if SHT_VERBOSE > 2
				printf(" *** dct_m0 plan :\n");		fftw_print_plan(shtns->dct_m0);	printf("\n");
#endif
		}
	}
	fftw_free(Sh);
#endif
}
#endif

/// SET MTR_DCT and updates fftw_plan for DCT's
void Set_MTR_DCT(shtns_cfg shtns, int m)
{
#ifndef SHT_NO_DCT
	if ((shtns->zlm_dct0 == NULL)||(m == MTR_DCT)) return;
	if ( m < 0 ) {	// don't use dct
		MTR_DCT = -1;
	} else {
		if (m>MMAX) m=MMAX;
		MTR_DCT = m;
		planDCT(shtns);
	}
#endif
}

int Get_MTR_DCT(shtns_cfg shtns) {
	return MTR_DCT;
}

/// Sets the value tm[im] used for polar optimiation on-the-fly.
void PolarOptimize(shtns_cfg shtns, double eps)
{
	int im, m, l, it;
	double v;
	double y[LMAX+1];

	for (im=0;im<=MMAX;im++)	shtns->tm[im] = 0;

	if (eps > 0.0) {
		for (im=1;im<=MMAX;im++) {
			m = im*MRES;
			it = -1;
			do {
				it++;
				legendre_sphPlm_array(shtns, LMAX, im, shtns->ct[it], y+m);
				v = 0.0;
				for (l=m; l<=LMAX; l++) {
					double ya = fabs(y[l]);
					if ( v < ya )	v = ya;
				}
			} while (v < eps);
			shtns->tm[im] = it;
		}
	#if SHT_VERBOSE > 0
		printf("        + polar optimization threshold = %.1e\n",eps);
	#endif
	#if SHT_VERBOSE > 1
		printf("          tm[im]=");
		for (im=0;im<=MMAX;im++)
			printf(" %d",shtns->tm[im]);
		printf("\n");
	#endif
	}
}

/// Perform some optimization on the SHT matrices.
void OptimizeMatrices(shtns_cfg shtns, double eps)
{
	unsigned short *tm;
	double **ylm, **zlm;
	struct DtDp** dylm;
	struct DtDp** dzlm;
	double *yg;
	int im,m,l,it;

	tm = shtns->tm;
	ylm = shtns->ylm;		dylm = shtns->dylm;
	zlm = shtns->zlm;		dzlm = shtns->dzlm;
	
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
void init_SH_synth(shtns_cfg shtns)
{
	double dtylm[LMAX+1];		// temp storage for derivative : d(P_l^m(x))/dx
	long int it,im,m,l;
	long int lstride;
	double *ct = shtns->ct;
	double *st = shtns->st;
	double **ylm = shtns->ylm;
	struct DtDp** dylm = shtns->dylm;

	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
		lstride = LMAX-m+1;
		if (im==0) lstride += (lstride&1);		// even stride for m=0.
		for (it=0; it<NLAT_2; it++) {
			legendre_sphPlm_deriv_array(shtns, LMAX, im, ct[it], st[it], ylm[im] + it*lstride, dtylm);	// fixed im legendre functions lookup table.
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
void init_SH_gauss(shtns_cfg shtns, int on_the_fly)
{
	double t,tmax;
	long int it,im,m,l;
	long double iylm_fft_norm;
	long double xg[NLAT], wgl[NLAT];	// gauss points and weights.

	shtns->wg = malloc((NLAT_2 +15) * sizeof(double));	// gauss weights, double precision.

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
		shtns->ct[it] = xg[it];
		shtns->st[it] = sqrtl((1.-xg[it])*(1.+xg[it]));
		shtns->st_1[it] = 1.0/sqrtl((1.-xg[it])*(1.+xg[it]));
	}
	for (it=0; it<NLAT_2; it++)
		shtns->wg[it] = wgl[it]*iylm_fft_norm;		// faster double-precision computations.
	if (NLAT & 1) {		// odd NLAT : adjust weigth of middle point.
		shtns->wg[NLAT_2-1] *= 0.5;
	}
	for (it=NLAT_2; it < NLAT_2 +15; it++) shtns->wg[it] = 0.0;		// padding for multi-way algorithm.

#if SHT_VERBOSE > 1
	printf(" NLAT=%d, NLAT_2=%d\n",NLAT,NLAT_2);
// TEST if gauss points are ok.
	tmax = 0.0;
	for (it = 0; it<NLAT_2; it++) {
		t = legendre_Pl(NLAT, shtns->ct[it]);
		if (t>tmax) tmax = t;
//		printf("i=%d, x=%12.12g, p=%12.12g\n",it,ct[it],t);
	}
	printf("          max zero at Gauss nodes for Pl[l=NLAT] : %g\n",tmax);
	if (NLAT_2 < 100) {
		printf("          Gauss nodes :");
		for (it=0;it<NLAT_2; it++)
			printf(" %g",shtns->ct[it]);
		printf("\n");
	}
#endif

	if (on_the_fly != 0) return;

	init_SH_synth(shtns);

// for analysis (decomposition, direct transform) : transpose and multiply by gauss weight and other normalizations.
// interleave l and l+1 : this stores data in the way it will be read.
	double **ylm = shtns->ylm;
	double **zlm = shtns->zlm;
	struct DtDp** dylm = shtns->dylm;
	struct DtDp** dzlm = shtns->dzlm;
	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
		long int talign = 0;
		long int lstride = LMAX-m+1;
		if (im==0) lstride += (lstride&1);		// even stride for m=0.

		for (it=0;it<NLAT_2;it++) {
			double nz0, nz1, norm;
			norm = shtns->wg[it];
			if ( (m>0) && (shtns->norm & SHT_REAL_NORM) )	norm *= 2;		// "Real" norm : zlm must be doubled for m>0
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


void free_SH_dct(shtns_cfg shtns)
{
	if (shtns->zlm_dct0 == NULL) return;

	if (shtns->dzlm_dct0 != NULL)	fftw_free(shtns->dzlm_dct0);
	shtns->dzlm_dct0 = NULL;
	if (shtns->zlm_dct0 != NULL)	fftw_free(shtns->zlm_dct0);
	shtns->zlm_dct0 = NULL;

	if ((shtns->dykm_dct != NULL)&&(shtns->dykm_dct[0] != NULL))	fftw_free(shtns->dykm_dct[0]);
	shtns->dykm_dct[0] = NULL;
	if ((shtns->ykm_dct != NULL)&&(shtns->ykm_dct[0] != NULL))	fftw_free(shtns->ykm_dct[0]);
	shtns->ykm_dct[0] = NULL;

	if (shtns->idct != NULL)	fftw_destroy_plan(shtns->idct);	// free unused dct plans
	if (shtns->dct_m0 != NULL)	fftw_destroy_plan(shtns->dct_m0);
	if (shtns->dct_r1 != NULL)	fftw_destroy_plan(shtns->dct_r1);
	if (shtns->idct_r1 != NULL)	fftw_destroy_plan(shtns->idct_r1);
	shtns->idct = NULL;		shtns->dct_m0 = NULL;		shtns->idct_r1 = NULL;		shtns->dct_r1 = NULL;
}

#ifndef SHT_NO_DCT
/// Computes the matrices required for SH transform on a regular grid (with or without DCT).
/// \param analysis : 0 => synthesis only.
void init_SH_dct(shtns_cfg shtns, int analysis)
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
	
	double *ct = shtns->ct;
	double *st = shtns->st;
	double *st_1 = shtns->st_1;

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
	shtns->ykm_dct[0] = (double *) fftw_malloc(sizeof(double)* sk);
#ifndef SHT_SCALAR_ONLY
	shtns->dykm_dct[0] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* dsk);
#endif
	shtns->zlm_dct0 = (double *) fftw_malloc( sizeof(double)* it );
#ifndef SHT_SCALAR_ONLY
	shtns->dzlm_dct0 = (double *) fftw_malloc( sizeof(double)* im );
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
		shtns->ykm_dct[im+1] = shtns->ykm_dct[im] + sk;
#ifndef SHT_SCALAR_ONLY
		shtns->dykm_dct[im+1] = shtns->dykm_dct[im] + dsk;
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

	init_SH_synth(shtns);

	for (im=0; im<=MMAX; im++) {
		double* yl = shtns->ylm[im];
		struct DtDp* dyl = shtns->dylm[im];
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
					Z[it] = yl[it*lstride + (l-m)] * st[it];	// P[l+1](x)	*st
#ifndef SHT_SCALAR_ONLY
					dZt[it] = dyl[it*(LMAX-m+1) + (l-m)].t;	// P[l](x)	*1
					dZp[it] = dyl[it*(LMAX-m+1) + (l-m)].p;		// P[l-1](x)	*1
#endif
				}
			} else {	// m even
				for (it=0; it<NLAT_2; it++) {
					Z[it] = yl[it*lstride + (l-m)];		// P[l](x)	*1
#ifndef SHT_SCALAR_ONLY
					dZt[it] = dyl[it*(LMAX-m+1) + (l-m)].t *st[it];	// P[l+1](x)	*st
					dZp[it] = yl[it*lstride + (l-m)] * m;	// P[l](x)	*st
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
			if ( (m>0) && (shtns->norm & SHT_REAL_NORM) )	lnorm *= 2;		// "real" norm : zlm must be doubled for m>0

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
					shtns->zlm[im][it] =  Z[it];
				}
			} else if ((sk == 0)&&(l == LMAX)) {
				for (it=0; it<NLAT_2; it++) {
					shtns->zlm[im][(l-m)*NLAT_2 + it + talign] =  Z[it];
	#ifndef SHT_SCALAR_ONLY
					shtns->dzlm[im][(l-l0)*NLAT_2 + it].p = dZp[it];
					shtns->dzlm[im][(l-l0)*NLAT_2 + it].t = dZt[it];
	#endif
				}
			} else {
				for (it=0; it<NLAT_2; it++) {
					shtns->zlm[im][(l-m-sk)*NLAT_2 + it*2 +sk + talign] = Z[it];
	#ifndef SHT_SCALAR_ONLY
					shtns->dzlm[im][(l-l0-sk)*NLAT_2 + it*2 +sk].p = dZp[it];
					shtns->dzlm[im][(l-l0-sk)*NLAT_2 + it*2 +sk].t = dZt[it];
	#endif
				}
			}
		}
		}

		// Compact the coefficients for improved cache efficiency.
		yg = shtns->ykm_dct[im];
		for (it=0; it<= KMAX; it+=2) {
			l = (it < m) ? m : it-(m&1);
			while (l<=LMAX) {
				yg[0] = yk[(it/2)*(LMAX+1-m) + (l-m)];
				l++;	yg++;
			}
			if ((m==0) && ((LMAX & 1) == 0)) {	yg[0] = 0;		yg++;	}		// SSE2 padding.
		}
#ifndef SHT_SCALAR_ONLY		
		dyg = shtns->dykm_dct[im];
		for (it=0; it<= KMAX; it+=2) {
			l = (it-2 < m) ? m : it-2+(m&1);
			while (l<=LMAX) {
				dyg[0].t = dyk[(it/2)*(LMAX+1-m) + (l-m)].t;
				dyg[0].p = dyk[(it/2)*(LMAX+1-m) + (l-m)].p;
				l++;	dyg++;
			}
		}
		if (im == 0) {		// compact and reorder m=0 dylm because .p = 0 :
			dyg = shtns->dykm_dct[im];
			yg = (double *) shtns->dykm_dct[im];
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
		shtns->klim = klim;		// store for use in codelets.
		yg = shtns->zlm_dct0;
		for (l=0; l<=LMAX; l+=2) {
			for (it=l; it<klim; it++) {	// for m=0, zl coeff with i<l are zeros.
				*yg = yk0[it];
				yg++;
			}
			yk0 += 2*NLAT_2;
		}
	#ifndef SHT_SCALAR_ONLY
		yg = shtns->dzlm_dct0;
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
double SHT_error(shtns_cfg shtns)
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
	Sh = (double *) fftw_malloc( NSPAT_ALLOC(shtns) * sizeof(double) );
	Th = (double *) fftw_malloc( NSPAT_ALLOC(shtns) * sizeof(double) );

// m = nphi/2 is also real if nphi is even.
	nlm_cplx = ( MMAX*2 == NPHI ) ? LiM(shtns, MRES*MMAX,MMAX) : NLM;
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

	SH_to_spat(shtns, Slm0,Sh);		// scalar SHT
	spat_to_SH(shtns, Sh, Slm);
	for (i=0, tmax=0., n2=0., jj=0; i<NLM; i++) {		// compute error
		t = cabs(Slm[i] - Slm0[i]);
		n2 += t*t;
		if (t>tmax) { tmax = t; jj = i; }
	}
	err = tmax;
#if SHT_VERBOSE > 1
	printf("        scalar SH - poloidal   rms error = %.3g  max error = %.3g for l=%hu,lm=%d\n",sqrt(n2/NLM),tmax,shtns->li[jj],jj);
#endif

#ifndef SHT_SCALAR_ONLY
	Slm0[0] = 0.0; 	Tlm0[0] = 0.0;		// l=0, m=0 n'a pas de signification sph/tor
	SHsphtor_to_spat(shtns, Slm0, Tlm0, Sh, Th);		// vector SHT
	spat_to_SHsphtor(shtns, Sh, Th, Slm, Tlm);
	for (i=0, tmax=0., n2=0., jj=0; i<NLM; i++) {		// compute error
		t = cabs(Slm[i] - Slm0[i]);
		n2 += t*t;
		if (t>tmax) { tmax = t; jj = i; }
	}
	if (tmax > err) err = tmax;
#if SHT_VERBOSE > 1
	printf("        vector SH - spheroidal rms error = %.3g  max error = %.3g for l=%hu,lm=%d\n",sqrt(n2/NLM),tmax,shtns->li[jj],jj);
#endif
	for (i=0, tmax=0., n2=0., jj=0; i<NLM; i++) {		// compute error
		t = cabs(Tlm[i] - Tlm0[i]);
		n2 += t*t;
		if (t>tmax) { tmax = t; jj = i; }
	}
	if (tmax > err) err = tmax;
#if SHT_VERBOSE > 1
	printf("                  - toroidal   rms error = %.3g  max error = %.3g for l=%hu,lm=%d\n",sqrt(n2/NLM),tmax,shtns->li[jj],jj);
#endif
#endif
	return(err);		// return max error.
}


char* sht_name[SHT_NALG] = {"hyb", "fly1", "fly2", "fly3", "fly4", "fly6", "fly8" };
char* sht_var[SHT_NVAR] = {"std", "ltr"};
char *sht_type[SHT_NTYP] = {"syn", "ana", "vsy", "van", "gsy", "gs2", "v3s", "v3a" };
int sht_npar[SHT_NTYP] = {2, 2, 4, 4, 3, 3, 6, 6};

extern void* sht_array[SHT_NTYP][SHT_NALG];
extern void* sht_array_l[SHT_NTYP][SHT_NALG];

// big array holding all sht functions, variants and algorithms
void* sht_func[SHT_NVAR][SHT_NTYP][SHT_NALG];

/// \internal use on-the-fly alogorithm (good guess without measuring)
void set_sht_fly(shtns_cfg shtns)
{
  #define SET_SHT_FUNC(ivar) \
	shtns->fptr[ivar][SHT_TYP_SSY] = sht_func[ivar][SHT_TYP_SSY][SHT_FLY2]; \
	shtns->fptr[ivar][SHT_TYP_GS1] = sht_func[ivar][SHT_TYP_GS1][SHT_FLY2]; \
	shtns->fptr[ivar][SHT_TYP_GS2] = sht_func[ivar][SHT_TYP_GS2][SHT_FLY2]; \
	shtns->fptr[ivar][SHT_TYP_VSY] = sht_func[ivar][SHT_TYP_VSY][SHT_FLY1];	\
	shtns->fptr[ivar][SHT_TYP_3SY] = sht_func[ivar][SHT_TYP_3SY][SHT_FLY1]; \
	if (shtns->wg != NULL) { \
		shtns->fptr[ivar][SHT_TYP_SAN] = sht_func[ivar][SHT_TYP_SAN][SHT_FLY4]; \
		shtns->fptr[ivar][SHT_TYP_VAN] = sht_func[ivar][SHT_TYP_VAN][SHT_FLY2]; \
		shtns->fptr[ivar][SHT_TYP_3AN] = sht_func[ivar][SHT_TYP_3AN][SHT_FLY2]; \
	}

	SET_SHT_FUNC(SHT_STD)
	SET_SHT_FUNC(SHT_LTR)
  #undef SET_SHT_FUNC
}

/// \internal set hyb alogorithm and copy all algos to sht_func array.
void set_sht_default(shtns_cfg shtns)
{
	for (int it=0; it<SHT_NTYP; it++) {
		for (int j=0; j<SHT_NALG; j++) {	// copy variants to global array.
			sht_func[SHT_STD][it][j] = sht_array[it][j];
			sht_func[SHT_LTR][it][j] = sht_array_l[it][j];
		}
		shtns->fptr[SHT_STD][it] = sht_func[SHT_STD][it][SHT_HYB];
		shtns->fptr[SHT_LTR][it] = sht_func[SHT_LTR][it][SHT_HYB];
	}
}


#include "cycle.h"
#include <time.h>

#if SHT_VERBOSE > 1
  #define PRINT_VERB(msg) printf(msg)
#else
  #define PRINT_VERB(msg) (0)
#endif

#if SHT_VERBOSE == 1
  #define PRINT_DOT 	{	printf(".");	fflush(stdout);	}
#else
  #define PRINT_DOT (0);
#endif

double get_time(shtns_cfg shtns, int nloop, int npar, char* name, void *fptr, void *i1, void *i2, void *i3, void *o1, void *o2, void *o3, int l)
{
	int i;
	ticks tik0, tik1;

	for (i=0; i<=nloop; i++) {
		switch(npar) {
/*			case 2: (*(pf2l)fptr)(shtns, i1,o1, l); break;			// l may be discarded.
			case 3: (*(pf3l)fptr)(shtns, i1,o1,o2, l); break;
			case 4: (*(pf4l)fptr)(shtns, i1,i2,o1,o2, l); break;
			default: (*(pf6l)fptr)(shtns, i1,i2,i3, o1,o2,o3, l); break;
*/			case 2: (*(pf2)fptr)(shtns, i1,o1); break;			// l may be discarded.
			case 3: (*(pf3)fptr)(shtns, i1,o1,o2); break;
			case 4: (*(pf4)fptr)(shtns, i1,i2,o1,o2); break;
			default: (*(pf6)fptr)(shtns, i1,i2,i3, o1,o2,o3); break;
		}
		if (i==0) tik0 = getticks();
	}
	tik1 = getticks();
	double t = elapsed(tik1,tik0)/nloop;
	#if SHT_VERBOSE > 1
		printf("  t(%s) = %.3g",name,t);
	#endif
	return t;
}


/// \internal choose fastest between on-the-fly and gauss algorithms.
/// *nlp is the number of loops. If zero, it is set to a good value.
/// on_the_fly : 1 = skip all memory algorithm. 0 = include memory and on-the-fly. -1 = test only DCT.
double choose_best_sht(shtns_cfg shtns, int* nlp, int on_the_fly, int l)
{
	complex double *Qlm, *Slm, *Tlm;
	double *Qh, *Sh, *Th;
	char *nb;
	int m, i, i0, minc, nloop;
	int dct = 0;
	int analys = 1;		// check also analysis.
	double t0, t, tt, r;
	double tdct, tnodct;
	ticks tik0, tik1;
	clock_t tcpu;
	char ndct[20];

	if (NLAT < 32) return(0.0);		// on-the-fly not possible for NLAT_2 < 2*NWAY (overflow) and DCT not efficient for low NLAT.
	if (on_the_fly == -1) {
		on_the_fly = 0;		dct = 1;		// choose mtr_dct.
	}
	if (shtns->wg == NULL)	analys = 0;		// on-the-fly analysis not supported.

	m = 2*(NPHI/2+1) * NLAT * sizeof(double);
	i = sizeof(complex double)* NLM;
	if (i>m) m=i;
	Qh = (double *) fftw_malloc(m);		Sh = (double *) fftw_malloc(m);		Th = (double *) fftw_malloc(m);
	Qlm = (complex double *) fftw_malloc(m);	Slm = (complex double *) fftw_malloc(m);	Tlm = (complex double *) fftw_malloc(m);

	for (i=0;i<NLM;i++) {
		int l = shtns->li[i];
		Slm[i] = shtns->l_2[l] + 0.5*I*shtns->l_2[l];
		Tlm[i] = 0.5*shtns->l_2[l] + I*shtns->l_2[l];
	}

	if (*nlp <= 0) {
		// find good nloop by requiring less than 3% difference between 2 consecutive timings.
		m=0;	nloop = 1;                     // number of loops to get timings.
		r = 0.0;	tt = 1.0;
		do {
			if ((r > 0.03)||(tt<0.1)) {
				m = 0;		nloop *= 3;
			} else 	m++;
			tcpu = clock();
			t0 = get_time(shtns, nloop, 2, "", shtns->fptr[SHT_STD][SHT_TYP_SSY], Slm, Tlm, Qlm, Sh, Th, Qh, l);
			t = get_time(shtns, nloop, 2, "", shtns->fptr[SHT_STD][SHT_TYP_SSY], Slm, Tlm, Qlm, Sh, Th, Qh, l);
			tcpu = clock() - tcpu;
			r = fabs(2.0*(t-t0)/(t+t0));
			tt = 1.e-6 * tcpu;		// real time should not exceed 1 sec.
			#if SHT_VERBOSE > 1
				printf(", nloop=%d, t0=%g, t=%g, r=%g, m=%d (real time = %g s)\n",nloop,t0,t,r,m,tt);
			#endif
			PRINT_DOT
		} while((nloop<10000)&&(m < 3)&&(tt<0.35));
		*nlp = nloop;
	} else {
		nloop = *nlp;
	}
	#if SHT_VERBOSE > 1
		printf("nloop=%d\n",nloop);
	#endif
	
	int ityp = 0;	do {
		if ((dct != 0) && (ityp >= 4)) break;		// dct !=0 : only scalar and vector.
		#if SHT_VERBOSE > 1
			printf("finding best %s ...",sht_type[ityp]);	fflush(stdout);
		#endif
		if (ityp == 2) nloop = (nloop+1)/2;		// scalar ar done.
		t0 = 1e100;
		i0 = -1;		i = -1;
		if (on_the_fly == 1) i = SHT_FLY1 -1;		// only on-the-fly
		while (++i < SHT_NALG) {
			void *pf = sht_func[0][ityp][i];
			if (pf != NULL) {
				if (ityp&1) {	// analysis
					t = get_time(shtns, nloop, sht_npar[ityp], sht_name[i], pf, Sh, Th, Qh, Slm, Tlm, Qlm, l);
				} else {
					t = get_time(shtns, nloop, sht_npar[ityp], sht_name[i], pf, Slm, Tlm, Qlm, Sh, Th, Qh, l);
				}
				if (i==0) t *= 1.0/MIN_PERF_IMPROVE_DCT;
				if (t < t0) {	i0 = i;		t0 = t;		PRINT_VERB("*");	}
			}
		}
		if (i0 >= 0) {
			for (int j=0; j<SHT_NVAR; j++) {
				shtns->fptr[j][ityp] = sht_func[j][ityp][i0];
				if (ityp == 4) shtns->fptr[j][ityp+1] = sht_func[j][ityp+1][i0];		// only one timing for both gradients variants.
			}
			PRINT_DOT
			#if SHT_VERBOSE > 1
				printf(" => %s\n",sht_name[i0]);
			#endif
		}
		if (ityp == 4) ityp++;		// skip second gradient
		if ( ((ityp&1) == 0) && (analys == 0) ) ityp++;		// skip analysis
	} while(++ityp < SHT_NTYP);

#ifndef SHT_NO_DCT
	if (dct > 0) {		// find the best DCT timings...
        minc = MMAX/20 + 1;             // don't test every single m.
	#if SHT_VERBOSE > 1
		printf("\nfinding best dct synthesis ...");
	#endif
		m = -1;		i = -1;		// reference = no dct.
			t0 = get_time(shtns, nloop*2, 2, "s", shtns->fptr[SHT_STD][SHT_TYP_SSY], Qlm, Slm, Tlm, Qh, Sh, Th, l);
			t0 += get_time(shtns, nloop, 4, "v", shtns->fptr[SHT_STD][SHT_TYP_VSY], Slm, Tlm, Qlm, Sh, Th, Qh, l);
			tnodct = t0;
		for (m=0; m<=MMAX; m+=minc) {
			#if SHT_VERBOSE > 1
				printf("\nm=%d  ",m);
			#endif
			Set_MTR_DCT(shtns, m);
			t = get_time(shtns, nloop*2, 2, "sdct", sht_array[SHT_TYP_SSY][SHT_HYB], Qlm, Slm, Tlm, Qh, Sh, Th, l);
			t += get_time(shtns, nloop, 4, "vdct", sht_array[SHT_TYP_VSY][SHT_HYB], Slm, Tlm, Qlm, Sh, Th, Qh, l);
			if (t < t0) {	t0 = t;		i = m;	PRINT_VERB("*"); }
			PRINT_DOT
		}
		tdct = t0;
		Set_MTR_DCT(shtns, i);		// the best DCT is chosen.
	}
#endif

	#if SHT_VERBOSE > 1
		printf("\n");
	#endif
	fftw_free(Tlm);	fftw_free(Slm);	fftw_free(Qlm);	fftw_free(Th);	fftw_free(Sh);	fftw_free(Qh);

	if (dct > 0) {
		return(tdct/tnodct);
	} else	return(0.0);
}



#ifndef _HGID_
  #define _HGID_ "unknown"
#endif

void print_shtns_version() {
	printf("[SHTns] build " __DATE__ ", " __TIME__ ", id: " _HGID_ "\n");
}

void print_shtns_cfg(shtns_cfg shtns, int opt)
{
	printf("        Lmax=%d, Mmax*Mres=%d, Mres=%d, Nlm=%d,  Nphi=%d, Nlat=%d  [",LMAX, MMAX*MRES, MRES, NLM, NPHI, NLAT);
	if (shtns->norm & SHT_REAL_NORM) printf("'real' norm, ");
	if (shtns->norm & SHT_NO_CS_PHASE) printf("no Condon-Shortley phase, ");
	if (SHT_NORM == sht_fourpi) printf("4.pi normalized]\n");
	else if (SHT_NORM == sht_schmidt) printf("Schmidt semi-normalized]\n");
	else printf("orthonormalized]\n");

	if (opt == 0) return;

	printf("            ");
	for (int it=0; it<SHT_NTYP; it++)
		printf("%5s ",sht_type[it]);
	for (int iv=0; iv<SHT_NVAR; iv++) {
		printf("\n        %4s:",sht_var[iv]);
		for (int it=0; it<SHT_NTYP; it++) {
			for (int ia=0; ia<SHT_NALG; ia++) {
				if (sht_func[iv][it][ia] == shtns->fptr[iv][it]) {
					printf("%5s ",sht_name[ia]);
					break;
				}
			}
		}
	}
	printf("\n");
}

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
shtns_cfg shtns_create(int lmax, int mmax, int mres, enum shtns_norm norm)
{
	shtns_cfg shtns, s2;
	int im, m, l, lm;
	int with_cs_phase = 1;		/// Condon-Shortley phase (-1)^m is used by default.
	double mpos_renorm = 1.0;	/// renormalization of m>0.
	int larrays_ok = 0;
	int legendre_ok = 0;
	int l_2_ok = 0;

//	if (lmax < 1) shtns_runerr("lmax must be larger than 1");
	if (lmax < 2) shtns_runerr("lmax must be at least 2");
	if (sizeof(shtns->lmax) < sizeof(int)) {
		int llim = 1 << (8*sizeof(shtns->lmax)-1);
		if (lmax >= llim) shtns_runerr("lmax too large");
	}
#ifdef SHT_AXISYM
	if (mmax != 0) shtns_runerr("axisymmetric version : only Mmax=0 allowed");
#endif
	if (mmax*mres > lmax) shtns_runerr("MMAX*MRES should not exceed LMAX");
	if (mres <= 0) shtns_runerr("MRES must be > 0");

	// allocate new setup and initialize some variables (used as flags) :
	shtns = malloc( sizeof(struct shtns_info) + (mmax+1)*( sizeof(int)+sizeof(unsigned short) ) );
	if (shtns == NULL) return shtns;	// FAIL
	{
		void **p0 = (void**) &shtns->tm;	// first pointer in struct.
		void **p1 = (void**) &shtns->Y00_1;	// first non-pointer.
		while(p0 < p1)	 *p0++ = NULL;		// write NULL to every pointer.
		set_sht_default(shtns);		// default SHT is hyb.
		shtns->lmidx = (int*) (shtns + 1);		// lmidx is stored at the end of the struct...
		shtns->tm = (unsigned short*) (shtns->lmidx + (mmax+1));		// and tm just after.
		shtns->ct = NULL;	shtns->st = NULL;
	}

	// copy sizes.
	shtns->norm = norm;
	if (norm & SHT_NO_CS_PHASE)
		with_cs_phase = 0;
	if (norm & SHT_REAL_NORM)
		mpos_renorm = 0.5;		// normalization for 'real' spherical harmonics.

#ifdef SHT_AXISYM
	shtns->mmax = 0;		shtns->mres = 2;		shtns->nphi = 1;
#else
	MMAX = mmax;	MRES = mres;
#endif
	LMAX = lmax;
	NLM = nlm_calc(LMAX, MMAX, MRES);
#if SHT_VERBOSE > 0
	printf("[SHTns] build " __DATE__ ", " __TIME__ ", id: " _HGID_ "\n");
	printf("        Lmax=%d, Mmax*Mres=%d, Mres=%d, Nlm=%d  [",LMAX, MMAX*MRES, MRES, NLM);
	if (norm & SHT_REAL_NORM) printf("'real' norm, ");
	if (!with_cs_phase) printf("no Condon-Shortley phase, ");
	if (SHT_NORM == sht_fourpi) printf("4.pi normalized]\n");
	else if (SHT_NORM == sht_schmidt) printf("Schmidt semi-normalized]\n");
	else printf("orthonormalized]\n");
  #ifdef SHT_SCALAR_ONLY
	printf("  *** Compiled with SCALAR support only (no vector support) ***\n");
	#warning "Compilation with SCALAR support only."
  #endif
#endif

	s2 = sht_data;		// check if some data can be shared ...
	while(s2 != NULL) {
		if ( (s2->lmax == lmax) && (s2->mmax == mmax) && (s2->mres == mres) )
		{	// we can reuse the l-related arrays (li + copy lmidx)
			shtns->li = s2->li;
			for (im=0; im<=mmax; im++)	shtns->lmidx[im] = s2->lmidx[im];
			larrays_ok = 1;
		}
		if ( (s2->lmax >= lmax) && (s2->mmax >= mmax) && (s2->mres == mres) && (s2->norm == norm) )
		{	// we can reuse the legendre tables.
			shtns->al0 = s2->al0;		shtns->alm = s2->alm;
			shtns->bl0 = s2->bl0;		shtns->blm = s2->blm;
			legendre_ok = 1;
		}
		if (s2->lmax >= lmax) {		// we can reuse l_2
			shtns->l_2 = s2->l_2;
			l_2_ok = 1;
		}
		s2 = s2->next;
	}
	if (larrays_ok == 0) {
		// alloc spectral arrays
		shtns->li = (unsigned short *) malloc( NLM*sizeof(unsigned short) );	// NLM defined at runtime.
		for (im=0, lm=0; im<=MMAX; im++) {		// init l-related arrays.
			m = im*MRES;
			shtns->lmidx[im] = lm -m;		// virtual pointer for l=0
			for (l=im*MRES;l<=LMAX;l++) {
				shtns->li[lm] = l;		lm++;
			}
		}
		if (lm != NLM) shtns_runerr("unexpected error");
	}
	if (legendre_ok == 0) {	// this quickly precomputes some values for the legendre recursion.
		legendre_precomp(shtns, SHT_NORM, with_cs_phase, mpos_renorm);
	}
	if (l_2_ok == 0) {
		shtns->l_2 = (double *) malloc( (LMAX+1)*sizeof(double) );
		shtns->l_2[0] = 0.0;	// undefined for l=0 => replace with 0.	
		for (l=1; l<=LMAX; l++)		shtns->l_2[l] = 1.0/(l*(l+1.0));
	}

	switch(SHT_NORM) {
		case sht_schmidt:
			shtns->Y00_1 = 1.0;		shtns->Y10_ct = 1.0;
			break;
		case sht_fourpi:
			shtns->Y00_1 = 1.0;		shtns->Y10_ct = sqrt(1./3.);
			break;
		case sht_orthonormal:
		default:
			shtns->Y00_1 = sqrt(4.*M_PI);		shtns->Y10_ct = sqrt(4.*M_PI/3.);
//			Y11_st = sqrt(2.*M_PI/3.);		// orthonormal :  \f$ \sin\theta\cos\phi/(Y_1^1 + Y_1^{-1}) = -\sqrt{2 \pi /3} \f$
	}
	shtns->Y11_st = shtns->Y10_ct * sqrt(0.5/mpos_renorm);
	if (with_cs_phase)	shtns->Y11_st *= -1.0;		// correct Condon-Shortley phase

// save a pointer to this setup and return.
	shtns->next = sht_data;		// reference of previous setup (may be NULL).
	sht_data = shtns;			// keep track of new setup.
	return(shtns);
}


int free_unused_array(shtns_cfg shtns, void* p)
{
	int i = 0;		// reference count.
	shtns_cfg s2 = sht_data;
	if (p==NULL) return i;

	while (s2 != NULL) {		// we must not free shared resources.
		if (s2 != shtns) {		// don't count the one we want to free.
			if (s2->alm == p) i++;
			if (s2->blm == p) i++;
			if (s2->l_2 == p) i++;
			if (s2->li == p) i++;
		}
		s2 = s2->next;
	}
	if (i == 0)  free(p);
	return i;
}

/// release all resources allocated by a given shtns_cfg
void shtns_destroy(shtns_cfg shtns)
{
	free_unused_array(shtns, shtns->l_2);
	if (shtns->blm != shtns->alm)
		free_unused_array(shtns, shtns->blm);
	free_unused_array(shtns, shtns->alm);
	free_unused_array(shtns, shtns->li);

	if (shtns->wg != NULL) free(shtns->wg);
	shtns->wg = NULL;
	free_SH_dct(shtns);
	free_SHTarrays(shtns);

	if (sht_data == shtns) {
		sht_data = shtns->next;		// forget shtns
	} else {
		shtns_cfg s2 = sht_data;
		while (s2 != NULL) {
			if (s2->next == shtns) {
				s2->next = shtns->next;		// forget shtns
				break;
			}
			s2 = s2->next;
		}
	}
	shtns->nlm = 0;		shtns->lmax = 0;	shtns->mmax = 0;		// security : mark as 0.
	free(shtns);
}

// clear all allocated memory (hopefully) and go back to 0 state.
void shtns_reset()
{
	while (sht_data != NULL) {
		shtns_destroy(sht_data);
	}
}

/*! Initialization of Spherical Harmonic transforms (backward and forward, vector and scalar, ...) of given size.
 * <b>This function must be called after \ref shtns_create and before any SH transform.</b> and sets all global variables and internal data.
 * returns the required number of doubles to be allocated for a spatial field.
 * \param nlat,nphi pointers to the number of latitudinal and longitudinal grid points respectively. If 0, they are set to optimal values.
 * \param nl_order defines the maximum SH degree to be resolved by analysis : lmax_analysis = lmax*nl_order. It is used to set an optimal and anti-aliasing nlat. If 0, the default SHT_DEFAULT_NL_ORDER is used.
 * \param flags allows to choose the type of transform (see \ref shtns_type) and the spatial data layout (see \ref spat)
 * \param eps polar optimization threshold : polar values of Legendre Polynomials below that threshold are neglected (for high m), leading to increased performance (a few percents)
 *  0 = no polar optimization;  1.e-14 = VERY safe;  1.e-10 = safe;  1.e-6 = aggresive, but still good accuracy.
*/
int shtns_set_grid_auto(shtns_cfg shtns, enum shtns_type flags, double eps, int nl_order, int *nlat, int *nphi)
{
	double t;
	int im,m;
	int layout, nspat;
	int nloop = 0;
	int n_gauss = 0;
	int on_the_fly = 0;
	int quick_init = 0;

	if (nl_order <= 0) nl_order = SHT_DEFAULT_NL_ORDER;
/*	shtns.lshift = 0;
	if (nl_order == 0) nl_order = SHT_DEFAULT_NL_ORDER;
	if (nl_order < 0) {	shtns.lshift = -nl_order;	nl_order = 1; }		// linear with a shift in l.
*/
	shtns->nlorder = nl_order;
	shtns->mtr_dct = -1;
	layout = flags & 0xFFFF00;
	flags = flags & 255;	// clear higher bits.

	shtns->fftw_plan_mode = FFTW_EXHAUSTIVE;		// defines the default FFTW planner mode.
	switch (flags) {
	  #ifdef SHT_NO_DCT
		case sht_auto :				flags = sht_gauss;		// auto means gauss if dct is disabled.
			break;
	  #endif
		case sht_gauss_fly :		flags = sht_gauss;		on_the_fly = 1;
			break;
		case sht_quick_init :	 	flags = sht_gauss;
		case sht_reg_poles :		quick_init = 1;		shtns->fftw_plan_mode = FFTW_ESTIMATE;		// quick fftw init.
			break;
		default : break;
	}

	if (*nphi == 0) {
		*nphi = fft_int((nl_order+1)*MMAX+1, 7);		// required fft nodes
	}
	if (*nlat == 0) {
		n_gauss = ((nl_order+1)*LMAX)/2 +1;		// required gauss nodes
		n_gauss += (n_gauss&1);		// even is better.
		if ((flags == sht_auto)||(flags == sht_reg_fast)) {
			m = fft_int(nl_order*LMAX+2, 7);		// required dct nodes
			*nlat = m + (m&1);		// even is better.
		} else *nlat = n_gauss;
	}

	t = sht_mem_size(shtns->lmax, shtns->mmax, shtns->mres, *nlat);
	if ( (t > SHTNS_MAX_MEMORY) && (on_the_fly == 0) ) {		// huge transform has been requested
		on_the_fly = 1;
		if ( (flags == sht_reg_dct) || (flags == sht_reg_fast) ) shtns_runerr("Memory limit exceeded, try using sht_gauss or increase SHTNS_MAX_MEMORY in sht_config.h");
		if (flags != sht_reg_poles) {
			flags = sht_gauss;
			if (n_gauss > 0) *nlat = n_gauss;
		}
		if (quick_init == 0) {
			if (*nphi > 256) shtns->fftw_plan_mode = FFTW_PATIENT;		// do not waste too much time finding optimal fftw.
			if (*nphi > 512) shtns->fftw_plan_mode = FFTW_MEASURE;
		}
		if (t > 10*SHTNS_MAX_MEMORY) quick_init =1;			// do not time such large transforms.
	}

	if (flags == sht_auto) {
		if ( ((nl_order>=2)&&(MMAX*MRES > LMAX/2)) || (n_gauss <= 32) ) {
			flags = sht_gauss;		// avoid computing DCT stuff when it is not expected to be faster.
			if (n_gauss > 0) *nlat = n_gauss;
		}
	}

	// copy to global variables.
#ifdef SHT_AXISYM
	shtns->nphi = 1;
	if (*nphi != 1) shtns_runerr("axisymmetric version : only Nphi=1 allowed");
#else
	NPHI = *nphi;
#endif
	NLAT_2 = (*nlat+1)/2;	NLAT = *nlat;
	if ((NLAT_2)*2 <= LMAX) shtns_runerr("NLAT_2*2 should be at least LMAX+1");

	alloc_SHTarrays(shtns, on_the_fly);		// allocate dynamic arrays
	nspat = planFFT(shtns, layout);		// initialize fftw
	shtns->zlm_dct0 = NULL;		// used as a flag.

	if (flags == sht_reg_dct) {		// pure dct.
		#ifdef SHT_NO_DCT
			runerr("DCT not supported. recompile without setting SHT_NO_DCT in sht_config.h");
		#endif
		init_SH_dct(shtns, 1);
		OptimizeMatrices(shtns, 0.0);
		Set_MTR_DCT(shtns, MMAX);
	}
	if ((flags == sht_auto)||(flags == sht_reg_fast))
	{
		init_SH_dct(shtns, 1);
		OptimizeMatrices(shtns, eps);
  #ifdef SHT_NO_DCT
			Set_MTR_DCT(shtns, -1);		// turn off DCT.
  #else
	#if SHT_VERBOSE > 0
		printf("        finding optimal MTR_DCT");	fflush(stdout);
	#endif
		t = choose_best_sht(shtns, &nloop, -1, 0);		// find optimal MTR_DCT.
	#if SHT_VERBOSE > 0
		printf("\n        + optimal MTR_DCT = %d  (%.1f%% performance gain)\n", MTR_DCT*MRES, 100.*(1/t-1));
	#endif
  		if ((n_gauss > 0)&&(flags == sht_auto)) t *= ((double) NLAT)/n_gauss;	// we can revert to gauss with a smaller nlat.
		if (t > MIN_PERF_IMPROVE_DCT) {
			Set_MTR_DCT(shtns, -1);		// turn off DCT.
		} else {
			t = SHT_error(shtns);
			if (t > MIN_ACCURACY_DCT) {
	#if SHT_VERBOSE > 0
				printf("     !! Not enough accuracy (%.3g) => turning off DCT.\n",t);
	#endif
				Set_MTR_DCT(shtns, -1);		// turn off DCT.
			}
		}
  #endif
		if (MTR_DCT == -1) {			// free memory used by DCT and disables DCT.
			free_SH_dct(shtns);			// free now useless arrays.
			if (flags == sht_auto) {
				flags = sht_gauss;		// switch to gauss grid, even better accuracy.
		#if SHT_VERBOSE > 0
				printf("        => switching back to Gauss Grid\n");
		#endif
				for (im=1; im<=MMAX; im++) {	//	im >= 1
					m = im*MRES;
					shtns->ylm[im]  -= shtns->tm[im]*(LMAX-m+1);		// restore pointers altered by OptimizeMatrices().
		#ifndef SHT_SCALAR_ONLY
					shtns->dylm[im] -= shtns->tm[im]*(LMAX-m+1);
		#endif
				}
				if (n_gauss > 0) {		// we should use the optimal size for gauss-legendre
					free_SHTarrays(shtns);
					*nlat = n_gauss;
					NLAT_2 = (*nlat+1)/2;	NLAT = *nlat;
					alloc_SHTarrays(shtns, on_the_fly);
					nspat = planFFT(shtns, layout);		// fft must be replanned because NLAT has changed.
				}
			}
		}
	}
	if (flags == sht_gauss)
	{
		MTR_DCT = -1;		// we do not use DCT !!!
		init_SH_gauss(shtns, on_the_fly);
		if (on_the_fly == 0) OptimizeMatrices(shtns, eps);
	}
	if (flags == sht_reg_poles)
	{
		MTR_DCT = -1;		// we do not use DCT !!!
	  #ifndef SHT_SCALAR_ONLY
		fftw_free(shtns->dzlm[0]);		shtns->dzlm[0] = NULL;
	  #endif
		fftw_free(shtns->zlm[0]);		shtns->zlm[0] = NULL;		// no inverse transform, mark as unused.
		EqualPolarGrid(shtns);
		if (on_the_fly == 0) {
			init_SH_synth(shtns);
			OptimizeMatrices(shtns, 0.0);
		}
	}
	
	if (on_the_fly == 1) {
  #if SHT_VERBOSE > 0
		printf("        + using on-the-fly transforms.\n");
  #endif
		if (NLAT < 32) shtns_runerr("on-the-fly only available for nlat>=32");		// avoid overflow with NLAT_2 < 2*NWAY
		PolarOptimize(shtns, eps);
		set_sht_fly(shtns);		// switch function pointers to "on-the-fly" functions.
	}

	if (quick_init == 0) {
		choose_best_sht(shtns, &nloop, on_the_fly, 2*LMAX/3);
		if (MMAX == 0) {		// use SHT_AXISYM version
/*			SH_to_spat_ptr = SH_to_spat_ptr_m0;		SHsphtor_to_spat_ptr = SHsphtor_to_spat_ptr_m0;		SHqst_to_spat_ptr = SHqst_to_spat_ptr_m0;
			spat_to_SH_ptr = spat_to_SH_ptr_m0;		spat_to_SHsphtor_ptr = spat_to_SHsphtor_ptr_m0;		spat_to_SHqst_ptr = spat_to_SHqst_ptr_m0;
			SH_to_spat_ptr_l = SH_to_spat_ptr_m0l;		SHsphtor_to_spat_ptr_l = SHsphtor_to_spat_ptr_m0l;		SHqst_to_spat_ptr_l = SHqst_to_spat_ptr_m0l;
			spat_to_SH_ptr_l = spat_to_SH_ptr_m0l;		spat_to_SHsphtor_ptr_l = spat_to_SHsphtor_ptr_m0l;		spat_to_SHqst_ptr_l = spat_to_SHqst_ptr_m0l;
*/		}
  #if SHT_VERBOSE > 0
		printf("\n");
  #endif
		t = SHT_error(shtns);		// compute SHT accuracy.
  #if SHT_VERBOSE > 0
		printf("        + SHT accuracy = %.3g\n",t);
  #endif
  #if SHT_VERBOSE < 2
		if (t > 1.e-3) shtns_runerr("bad SHT accuracy");		// stop if something went wrong (but not in debug mode)
  #endif
	}
  #if SHT_VERBOSE > 0
	printf("        => SHTns is ready.\n");
  #endif
	return(nspat);	// returns the number of doubles to be allocated for a spatial field.
}


/*! Initialization of Spherical Harmonic transforms (backward and forward, vector and scalar, ...) of given size.
 * <b>This function must be called after \ref shtns_create and before any SH transform.</b> and sets all global variables.
 * returns the required number of doubles to be allocated for a spatial field.
 * \param nlat,nphi respectively the number of latitudinal and longitudinal grid points.
 * \param flags allows to choose the type of transform (see \ref shtns_type) and the spatial data layout (see \ref spat)
 * \param eps polar optimization threshold : polar values of Legendre Polynomials below that threshold are neglected (for high m), leading to increased performance (a few percents)
 *  0 = no polar optimization;  1.e-14 = VERY safe;  1.e-10 = safe;  1.e-6 = aggresive, but still good accuracy.
*/
int shtns_set_grid(shtns_cfg shtns, enum shtns_type flags, double eps, int nlat, int nphi)
{
	if ((nlat == 0)||(nphi == 0)) shtns_runerr("nlat or nphi is zero !");
	return( shtns_set_grid_auto(shtns, flags, eps, 0, &nlat, &nphi) );
}


/*! Simple initialization of Spherical Harmonic transforms (backward and forward, vector and scalar, ...) of given size.
 * This function sets all global variables by calling \ref shtns_create followed by \ref shtns_set_grid, with the
 * default normalization and the default polar optimization (see \ref sht_config.h).
 * returns the number of modes to describe a scalar field.
 * \param lmax : maximum SH degree that we want to describe.
 * \param mmax : number of azimutal wave numbers.
 * \param mres : \c 2.pi/mres is the azimutal periodicity. \c mmax*mres is the maximum SH order.
 * \param nlat,nphi : respectively the number of latitudinal and longitudinal grid points.
 * \param flags allows to choose the type of transform (see \ref shtns_type) and the spatial data layout (see \ref spat)
*/
shtns_cfg shtns_init(enum shtns_type flags, int lmax, int mmax, int mres, int nlat, int nphi)
{
	shtns_cfg shtns = shtns_create(lmax, mmax, mres, SHT_DEFAULT_NORM);
	if (shtns != NULL)
		shtns_set_grid(shtns, flags, SHT_DEFAULT_POLAR_OPT, nlat, nphi);
	return shtns;
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
