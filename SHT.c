///////////////////////////////////////////////
// SHT : Spherical Harmonic Transform
//   requires SHT.h for size parameters.
//////////////////////////////////////////////

#include <complex.h>
#include <math.h>
// FFTW la derivee d/dx = ik	(pas de moins !)
#include <fftw3.h>
// GSL for Legendre functions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>


// SHT.h : parameter for SHT (sizes : LMAX, NLAT, MMAX, MRES, NPHI)
#include "SHT.h"
// NLM : total number of Spherical Harmonic coefficients.
//#define NLM ( ((MMAX+1)*(2*LMAX+2-MMAX))/2 )
#define NLM ( (MMAX+1)*(LMAX+1) - MRES* MMAX*(MMAX+1)/2 )
// LM(l,m) : index in the Spherical Harmonic coefficient array [ (l,m) space ]
//#define LM(l,m) ((m*(2*LMAX +1-m))/2 + l)
#define LM(l,m) ( (m/MRES)*(2*LMAX+3 -m)/2 + l-m )
#define LiM(l,im) ( im*(2*LMAX+3 -(im+2)*MRES)/2 + l )

struct DtDp {		// theta and phi derivatives stored together.
	double t, p;
};

/*
union cplx {		// to access easily the real part of a complex number.
	complex double c;
	double r;
};
*/

double pi = atan(1.0)*4.0;
double l2[NLM], l_2[NLM];			// l(l+1) and 1/(l(l+1))
double ct[NLAT/2], st[NLAT/2], st_1[NLAT/2];	// cos(theta), sin(theta), 1/sin(theta);

long int tm[MMAX+1];	// start theta value for SH (polar optimization : near the poles the legendre polynomials go to zero for high m's)

double* ylm[MMAX+1];		// matrix for direct transform
struct DtDp* dylm[MMAX+1];	// theta and phi derivative of Ylm matrix
double* iylm[MMAX+1];		// matrix for inverse transform
struct DtDp* idylm[MMAX+1];

fftw_plan ifft, fft;	// plans for FFTW.
unsigned fftw_plan_mode = FFTW_PATIENT;		// defines the default FFTW planner mode.


// compute non-linear terms in m-spectral space and l-physical space, with FFTW conventions, ie :
// y = x0 + sum(m=1..MMAX) [ xm.exp(i.m.phi) + xm*.exp(-i.m.phi) ]
// this should be faster for small m's, and requires a legendre transform before.
// No aliasing problems.
void NLspec(complex double *x, complex double *y, complex double *nl)
{
	long int m,i,l;

	m=0;	// 2MMAX+1 terms
		for (l=0;l<NLAT;l++)
			nl[m*NLAT+l] = creal(x[0*NLAT+l]) * creal(y[0*NLAT+l]);	// 1 term
		for (i=1; i<=MMAX; i++) {		// 2(MMAX)
			for (l=0;l<NLAT;l++)
//				(double) nl[m*NLAT+l] += 2* creal(x[i*NLAT+l]*conj(y[i*NLAT+l]));	//xi*conj(yi) + conj(xi)*yi
				(double) nl[m*NLAT+l] += 2*( creal(x[i*NLAT+l])*creal(y[i*NLAT+l]) + cimag(x[i*NLAT+l])*cimag(y[i*NLAT+l]) );
		}
	for (m=1; m<=MMAX; m++) {	// total 2(MMAX)+1-m terms
		for (l=0;l<NLAT;l++)
			nl[m*NLAT+l] = creal(x[0*NLAT+l])*y[m*NLAT+l] + x[m*NLAT+l]*creal(y[0*NLAT+l]);	// 2 terms
		for (i=1; i<m; i++) {		// 3(m-1) terms
			for (l=0;l<NLAT;l++)
				nl[m*NLAT+l] += conj(x[i*NLAT+l])*y[(i+m)*NLAT+l] + x[(i+m)*NLAT+l]*conj(y[i*NLAT+l]) + x[i*NLAT+l]*y[(m-i)*NLAT+l];
		}
		for (i=m; i<=MMAX-m; i++) {	// 2(MMAX-2m+1) terms
			for (l=0;l<NLAT;l++)
				nl[m*NLAT+l] += conj(x[i*NLAT+l])*y[(i+m)*NLAT+l] + x[(i+m)*NLAT+l]*conj(y[i*NLAT+l]);
		}
	}
}

/////////////////////////////////////////////////////
//   Scalar Spherical Harmonics Transform
// input  : ShF = spatial/fourrier data : complex double array of size NLAT*(NPHI/2+1) or double array of size NLAT*(NPHI/2+1)*2
// output : Slm = spherical harmonics coefficients : complex double array of size NLM
void spat_to_SH(complex double *ShF, complex double *Slm)
{
	complex double fp[NLAT/2];	// symmetric (even) part
	complex double fm[NLAT/2];	// anti-symmetric (odd) part
	complex double *Sl;		// virtual pointers for given im
	double *iyl;
	long int i,im,m,l;

	fftw_execute_dft_r2c(fft,(double *) ShF, ShF);

	im = 0;
		m=im*MRES;
		for (i=0;i<NLAT/2;i++) {	// compute symmetric and antisymmetric parts. m=0 : everything is REAL
			(double) fp[i] = (double) ShF[i] + (double) ShF[NLAT-(i+1)];
			(double) fm[i] = (double) ShF[i] - (double) ShF[NLAT-(i+1)];
		}
		l=m;
		Sl = &Slm[LiM(0,im)];		// virtual pointer for l=0 and im
		iyl = iylm[im];
		ShF += NLAT;
		while (l<LMAX) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
			Sl[l] = 0.0;	Sl[l+1] = 0.0;		// Slm[LiM(l,im)] = 0.0;	Slm[LiM(l+1,im)] = 0.0;
			for (i=0;i<NLAT/2;i++) {
				(double) Sl[l] += (double) fp[i] * iyl[i];		// Slm[LiM(l,im)] += iylm[im][(l-m)*NLAT/2 + i] * fp[i];
				(double) Sl[l+1] += (double) fm[i] * iyl[NLAT/2 + i];	// Slm[LiM(l+1,im)] += iylm[im][(l+1-m)*NLAT/2 + i] * fm[i];
			}
			l+=2;
			iyl += NLAT;
		}
		if (l==LMAX) {
			Sl[l] = 0.0;	// Slm[LiM(l,im)] = 0.0;
			for (i=tm[im];i<NLAT/2;i++) {	// polar optimization
				(double) Sl[l] += iyl[i] * (double) fp[i];	// Slm[LiM(l,im)] += iylm[im][(l-m)*NLAT/2 + i] * fp[i];
			}
		}
	for (im=1;im<=MMAX;im++) {
		m=im*MRES;
		for (i=tm[im];i<NLAT/2;i++) {	// compute symmetric and antisymmetric parts.
			fp[i] = ShF[i] + ShF[NLAT-(i+1)];
			fm[i] = ShF[i] - ShF[NLAT-(i+1)];
		}
		l=m;
		Sl = &Slm[LiM(0,im)];		// virtual pointer for l=0 and im
		iyl = iylm[im];
		ShF += NLAT;
		while (l<LMAX) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
			Sl[l] = 0.0;	Sl[l+1] = 0.0;		// Slm[LiM(l,im)] = 0.0;	Slm[LiM(l+1,im)] = 0.0;
			for (i=tm[im];i<NLAT/2;i++) {	// tm[im] : polar optimization
				Sl[l] += fp[i] * iyl[i];		// Slm[LiM(l,im)] += iylm[im][(l-m)*NLAT/2 + i] * fp[i];
				Sl[l+1] += fm[i] * iyl[NLAT/2 + i];	// Slm[LiM(l+1,im)] += iylm[im][(l+1-m)*NLAT/2 + i] * fm[i];
			}
			l+=2;
			iyl += NLAT;
		}
		if (l==LMAX) {
			Sl[l] = 0.0;	// Slm[LiM(l,im)] = 0.0;
			for (i=tm[im];i<NLAT/2;i++) {	// polar optimization
				Sl[l] += iyl[i] * fp[i];	// Slm[LiM(l,im)] += iylm[im][(l-m)*NLAT/2 + i] * fp[i];
			}
		}
	}
}

/////////////////////////////////////////////////////
//   Scalar inverse Spherical Harmonics Transform
// input  : Slm = spherical harmonics coefficients : complex double array of size NLM [unmodified]
// output : ShF = spatial/fourrier data : complex double array of size NLAT*(NPHI/2+1) or double array of size NLAT*(NPHI/2+1)*2
void SH_to_spat(complex double *Slm, complex double *ShF)
{
	complex double fe, fo;		// even and odd parts
	complex double *Sl;
	double *yl;
	long int i,im,m,l;

	im = 0;
		m = im*MRES;
		Sl = &Slm[LiM(0,im)];	// virtual pointer for l=0 and im
		i=0;
		yl = ylm[im] + i*(LMAX-m+1) -m;
		while (i<NLAT/2) {	// ops : NLAT/2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
			fe = 0.0;	fo = 0.0;
			while (l<LMAX) {	// compute even and odd parts
				(double) fe += yl[l] * (double) Sl[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Slm[LiM(l,im)];
				(double) fo += yl[l+1] * (double) Sl[l+1];	// fo += ylm[im][i*(LMAX-m+1) + (l+1-m)] * Slm[LiM(l+1,im)];
				l+=2;
			}
			if (l==LMAX) {
				(double) fe += yl[l] * (double) Sl[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Slm[LiM(l,im)];
			}
			ShF[i] = fe + fo;
			ShF[NLAT-(i+1)] = fe - fo;
			i++;
			yl += (LMAX-m+1);
		}
		ShF += NLAT;
	for (im=1; im<=MMAX; im++) {
		m = im*MRES;
		Sl = &Slm[LiM(0,im)];	// virtual pointer for l=0 and im
		i=0;
		while (i<tm[im]) {	// polar optimization
			ShF[i] = 0.0;
			ShF[NLAT-tm[im] + i] = 0.0;	// south pole zeroes <=> ShF[im*NLAT + NLAT-(i+1)] = 0.0;
			i++;
		}
		yl = ylm[im] + i*(LMAX-m+1) -m;
		while (i<NLAT/2) {	// ops : NLAT/2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
			fe = 0.0;	fo = 0.0;
			while (l<LMAX) {	// compute even and odd parts
				fe += yl[l] * Sl[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Slm[LiM(l,im)];
				fo += yl[l+1] * Sl[l+1];	// fo += ylm[im][i*(LMAX-m+1) + (l+1-m)] * Slm[LiM(l+1,im)];
				l+=2;
			}
			if (l==LMAX) {
				fe += yl[l] * Sl[l];		// fe += ylm[im][i*(LMAX-m+1) + (l-m)] * Slm[LiM(l,im)];
			}
			ShF[i] = fe + fo;
			ShF[NLAT-(i+1)] = fe - fo;
			i++;
			yl += (LMAX-m+1);
		}
		ShF += NLAT;
	}
	for(im=MMAX+1; im<=NPHI/2; im++) {	// padding for high m's
		for (i=0;i<NLAT;i++)
			ShF[i] = 0.0;
		ShF += NLAT;
	}

	ShF -= NLAT*(NPHI/2+1);		// restore original pointer
	fftw_execute_dft_c2r(ifft, ShF, (double *) ShF);
}

void SH_to_grad_spat(complex double *Slm, complex double *GtF, complex double *GpF)
{
	complex double gte, gto, gpe, gpo;		// even and odd parts
	complex double *Sl;
	struct DtDp *dyl;
	long int i,im,m,l;

	im=0;		// zonal part : Gp = 0.0
		m = im*MRES;
		Sl = &Slm[LiM(0,im)];	// virtual pointer for l=0 and im
		i=0;
		dyl = dylm[im] + i*(LMAX-m+1) -m;
		while (i<NLAT/2) {	// ops : NLAT/2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
			gte = 0.0;  gto = 0.0;
			while (l<LMAX) {	// compute even and odd parts
				(double) gto += dyl[l].t * (double) Sl[l];	// m=0 : everything is REAL
				(double) gte += dyl[l+1].t * (double) Sl[l+1];
				l+=2;
			}
			if (l==LMAX) {
				(double) gto += dyl[l].t * Sl[l];
			}
			GtF[i] = (gte + gto);
			GtF[NLAT-(i+1)] = (gte - gto);
			GpF[i] = 0.0;
			GpF[NLAT-(i+1)] = 0.0;
			i++;
			dyl += (LMAX-m+1);
		}
		GtF += NLAT;	GpF += NLAT;
	for (im=1; im<=MMAX; im++) {
		m = im*MRES;
		Sl = &Slm[LiM(0,im)];	// virtual pointer for l=0 and im
		i=0;
		while (i<tm[im]) {	// polar optimization
			GtF[i] = 0.0;
			GtF[NLAT-tm[im] + i] = 0.0;	// south pole zeroes <=> ShF[im*NLAT + NLAT-(i+1)] = 0.0;
			GpF[i] = 0.0;
			GpF[NLAT-tm[im] + i] = 0.0;
			i++;
		}
		dyl = dylm[im] + i*(LMAX-m+1) -m;
		while (i<NLAT/2) {	// ops : NLAT/2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
			gte = 0.0;  gto = 0.0;  gpe = 0.0;  gpo = 0.0;
			while (l<LMAX) {	// compute even and odd parts
				gto += dyl[l].t * Sl[l];
				gpe += dyl[l].p * Sl[l];
				gte += dyl[l+1].t * Sl[l+1];
				gpo += dyl[l+1].p * Sl[l+1];
				l+=2;
			}
			if (l==LMAX) {
				gto += dyl[l].t * Sl[l];
				gpe += dyl[l].p * Sl[l];
			}
			GtF[i] = gte + gto;
			GtF[NLAT-(i+1)] = gte - gto;
			GpF[i] = I*(gpe + gpo);
			GpF[NLAT-(i+1)] = I*(gpe - gpo);
			i++;
			dyl += (LMAX-m+1);
		}
		GtF += NLAT;	GpF += NLAT;
	}
	for(im=MMAX+1; im<=NPHI/2; im++) {	// padding for high m's
		for (i=0;i<NLAT;i++) {
			GtF[i] = 0.0;	GpF[i] = 0.0;
		}
		GtF += NLAT;	GpF[i] += NLAT;
	}

	GtF -= NLAT*(NPHI/2+1);		// restore original pointer
	GpF -= NLAT*(NPHI/2+1);
	fftw_execute_dft_c2r(ifft, GtF, (double *) GtF);
	fftw_execute_dft_c2r(ifft, GpF, (double *) GpF);
}


/////////////////////////////////////////////////////
//   Spheroidal/Toroidal to (theta,phi) components inverse Spherical Harmonics Transform
// input  : Slm,Tlm = spherical harmonics coefficients of Spheroidal and Toroidal scalars : 
//          complex double array of size NLM [unmodified]
// output : BtF, BpF = theta, and phi vector components, spatial/fourrier data : 
//          complex double array of size NLAT*(NPHI/2+1) or double array of size NLAT*(NPHI/2+1)*2
void SHsphtor_to_spat(complex double *Slm, complex double *Tlm, complex double *BtF, complex double *BpF)
{
	complex double se, so, dse, dso;	// spheroidal even and odd parts
	complex double te, to, dte, dto;	// toroidal ...
	complex double *Sl, *Tl;
	struct DtDp *dyl;
	long int i,im,m,l;

	im = 0;		// zonal part : d/dphi = 0;
		m = im*MRES;
		Sl = &Slm[LiM(0,im)];	// virtual pointer for l=0 and im
		Tl = &Tlm[LiM(0,im)];
		i=0;
		dyl = dylm[im] + i*(LMAX-m+1) -m;
		while (i<NLAT/2) {	// ops : NLAT/2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
			dse = 0.0;	dso = 0.0;	dte = 0.0;	dto = 0.0;
			while (l<LMAX) {	// compute even and odd parts
				(double) dto += dyl[l].t * (double) Tl[l];	// m=0 : everything is real.
				(double) dso += dyl[l].t * (double) Sl[l];
				(double) dte += dyl[l+1].t * (double) Tl[l+1];
				(double) dse += dyl[l+1].t * (double) Sl[l+1];
				l+=2;
			}
			if (l==LMAX) {
				(double) dto += dyl[l].t * Tl[l];
				(double) dso += dyl[l].t * Sl[l];
			}
			BtF[i] = (dse+dso);			// Bt = dS/dt
			BtF[NLAT-(i+1)] = (dse-dso);
			BpF[i] = - (dte+dto);			// Bp = - dT/dt
			BpF[NLAT-(i+1)] = - (dte-dto);
			i++;
			dyl += (LMAX-m+1);
		}
		BpF += NLAT;	BtF += NLAT;
	for (im=1; im<=MMAX; im++) {
		m = im*MRES;
		Sl = &Slm[LiM(0,im)];	// virtual pointer for l=0 and im
		Tl = &Tlm[LiM(0,im)];
		i=0;
		while (i<tm[im]) {	// polar optimization
			BtF[i] = 0.0;
			BtF[NLAT-tm[im] + i] = 0.0;	// south pole zeroes
			BpF[i] = 0.0;
			BpF[NLAT-tm[im] + i] = 0.0;	// south pole zeroes
			i++;
		}
		dyl = dylm[im] + i*(LMAX-m+1) -m;
		while (i<NLAT/2) {	// ops : NLAT/2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
			dse = 0.0;	dso = 0.0;	dte = 0.0;	dto = 0.0;
			se = 0.0;	so = 0.0;	te = 0.0;	to = 0.0;
			while (l<LMAX) {	// compute even and odd parts
				dto += dyl[l].t * Tl[l];
				te += dyl[l].p * Tl[l];
				dso += dyl[l].t * Sl[l];
				se += dyl[l].p * Sl[l];
				dte += dyl[l+1].t * Tl[l+1];
				to += dyl[l+1].p * Tl[l+1];
				dse += dyl[l+1].t * Sl[l+1];
				so += dyl[l+1].p * Sl[l+1];
				l+=2;
			}
			if (l==LMAX) {
				dto += dyl[l].t * Tl[l];
				te += dyl[l].p * Tl[l];
				dso += dyl[l].t * Sl[l];
				se += dyl[l].p * Sl[l];
			}
			BtF[i] = (dse+dso) + I*(te+to);			// Bt = dS/dt       + I.m/sint *T
			BtF[NLAT-(i+1)] = (dse-dso) + I*(te-to);
			BpF[i] = I*(se+so) - (dte+dto);			// Bp = I.m/sint * S - dT/dt
			BpF[NLAT-(i+1)] = I*(se-so) - (dte-dto);
			i++;
			dyl += (LMAX-m+1);
		}
		BpF += NLAT;	BtF += NLAT;
	}
	for(im=MMAX+1; im<=NPHI/2; im++) {	// padding for high m's
		for (i=0;i<NLAT;i++) {
			BpF[i] = 0.0;	BtF[i] = 0.0;
		}
		BpF += NLAT;	BtF += NLAT;
	}

	BpF -= NLAT*(NPHI/2+1);		// restore original pointers
	BtF -= NLAT*(NPHI/2+1);
	fftw_execute_dft_c2r(ifft, BpF, (double *) BpF);
	fftw_execute_dft_c2r(ifft, BtF, (double *) BtF);
}

void spat_to_SHsphtor(complex double *BtF, complex double *BpF, complex double *Slm, complex double *Tlm)
{
	complex double te[NLAT/2], to[NLAT/2];	// theta even and odd parts
	complex double pe[NLAT/2], po[NLAT/2];	// phi even and odd parts.
	complex double *Sl, *Tl;		// virtual pointers for given im
	struct DtDp *idyl;
	long int i,im,m,l;

	fftw_execute_dft_r2c(fft,(double *) BtF, BtF);
	fftw_execute_dft_r2c(fft,(double *) BpF, BpF);

	im = 0;		// idyl.p = 0.0 : and evrything is REAL
		m=im*MRES;
		for (i=0;i<NLAT/2;i++) {	// compute symmetric and antisymmetric parts.
			(double) te[i] = (double) BtF[i] + (double) BtF[NLAT-(i+1)];
			(double) to[i] = (double) BtF[i] - (double) BtF[NLAT-(i+1)];
			(double) pe[i] = (double) BpF[i] + (double) BpF[NLAT-(i+1)];
			(double) po[i] = (double) BpF[i] - (double) BpF[NLAT-(i+1)];
		}
		l=m;
		Sl = &Slm[LiM(0,im)];		// virtual pointer for l=0 and im
		Tl = &Tlm[LiM(0,im)];
		idyl = idylm[im];
		BtF += NLAT;	BpF += NLAT;
		while (l<LMAX) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
			Sl[l] = 0.0;	Sl[l+1] = 0.0;		// Slm[LiM(l,im)] = 0.0;	Slm[LiM(l+1,im)] = 0.0;
			Tl[l] = 0.0;	Tl[l+1] = 0.0;
			for (i=0;i<NLAT/2;i++) {
				(double) Sl[l] += idyl[i].t * (double) to[i];
				(double) Tl[l] -= idyl[i].t * (double) po[i];
				
				(double) Sl[l+1] += idyl[NLAT/2 +i].t * (double) te[i];
				(double) Tl[l+1] -= idyl[NLAT/2 +i].t * (double) pe[i];
			}
			l+=2;
			idyl += NLAT;
		}
		if (l==LMAX) {
			Sl[l] = 0.0;	Tl[l] = 0.0;
			for (i=0;i<NLAT/2;i++) {
				(double) Sl[l] += idyl[i].t * (double) to[i];
				(double) Tl[l] -= idyl[i].t * (double) po[i];
			}
		}
	for (im=1;im<=MMAX;im++) {
		m=im*MRES;
		for (i=tm[im];i<NLAT/2;i++) {	// compute symmetric and antisymmetric parts.
			te[i] = BtF[i] + BtF[NLAT-(i+1)];
			to[i] = BtF[i] - BtF[NLAT-(i+1)];
			pe[i] = BpF[i] + BpF[NLAT-(i+1)];
			po[i] = BpF[i] - BpF[NLAT-(i+1)];
		}
		l=m;
		Sl = &Slm[LiM(0,im)];		// virtual pointer for l=0 and im
		Tl = &Tlm[LiM(0,im)];
		idyl = idylm[im];
		BtF += NLAT;	BpF += NLAT;
		while (l<LMAX) {		// ops : NLAT/2 * (2*(LMAX-m+1) + 4) : almost twice as fast.
			Sl[l] = 0.0;	Sl[l+1] = 0.0;		// Slm[LiM(l,im)] = 0.0;	Slm[LiM(l+1,im)] = 0.0;
			Tl[l] = 0.0;	Tl[l+1] = 0.0;
			for (i=tm[im];i<NLAT/2;i++) {	// tm[im] : polar optimization
				Sl[l] += idyl[i].t *to[i] - idyl[i].p *pe[i]*I;		// ref: these E. Dormy p 72.
				Tl[l] -= idyl[i].t *po[i] + idyl[i].p *te[i]*I;
				
				Sl[l+1] += idyl[NLAT/2 +i].t *te[i] - idyl[NLAT/2 +i].p *po[i]*I;
				Tl[l+1] -= idyl[NLAT/2 +i].t *pe[i] - idyl[NLAT/2 +i].p *to[i]*I;
			}
			l+=2;
			idyl += NLAT;
		}
		if (l==LMAX) {
			Sl[l] = 0.0;	Tl[l] = 0.0;
			for (i=tm[im];i<NLAT/2;i++) {	// polar optimization
				Sl[l] += idyl[i].t *to[i] - idyl[i].p *pe[i]*I;
				Tl[l] -= idyl[i].t *po[i] + idyl[i].p *te[i]*I;
			}
		}
	}
}


/*
	INITIALIZATION FUNCTIONS
*/

void runerr(const char * error_text)
{
	printf("*** Run-time error : %s\n",error_text);
	exit(1);
}

// Generates the abscissa and weights for a Gauss-Legendre quadrature.
// Newton method from initial Guess to find the zeros of the Legendre Polynome
// x = abscissa, w = weights, n points.
// Reference:  Numerical Recipes, Cornell press.
void Gauss(double *x, double *w, int n)
{
	double z, z1, p1, p2, p3, pp, eps;
	long int i,j,m;

	eps = 1.0e-15;	// desired precision, minimum = 2.2204e-16 (double)

	m = (n+1)/2;
	for (i=1;i<=m;i++) {
		z = cos(pi*((double)i-0.25)/((double)n+0.5));
		z1 = z+1;
		while ( fabs(z-z1) > eps )
		{
			p1 = 1.0;
			p2 = 0.0;
			for(j=1;j<=n;j++) {
        			p3 = p2;
        			p2 = p1;
        			p1 = ((2*j-1)*z*p2-(j-1)*p3)/j;	// The Legendre polynomial...
			}
			pp = ((double)n)*(z*p1-p2)/(z*z-1.0);                       // ... and its derivative.
			z1 = z;
			z = z1-p1/pp;
		}
		x[i-1] = -z;		// Build up the abscissas.
		x[n-i] = z;
		w[i-1] = 2.0/((1-z*z)*(pp*pp));		// Build up the weights.
		w[n-i] = w[i-1];
	}

// as we started with initial guesses, we should check if the gauss points are actually unique.
	for (i=n; i>0; i--) {
		if (x[i] == x[i-1]) runerr("bad gauss points\n");
	}
}

// initialize FFTs using FFTW. stride = NLAT, (contiguous l)
void planFFT()
{
	complex double *ShF;
	double *Sh;
	int nfft = NPHI;
	int ncplx = NPHI/2 +1;
	int nreal;
	
	nreal = 2*ncplx;
	
// Allocate dummy Spatial Fields.
	ShF = (complex double *) fftw_malloc(ncplx * NLAT * sizeof(complex double));
	Sh = (double *) ShF;

	printf("[FFTW] Mmax=%d, Nphi=%d\n",MMAX,NPHI);

	if (NPHI < 2*MMAX) runerr("[FFTW] the condition Nphi >= 2*Mmax is not met.");
	if (NPHI < 3*MMAX) printf("       ! Warning : 2/3 rule for anti-aliasing not met !\n");
	
// IFFT : unnormalized.
	ifft = fftw_plan_many_dft_c2r(1, &nfft, NLAT, ShF, &ncplx, NLAT, 1, Sh, &nreal, NLAT, 1, fftw_plan_mode);
	if (ifft == NULL)
		runerr("[FFTW] ifft planning failed !");

// FFT : must be normalized.
	fft = fftw_plan_many_dft_r2c(1, &nfft, NLAT, Sh, &nreal, NLAT, 1, ShF, &ncplx, NLAT, 1, fftw_plan_mode);
	if (fft == NULL)
		runerr("[FFTW] fft planning failed !");

//	fft_norm = 1.0/nfft;
	fftw_free(ShF);
	printf("       done.\n");
}

// initialize SH transform.
void init_SH()
{
	double xg[NLAT], wg[NLAT];	// gauss points and weights.
	double dtylm[LMAX+1];		// temp storage for derivative : d(P_l^m(x))/dx
	double eps = POLAR_OPT_THRESHOLD;	// eps : polar coefficients below that threshold are neglected (for high ms)
	double iylm_fft_norm = 2.0*pi/NPHI;	// normation FFT pour iylm
	double t,tmax;
	long int it,im,m,l;

	printf("[init_SH] Lmax=%d, Nlat=%d, Mres=%d, Mmax*Mres=%d, LMmax=%d\n",LMAX,NLAT,MRES,MMAX*MRES,NLM);
	if (MMAX*MRES > LMAX) runerr("[init_SH] MMAX*MRES should not exceed LMAX");
	if (NLAT <= LMAX) runerr("[init_SH] NLAT should be at least LMAX+1");
	
	Gauss(xg,wg,NLAT);	// generate gauss nodes and weights : xg = ]-1,1[ = cos(theta) 
	for (it=0; it<NLAT/2; it++) {
		ct[it] = xg[NLAT-1-it];			// on prend theta = ]0,pi/2[ => cos(theta) = ]1,0[
		st[it] = sqrt(1.0 - ct[it]*ct[it]);
		st_1[it] = 1.0/sqrt(1.0 - ct[it]*ct[it]);
	}

#ifdef _SH_DEBUG_
// TEST if gauss points are ok.
	tmax = 0.0;
	for (it = 0; it<NLAT/2; it++) {
		t = gsl_sf_legendre_sphPlm(NLAT, 0, ct[it]);
		if (t>tmax) tmax = t;
//		printf("i=%d, x=%12.12g, p=%12.12g\n",i,ct[i],t);
	}
	printf("          max zero at Gauss node for Plm[l=LMAX+1,m=0] : %g\n",tmax);
#endif

// Allocate legendre functions lookup tables.
	ylm[0] = (double *) fftw_malloc(sizeof(double)* NLM*NLAT/2);
	dylm[0] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* NLM*NLAT/2);
	iylm[0] = (double *) fftw_malloc(sizeof(double)* NLM*NLAT/2);
	idylm[0] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* NLM*NLAT/2);
	for (im=0; im<MMAX; im++) {
		m = im*MRES;
		ylm[im+1] = ylm[im] + NLAT/2*(LMAX+1 -m);
		dylm[im+1] = dylm[im] + NLAT/2*(LMAX+1 -m);
		iylm[im+1] = iylm[im] + NLAT/2*(LMAX+1 -m);
		idylm[im+1] = idylm[im] + NLAT/2*(LMAX+1 -m);
	}

// Even/Odd symmetry : ylm is even or odd across equator, as l-m is even or odd => only NLAT/2 points required.
// for synthesis (inverse transform)
	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
//		ylm[im] = (double *) fftw_malloc(sizeof(double)* (LMAX+1-m)*NLAT/2);
//		dylm[im] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* (LMAX+1-m)*NLAT/2);
		for (it=0;it<NLAT/2;it++) {
			gsl_sf_legendre_sphPlm_deriv_array(LMAX, m, ct[it], ylm[im] + it*(LMAX-m+1), dtylm);	// fixed im legendre functions lookup table.
			for (l=m; l<=LMAX; l++) {
				dylm[im][it*(LMAX-m+1) + (l-m)].t = -st[it] * dtylm[l-m];	// d(Plm(cos(t)))/dt = -sin(t) d(Plm(x))/dx
				dylm[im][it*(LMAX-m+1) + (l-m)].p = ylm[im][it*(LMAX-m+1) + (l-m)] *m/st[it];	// 1/sint(t) dYlm/dphi
			}
		}
	}
	
// for analysis (decomposition, direct transform) : transpose and multiply by gauss weight and other normalizations.
	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
//		iylm[im] = (double *) fftw_malloc(sizeof(double)* (LMAX+1-m)*NLAT/2);
//		idylm[im] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* (LMAX+1-m)*NLAT/2);
		for (it=0;it<NLAT/2;it++) {
			for (l=m;l<=LMAX;l++) {
				iylm[im][(l-m)*NLAT/2 + it] = ylm[im][it*(LMAX-m+1) + (l-m)] * wg[it] *iylm_fft_norm;
				if (l!=0) {
					idylm[im][(l-m)*NLAT/2 + it].t = dylm[im][it*(LMAX-m+1) + (l-m)].t * wg[it] *iylm_fft_norm /(l*(l+1));
					idylm[im][(l-m)*NLAT/2 + it].p = dylm[im][it*(LMAX-m+1) + (l-m)].p * wg[it] *iylm_fft_norm /(l*(l+1));
				} else {	// les derivees sont nulles pour l=0 (=> m=0)
					idylm[im][(l-m)*NLAT/2 + it].t = 0.0;
					idylm[im][(l-m)*NLAT/2 + it].p = 0.0;
				}
			}
		}
	}

// POLAR OPTIMIZATION : analyzing coefficients, some can be safely neglected.
	for (im=0;im<=MMAX;im++) {
		m = im*MRES;
		tm[im] = 1.0*NLAT;
		for (l=m;l<=LMAX;l++) {
			it=0;
			while( fabs(ylm[im][it*(LMAX-m+1) + (l-m)]) < eps ) { it++; }
			if (tm[im] > it) tm[im] = it;
		}
	}
	if (eps > 0.0) {
		printf("          polar optimization threshold = %e\n",eps);
#ifdef _SH_DEBUG_
		printf("          tm[im]=");
		for (im=0;im<=MMAX;im++)
			printf(" %d",tm[im]);
		printf("\n");
#endif
	}

	planFFT();		// initialize fftw

// Additional arrays :
	it = 0;
	for (im=0;im<=MMAX;im++) {
		for (l=m;l<=LMAX;l++) {
			l2[it] = l*(l+1);	l_2[it] = 1.0/(l*(l+1));
			it++;
		}
	}
	l_2[0] = 0.0;	// undefined for l=0 => replace by 0.
}
