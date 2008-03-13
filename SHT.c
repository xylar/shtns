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

double pi = atan(1.0)*4.0;
double l2[NLM], l_2[NLM];			// l(l+1) and 1/(l(l+1))
double ct[NLAT/2], st[NLAT/2], st_1[NLAT/2];	// cos(theta), sin(theta), 1/sin(theta);

long int tm[MMAX+1];	// start theta value for SH (polar optimization : near the poles the legendre polynomials go to zero for high m's)
double* ylm[MMAX+1];	// matrix for direct transform
double* dtylm[MMAX+1];	// theta derivative matrix for direct transform
double* iylm[MMAX+1];	// matrix for inverse transform

fftw_plan ifft, fft;	// plans for FFTW.
unsigned fftw_plan_mode = FFTW_PATIENT;		// defines the default FFTW planner mode.


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

	for (im=0;im<=MMAX;im++) {
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
				Sl[l] += iyl[i] * fp[i];		// Slm[LiM(l,im)] += iylm[im][(l-m)*NLAT/2 + i] * fp[i];
				Sl[l+1] += iyl[NLAT/2 + i] * fm[i];	// Slm[LiM(l+1,im)] += iylm[im][(l+1-m)*NLAT/2 + i] * fm[i];
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

	for (im=0; im<=MMAX; im++) {
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
			fe = 0.0;	fo = 0.0;
			l=m;
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
	fftw_execute_dft_c2r(ifft, ShF, (double *) ShF);
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
	double *yl, *dtyl;
	long int i,im,m,l;

	for (im=0; im<=MMAX; im++) {
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
		yl = ylm[im] + i*(LMAX-m+1) -m;
		dtyl = dtylm[im] + i*(LMAX-m+1) -m;
		while (i<NLAT/2) {	// ops : NLAT/2 * [ (lmax-m+1)*2 + 4]	: almost twice as fast.
			l=m;
			se = 0.0;	so = 0.0;	dse = 0.0;	dso = 0.0;
			te = 0.0;	to = 0.0;	dte = 0.0;	dto = 0.0;
			while (l<LMAX) {	// compute even and odd parts
				se += yl[l] * Sl[l];
				so += yl[l+1] * Sl[l+1];
				dso += dtyl[l] * Sl[l];
				dse += dtyl[l+1] * Sl[l+1];
				te += yl[l] * Tl[l];
				to += yl[l+1] * Tl[l+1];
				dto += dtyl[l] * Tl[l];
				dte += dtyl[l+1] * Tl[l+1];
				l+=2;
			}
			if (l==LMAX) {
				se += yl[l] * Sl[l];
				dso += dtyl[l] * Sl[l];
				te += yl[l] * Tl[l];
				dto += dtyl[l] * Tl[l];
			}
			se *= I*st_1[i]*m;	so *= I*st_1[i]*m;
			te *= I*st_1[i]*m;	to *= I*st_1[i]*m;

			BtF[i] = (dse+dso) + (te+to);			// Bt = dS/dt       + Im/sint *T
			BtF[NLAT-(i+1)] = (dse-dso) + (te-to);
			BpF[i] = (se+so) - (dte+dto);			// Bp = Im/sint * S - dT/dt
			BpF[NLAT-(i+1)] = (se-so) - (dte-dto);

			i++;
			yl += (LMAX-m+1);	dtyl += (LMAX-m+1);
		}
		BpF += NLAT;	BtF += NLAT;
	}
	for(im=MMAX+1; im<=NPHI/2; im++) {	// padding for high m's
		for (i=0;i<NLAT;i++) {
			BpF[i] = 0.0;	BtF[i] = 0.0;
		}
		BpF += NLAT;	BtF += NLAT;
	}

	fftw_execute_dft_c2r(ifft, BpF, (double *) BpF);
	fftw_execute_dft_c2r(ifft, BtF, (double *) BtF);
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
	complex double * ShF;
	int nfft = NPHI;
	int ncplx = NPHI/2 +1;
	int nreal;
	
	nreal = 2*ncplx;
	
// Allocate dummy Spatial Fields.
	ShF = (complex double *) fftw_malloc(ncplx * NLAT * sizeof(complex double));

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
	double eps = POLAR_OPT_THRESHOLD;	// eps : polar coefficients below that threshold are neglected (for high ms)
	double iylm_fft_norm = 2.0*pi/NPHI;	// normation FFT pour iylm
	double t,tmax;
	long int it,im,m,l;

	printf("[init_SH] Lmax=%d, Nlat=%d, Mres=%d, Mmax*Mres=%d, LMmax=%d\n",LMAX,NLAT,MRES,MMAX*MRES,NLM);
	if (MMAX*MRES > LMAX) runerr("[init_SH] MMAX*MRES should not exceed LMAX");
	if (NLAT <= LMAX) runerr("[init_SH] NLAT should be at least LMAX+1");
	
	Gauss(xg,wg,NLAT);	// generate gauss nodes and weights [ x = cos(theta) ]
	for (it=0; it<NLAT/2; it++) {
		ct[it] = xg[it];
		st[it] = sqrt(1.0 - xg[it]*xg[it]);
		st_1[it] = 1.0/sqrt(1.0 - xg[it]*xg[it]);
	}

#ifdef _SH_DEBUG_
// TEST if gauss points are ok.
	tmax = 0.0;
	for (it = 0; it<NLAT/2; it++) {
		t = gsl_sf_legendre_sphPlm(NLAT, 0, xg[it]);
		if (t>tmax) tmax = t;
//		printf("i=%d, x=%12.12g, p=%12.12g\n",i,xg[i],t);
	}
	printf("          max zero at Gauss node for Plm[l=LMAX+1,m=0] : %g\n",tmax);
#endif

// Even/Odd symmetry : ylm is even or odd across equator, as l-m is even or odd => only NLAT/2 points required.
// for synthesis (inverse transform)
	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
		ylm[im] = (double *) malloc(sizeof(double)* (LMAX-m+1)*NLAT);		// P_l^m(x)   |x| <= 1.0
		dtylm[im] = (double *) malloc(sizeof(double)* (LMAX-m+1)*NLAT/2);	// d(P_l^m(x))/dx
		for (it=0;it<NLAT/2;it++) {
			gsl_sf_legendre_sphPlm_deriv_array(LMAX, m, xg[it], ylm[im] + it*(LMAX-m+1), dtylm[im] + it*(LMAX-m+1));	// fixed im legendre functions lookup table.
			for (l=m; l<=LMAX; l++) {
				dtylm[im][it*(LMAX-m+1) + (l-m)] *= -st[it];	// d(Plm(cos(t)))/dt = -sin(t) d(Plm(x))/dx
			}
		}
	}
	
// for analysis (decomposition, direct transform) : transpose and multiply by gauss weight.
	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
		iylm[im] = (double *) malloc(sizeof(double)* (LMAX-m+1)*NLAT/2);
		for (it=0;it<NLAT/2;it++) {
			for (l=m;l<=LMAX;l++) {
				iylm[im][(l-m)*NLAT/2 + it] = ylm[im][it*(LMAX-m+1) + (l-m)] * wg[it] *iylm_fft_norm;
			}
		}
	}

// POLAR OPTIMIZATION : analysing coefficients, some can be safely neglected.
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
			l2[it] = l*(l+1);
			it++;
		}
	}
}
