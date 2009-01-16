///////////////////////////////////////////////
// SHT : Spherical Harmonic Transform, with DCT.
//   uses function from S2kit
//   requires SHT.h for size parameters.
//////////////////////////////////////////////

#include <complex.h>
#include <math.h>
// FFTW la derivee d/dx = ik	(pas de moins !)
#include <fftw3.h>
// GSL for Legendre functions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>


#ifndef LMAX

//  SIZES  //
// LMAX : maximum degree of Spherical Harmonic
#define LMAX 8
// NLAT : number of latitudinal (theta) gauss points, at least (LMAX+1), must be EVEN (and (LMAX+1)*2 for dealias)
#define NLAT (LMAX+2)

// MMAX : max number of fourrier decomposition (degree = MMAX * MRES)
// power of two is better (but not mandatory), MMAX*MRES <= LMAX
#define MMAX 4
// MRES : azimutal symmetry
#define MRES 2
// NPHI : number of azimutal grid points, at least 2*MMAX, 3*MMAX for antialiasing and 4*MMAX for full dealias
#define NPHI (3*MMAX)
// compute and print some debugging information...
#define _SH_DEBUG_

#else
 #define NLAT (2*NLAT_2)
#endif


// NLM : total number of Spherical Harmonic coefficients.
//#define NLM ( ((MMAX+1)*(2*LMAX+2-MMAX))/2 )
#define NLM ( (MMAX+1)*(LMAX+1) - MRES* MMAX*(MMAX+1)/2 )
// LM(l,m) : index in the Spherical Harmonic coefficient array [ (l,m) space ]
//#define LM(l,m) ((m*(2*LMAX +1-m))/2 + l)
#define LiM(l,im) ( (im*(2*LMAX+2 -MRES*(im+1)))/2 + l )
#define LM(l,m) ( (m*(2*LMAX+2 -(m+MRES)))/(2*MRES) + l )

// useful values for some basic spherical harmonic representations
// Y00_1 = 1/Y00 = spherical harmonic representation of 1 (l=0,m=0)
#define Y00_1 sqrt(4*pi)
// Y10_ct = spherical harmonic representation of cos(theta) (l=1,m=0)
#define Y10_ct sqrt(4.0*pi/3.0)


double pi = atan(1.0)*4.0;
double l2[NLM], l_2[NLM];			// l(l+1) and 1/(l(l+1))

double wc[2*NLAT];	// chebychev weights. for even and odd m's.
double ct[NLAT],st[NLAT],st_1[NLAT];	// cos(theta), sin(theta), and 1/sin(theta) arrays.

double* ylm[MMAX+1];		// matrix for inverse transform (synthesis)
double* zlm[MMAX+1];		// matrix for direct transform (analysis)

fftw_plan ifft, fft;	// plans for FFTW.
fftw_plan idct, dct;
//fftw_plan idct2, dct2;
unsigned fftw_plan_mode = FFTW_PATIENT;		// defines the default FFTW planner mode.

/*
/////////////////////////////////////////////////////
//   Scalar Spherical Harmonics Transform
// input  : ShF = spatial/fourrier data : complex double array of size NLAT*(NPHI/2+1) or double array of size NLAT*(NPHI/2+1)*2
// output : Slm = spherical harmonics coefficients : complex double array of size NLM
void spat_to_SH(complex double *ShF, complex double *Slm)
{
	complex double *Sl;		// virtual pointers for given im
	double *yl;
	long int it,im,m,l;

	fftw_execute_dft_r2c(fft,(double *) ShF, ShF);
	if (MRES % 2) {		// odd m's are present
		for (im=0; im<MMAX; im+=2) {
			for (it=0; it<NLAT*2; it++)				// span two im's because
				ShF[im*NLAT + it] *= wc[it];		// even and odd weights are different.
		}
		if (im == MMAX) {
			for (it=0; it<NLAT; it++)
				ShF[im*NLAT + it] *= wc[it];
		}
	} else {			// m is always even
		for (im=0; im<=MMAX; im++) {
			for (it=0; it<NLAT; it++)
				ShF[im*NLAT + it] *= wc[it];
		}
	}
	fftw_execute_r2r(dct,(double *) ShF, (double *) ShF);		// DCT of weighted data

	for (im=0; im<=MMAX; im++) {
		m=im*MRES;
		Sl = &Slm[LiM(0,im)];		// virtual pointer for l=0 and im
		yl = ylm[im];
		ShF[0] *= 1.0/sqrt(2.);
		for (l=m; l<LMAX; l+=2) {		// l has parity of m
			Sl[l] = 0.0;	Sl[l+1] = 0.0;
			for (it=0; it<=l; it+=2) {		// it : DCT index
				Sl[l]   += ShF[it]   * yl[it];		// * P(l)(it)
				Sl[l+1] += ShF[it+1] * yl[it+1];	// * P(l+1)(it+1)
			}
			yl += (l+2 - (m&1));
		}
		if (l==LMAX) {
			Sl[l] = 0.0;
			for (it=0; it<=l; it+=2)
				Sl[l]   += ShF[it] * yl[it];
		}
		ShF += NLAT;
	}
}
*/
/////////////////////////////////////////////////////
//   Scalar Spherical Harmonics Transform
// input  : ShF = spatial/fourrier data : complex double array of size NLAT*(NPHI/2+1) or double array of size NLAT*(NPHI/2+1)*2
// output : Slm = spherical harmonics coefficients : complex double array of size NLM
void spat_to_SH(complex double *ShF, complex double *Slm)
{
	complex double *Sl;		// virtual pointers for given im
	double *zl;
	long int k,im,m,l;

	fftw_execute_dft_r2c(fft,(double *) ShF, ShF);
	fftw_execute_r2r(dct,(double *) ShF, (double *) ShF);		// DCT

	for (im=0; im<=MMAX; im++) {
		m=im*MRES;
		Sl = &Slm[LiM(0,im)];		// virtual pointer for l=0 and im
		zl = zlm[im];
		ShF[0] *= 1.0/sqrt(2.);
		for (l=m; l<=LMAX; l++) {		// l has parity of m
			Sl[l] = 0.0;
			for (k=0; k<=l; k++) {
				Sl[l]   += ShF[k]   * zl[k];
			}
			zl += LMAX+1;
		}
		ShF += NLAT;
	}
}


/////////////////////////////////////////////////////
//   Scalar inverse Spherical Harmonics Transform
// input  : Slm = spherical harmonics coefficients : complex double array of size NLM [unmodified]
// output : ShF = spatial/fourrier data : complex double array of size NLAT*(NPHI/2+1) or double array of size NLAT*(NPHI/2+1)*2
void SH_to_spat(complex double *Slm, complex double *ShF)
{
	complex double *Sl;
	double *yl;
	long int it,im,m,l;

	for (im=0;im<=MMAX;im++) {
		m=im*MRES;
		Sl = &Slm[LiM(0,im)];		// virtual pointer for l=0 and im
		yl = ylm[im];
		for (it=0; it<NLAT; it++)
			ShF[it] = 0.0;		// zero out array (includes DCT padding)
		for (l=m; l<LMAX; l+=2) {
			for (it=0; it<=l; it+=2) {
				ShF[it]   += Sl[l]   * yl[it];
				ShF[it+1] += Sl[l+1] * yl[it+1];
			}
			yl += (l+2 - (m&1));
		}
		if (l==LMAX) {
			for (it=0; it<=l; it+=2)
				ShF[it] += Sl[l] * yl[it];
		}
		ShF[0] *= sqrt(2.);
		ShF += NLAT;
	}
	for(im=MMAX+1; im<=NPHI/2; im++) {	// FFT padding for high m's
		for (it=0; it<NLAT; it++)
			ShF[it] = 0.0;
		ShF += NLAT;
	}

	ShF -= NLAT*(NPHI/2+1);		// restore original pointer
	fftw_execute_r2r(idct,(double *) ShF, (double *) ShF);		// iDCT

	if (MRES % 2) {		// odd m's must be multiplied by sin(theta) which was removed from ylm's
		for (im=1; im<=MMAX; im+=2) {
			for (it=0; it<NLAT; it++)
				ShF[im*NLAT + it] *= st[it];
		}
	}
	fftw_execute_dft_c2r(ifft, ShF, (double *) ShF);
}



/*
	INITIALIZATION FUNCTIONS
*/

void runerr(const char * error_text)
{
	printf("*** Run-time error : %s\n",error_text);
	exit(1);
}


// initialize FFTs using FFTW. stride = NLAT, (contiguous l)
void planFFT()
{
	complex double *ShF;
	double *Sh;
	int nfft = NPHI;
	int ncplx = NPHI/2 +1;
	int nreal;
	int ndct = NLAT;
	fftw_r2r_kind r2r_kind;
	fftw_iodim dims, hdims[2];
	
	nreal = 2*ncplx;
	
/* AZIMUTAL FFT (PHI)*/
// Allocate dummy Spatial Fields.
	ShF = (complex double *) fftw_malloc(ncplx * NLAT * sizeof(complex double));
	Sh = (double *) ShF;

	printf("[FFTW] Mmax=%d, Nphi=%d\n",MMAX,NPHI);

	if (NPHI < 2*MMAX) runerr("[FFTW] the condition Nphi >= 2*Mmax is not met.");
	if (NPHI < 3*MMAX) printf("       ! Warning : 2/3 rule for anti-aliasing not met !\n");
	
// IFFT : unnormalized.
	ifft = fftw_plan_many_dft_c2r(1, &nfft, NLAT, ShF, &ncplx, NLAT, 1, Sh, &nreal, NLAT, 1, fftw_plan_mode);
// FFT : must be normalized.
	fft = fftw_plan_many_dft_r2c(1, &nfft, NLAT, Sh, &nreal, NLAT, 1, ShF, &ncplx, NLAT, 1, fftw_plan_mode);
	if ((fft == NULL)||(ifft == NULL))
		runerr("[FFTW] fft planning failed !");

//	fft_norm = 1.0/nfft;

/* LATITUDINAL DCT (THETA) */
/*	r2r_kind = FFTW_REDFT10;
	dct = fftw_plan_many_r2r(1, &ndct, MMAX+1, Sh, &ndct, 2, 2*NLAT, Sh, &ndct, 2, 2*NLAT, &r2r_kind, fftw_plan_mode );
	dct2 = fftw_plan_many_r2r(1, &ndct, MMAX+1, Sh+1, &ndct, 2, 2*NLAT, Sh+1, &ndct, 2, 2*NLAT, &r2r_kind, fftw_plan_mode );
	r2r_kind = FFTW_REDFT01;
	idct = fftw_plan_many_r2r(1, &ndct, MMAX+1, Sh, &ndct, 2, 2*NLAT, Sh, &ndct, 2, 2*NLAT, &r2r_kind, fftw_plan_mode );
	idct2 = fftw_plan_many_r2r(1, &ndct, MMAX+1, Sh+1, &ndct, 2, 2*NLAT, Sh+1, &ndct, 2, 2*NLAT, &r2r_kind, fftw_plan_mode );
*/
	dims.n = NLAT;	dims.is = 2;	dims.os = 2;
	hdims[0].n = MMAX+1;	hdims[0].is = 2*NLAT; 	hdims[0].os = 2*NLAT;
	hdims[1].n = 2;			hdims[1].is = 1; 	hdims[1].os = 1;
	r2r_kind = FFTW_REDFT10;
	dct = fftw_plan_guru_r2r(1, &dims, 2, hdims, Sh, Sh, &r2r_kind, fftw_plan_mode );
	r2r_kind = FFTW_REDFT01;
	idct = fftw_plan_guru_r2r(1, &dims, 2, hdims, Sh, Sh, &r2r_kind, fftw_plan_mode );
	if ((dct == NULL)||(idct == NULL))
		runerr("[FFTW] dct planning failed !");

//	dct_norm = 1.0/(2*NLAT);

	fftw_free(ShF);
	printf("       done.\n");
}

/*	FROM S2KIT :
  makeweights: given a bandwidth bw, make weights for
  both even *and* odd order Legendre transforms.

  bw -> bandwidth of transform
  weights -> pointer to double array of length 4*bw; this
             array will contain the even and odd weights;
	     even weights start at weight[0], and odd weights
	     start at weights[2*bw]

*/
void ChebychevNodes( double *nodes, double *weights, int bw )
{
	int j, k;
	double fudge ;
	double tmpsum ;
	double twoN;
	
// NODES
	twoN = (double) (2*bw);

	for (j=0; j<bw; j++)
		nodes[j] = cos((( 2.0*((double)j)+1.0 ) * pi) / twoN);

// WEIGHTS
	bw = bw/2;

	fudge = M_PI/((double)(4*bw)) ;

	for ( j = 0 ; j < 2*bw ; j ++ )
	{
		tmpsum = 0.0 ;
		for ( k = 0 ; k < bw ; k ++ )
			tmpsum += 1./((double)(2*k+1)) *
					sin((double)((2*j+1)*(2*k+1))*fudge);
		tmpsum *= sin((double)(2*j+1)*fudge);
		tmpsum *= 2./((double) bw) ;

		weights[j] = tmpsum ;
		weights[j + 2*bw] = tmpsum * sin((double)(2*j+1)*fudge);
	}
}

// Generates the abscissa and weights for a Gauss-Legendre quadrature.
// Newton method from initial Guess to find the zeros of the Legendre Polynome
// x = abscissa, w = weights, n points.
// Reference:  Numerical Recipes, Cornell press.
void GaussNodes(double *x, double *w, int n)
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


//#define NG (LMAX+1)
#define NG ((LMAX+1)*2)
// initialize SH transform.
void init_SH(double eps)
{
	double tg[NG],xg[NG],wg[NG];	// Gauss quadrature weights.
	double *ft,*Zlm;		// for DCT.
	double *ylmt;				// temp storage for Plm's
	double *yl;			// virtual pointer.
	fftw_plan dct, idct;
	double iylm_fft_norm = 2.0*pi/(NPHI);	// normation FFT +DCT pour iylm
	double t,tmax;
	long int it,im,m,l,c,k;
#ifdef _SH_DEBUG_
	FILE *fp;
	fp = fopen("ylm_dct","w");
#endif

	printf("[init_SH] Lmax=%d, Nlat=%d, Mres=%d, Mmax*Mres=%d, LMmax=%d\n",LMAX,NLAT,MRES,MMAX*MRES,NLM);
	printf("          => using regular theta grid and DCT\n");
	if (MMAX*MRES > LMAX) runerr("[init_SH] MMAX*MRES should not exceed LMAX");
	if (NLAT <= LMAX) runerr("[init_SH] NLAT should be at least LMAX+1");
	
	ChebychevNodes(ct,wc,NLAT);		// chebychev (equaly-spaced) points and weights.
	for (it=0;it<NLAT; it++) {
		st[it] = sqrt(1.0 - ct[it]*ct[it]);
		st_1[it] = 1.0/st[it];
		wc[it] *= iylm_fft_norm;		// include FFT normation in weigths
		wc[NLAT+it] *= iylm_fft_norm;	// odd m's weight
	}

// Allocate legendre functions lookup tables.
	c = 0;			// how much memory to allocate for ylm dct
	for(im=0; im<=MMAX; im++) {
		m = im*MRES;
		for (l=m;l<=LMAX;l+=2)
			c += (l+2 - m%2);
	}
#ifdef _SH_DEBUG_
	printf("  ylm dct storage = %d  (w/o DCT : NLM*NLAT/4 = %d, ratio=%f)\n",c,NLM*NLAT/4,(double)c /(double)(NLM*NLAT/4));
#endif

	ylm[0] = (double *) fftw_malloc(sizeof(double)* c);
	for (im=0; im<MMAX; im++) {
		m = im*MRES;
		c = 0;
		for (l=m;l<=LMAX;l+=2)
			c += (l+2 - m%2);
		ylm[im+1] = ylm[im] + c;
	}

	ft = (double *) fftw_malloc(sizeof(double)* NLAT);
	dct = fftw_plan_r2r_1d( NLAT, ft, ft, FFTW_REDFT10, FFTW_ESTIMATE );	// quick and dirty dfts.
	if (dct == NULL) runerr("FFTW : dct could not be created...");

	ylmt = (double *) fftw_malloc(sizeof(double)* (LMAX+1)*NLAT);
	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
//		ylm[im] = (double *) fftw_malloc(sizeof(double)* (LMAX+1-m)*NLAT/2);
//		dylm[im] = (struct DtDp *) fftw_malloc(sizeof(struct DtDp)* (LMAX+1-m)*NLAT/2);
		for (it=0;it<NLAT;it++) {
			gsl_sf_legendre_sphPlm_array(LMAX, m, ct[it], ylmt + it*(LMAX+1));	// fixed im legendre functions lookup table.
		}
		yl = ylm[im];				// starting point for Plm storage.
		for (l=m;l<=LMAX;l++) {
			for (it=0;it<NLAT;it++) {
				if (m%2 == 0) {		// even m's
					ft[it] = ylmt[it*(LMAX+1) + (l-m)];
				} else {			// odd m's : divide by sin(theta) to get a cos poly.
					ft[it] = ylmt[it*(LMAX+1) + (l-m)] / st[it];
				}
			}
			fftw_execute(dct);		// take the DCT.
			for (it=(l-m)%2; it<=l; it+=2) {	// store non-zero coeffs : parity of (l-m), not larger than l.
				yl[it] = ft[it]/(2.0*NLAT);		// odd and even (l-m) are stored interleaved...
			}
			if ((l-m)%2 == 0) {
				yl[0] = yl[0] / sqrt(2.0);
			}
#ifdef _SH_DEBUG_
			fprintf(fp,"*** m=%d, l=%d ::",m,l);
			c = 0;
			for(it=0;it<NLAT;it++) {
				if (fabs(ft[it]) > 1.e-10) {
					fprintf(fp," %d %.3f",it,ft[it]); c++;
				}
			}
			fprintf(fp,"  C=%d\n",c);
			if (((l-m)%2)||(l==LMAX)) {
				fprintf(fp,"  ylm ::");
				for(it=0;it<=(l-m%2);it++) {
					fprintf(fp," %g",yl[it]);
				}
				fprintf(fp,"\n");
			}
#endif
			if ((l-m)%2) {		// l-m odd : go to next line of storage.
				yl += (l+1 - m%2);
			}
		}
	}
	fftw_free(ylmt);
	fftw_destroy_plan(dct);
	fftw_free(ft);
#ifdef _SH_DEBUG_
	fclose(fp);
#endif

	zlm[0] = (double *) fftw_malloc(sizeof(double)* NLM*(LMAX+1));
	for (im=0; im<MMAX; im++) {
		zlm[im+1] = zlm[im] + (LMAX+1-m)*(LMAX+1);
	}
	printf("     Gauss quadrature for equaly-spaced grid ...\n");
	GaussNodes(xg,wg,NG);	// for quadrature, gauss nodes and weights : xg = ]-1,1[ = cos(theta) 
	for (it=0; it<NG; it++) {
		tg[it] = acos(xg[it]);		// theta at gauss points.
	}
/* GAUSS QUADRATURE TO COMPUTE DIRECT TRANSFORM MATRIX (analysis) */
	Zlm = (double *) fftw_malloc(sizeof(double)* NLAT);
//	idct = fftw_plan_r2r_1d( NLAT, Zlm, Zlm, FFTW_REDFT01, FFTW_ESTIMATE );	// quick and dirty dfts.
	ylmt = (double *) fftw_malloc(sizeof(double)* NG*(LMAX+1));
	for (im=0; im<=MMAX; im++) {
		m = im*MRES;
		for (it=0;it<NG;it++) {
			gsl_sf_legendre_sphPlm_array(LMAX, m, xg[it], ylmt + it*(LMAX+1));	// fixed im legendre functions lookup table.
		}
		for (l=m;l<=LMAX;l++) {
			printf("* m=%d, l=%d ::",m,l);
			for (k=0;k<=LMAX;k++) {
				Zlm[k] = 0.0;
				for(it=0;it<NG;it++) {
					Zlm[k] += cos(tg[it]*k)*wg[it]*ylmt[it*(LMAX+1) + (l-m)];	// Gauss Quadrature of cos(k.theta)
				}
				printf(" %.3f",Zlm[k]);
				zlm[im][(LMAX+1)*(l-m) +k] = Zlm[k] /(2.0*NLAT);
				if (k==0) {
					zlm[im][(LMAX+1)*(l-m) +k] = Zlm[k] /(sqrt(2.0)*NLAT);
				}
			}
			for (k=LMAX+1;k<NLAT;k++)
				Zlm[k] = 0.0;
/*
			printf("\n      idct ::");
			fftw_execute(idct);
			// Zlm(theta) = iDCT( Zlm[k] )
			for (it=0;it<NLAT;it++)
				printf(" %.3f",Zlm[it]);
			printf("\n");
*/
		}
	}
	fftw_free(ylmt);
//	fftw_destroy_plan(idct);
	fftw_free(Zlm);

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
