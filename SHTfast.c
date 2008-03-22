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


// SHT.h : parameter for SHT (sizes : LMAX, NLAT, MMAX, MRES, NPHI)

//  SIZES  //
// LMAX : maximum degree of Spherical Harmonic
#define LMAX 96
// NLAT : number of latitudinal (theta) gauss points, at least (LMAX+1), must be EVEN (and (LMAX+1)*2 for dealias)
#define NLAT (3*LMAX)

// MMAX : max number of fourrier decomposition (degree = MMAX * MRES)
// power of two is better (but not mandatory), MMAX*MRES <= LMAX
#define MMAX 32
// MRES : azimutal symmetry
#define MRES 1
// NPHI : number of azimutal grid points, at least 2*MMAX, 3*MMAX for antialiasing and 4*MMAX for full dealias
#define NPHI (3*MMAX)
// compute and print some debugging information...
#define _SH_DEBUG_


// NLM : total number of Spherical Harmonic coefficients.
//#define NLM ( ((MMAX+1)*(2*LMAX+2-MMAX))/2 )
#define NLM ( (MMAX+1)*(LMAX+1) - MRES* MMAX*(MMAX+1)/2 )
// LM(l,m) : index in the Spherical Harmonic coefficient array [ (l,m) space ]
//#define LM(l,m) ((m*(2*LMAX +1-m))/2 + l)
#define LM(l,m) ( (m/MRES)*(2*LMAX+3 -m)/2 + l-m )
#define LiM(l,im) ( im*(2*LMAX+3 -(im+2)*MRES)/2 + l )

double pi = atan(1.0)*4.0;
double l2[NLM], l_2[NLM];			// l(l+1) and 1/(l(l+1))

double wc[2*NLAT];	// chebychev weights. for even and odd m's.
double ct[NLAT],st[NLAT],st_1[NLAT];	// cos(theta), sin(theta), and 1/sin(theta) arrays.

double* ylm[MMAX+1];		// matrix for direct transform

fftw_plan ifft, fft;	// plans for FFTW.
fftw_plan idct, dct;
fftw_plan idct2, dct2;
unsigned fftw_plan_mode = FFTW_PATIENT;		// defines the default FFTW planner mode.


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
	for (im=0; im<MMAX; im+=2) {
		for (it=0; it<NLAT*2; it++)				// span two im's because
			ShF[im*NLAT + it] *= wc[it];		// even and odd weights are different.
	}
	if (im == MMAX) {
		for (it=0; it<NLAT; it++)
			ShF[im*NLAT + it] *= wc[it];
	}
	fftw_execute_r2r(dct,(double *) ShF, (double *) ShF);		// DCT of weighted data
	fftw_execute_r2r(dct2,((double *) ShF)+1, ((double *) ShF)+1);		// DCT of weighted data

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
	fftw_execute_r2r(idct2,((double *) ShF)+1, ((double *) ShF)+1);		// iDCT

	for (im=1; im<=MMAX; im+=2) {	// odd m's must be multiplied by sin(theta) which was removed from ylm's
		for (it=0; it<NLAT; it++)
			ShF[im*NLAT + it] *= st[it];
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
	if (ifft == NULL)
		runerr("[FFTW] ifft planning failed !");

// FFT : must be normalized.
	fft = fftw_plan_many_dft_r2c(1, &nfft, NLAT, Sh, &nreal, NLAT, 1, ShF, &ncplx, NLAT, 1, fftw_plan_mode);
	if (fft == NULL)
		runerr("[FFTW] fft planning failed !");

//	fft_norm = 1.0/nfft;

/* LATITUDINAL DCT (THETA) */
	r2r_kind = FFTW_REDFT10;
	dct = fftw_plan_many_r2r(1, &ndct, MMAX+1, Sh, &ndct, 2, 2*NLAT, Sh, &ndct, 2, 2*NLAT, &r2r_kind, fftw_plan_mode );
	dct2 = fftw_plan_many_r2r(1, &ndct, MMAX+1, Sh+1, &ndct, 2, 2*NLAT, Sh+1, &ndct, 2, 2*NLAT, &r2r_kind, fftw_plan_mode );
	r2r_kind = FFTW_REDFT01;
	idct = fftw_plan_many_r2r(1, &ndct, MMAX+1, Sh, &ndct, 2, 2*NLAT, Sh, &ndct, 2, 2*NLAT, &r2r_kind, fftw_plan_mode );
	idct2 = fftw_plan_many_r2r(1, &ndct, MMAX+1, Sh+1, &ndct, 2, 2*NLAT, Sh+1, &ndct, 2, 2*NLAT, &r2r_kind, fftw_plan_mode );

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

// initialize SH transform.
void init_SH()
{
	double *ft;					// for DCT.
	double *ylmt;				// temp storage for Plm's
	double *yl;			// virtual pointer.
	fftw_plan dct;
	double iylm_fft_norm = 2.0*pi/(NPHI);	// normation FFT +DCT pour iylm
	double t,tmax;
	long int it,im,m,l,c;
#ifdef _SH_DEBUG_
	FILE *fp;
	fp = fopen("ylm_dct","w");
#endif

	printf("[init_SH] Lmax=%d, Nlat=%d, Mres=%d, Mmax*Mres=%d, LMmax=%d\n",LMAX,NLAT,MRES,MMAX*MRES,NLM);
	if (MMAX*MRES > LMAX) runerr("[init_SH] MMAX*MRES should not exceed LMAX");
	if (NLAT <= LMAX) runerr("[init_SH] NLAT should be at least LMAX+1");
	
	ChebychevNodes(ct,wc,NLAT);		// chebychev points and weights.
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
	ylmt = (double *) fftw_malloc(sizeof(double)* (LMAX+1)*NLAT);
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
			fprintf(fp,"*** im=%d, l=%d ::",im,l);
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
					fprintf(fp," %.3f",yl[it]);
				}
				fprintf(fp,"\n");
			}
#endif
			if ((l-m)%2) {		// l-m odd : go to next line of storage.
				yl += (l+1 - m%2);
			}
		}
	}
	
	fftw_destroy_plan(dct);
	fftw_free(ft);
	fftw_free(ylmt);
	
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
