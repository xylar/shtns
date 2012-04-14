/* time different variants of fft */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
// FFTW la derivee d/dx = ik	(pas de moins !)
#include <fftw3.h>


int main(int argc, char *argv[])
{
	double *Sh, *ShF;
	fftw_iodim dims[3];
	fftw_iodim many[3];
	fftw_plan fft, ifft;
	double t,tref;

	int fftw_plan_mode = FFTW_EXHAUSTIVE;
	int NPHI = 128;
	int NLAT = 142;

	if (argc > 1) {
		sscanf(argv[1], "%dx%d", &NLAT, &NPHI);
	}
	printf("planning %d FFT of size %d ...\n",NLAT, NPHI);
	tref = 1e200;

	int ncplx = NPHI/2 +1;
	ncplx += ncplx & 1;		// even.
	// Allocate dummy Spatial Fields.
	ShF = (double *) fftw_malloc(ncplx * NLAT * sizeof(complex double));
	Sh = (double *) fftw_malloc(ncplx * NLAT * sizeof(complex double));

	dims[0].n = NPHI;		many[0].n = NLAT;

#ifdef DONT
/* split */
	dims[0].is = 2*NLAT;	dims[0].os = NLAT;
	many[0].is = 1;	many[0].os = 1;
	ifft = fftw_plan_guru_split_dft_c2r( 1, dims, 1, many, (double *) ShF, ((double*)ShF)+ncplx, Sh, fftw_plan_mode);
	printf("\n split (native) : ifft=%g", fftw_cost(ifft));
	dims[0].is = NLAT;	dims[0].os = 2*NLAT;
	many[0].is = 1;	many[0].os = 1;
	fft = fftw_plan_guru_split_dft_r2c( 1, dims, 1, many, Sh, (double *) ShF, ((double*)ShF)+ncplx, fftw_plan_mode);
	printf(" fft=%g", fftw_cost(fft));
	t = fftw_cost(fft) + fftw_cost(ifft);
	if (t < tref) {
		tref = t;		printf("  ***");
	}

	dims[0].is = 2*NLAT;	dims[0].os = 1;
	many[0].is = 1;	many[0].os = NPHI;
	ifft = fftw_plan_guru_split_dft_c2r( 1, dims, 1, many, (double *) ShF, ((double*)ShF)+ncplx, Sh, fftw_plan_mode);
	printf("\n split (transpose) : ifft=%g", fftw_cost(ifft));
	dims[0].is = 1;	dims[0].os = 2*NLAT;
	many[0].is = NPHI;	many[0].os = 1;
	fft = fftw_plan_guru_split_dft_r2c( 1, dims, 1, many, Sh, (double *) ShF, ((double*)ShF)+ncplx, fftw_plan_mode);
	printf(" fft=%g", fftw_cost(fft));
	t = fftw_cost(fft) + fftw_cost(ifft);
	if (t < tref) {
		tref = t;		printf("  ***");
	}

	dims[0].is = 1;	dims[0].os = 1;
	many[0].is = 2*ncplx;	many[0].os = NPHI;
	ifft = fftw_plan_guru_split_dft_c2r( 1, dims, 1, many, (double *) ShF, ((double*)ShF)+ncplx, Sh, fftw_plan_mode);
	printf("\n split (phi-conti) : ifft=%g", fftw_cost(ifft));
	dims[0].is = 1;	dims[0].os = 1;
	many[0].is = NPHI;	many[0].os = 2*ncplx;
	fft = fftw_plan_guru_split_dft_r2c( 1, dims, 1, many, Sh, (double *) ShF, ((double*)ShF)+ncplx, fftw_plan_mode);
	printf(" fft=%g", fftw_cost(fft));
	t = fftw_cost(fft) + fftw_cost(ifft);
	if (t < tref) {
		tref = t;		printf("  ***");
	}

/* interleaved */
	dims[0].is = NLAT;	dims[0].os = NLAT;
	many[0].is = 1;	many[0].os = 1;
	ifft = fftw_plan_guru_dft_c2r( 1, dims, 1, many, (complex double *) ShF, Sh, fftw_plan_mode);
	printf("\n inter (native) : ifft=%g", fftw_cost(ifft));
	dims[0].is = NLAT;	dims[0].os = NLAT;
	many[0].is = 1;	many[0].os = 1;
	fft = fftw_plan_guru_dft_r2c( 1, dims, 1, many, Sh, (complex double *) ShF, fftw_plan_mode);
	printf(" fft=%g", fftw_cost(fft));
	t = fftw_cost(fft) + fftw_cost(ifft);
	if (t < tref) {
		tref = t;		printf("  ***");
	}

	dims[0].is = NLAT;	dims[0].os = 1;
	many[0].is = 1;	many[0].os = NPHI;
	ifft = fftw_plan_guru_dft_c2r( 1, dims, 1, many, (complex double *) ShF, Sh, fftw_plan_mode);
	printf("\n inter (transpose) : ifft=%g", fftw_cost(ifft));
	dims[0].is = 1;	dims[0].os = NLAT;
	many[0].is = NPHI;	many[0].os = 1;
	fft = fftw_plan_guru_dft_r2c( 1, dims, 1, many, Sh, (complex double *) ShF, fftw_plan_mode);
	printf(" fft=%g", fftw_cost(fft));
	t = fftw_cost(fft) + fftw_cost(ifft);
	if (t < tref) {
		tref = t;		printf("  ***");
	}

	dims[0].is = 1;	dims[0].os = 1;
	many[0].is = ncplx;	many[0].os = NPHI;
	ifft = fftw_plan_guru_dft_c2r( 1, dims, 1, many, (complex double *) ShF, Sh, fftw_plan_mode);
	printf("\n inter (phi-conti) : ifft=%g", fftw_cost(ifft));
	dims[0].is = 1;	dims[0].os = 1;
	many[0].is = NPHI;	many[0].os = ncplx;
	fft = fftw_plan_guru_dft_r2c( 1, dims, 1, many, Sh, (complex double *) ShF, fftw_plan_mode);
	printf(" fft=%g", fftw_cost(fft));
	t = fftw_cost(fft) + fftw_cost(ifft);
	if (t < tref) {
		tref = t;		printf("  ***");
	}

/* in-place */
	dims[0].is = NLAT;	dims[0].os = NLAT;
	many[0].is = 1;	many[0].os = 1;
	ifft = fftw_plan_guru_dft_c2r( 1, dims, 1, many, (complex double *) Sh, Sh, fftw_plan_mode);
	printf("\n cost (in-place) : ifft=%g", fftw_cost(ifft));
	dims[0].is = NLAT;	dims[0].os = NLAT;
	many[0].is = 1;	many[0].os = 1;
	fft = fftw_plan_guru_dft_r2c( 1, dims, 1, many, Sh, (complex double *) Sh, fftw_plan_mode);
	printf(" fft=%g", fftw_cost(fft));
	t = fftw_cost(fft) + fftw_cost(ifft);
	if (t < tref) {
		tref = t;		printf("  ***");
	}
#endif

/* complex in-place */
	dims[0].n = NPHI;		many[0].n = NLAT/2;
	dims[0].is = NLAT/2;	dims[0].os = NLAT/2;
	many[0].is = 1;			many[0].os = 1;
	ifft = fftw_plan_guru_dft( 1, dims, 1, many, (complex double *) Sh, (complex double *) Sh, FFTW_BACKWARD, fftw_plan_mode);
	printf("\n complex (in-place) : ifft=%g", fftw_cost(ifft));
	fft = fftw_plan_guru_dft( 1, dims, 1, many, (complex double *) Sh, (complex double *) Sh, FFTW_FORWARD, fftw_plan_mode);
	printf(" fft=%g", fftw_cost(fft));
	t = fftw_cost(fft) + fftw_cost(ifft);
	if (t < tref) {
		tref = t;		printf("  ***");
	}
	ifft = fftw_plan_guru_dft( 1, dims, 1, many, (complex double *) ShF, (complex double *) Sh, FFTW_BACKWARD, fftw_plan_mode);
	printf("\n complex (oop) : ifft=%g", fftw_cost(ifft));
	fft = fftw_plan_guru_dft( 1, dims, 1, many, (complex double *) Sh, (complex double *) ShF, FFTW_FORWARD, fftw_plan_mode);
	printf(" fft=%g", fftw_cost(fft));
	t = fftw_cost(fft) + fftw_cost(ifft);
	if (t < tref) {
		tref = t;		printf("  ***");
	}

/* phi-first in-place */
	dims[0].n = NPHI;		many[0].n = NLAT/2;
	dims[0].is = 1;			dims[0].os = 1;
	many[0].is = NPHI;		many[0].os = NPHI;
	ifft = fftw_plan_guru_dft( 1, dims, 1, many, (complex double *) Sh, (complex double *) Sh, FFTW_BACKWARD, fftw_plan_mode);
	printf("\n ref complex phi-first : ifft=%g", fftw_cost(ifft));
	fft = fftw_plan_guru_dft( 1, dims, 1, many, (complex double *) Sh, (complex double *) Sh, FFTW_FORWARD, fftw_plan_mode);
	printf(" fft=%g", fftw_cost(fft));
	t = fftw_cost(fft) + fftw_cost(ifft);
	if (t < tref) {
		tref = t;		printf("  ***");
	}
	ifft = fftw_plan_guru_dft( 1, dims, 1, many, (complex double *) ShF, (complex double *) Sh, FFTW_BACKWARD, fftw_plan_mode);
	printf("\n ref complex phi-first oop : ifft=%g", fftw_cost(ifft));
	fft = fftw_plan_guru_dft( 1, dims, 1, many, (complex double *) Sh, (complex double *) ShF, FFTW_FORWARD, fftw_plan_mode);
	printf(" fft=%g", fftw_cost(fft));
	t = fftw_cost(fft) + fftw_cost(ifft);
	if (t < tref) {
		tref = t;		printf("  ***");
	}

	printf("\n");
}
