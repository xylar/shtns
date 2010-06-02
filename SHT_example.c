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

/** \example SHT_example.c
 \brief An example program that performs backward and forward Spherical Harmonic Transforms using SHTns.

  Compile using : \code make SHT_example \endcode
  \see \ref usage
**/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <fftw3.h>

#include <shtns.h>

/// a simple function that writes a vector to a file
void write_vect(char *fn, double *vec, int N);

/// a simple function that writes a matrix to a file.
void write_mx(char *fn, double *mx, int N1, int N2);


int main()
{
	long int lmax,mmax,nlat,nphi,mres, NLM;
	complex double *Slm, *Tlm;	// spherical harmonics l,m space
	double *Sh, *Th;		// real space : theta,phi
	long int i,im,lm;
	double t;

	lmax = 5;	nlat = 32;
	mmax = 3;	nphi = 10;
	mres = 1;
	NLM = shtns_init( sht_gauss, lmax, mmax, mres, nlat, nphi );
//	NLM = shtns_set_size(lmax, mmax, mres, sht_orthonormal | SHT_REAL_NORM);
//	shtns_precompute( sht_gauss, 0.0, nlat, nphi);

// Memory allocation : the use of fftw_malloc is required because we need proper 16-byte alignement.
// allocate spatial fields.
	Sh = (double *) fftw_malloc( NSPAT_ALLOC * sizeof(double));
	Th = (double *) fftw_malloc( NSPAT_ALLOC * sizeof(double));

// allocate SH representations.
	Slm = (complex double *) fftw_malloc( NLM * sizeof(complex double));
	Tlm = (complex double *) fftw_malloc( NLM * sizeof(complex double));

// SH_to_spat
	LM_LOOP( Slm[lm]=0.0;  Tlm[lm] = 0.0; )		/* this is the same as :
						for (lm=0;lm<NLM;lm++) {
							Slm[lm] = 0.0;	Tlm[lm] = 0.0;
						} */

//	Slm[LM(1,1)] = sh11_st();				// access to SH coefficient
	Slm[LM(2,0)] = 1.0;
// 	Slm[LiM(1,0)] = sh10_ct();
//	Slm[LiM(0,0)] = 0.5*sh00_1();
//	SH_to_spat(Slm,Sh);
	SHtor_to_spat_m0(Slm,Sh);
	write_vect("ylm",(double *) Slm,NLM*2);
	write_mx("spat",Sh,nphi,nlat);

// compute value of SH expansion at a given physical point.
 	double t2;
	SHqst_to_point(Tlm, Tlm, Slm, ct[nlat/3], 2.*M_PI/(mres*nphi),&t2,&t2,&t);
	printf("check if SH_to_point coincides with SH_to_spat : %f = %f\n",t,Sh[nlat/3]);
	printf("ct*st = %f\n",ct[nlat/3]*st[nlat/3]);

// check non-linear behaviour
	for (im=0;im<nphi;im++) {
		for (i=0;i<nlat;i++) {
			Sh[im*nlat+i] *= Sh[im*nlat+i];
		}
	}
	spat_to_SH(Sh,Tlm);
	write_vect("ylm_nl",(double *) Tlm,NLM*2);

// vector transform
	LM_LOOP( Slm[lm]=0.0;  Tlm[lm] = 0.0; )
	SHsphtor_to_spat(Slm,Tlm,Sh,Th);		// vector transform
	write_mx("spatt",Sh,nphi,nlat);
	write_mx("spatp",Th,nphi,nlat);

// spat_to_SH
	for (im=0;im<nphi;im++) {
		for (i=0;i<nlat;i++) {
			Sh[im*nlat+i] = ct[i];			// use cos(theta) array
		}
	}
	spat_to_SH(Sh,Slm);
	write_vect("ylm_v",(double *) Slm,NLM*2);
}


void write_vect(char *fn, double *vec, int N)
{
	FILE *fp;
	int i;
	
	fp = fopen(fn,"w");
	for (i=0;i<N;i++) {
		fprintf(fp,"%.6g ",vec[i]);
	}
	fclose(fp);
}

void write_mx(char *fn, double *mx, int N1, int N2)
{
	FILE *fp;
	int i,j;
	
	fp = fopen(fn,"w");
	for (i=0;i<N1;i++) {
		for(j=0;j<N2;j++) {
			fprintf(fp,"%.6g ",mx[i*N2+j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

