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

/// \file time_SHT.c This program performs some spherical harmonic transforms, and displays timings and accuracy.
/// \c make \c time_SHT to compile, and then run.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>

#include "shtns.h"

shtns_cfg shtns;

complex double *Slm, *Slm0, *Tlm, *Tlm0, *Qlm;	// spherical harmonics l,m space
complex double *ShF, *ThF, *NLF;	// Fourier space : theta,m
double *Sh, *Th, *NL;		// real space : theta,phi (alias of ShF)

int LMAX,MMAX,MRES,NLM;
int NLAT = 0;
int NPHI = 0;

int main(int argc, char *argv[])
{
	complex double t1, t2;
	double t,tmax,n2;
	int i,im,m,l,jj;
	double e0,e1;
	double polaropt = 1.e-8;		// default for polar optimization.
	enum shtns_type shtmode = sht_auto;		// default to "auto" (fastest) mode.
	enum shtns_norm shtnorm = sht_orthonormal;		// default to "orthonormal" SH.
	int layout = SHT_NATIVE_LAYOUT;
	int nlorder = 0;
	int vector = 1;
	char name[20];

	srand( time(NULL) );	// initialise les nombres.
	
	MMAX=LMAX=2;
	MRES=1;
	NLAT=32;

	shtns = shtns_create(LMAX, MMAX, MRES, shtnorm);
	NLM = shtns->nlm;
	shtns_set_grid_auto(shtns, sht_gauss_fly, 0.0, 1e-10, &NLAT, &NPHI);
	shtns_print_cfg(shtns);

	ShF = (complex double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Sh = (double *) ShF;
	ThF = (complex double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Th = (double *) ThF;
	NLF = (complex double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	NL = (double *) NLF;

	Tlm0 = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Slm0 = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Slm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Tlm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Qlm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);


	for (l=0; l<NLM; l++) {
		Qlm[l] = 0.0;
	}
	Qlm[0] = 1.0e3;
//	Qlm[LiM(shtns, 1, 0)] = 10.0;
	Qlm[LiM(shtns, 1, 1)] = 1.0;	//0.1 + I*0.05;
	
	LMAX=1;
	SH_Yrotate90(shtns, Qlm, Slm, LMAX);

	for (l=0; l<=LMAX; l++) {
		double norm0=0;		double normR=0;
		for (m=0; m<=l; m++) {
			printf("l=%d, m=%d,  Q=%f,%f,  \t  S=%f,%f\n",l,m, creal(Qlm[LiM(shtns,l,m)]), cimag(Qlm[LiM(shtns,l,m)]),
			creal(Slm[LiM(shtns,l,m)]), cimag(Slm[LiM(shtns,l,m)]));
			double e0 = creal(Qlm[LiM(shtns,l,m)])*creal(Qlm[LiM(shtns,l,m)]) + cimag(Qlm[LiM(shtns,l,m)])*cimag(Qlm[LiM(shtns,l,m)]);
			double eR = creal(Slm[LiM(shtns,l,m)])*creal(Slm[LiM(shtns,l,m)]) + cimag(Slm[LiM(shtns,l,m)])*cimag(Slm[LiM(shtns,l,m)]);			
			if (m>0) {
				e0 *=2;		eR *= 2;
			}
			norm0 += e0;
			normR += eR;
		}
		printf("norm0 = %f, norm1 = %f\n", sqrt(norm0), sqrt(normR));
	}
		

	SH_Yrotate90(shtns, Slm, Tlm0, LMAX);
	SH_Yrotate90(shtns, Tlm0, Slm, LMAX);
	SH_Yrotate90(shtns, Slm, Tlm, LMAX);

	for (l=0; l<=LMAX; l++) {
		double norm0=0;		double normR=0;
		for (m=0; m<=l; m++) {
			printf("l=%d, m=%d,  Q=%f,%f,  \t  S=%f,%f\n",l,m, creal(Qlm[LiM(shtns,l,m)]), cimag(Qlm[LiM(shtns,l,m)]),
			creal(Tlm[LiM(shtns,l,m)]), cimag(Tlm[LiM(shtns,l,m)]));
			double e0 = creal(Qlm[LiM(shtns,l,m)])*creal(Qlm[LiM(shtns,l,m)]) + cimag(Qlm[LiM(shtns,l,m)])*cimag(Qlm[LiM(shtns,l,m)]);
			double eR = creal(Tlm[LiM(shtns,l,m)])*creal(Tlm[LiM(shtns,l,m)]) + cimag(Tlm[LiM(shtns,l,m)])*cimag(Tlm[LiM(shtns,l,m)]);			
			if (m>0) {
				e0 *=2;		eR *= 2;
			}
			norm0 += e0;
			normR += eR;
		}
		printf("norm0 = %f, norm1 = %f\n", sqrt(norm0), sqrt(normR));
	}

}

