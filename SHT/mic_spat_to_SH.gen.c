/*
 * Copyright (c) 2010-2019 Centre National de la Recherche Scientifique.
 * written by Nathanael Schaeffer (CNRS, ISTerre, Grenoble, France).
 * 
 * nathanael.schaeffer@univ-grenoble-alpes.fr
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

# This file is meta-code for SHT.c (spherical harmonic transform).
# it is intended for "make" to generate C code for 3 similar SHT functions,
# (namely spat_to_SH [Q tag]), spat_to_SHsphtor [V tag], spat_to_SH3 [both Q&V tags])
# from one generic function + tags.
# Basically, there are tags at the beginning of lines (Q,V) that are information
# to keep or remove the line depending on the function to build. (Q for scalar, V for vector, # for comment)
#
//////////////////////////////////////////////////

	#define HI_LLIM
QX	#define BASE _an1_hi
QX	#include "spat_to_SH_kernel.c"
VX	#define BASE _an2_hi
VX	#include "spat_to_SHst_kernel.c"
3	#define BASE _an3_hi
3	#include "spat_to_SHqst_kernel.c"
	#undef BASE
	
	#undef HI_LLIM
QX	#define BASE _an1
QX	#include "spat_to_SH_kernel.c"
VX	#define BASE _an2
VX	#include "spat_to_SHst_kernel.c"
3	#define BASE _an3
3	#include "spat_to_SHqst_kernel.c"
	#undef BASE


	static
QX	void GEN3(spat_to_SH_mic,NWAY,SUFFIX)(shtns_cfg shtns, double *Vr, cplx *Qlm, long int llim) {
VX	void GEN3(spat_to_SHsphtor_mic,NWAY,SUFFIX)(shtns_cfg shtns, double *Vt, double *Vp, cplx *Slm, cplx *Tlm, long int llim) {
3	void GEN3(spat_to_SHqst_mic,NWAY,SUFFIX)(shtns_cfg shtns, double *Vr, double *Vt, double *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm, long int llim) {

Q	double *BrF;		// contains the Fourier transformed data
V	double *BtF, *BpF;	// contains the Fourier transformed data
	unsigned imlim=0;

Q	BrF = Vr;
V	BtF = Vt;	BpF = Vp;
  #ifndef SHT_AXISYM
	imlim = MTR;
	#ifdef SHT_VAR_LTR
		if (imlim*MRES > (unsigned) llim) imlim = ((unsigned) llim)/MRES;		// 32bit mul and div should be faster
	#endif
	if (shtns->fftc_mode >= 0) {
		if (shtns->fftc_mode > 0) {		// alloc memory for out-of-place FFT
			unsigned long nv = shtns->nspat;
QX			BrF = (double*) VMALLOC( nv * sizeof(double) );
VX			BtF = (double*) VMALLOC( 2*nv * sizeof(double) );
VX			BpF = BtF + nv;
3			BrF = (double*) VMALLOC( 3*nv * sizeof(double) );
3			BtF = BrF + nv;		BpF = BtF + nv;
		}
	    if (shtns->fftc_mode != 1) {	// regular FFT
Q			fftw_execute_dft(shtns->fftc,(cplx*)Vr, (cplx*)BrF);
V			fftw_execute_dft(shtns->fftc,(cplx*)Vt, (cplx*)BtF);
V			fftw_execute_dft(shtns->fftc,(cplx*)Vp, (cplx*)BpF);
		} else {	// split FFT
Q			fftw_execute_split_dft(shtns->fftc, Vr+NPHI, Vr, BrF+1, BrF);
V			fftw_execute_split_dft(shtns->fftc, Vt+NPHI, Vt, BtF+1, BtF);
V			fftw_execute_split_dft(shtns->fftc, Vp+NPHI, Vp, BpF+1, BpF);
	    }
	}
  #endif
	imlim += 1;

	if (llim < SHT_L_RESCALE_FLY) {
		#pragma omp parallel num_threads(shtns->nthreads)
		{
QX			GEN3(_an1,NWAY,SUFFIX)(shtns, BrF, Qlm, llim, imlim);
VX			GEN3(_an2,NWAY,SUFFIX)(shtns, BtF, BpF, Slm, Tlm, llim, imlim);
3			GEN3(_an3,NWAY,SUFFIX)(shtns, BrF, BtF, BpF, Qlm, Slm, Tlm, llim, imlim);
		}
	} else {
		#pragma omp parallel num_threads(shtns->nthreads)
		{
QX			GEN3(_an1_hi,NWAY,SUFFIX)(shtns, BrF, Qlm, llim, imlim);
VX			GEN3(_an2_hi,NWAY,SUFFIX)(shtns, BtF, BpF, Slm, Tlm, llim, imlim);
3			GEN3(_an3_hi,NWAY,SUFFIX)(shtns, BrF, BtF, BpF, Qlm, Slm, Tlm, llim, imlim);
		}
	}

  #ifndef SHT_AXISYM
  	if (shtns->fftc_mode > 0) {		// free memory
Q	    VFREE(BrF);
VX	    VFREE(BtF);	// this frees also BpF.
	}
  #endif

  }

	#undef LSPAN
