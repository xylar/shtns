/*
 * Copyright (c) 2010 Centre National de la Recherche Scientifique.
 * written by Nathanael Schaeffer (CNRS, LGIT, Grenoble, France).
 * 
 * nathanael.schaeffer@ujf-grenoble.fr
 * 
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 * 
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * 
 */

/* shtns.i : SWIG interface to Python*/

%module (docstring="Python interface to the SHTns spherical harmonic library") shtns

%{

#include "sht_private.h"

%}


/* Py_complex is the equivalent of "complex double" */
%typemap(in) Py_complex {
        $1 = PyComplex_AsCComplex($input);
}
%typemap(out) Py_complex {
        $result = PyComplex_FromCComplex($1);
}

/* Helper functions to allocate and access sh coefficients arrays */
%inline %{
  Py_complex ylm_get(Py_complex *y, int l, int m) {
    return( y[LM(l,m)] );
  }
  void ylm_set(Py_complex *y, int l, int m, Py_complex val) {
    y[LM(l,m)] = val;
  }
  Py_complex *ylm_alloc() {
     int i;
     complex double *y = (complex double *) malloc(sizeof(complex double) * shtns.nlm);
     for (i=0;i<shtns.nlm;i++) y[i] = 0.0;
     return (Py_complex *) y;
  }
%}

int shtns_set_size(int lmax, int mmax, int mres, int norm);
int shtns_precompute(int flags, double eps, int nlat, int nphi);
int shtns_init(int flags, int lmax, int mmax, int mres, int nlat, int nphi);
long int nlm_calc(long int lmax, long int mmax, long int mres);
void Set_MTR_DCT(int m);
int Get_MTR_DCT();

double sh00_1();	///< returns the spherical harmonic representation of 1 (l=0,m=0)
double sh10_ct();	///< returns the spherical harmonic representation of cos(theta) (l=1,m=0)
double sh11_st();	///< returns the spherical harmonic representation of sin(theta)*cos(phi) (l=1,m=1)

void spat_to_SH(double *Vr, complex double *Qlm);
void SH_to_spat(complex double *Qlm, double *Vr);

void spat_to_SHsphtor(double *Vt, double *Vp, complex double *Slm, complex double *Tlm);
void SHsphtor_to_spat(complex double *Slm, complex double *Tlm, double *Vt, double *Vp);
void SHsph_to_spat(complex double *Slm, double *Vt, double *Vp);
void SHtor_to_spat(complex double *Tlm, double *Vt, double *Vp);

double SH_to_point(complex double *Qlm, double cost, double phi);
void SHqst_to_point(complex double *Qlm, complex double *Slm, complex double *Tlm,
					double cost, double phi, double *vr, double *vt, double *vp);
void SHqst_to_lat(complex double *Qlm, complex double *Slm, complex double *Tlm, double cost,
					double *vr, double *vt, double *vp, int nphi, int ltr, int mtr);

void spat_to_SH_l(double *Vr, complex double *Qlm, int LTR);
void SH_to_spat_l(complex double *Qlm, double *Vr, int LTR);

void SHsphtor_to_spat_l(complex double *Slm, complex double *Tlm, double *Vt, double *Vp, int LTR);
void SHsph_to_spat_l(complex double *Slm, double *Vt, double *Vp, int LTR);
void SHtor_to_spat_l(complex double *Tlm, double *Vt, double *Vp, int LTR);
void spat_to_SHsphtor_l(double *Vt, double *Vp, complex double *Slm, complex double *Tlm, int LTR);

void spat_to_SH_m0(double *Vr, complex double *Qlm);
void SH_to_spat_m0(complex double *Qlm, double *Vr);
void SHsphtor_to_spat_m0(complex double *Slm, complex double *Tlm, double *Vt, double *Vp);
void SHsph_to_spat_m0(complex double *Slm, double *Vt);
void SHtor_to_spat_m0(complex double *Tlm, double *Vp);
void spat_to_SHsphtor_m0(double *Vt, double *Vp, complex double *Slm, complex double *Tlm);

void SHeo_to_spat(complex double *Qlm, double *Vr, int parity);
void spat_to_SHeo(double *Vr, complex double *Qlm, int parity);
void SHeo_sphtor_to_spat(complex double *Slm, complex double *Tlm, double *Vt, double *Vp, int parity);
void spat_to_SHeo_sphtor(double *Vt, double *Vp, complex double *Slm, complex double *Tlm, int parity);

