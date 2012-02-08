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

/* shtns.i : SWIG interface to Python */

%module (docstring="Python interface to the SHTns spherical harmonic transform library") shtns

%{

#include "sht_private.h"

%}

%include "shtns.h"

%feature("autodoc");

long nlm_calc(long lmax, long mmax, long mres);


%extend shtns_info {
	shtns_info(int lmax, int mmax, int mres=1, int norm=sht_orthonormal) {	// default arguments : mres and norm
		return shtns_create(lmax, mmax, mres, norm);
	}
	~shtns_info() {
		shtns_destroy($self);		// free memory.
	}
	int set_grid(int nlat, int nphi, enum shtns_type flags=sht_quick_init, double eps=1.0e-8) {	// default arguments
		return shtns_set_grid($self, flags, eps, nlat, nphi);
	}
	int set_grid_auto(int nlat=0, int nphi=0, int nl_order=1, enum shtns_type flags=sht_quick_init, double eps=1.0e-8) {
		return shtns_set_grid_auto($self, flags, eps, nl_order, &nlat, &nphi);
	}

	void print() {
		shtns_print_cfg($self);
	}
	double sh00_1() {
		return sh00_1($self);
	}
	double sh10_ct() {
		return sh10_ct($self);
	}
	double sh11_st() {
		return sh11_st($self);
	}
	double shlm_e1(int l, int m) {
		return shlm_e1($self, l, m);
	}

/* Helper functions to allocate and access sh coefficients arrays */
	PyObject* ylm_get(Py_complex *y, int l, int m) {
		return PyComplex_FromCComplex(y[LM($self, l,m)]);
	}
	void ylm_set(Py_complex *y, int l, int m, PyObject *val) {
		y[LM($self, l,m)] = PyComplex_AsCComplex(val);
	}
	Py_complex *ylm_alloc() {
		int i;
		complex double *y = (complex double *) malloc(sizeof(complex double) * $self->nlm);
		for (i=0; i<$self->nlm; i++) y[i] = 0.0;
		return (Py_complex *) y;
	}

	/* scalar transforms */
	void spat_to_SH(double *Vr, complex double *Qlm) {
		spat_to_SH($self, Vr, Qlm);	}

	void SH_to_spat(complex double *Qlm, double *Vr) {
		SH_to_spat($self, Qlm, Vr);	}

	/* 2D vectors */
	void spat_to_SHsphtor(double *Vt, double *Vp, complex double *Slm, complex double *Tlm) {
		spat_to_SHsphtor($self, Vt, Vp, Slm, Tlm);	}
	void SHsphtor_to_spat(complex double *Slm, complex double *Tlm, double *Vt, double *Vp) {
		SHsphtor_to_spat($self, Slm, Tlm, Vt, Vp);	}
	void SHsph_to_spat(complex double *Slm, double *Vt, double *Vp) {
		SHsph_to_spat($self, Slm, Vt, Vp);	}
	void SHtor_to_spat(complex double *Tlm, double *Vt, double *Vp) {
		SHtor_to_spat($self, Tlm, Vt, Vp);	}

	/* 3D vectors */
	void spat_to_SHqst(double *Vr, double *Vt, double *Vp, complex double *Qlm, complex double *Slm, complex double *Tlm) {
		spat_to_SHqst($self, Vr, Vt, Vp, Qlm, Slm, Tlm);	}
	void SHqst_to_spat(complex double *Qlm, complex double *Slm, complex double *Tlm, double *Vr, double *Vt, double *Vp) {
		SHqst_to_spat($self, Qlm, Slm, Tlm, Vr, Vt, Vp);	}

	/* local evaluations */
	double SH_to_point(complex double *Qlm, double cost, double phi) {
		SH_to_point($self, Qlm, cost, phi);	}
	void SHqst_to_point(complex double *Qlm, complex double *Slm, complex double *Tlm,
					double cost, double phi, double *vr, double *vt, double *vp) {
		SHqst_to_point($self, Qlm, Slm, Tlm, cost, phi, vr, vt, vp);	}

	/* rotation of SH representations (experimental) */
	void SH_Zrotate(complex double *Qlm, double alpha, complex double *Rlm) {
		SH_Zrotate($self, Qlm, alpha, Rlm);	}
	void SH_Yrotate(complex double *Qlm, double alpha, complex double *Rlm) {
		SH_Yrotate($self, Qlm, alpha, Rlm);	}
	void SH_Yrotate90(complex double *Qlm, complex double *Rlm) {
		SH_Yrotate90($self, Qlm, Rlm);	}
	void SH_Xrotate90(complex double *Qlm, complex double *Rlm) {
		SH_Xrotate90($self, Qlm, Rlm);
	}

};


