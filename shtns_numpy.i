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

/* shtns.i : SWIG interface to Python using NumPy */

/* TODO and known problems :
 * - exception handling does not work well : some functions will still return incorect values when parameters are wrong
 * - alignement on 16 bytes of NumPy arrays is not guaranteed. It should however work on 64bit systems or on modern 32bit systems.
 */

%module (docstring="Python/NumPy interface to the SHTns spherical harmonic transform library") shtns

%{
	
#include "sht_private.h"
#include "/usr/lib/python2.7/site-packages/numpy/core/include/numpy/arrayobject.h"

int check_spatial(PyObject *a, int size) {
	if (size == 0) {
		PyErr_SetString(PyExc_RuntimeError,"grid not set");
		return 0;
	}
	if (!PyArray_Check(a)) {
		PyErr_SetString(PyExc_RuntimeError,"NumPy array expected");
		return 0;
	}
	if (PyArray_TYPE(a) != PyArray_DOUBLE) {
		PyErr_SetString(PyExc_RuntimeError,"array must consist of float");
		return 0;
	}
	if (!PyArray_ISCONTIGUOUS(a)) {
		PyErr_SetString(PyExc_RuntimeError,"array not contiguous");
		return 0;
	}
	if (PyArray_SIZE(a) < size) {
		PyErr_SetString(PyExc_RuntimeError,"array too small");
		return 0;
	}
	return 1;
}

int check_spectral(PyObject *a, int size) {
	if (!PyArray_Check(a)) {
		PyErr_SetString(PyExc_RuntimeError,"NumPy array expected");
		return 0;
	}
	if (PyArray_TYPE(a) != PyArray_CDOUBLE) {
		PyErr_SetString(PyExc_RuntimeError,"array must consist of complex float");
		return 0;
	}
	if (!PyArray_ISCONTIGUOUS(a)) {
		PyErr_SetString(PyExc_RuntimeError,"array not contiguous");
		return 0;
	}
	if (PyArray_SIZE(a) < size) {
		PyErr_SetString(PyExc_RuntimeError,"array too small");
		return 0;
	}
	return 1;
}

%}

// main object is renamed to sht.
%rename("sht") shtns_info;

%include "shtns.h"

%feature("autodoc");

long nlm_calc(long lmax, long mmax, long mres);


%extend shtns_info {
	shtns_info(int lmax, int mmax, int mres=1, int norm=sht_orthonormal) {	// default arguments : mres and norm
		import_array();
		return shtns_create(lmax, mmax, mres, norm);
	}
	~shtns_info() {
		shtns_destroy($self);		// free memory.
	}
	int set_grid(int nlat, int nphi, enum shtns_type flags=sht_quick_init|SHT_THETA_CONTIGUOUS, double eps=1.0e-8) {	// default arguments
		return shtns_set_grid($self, flags, eps, nlat, nphi);
	}
	int set_grid_auto(int nlat=0, int nphi=0, int nl_order=1, enum shtns_type flags=sht_quick_init|SHT_THETA_CONTIGUOUS, double eps=1.0e-8) {
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

	/* returns useful data */
	PyObject* l() {
		int i;
		npy_intp dims = $self->nlm;
		npy_intp strides = sizeof(double);
		PyObject *obj = PyArray_New(&PyArray_Type, 1, &dims, PyArray_DOUBLE, &strides, NULL, strides, 0, NULL);
		double *el = (double*) PyArray_DATA(obj);
		for (i=0; i<$self->nlm; i++)		el[i] = $self->li[i];		// convert and copy
		return obj;
	}
	PyObject* cos_theta() {		// grid must have been initialized.
		int i;
		npy_intp dims = $self->nlat;
		npy_intp strides = sizeof(double);
		if (dims == 0) {	// no grid
			PyErr_SetString(PyExc_RuntimeError,"grid not set");
			return NULL;
		}
		PyObject *obj = PyArray_New(&PyArray_Type, 1, &dims, PyArray_DOUBLE, &strides, NULL, strides, 0, NULL);
		double *ct = (double*) PyArray_DATA(obj);
		for (i=0; i<$self->nlat; i++)		ct[i] = $self->ct[i];		// copy
		return obj;
	}

	// returns the index in a spectral array of (l,m) coefficient.
	int idx(int l, int m) {
		if ( (l < 0) || (l > $self->lmax) || (m < 0) || (m > $self->mmax * $self->mres) || (m % $self->mres != 0) ) {
			PyErr_SetString(PyExc_IndexError,"l or m index out-of-bounds");
			return 0;
		}
		return LM($self, l, m);
	}

	/* scalar transforms */
	void spat_to_SH(PyObject *Vr, PyObject *Qlm) {
		if (check_spatial(Vr, $self->nspat) && check_spectral(Qlm, $self->nlm))
			spat_to_SH($self, PyArray_DATA(Vr), PyArray_DATA(Qlm));
	}
	void SH_to_spat(PyObject *Qlm, PyObject *Vr) {
		if (check_spatial(Vr, $self->nspat) && check_spectral(Qlm, $self->nlm))
			SH_to_spat($self, PyArray_DATA(Qlm), PyArray_DATA(Vr));
	}

	/* 2D vectors */
	void spat_to_SHsphtor(PyObject *Vt, PyObject *Vp, PyObject *Slm, PyObject *Tlm) {
		if (check_spatial(Vt, $self->nspat) && check_spatial(Vp, $self->nspat) && check_spectral(Slm, $self->nlm) && check_spectral(Tlm, $self->nlm))
			spat_to_SHsphtor($self, PyArray_DATA(Vt), PyArray_DATA(Vp), PyArray_DATA(Slm), PyArray_DATA(Tlm));
	}
	void SHsphtor_to_spat(PyObject *Slm, PyObject *Tlm, PyObject *Vt, PyObject *Vp) {
		if (check_spatial(Vt, $self->nspat) && check_spatial(Vp, $self->nspat) && check_spectral(Slm, $self->nlm) && check_spectral(Tlm, $self->nlm))
			SHsphtor_to_spat($self, PyArray_DATA(Slm), PyArray_DATA(Tlm), PyArray_DATA(Vt), PyArray_DATA(Vp));
	}
	void SHsph_to_spat(PyObject *Slm, PyObject *Vt, PyObject *Vp) {
		if (check_spatial(Vt, $self->nspat) && check_spatial(Vp, $self->nspat) && check_spectral(Slm, $self->nlm))
		SHsph_to_spat($self, PyArray_DATA(Slm), PyArray_DATA(Vt), PyArray_DATA(Vp));
	}
	void SHtor_to_spat(PyObject *Tlm, PyObject *Vt, PyObject *Vp) {
		if (check_spatial(Vt, $self->nspat) && check_spatial(Vp, $self->nspat) && check_spectral(Tlm, $self->nlm))
		SHtor_to_spat($self, PyArray_DATA(Tlm), PyArray_DATA(Vt), PyArray_DATA(Vp));
	}

	/* 3D vectors */
	void spat_to_SHqst(PyObject *Vr, PyObject *Vt, PyObject *Vp, PyObject *Qlm, PyObject *Slm, PyObject *Tlm) {
		if (check_spatial(Vr, $self->nspat) && check_spatial(Vt, $self->nspat) && check_spatial(Vp, $self->nspat)
			&& check_spectral(Qlm, $self->nlm) && check_spectral(Slm, $self->nlm) && check_spectral(Tlm, $self->nlm))
		spat_to_SHqst($self, PyArray_DATA(Vr), PyArray_DATA(Vt), PyArray_DATA(Vp), PyArray_DATA(Qlm), PyArray_DATA(Slm), PyArray_DATA(Tlm));
	}
	void SHqst_to_spat(PyObject *Qlm, PyObject *Slm, PyObject *Tlm, PyObject *Vr, PyObject *Vt, PyObject *Vp) {
		if (check_spatial(Vr, $self->nspat) && check_spatial(Vt, $self->nspat) && check_spatial(Vp, $self->nspat)
			&& check_spectral(Qlm, $self->nlm) && check_spectral(Slm, $self->nlm) && check_spectral(Tlm, $self->nlm))
		SHqst_to_spat($self, PyArray_DATA(Qlm), PyArray_DATA(Slm), PyArray_DATA(Tlm), PyArray_DATA(Vr), PyArray_DATA(Vt), PyArray_DATA(Vp));
	}

	/* local evaluations */
	double SH_to_point(PyObject *Qlm, double cost, double phi) {
		if (check_spectral(Qlm, $self->nlm))	return SH_to_point($self, PyArray_DATA(Qlm), cost, phi);
	}
	%apply double *OUTPUT { double *vr };
	%apply double *OUTPUT { double *vt };
	%apply double *OUTPUT { double *vp };
	void SHqst_to_point(PyObject *Qlm, PyObject *Slm, PyObject *Tlm,
					double cost, double phi, double *vr, double *vt, double *vp) {
		if (check_spectral(Qlm, $self->nlm) && check_spectral(Slm, $self->nlm) && check_spectral(Tlm, $self->nlm))
			SHqst_to_point($self, PyArray_DATA(Qlm), PyArray_DATA(Slm), PyArray_DATA(Tlm), cost, phi, vr, vt, vp);
	}
	%clear double *vr;
	%clear double *vt;
	%clear double *vp;

	/* rotation of SH representations (experimental) */
	void SH_Zrotate(PyObject *Qlm, double alpha, PyObject *Rlm) {
		if (check_spectral(Qlm, $self->nlm) && check_spectral(Rlm, $self->nlm))
			SH_Zrotate($self, PyArray_DATA(Qlm), alpha, PyArray_DATA(Rlm));
	}
	void SH_Yrotate(PyObject *Qlm, double alpha, PyObject *Rlm) {
		if (check_spectral(Qlm, $self->nlm) && check_spectral(Rlm, $self->nlm))
			SH_Yrotate($self, PyArray_DATA(Qlm), alpha, PyArray_DATA(Rlm));
	}
	void SH_Yrotate90(PyObject *Qlm, PyObject *Rlm) {
		if (check_spectral(Qlm, $self->nlm) && check_spectral(Rlm, $self->nlm))
			SH_Yrotate90($self, PyArray_DATA(Qlm), PyArray_DATA(Rlm));
	}
	void SH_Xrotate90(PyObject *Qlm, PyObject *Rlm) {
		if (check_spectral(Qlm, $self->nlm) && check_spectral(Rlm, $self->nlm))
			SH_Xrotate90($self, PyArray_DATA(Qlm), PyArray_DATA(Rlm));
	}

};
