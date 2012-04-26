/*
 * Copyright (c) 2010-2012 Centre National de la Recherche Scientifique.
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

/* shtns_numpy.i : SWIG interface to Python using NumPy */

/* TODO and known problems :
 * - alignement on 16 bytes of NumPy arrays is not guaranteed. It should however work on 64bit systems or on modern 32bit systems.
 * - you may have to adjust the path below to include the header file "arrayobject.h" from the NumPy package.
 */

%module (docstring="Python/NumPy interface to the SHTns spherical harmonic transform library") shtns

%{
	
#include "sht_private.h"
#include "/usr/lib/python2.7/site-packages/numpy/core/include/numpy/arrayobject.h"

// variables used for exception handling.
static int shtns_error = 0;
static char* shtns_err_msg;
static char msg_buffer[128];
static char msg_grid_err[] = "Grid not set. Call .set_grid_auto() or .set_grid() mehtod.";
static char msg_numpy_arr[] = "Numpy array expected.";

static void throw_exception(int error, int iarg, char* msg)
{
	shtns_error = error;
	shtns_err_msg = msg;
	if (iarg > 0) {
		sprintf(msg_buffer, "arg #%d : %.100s", iarg, msg);
		shtns_err_msg = msg_buffer;
	}
}

static int check_spatial(int i, PyObject *a, int size) {
	if (size == 0) {
		throw_exception(SWIG_RuntimeError,0,msg_grid_err);	return 0;
	}
	if (!PyArray_Check(a)) {
		throw_exception(SWIG_TypeError,i,msg_numpy_arr);		return 0;
	}
	if (PyArray_TYPE(a) != PyArray_DOUBLE) {
		throw_exception(SWIG_TypeError,i,"spatial array must consist of float");		return 0;
	}
	if (!PyArray_ISCONTIGUOUS(a)) {
		throw_exception(SWIG_RuntimeError,i,"spatial array not contiguous. Use 'b=a.copy()' to copy a to a contiguous array b.");		return 0;
	}
	if (PyArray_SIZE(a) != size) {
		throw_exception(SWIG_RuntimeError,i,"spatial array has wrong size");		return 0;
	}
	return 1;
}

static int check_spectral(int i, PyObject *a, int size) {
	if (!PyArray_Check(a)) {
		throw_exception(SWIG_RuntimeError,i,msg_numpy_arr);		return 0;
	}
	if (PyArray_TYPE(a) != PyArray_CDOUBLE) {
		throw_exception(SWIG_RuntimeError,i,"spectral array must consist of complex float. Create with: 'numpy.zeros(sh.nlm, dtype=complex)'");		return 0;
	}
	if (!PyArray_ISCONTIGUOUS(a)) {
		throw_exception(SWIG_RuntimeError,i,"spactral array not contiguous. Use .copy() to copy to a contiguous array.");		return 0;
	}
	if (PyArray_SIZE(a) != size) {
		throw_exception(SWIG_RuntimeError,i,"spectral array has wrong size");		return 0;
	}
	return 1;
}

%}

// main object is renamed to sht.
%rename("sht") shtns_info;
%ignore SHT_NATIVE_LAYOUT;

%include "shtns.h"

%include "exception.i"
%feature("autodoc");

long nlm_calc(long lmax, long mmax, long mres);

%extend shtns_info {

	%exception {
		shtns_error = 0;	// clear exception
		$function
		if (shtns_error) {	// test for exception
			SWIG_exception(shtns_error, shtns_err_msg);		return NULL;
		}
	}

	shtns_info(int lmax, int mmax=-1, int mres=1, int norm=sht_orthonormal) {	// default arguments : mmax, mres and norm
		if (lmax < 2) {
			throw_exception(SWIG_ValueError,1,"lmax < 2 not allowed");	return NULL;
		}
		import_array();		// required by NumPy
		shtns_use_threads(0);		// use openmp threads if available
		if (mmax < 0) mmax = lmax;
		return shtns_create(lmax, mmax, mres, norm);
	}
	~shtns_info() {
		shtns_destroy($self);		// free memory.
	}
	%apply int *OUTPUT { int *nlat_out };
	%apply int *OUTPUT { int *nphi_out };
	void set_grid(int nlat=0, int nphi=0, enum shtns_type flags=sht_quick_init, double eps=1.0e-8, int nl_order=1, int *nlat_out, int *nphi_out) {	// default arguments
		if (!(flags & SHT_THETA_CONTIGUOUS))  flags |= SHT_PHI_CONTIGUOUS;	// default to SHT_PHI_CONTIGUOUS.
		*nlat_out = nlat;		*nphi_out = nphi;
		shtns_set_grid_auto($self, flags, eps, nl_order, nlat_out, nphi_out);
	}
	void set_grid_auto(int nl_order=1, enum shtns_type flags=sht_quick_init, double eps=1.0e-8, int *nlat_out, int *nphi_out) {
		if (!(flags & SHT_THETA_CONTIGUOUS))  flags |= SHT_PHI_CONTIGUOUS;	// default to SHT_PHI_CONTIGUOUS.
		*nlat_out = 0;		*nphi_out = 0;
		shtns_set_grid_auto($self, flags, eps, nl_order, nlat_out, nphi_out);
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
			throw_exception(SWIG_RuntimeError,0,msg_grid_err);
			return NULL;
		}
		PyObject *obj = PyArray_New(&PyArray_Type, 1, &dims, PyArray_DOUBLE, &strides, NULL, strides, 0, NULL);
		double *ct = (double*) PyArray_DATA(obj);
		for (i=0; i<$self->nlat; i++)		ct[i] = $self->ct[i];		// copy
		return obj;
	}

	// returns the index in a spectral array of (l,m) coefficient.
	int idx(int l, int m) {
		if ( (l < 0) || (l > $self->lmax) ) {
			throw_exception(SWIG_ValueError,1,"l invalid");	return 0;
		}
		if ( (m < 0) || (m > $self->mmax * $self->mres) || (m % $self->mres != 0) ) {
			throw_exception(SWIG_ValueError,2,"m invalid");	return 0;
		}
		return LM($self, l, m);
	}

	/* scalar transforms */
	void spat_to_SH(PyObject *Vr, PyObject *Qlm) {
		if (check_spatial(1,Vr, $self->nspat) && check_spectral(2,Qlm, $self->nlm))
			spat_to_SH($self, PyArray_DATA(Vr), PyArray_DATA(Qlm));
	}
	void SH_to_spat(PyObject *Qlm, PyObject *Vr) {
		if (check_spatial(2,Vr, $self->nspat) && check_spectral(1,Qlm, $self->nlm))
			SH_to_spat($self, PyArray_DATA(Qlm), PyArray_DATA(Vr));
	}

	/* 2D vectors */
	void spat_to_SHsphtor(PyObject *Vt, PyObject *Vp, PyObject *Slm, PyObject *Tlm) {
		if (check_spatial(1,Vt, $self->nspat) && check_spatial(2,Vp, $self->nspat) && check_spectral(3,Slm, $self->nlm) && check_spectral(4,Tlm, $self->nlm))
			spat_to_SHsphtor($self, PyArray_DATA(Vt), PyArray_DATA(Vp), PyArray_DATA(Slm), PyArray_DATA(Tlm));
	}
	void SHsphtor_to_spat(PyObject *Slm, PyObject *Tlm, PyObject *Vt, PyObject *Vp) {
		if (check_spatial(3,Vt, $self->nspat) && check_spatial(4,Vp, $self->nspat) && check_spectral(1,Slm, $self->nlm) && check_spectral(2,Tlm, $self->nlm))
			SHsphtor_to_spat($self, PyArray_DATA(Slm), PyArray_DATA(Tlm), PyArray_DATA(Vt), PyArray_DATA(Vp));
	}
	void SHsph_to_spat(PyObject *Slm, PyObject *Vt, PyObject *Vp) {
		if (check_spatial(2,Vt, $self->nspat) && check_spatial(3,Vp, $self->nspat) && check_spectral(1,Slm, $self->nlm))
		SHsph_to_spat($self, PyArray_DATA(Slm), PyArray_DATA(Vt), PyArray_DATA(Vp));
	}
	void SHtor_to_spat(PyObject *Tlm, PyObject *Vt, PyObject *Vp) {
		if (check_spatial(2,Vt, $self->nspat) && check_spatial(3,Vp, $self->nspat) && check_spectral(1,Tlm, $self->nlm))
		SHtor_to_spat($self, PyArray_DATA(Tlm), PyArray_DATA(Vt), PyArray_DATA(Vp));
	}

	/* 3D vectors */
	void spat_to_SHqst(PyObject *Vr, PyObject *Vt, PyObject *Vp, PyObject *Qlm, PyObject *Slm, PyObject *Tlm) {
		if (check_spatial(1,Vr, $self->nspat) && check_spatial(2,Vt, $self->nspat) && check_spatial(3,Vp, $self->nspat)
			&& check_spectral(4,Qlm, $self->nlm) && check_spectral(5,Slm, $self->nlm) && check_spectral(6,Tlm, $self->nlm))
		spat_to_SHqst($self, PyArray_DATA(Vr), PyArray_DATA(Vt), PyArray_DATA(Vp), PyArray_DATA(Qlm), PyArray_DATA(Slm), PyArray_DATA(Tlm));
	}
	void SHqst_to_spat(PyObject *Qlm, PyObject *Slm, PyObject *Tlm, PyObject *Vr, PyObject *Vt, PyObject *Vp) {
		if (check_spatial(4,Vr, $self->nspat) && check_spatial(5,Vt, $self->nspat) && check_spatial(6,Vp, $self->nspat)
			&& check_spectral(1,Qlm, $self->nlm) && check_spectral(2,Slm, $self->nlm) && check_spectral(3,Tlm, $self->nlm))
		SHqst_to_spat($self, PyArray_DATA(Qlm), PyArray_DATA(Slm), PyArray_DATA(Tlm), PyArray_DATA(Vr), PyArray_DATA(Vt), PyArray_DATA(Vp));
	}

	/* local evaluations */
	double SH_to_point(PyObject *Qlm, double cost, double phi) {
		if (check_spectral(1,Qlm, $self->nlm))	return SH_to_point($self, PyArray_DATA(Qlm), cost, phi);
	}
	%apply double *OUTPUT { double *vr };
	%apply double *OUTPUT { double *vt };
	%apply double *OUTPUT { double *vp };
	void SHqst_to_point(PyObject *Qlm, PyObject *Slm, PyObject *Tlm,
					double cost, double phi, double *vr, double *vt, double *vp) {
		if (check_spectral(1,Qlm, $self->nlm) && check_spectral(2,Slm, $self->nlm) && check_spectral(3,Tlm, $self->nlm))
			SHqst_to_point($self, PyArray_DATA(Qlm), PyArray_DATA(Slm), PyArray_DATA(Tlm), cost, phi, vr, vt, vp);
	}
	%clear double *vr;
	%clear double *vt;
	%clear double *vp;

	/* rotation of SH representations (experimental) */
	void SH_Zrotate(PyObject *Qlm, double alpha, PyObject *Rlm) {
		if (check_spectral(1,Qlm, $self->nlm) && check_spectral(3,Rlm, $self->nlm))
			SH_Zrotate($self, PyArray_DATA(Qlm), alpha, PyArray_DATA(Rlm));
	}
	void SH_Yrotate(PyObject *Qlm, double alpha, PyObject *Rlm) {
		if (check_spectral(1,Qlm, $self->nlm) && check_spectral(3,Rlm, $self->nlm))
			SH_Yrotate($self, PyArray_DATA(Qlm), alpha, PyArray_DATA(Rlm));
	}
	void SH_Yrotate90(PyObject *Qlm, PyObject *Rlm) {
		if (check_spectral(1,Qlm, $self->nlm) && check_spectral(2,Rlm, $self->nlm))
			SH_Yrotate90($self, PyArray_DATA(Qlm), PyArray_DATA(Rlm));
	}
	void SH_Xrotate90(PyObject *Qlm, PyObject *Rlm) {
		if (check_spectral(1,Qlm, $self->nlm) && check_spectral(2,Rlm, $self->nlm))
			SH_Xrotate90($self, PyArray_DATA(Qlm), PyArray_DATA(Rlm));
	}

};
