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

%pythoncode{
	import numpy as np
}

%{
	
#include <numpy/arrayobject.h>
#include "sht_private.h"

// variables used for exception handling.
static int shtns_error = 0;
static char* shtns_err_msg;
static char msg_buffer[128];
static char msg_grid_err[] = "Grid not set. Call .set_grid() mehtod.";
static char msg_numpy_arr[] = "Numpy array expected.";
static char msg_rot_err[] = "truncation must be triangular (lmax=mmax, mres=1)";

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
%ignore nlat_2;
%ignore lmidx;
%ignore li;
%ignore ct;
%ignore st;

%feature("autodoc");
%include "shtns.h"
%include "exception.i"


%extend shtns_info {
	%exception {
		shtns_error = 0;	// clear exception
		$function
		if (shtns_error) {	// test for exception
			SWIG_exception(shtns_error, shtns_err_msg);		return NULL;
		}
	}

	%pythonappend shtns_info %{
		self.m = np.zeros(self.nlm, dtype=np.int32)
		self.l = np.zeros(self.nlm, dtype=np.int32)
		for mloop in range(0, self.mmax*self.mres+1, self.mres):
			for lloop in range(mloop, self.lmax+1):
				ii = self.idx(lloop,mloop)
				self.m[ii] = mloop
				self.l[ii] = lloop
		self.m.flags.writeable = False		# prevent writing in m and l arrays
		self.l.flags.writeable = False
	%}
	%feature("kwargs") shtns_info;
	shtns_info(int lmax, int mmax=-1, int mres=1, int norm=sht_orthonormal, int nthreads=0) {	// default arguments : mmax, mres and norm
		if (lmax < 2) {
			throw_exception(SWIG_ValueError,1,"lmax < 2 not allowed");	return NULL;
		}
		if (mres <= 0) {
			throw_exception(SWIG_ValueError,3,"mres <= 0 invalid");	return NULL;
		}
		if (mmax < 0) mmax = lmax/mres;		// default mmax
		if (mmax*mres > lmax) {
			throw_exception(SWIG_ValueError,1,"lmax < mmax*mres invalid");	return NULL;
		}
		import_array();		// required by NumPy
		shtns_use_threads(nthreads);		// use nthreads openmp threads if available (0 means auto)
		return shtns_create(lmax, mmax, mres, norm);
	}

	~shtns_info() {
		shtns_destroy($self);		// free memory.
	}
	
	%pythonappend set_grid %{
		self.cos_theta = self.__ct()
		self.cos_theta.flags.writeable = False
	%}
	%apply int *OUTPUT { int *nlat_out };
	%apply int *OUTPUT { int *nphi_out };
	%feature("kwargs") set_grid;
	void set_grid(int nlat=0, int nphi=0, int flags=sht_quick_init, double eps=1.0e-8, int nl_order=1, int *nlat_out, int *nphi_out) {	// default arguments
		if (nlat != 0) {
			if (nlat <= $self->lmax) {	// nlat too small
				throw_exception(SWIG_ValueError,1,"nlat <= lmax");		return;
			}
			if (nlat & 1) {		// nlat must be even
				throw_exception(SWIG_ValueError,1,"nlat must be even");		return;
			}
		}
		if ((nphi != 0) && (nphi <= $self->mmax *2)) {		// nphi too small
			throw_exception(SWIG_ValueError,2,"nphi <= 2*mmax");	return;
		}
		if (!(flags & SHT_THETA_CONTIGUOUS))  flags |= SHT_PHI_CONTIGUOUS;	// default to SHT_PHI_CONTIGUOUS.
		*nlat_out = nlat;		*nphi_out = nphi;
		shtns_set_grid_auto($self, flags, eps, nl_order, nlat_out, nphi_out);
	}

	void print_info() {
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
	double shlm_e1(unsigned l, unsigned m) {
		return shlm_e1($self, l, m);
	}

	/* returns useful data */
	PyObject* __ct() {		// grid must have been initialized.
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

	%pythoncode %{
		def spec_array(self):
			return np.zeros(self.nlm, dtype=complex)
	%}

	PyObject* spat_array() {
		npy_intp dims[2];
		if ($self->nlat == 0) {	// no grid
			throw_exception(SWIG_RuntimeError,0,msg_grid_err);
			return NULL;
		}
		dims[0] = $self->nphi;	dims[1] = $self->nlat;
		if ($self->fftc_mode == 1) {	// phi-contiguous
			dims[0] = $self->nlat;		dims[1] = $self->nphi;
		}
		return PyArray_ZEROS(2, dims, PyArray_DOUBLE, 0);
	}

	// returns the index in a spectral array of (l,m) coefficient.
	int idx(unsigned l, unsigned m) {
		if (l > $self->lmax) {
			throw_exception(SWIG_ValueError,1,"l invalid");	return 0;
		}
		if ( (m > l) || (m > $self->mmax * $self->mres) || (m % $self->mres != 0) ) {
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
		if (($self->mres != 1)||($self->mmax != $self->lmax)) {
			throw_exception(SWIG_RuntimeError,0,msg_rot_err);	return;
		}
		if (check_spectral(1,Qlm, $self->nlm) && check_spectral(3,Rlm, $self->nlm))
			SH_Yrotate($self, PyArray_DATA(Qlm), alpha, PyArray_DATA(Rlm));
	}
	void SH_Yrotate90(PyObject *Qlm, PyObject *Rlm) {
		if (($self->mres != 1)||($self->mmax != $self->lmax)) {
			throw_exception(SWIG_RuntimeError,0,msg_rot_err);	return;
		}
		if (check_spectral(1,Qlm, $self->nlm) && check_spectral(2,Rlm, $self->nlm))
			SH_Yrotate90($self, PyArray_DATA(Qlm), PyArray_DATA(Rlm));
	}
	void SH_Xrotate90(PyObject *Qlm, PyObject *Rlm) {
		if (($self->mres != 1)||($self->mmax != $self->lmax)) {
			throw_exception(SWIG_RuntimeError,0,msg_rot_err);	return;
		}
		if (check_spectral(1,Qlm, $self->nlm) && check_spectral(2,Rlm, $self->nlm))
			SH_Xrotate90($self, PyArray_DATA(Qlm), PyArray_DATA(Rlm));
	}

};
