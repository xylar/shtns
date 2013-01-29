/*
 * Copyright (c) 2010-2013 Centre National de la Recherche Scientifique.
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

%init{
	import_array();		// required by NumPy
}

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
		throw_exception(SWIG_TypeError,i,"spatial array must consist of float.");		return 0;
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
		throw_exception(SWIG_RuntimeError,i,"spectral array must consist of complex float. Create with: 'sh.spec_array()'");		return 0;
	}
	if (!PyArray_ISCONTIGUOUS(a)) {
		throw_exception(SWIG_RuntimeError,i,"spactral array not contiguous. Use .copy() to copy to a contiguous array.");		return 0;
	}
	if (PyArray_SIZE(a) != size) {
		throw_exception(SWIG_RuntimeError,i,"spectral array has wrong size");		return 0;
	}
	return 1;
}

inline PyObject* SpecArray_New(int size) {
	npy_intp dims = size;
	npy_intp strides = sizeof(complex double);
	return PyArray_New(&PyArray_Type, 1, &dims, PyArray_CDOUBLE, &strides, NULL, strides, 0, NULL);	
}

inline PyObject* SpatArray_New(int size) {
	npy_intp dims = size;
	npy_intp strides = sizeof(double);
	return PyArray_New(&PyArray_Type, 1, &dims, PyArray_DOUBLE, &strides, NULL, strides, 0, NULL);	
}

%}

// main object is renamed to sht.
%rename("sht") shtns_info;
%ignore SHT_NATIVE_LAYOUT;
%ignore nlat_2;
%ignore lmidx;
%ignore li;
%ignore mi;
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
		## array giving the degree of spherical harmonic coefficients.
		self.l = np.zeros(self.nlm, dtype=np.int32)
		## array giving the order of spherical harmonic coefficients.
		self.m = np.zeros(self.nlm, dtype=np.int32)
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
		shtns_use_threads(nthreads);		// use nthreads openmp threads if available (0 means auto)
		return shtns_create(lmax, mmax, mres, norm);
	}

	~shtns_info() {
		shtns_destroy($self);		// free memory.
	}

	%pythonappend set_grid %{
		## array giving the cosine of the colatitude for the grid.
		self.cos_theta = self.__ct()
		self.cos_theta.flags.writeable = False
		## shape of a spatial array for the grid (tuple of 2 values).
		self.spat_shape = tuple(self.__spat_shape())
	%}
	%apply int *OUTPUT { int *nlat_out };
	%apply int *OUTPUT { int *nphi_out };
	%feature("kwargs") set_grid;
	void set_grid(int nlat=0, int nphi=0, int flags=sht_quick_init, double polar_opt=1.0e-8, int nl_order=1, int *nlat_out, int *nphi_out) {	// default arguments
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
		shtns_set_grid_auto($self, flags, polar_opt, nl_order, nlat_out, nphi_out);
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
		if ($self->nlat == 0) {	// no grid
			throw_exception(SWIG_RuntimeError,0,msg_grid_err);
			return NULL;
		}
		PyObject *obj = SpatArray_New($self->nlat);
		double *ct = (double*) PyArray_DATA(obj);
		for (i=0; i<$self->nlat; i++)		ct[i] = $self->ct[i];		// copy
		return obj;
	}
	PyObject* gauss_wts() {		// gauss grid must have been initialized.
		if ($self->nlat == 0) {	// no grid
			throw_exception(SWIG_RuntimeError,0,msg_grid_err);
			return NULL;
		}
		if ($self->wg == NULL) {
			throw_exception(SWIG_RuntimeError,0,"not a gauss grid");
			return NULL;
		}
		PyObject *obj = SpatArray_New($self->nlat_2);
		shtns_gauss_wts($self, PyArray_DATA(obj));
		return obj;
	}
	PyObject* mul_ct_matrix() {
		PyObject *mx = SpatArray_New(2*$self->nlm);
		mul_ct_matrix($self, PyArray_DATA(mx));
		return mx;
	}
	PyObject* st_dt_matrix() {
		PyObject *mx = SpatArray_New(2*$self->nlm);
		st_dt_matrix($self, PyArray_DATA(mx));
		return mx;
	}

	%apply int *OUTPUT { int *dim0 };
	%apply int *OUTPUT { int *dim1 };
	void __spat_shape(int *dim0, int *dim1) {
		*dim0 = $self->nphi;	*dim1 = $self->nlat;
		if ($self->fftc_mode == 1) {	// phi-contiguous
			*dim0 = $self->nlat;		*dim1 = $self->nphi;
		}
	}

	%pythoncode %{
		def spec_array(self):
			"""return a numpy array of spherical harmonic coefficients (complex). Adress coefficients with index sh.idx(l,m)"""
			return np.zeros(self.nlm, dtype=complex)
		
		def spat_array(self):
			"""return a numpy array of 2D spatial field."""
			if self.nlat == 0: raise RuntimeError("Grid not set. Call .set_grid() mehtod.")
			return np.zeros(self.spat_shape)
	%}

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

	%pythoncode %{
		def synth(self,*arg):
			"""
			spectral to spatial transform, for scalar or vector data.
			v = synth(qlm) : compute the spatial representation of the scalar qlm
			vtheta,vphi = synth(slm,tlm) : compute the 2D spatial vector from its spectral spheroidal/toroidal scalars (slm,tlm)
			vr,vtheta,vphi = synth(qlm,slm,tlm) : compute the 3D spatial vector from its spectral radial/spheroidal/toroidal scalars (qlm,slm,tlm)
			"""
			if self.nlat == 0: raise RuntimeError("Grid not set. Call .set_grid() mehtod.")
			n = len(arg)
			if (n>3) or (n<1): raise RuntimeError("1,2 or 3 arguments required.")
			q = list(arg)
			for i in range(0,n):
				if q[i].size != self.nlm: raise RuntimeError("spectral array has wrong size.")
				if q[i].dtype.num != np.dtype('complex128').num: raise RuntimeError("spectral array should be dtype=complex.")
				if q[i].flags.contiguous == False: q[i] = q[i].copy()		# contiguous array required.
			if n==1:	#scalar transform
				vr = np.empty(self.spat_shape)
				self.SH_to_spat(q[0],vr)
				return vr
			elif n==2:	# 2D vector transform
				vt = np.empty(self.spat_shape)		# v_theta
				vp = np.empty(self.spat_shape)		# v_phi
				self.SHsphtor_to_spat(q[0],q[1],vt,vp)
				return vt,vp
			else:		# 3D vector transform
				vr = np.empty(self.spat_shape)		# v_r
				vt = np.empty(self.spat_shape)		# v_theta
				vp = np.empty(self.spat_shape)		# v_phi
				self.SHqst_to_spat(q[0],q[1],q[2],vr,vt,vp)
				return vr,vt,vp

		def analys(self,*arg):
			"""
			spatial to spectral transform, for scalar or vector data.
			qlm = analys(q) : compute the spherical harmonic representation of the scalar q
			slm,tlm = analys(vtheta,vphi) : compute the spectral spheroidal/toroidal scalars (slm,tlm) from 2D vector components (vtheta, vphi)
			qlm,slm,tlm = synth(vr,vtheta,vphi) : compute the spectral radial/spheroidal/toroidal scalars (qlm,slm,tlm) from 3D vector components (vr,vtheta,vphi)
			"""
			if self.nlat == 0: raise RuntimeError("Grid not set. Call .set_grid() mehtod.")
			n = len(arg)
			if (n>3) or (n<1): raise RuntimeError("1,2 or 3 arguments required.")
			v = list(arg)
			for i in range(0,n):
				if v[i].shape != self.spat_shape: raise RuntimeError("spatial array has wrong shape.")
				if v[i].dtype.num != np.dtype('float64').num: raise RuntimeError("spatial array should be dtype=float64.")
				if v[i].flags.contiguous == False: v[i] = v[i].copy()		# contiguous array required.
			if n==1:
				q = np.empty(self.nlm, dtype=complex)
				self.spat_to_SH(v[0],q)
				return q
			elif n==2:
				s = np.empty(self.nlm, dtype=complex)
				t = np.empty(self.nlm, dtype=complex)
				self.spat_to_SHsphtor(v[0],v[1],s,t)
				return s,t
			else:
				q = np.empty(self.nlm, dtype=complex)
				s = np.empty(self.nlm, dtype=complex)
				t = np.empty(self.nlm, dtype=complex)
				self.spat_to_SHqst(v[0],v[1],v[2],q,s,t)
				return q,s,t

		def synth_grad(self,slm):
			"""(vtheta,vphi) = synth_grad(sht self, slm) : compute the spatial representation of the gradient of slm"""
			if self.nlat == 0: raise RuntimeError("Grid not set. Call .set_grid() mehtod.")
			if slm.size != self.nlm: raise RuntimeError("spectral array has wrong size.")
			if slm.dtype.num != np.dtype('complex128').num: raise RuntimeError("spectral array should be dtype=complex.")
			if slm.flags.contiguous == False: slm = slm.copy()		# contiguous array required.
			vt = np.empty(self.spat_shape)
			vp = np.empty(self.spat_shape)
			self.SHsph_to_spat(slm,vt,vp)
			return vt,vp
	%}

	/* local evaluations */
	double SH_to_point(PyObject *Qlm, double cost, double phi) {
		if (check_spectral(1,Qlm, $self->nlm))	return SH_to_point($self, PyArray_DATA(Qlm), cost, phi);
		return 0.0;
	}
	%apply double *OUTPUT { double *vr };
	%apply double *OUTPUT { double *vt };
	%apply double *OUTPUT { double *vp };
	void SH_to_grad_point(PyObject *DrSlm, PyObject *Slm, double cost, double phi, double *vr, double *vt, double *vp) {
		if (check_spectral(1,DrSlm, $self->nlm) && check_spectral(2,Slm, $self->nlm))
			SH_to_grad_point($self, PyArray_DATA(DrSlm), PyArray_DATA(Slm), cost, phi, vr, vt, vp);
	}
	void SHqst_to_point(PyObject *Qlm, PyObject *Slm, PyObject *Tlm,
					double cost, double phi, double *vr, double *vt, double *vp) {
		if (check_spectral(1,Qlm, $self->nlm) && check_spectral(2,Slm, $self->nlm) && check_spectral(3,Tlm, $self->nlm))
			SHqst_to_point($self, PyArray_DATA(Qlm), PyArray_DATA(Slm), PyArray_DATA(Tlm), cost, phi, vr, vt, vp);
	}
	%clear double *vr;
	%clear double *vt;
	%clear double *vp;

	/* rotation of SH representations (experimental) */
	PyObject* Zrotate(PyObject *Qlm, double alpha) {
		if (check_spectral(1,Qlm, $self->nlm)) {
			PyObject *Rlm = SpecArray_New($self->nlm);
			SH_Zrotate($self, PyArray_DATA(Qlm), alpha, PyArray_DATA(Rlm));
			return Rlm;
		}
		return NULL;
	}
	PyObject* Yrotate(PyObject *Qlm, double alpha) {
		if (($self->mres != 1)||($self->mmax != $self->lmax)) {
			throw_exception(SWIG_RuntimeError,0,msg_rot_err);	return NULL;
		}
		if (check_spectral(1,Qlm, $self->nlm)) {
			PyObject *Rlm = SpecArray_New($self->nlm);
			SH_Yrotate($self, PyArray_DATA(Qlm), alpha, PyArray_DATA(Rlm));
			return Rlm;
		}
		return NULL;
	}
	PyObject* Yrotate90(PyObject *Qlm) {
		if (($self->mres != 1)||($self->mmax != $self->lmax)) {
			throw_exception(SWIG_RuntimeError,0,msg_rot_err);	return NULL;
		}
		if (check_spectral(1,Qlm, $self->nlm)) {
			PyObject *Rlm = SpecArray_New($self->nlm);
			SH_Yrotate90($self, PyArray_DATA(Qlm), PyArray_DATA(Rlm));
			return Rlm;
		}
		return NULL;
	}
	PyObject* Xrotate90(PyObject *Qlm) {
		if (($self->mres != 1)||($self->mmax != $self->lmax)) {
			throw_exception(SWIG_RuntimeError,0,msg_rot_err);	return NULL;
		}
		if (check_spectral(1,Qlm, $self->nlm)) {
			PyObject *Rlm = SpecArray_New($self->nlm);
			SH_Xrotate90($self, PyArray_DATA(Qlm), PyArray_DATA(Rlm));
			return Rlm;
		}
		return NULL;
	}

	/* multiplication by l+1 l-1 matrix (mul_ct_matrix or st_dt_matrix) */
	PyObject* SH_mul_mx(PyObject *mx, PyObject *Qlm) {
		if (check_spectral(2,Qlm, $self->nlm) && check_spatial(1, mx, 2* $self->nlm)) {
			PyObject *Rlm = SpecArray_New($self->nlm);
			SH_mul_mx($self, PyArray_DATA(mx), PyArray_DATA(Qlm), PyArray_DATA(Rlm));
			return Rlm;
		}
		return NULL;
	}

};
