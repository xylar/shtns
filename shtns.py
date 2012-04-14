# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.4
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.


"""
Python/NumPy interface to the SHTns spherical harmonic transform library
"""


from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_shtns', [dirname(__file__)])
        except ImportError:
            import _shtns
            return _shtns
        if fp is not None:
            try:
                _mod = imp.load_module('_shtns', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _shtns = swig_import_helper()
    del swig_import_helper
else:
    import _shtns
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


sht_orthonormal = _shtns.sht_orthonormal
sht_fourpi = _shtns.sht_fourpi
sht_schmidt = _shtns.sht_schmidt
SHT_NO_CS_PHASE = _shtns.SHT_NO_CS_PHASE
SHT_REAL_NORM = _shtns.SHT_REAL_NORM
sht_gauss = _shtns.sht_gauss
sht_auto = _shtns.sht_auto
sht_reg_fast = _shtns.sht_reg_fast
sht_reg_dct = _shtns.sht_reg_dct
sht_quick_init = _shtns.sht_quick_init
sht_reg_poles = _shtns.sht_reg_poles
sht_gauss_fly = _shtns.sht_gauss_fly
SHT_NATIVE_LAYOUT = _shtns.SHT_NATIVE_LAYOUT
SHT_THETA_CONTIGUOUS = _shtns.SHT_THETA_CONTIGUOUS
SHT_PHI_CONTIGUOUS = _shtns.SHT_PHI_CONTIGUOUS
SHT_SOUTH_POLE_FIRST = _shtns.SHT_SOUTH_POLE_FIRST
SHT_SCALAR_ONLY = _shtns.SHT_SCALAR_ONLY
class sht(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, sht, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, sht, name)
    __repr__ = _swig_repr
    __swig_getmethods__["nlm"] = _shtns.sht_nlm_get
    if _newclass:nlm = _swig_property(_shtns.sht_nlm_get)
    __swig_getmethods__["lmax"] = _shtns.sht_lmax_get
    if _newclass:lmax = _swig_property(_shtns.sht_lmax_get)
    __swig_getmethods__["mmax"] = _shtns.sht_mmax_get
    if _newclass:mmax = _swig_property(_shtns.sht_mmax_get)
    __swig_getmethods__["mres"] = _shtns.sht_mres_get
    if _newclass:mres = _swig_property(_shtns.sht_mres_get)
    __swig_getmethods__["nphi"] = _shtns.sht_nphi_get
    if _newclass:nphi = _swig_property(_shtns.sht_nphi_get)
    __swig_getmethods__["nlat"] = _shtns.sht_nlat_get
    if _newclass:nlat = _swig_property(_shtns.sht_nlat_get)
    __swig_getmethods__["nlat_2"] = _shtns.sht_nlat_2_get
    if _newclass:nlat_2 = _swig_property(_shtns.sht_nlat_2_get)
    __swig_getmethods__["lmidx"] = _shtns.sht_lmidx_get
    if _newclass:lmidx = _swig_property(_shtns.sht_lmidx_get)
    __swig_getmethods__["li"] = _shtns.sht_li_get
    if _newclass:li = _swig_property(_shtns.sht_li_get)
    __swig_getmethods__["ct"] = _shtns.sht_ct_get
    if _newclass:ct = _swig_property(_shtns.sht_ct_get)
    __swig_getmethods__["st"] = _shtns.sht_st_get
    if _newclass:st = _swig_property(_shtns.sht_st_get)
    __swig_getmethods__["nspat"] = _shtns.sht_nspat_get
    if _newclass:nspat = _swig_property(_shtns.sht_nspat_get)
    def __init__(self, *args): 
        """__init__(self, int lmax, int mmax = -1, int mres = 1, int norm = sht_orthonormal) -> sht"""
        this = _shtns.new_sht(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _shtns.delete_sht
    __del__ = lambda self : None;
    def set_grid(self, *args):
        """
        set_grid(self, int nlat, int nphi, enum shtns_type flags = sht_quick_init|256, 
            double eps = 1.0e-8) -> int
        """
        return _shtns.sht_set_grid(self, *args)

    def set_grid_auto(self, *args):
        """
        set_grid_auto(self, int nlat = 0, int nphi = 0, int nl_order = 1, enum shtns_type flags = sht_quick_init|256, 
            double eps = 1.0e-8) -> int
        """
        return _shtns.sht_set_grid_auto(self, *args)

    def _print(self):
        """_print(self)"""
        return _shtns.sht__print(self)

    def sh00_1(self):
        """sh00_1(self) -> double"""
        return _shtns.sht_sh00_1(self)

    def sh10_ct(self):
        """sh10_ct(self) -> double"""
        return _shtns.sht_sh10_ct(self)

    def sh11_st(self):
        """sh11_st(self) -> double"""
        return _shtns.sht_sh11_st(self)

    def shlm_e1(self, *args):
        """shlm_e1(self, int l, int m) -> double"""
        return _shtns.sht_shlm_e1(self, *args)

    def l(self):
        """l(self) -> PyObject"""
        return _shtns.sht_l(self)

    def cos_theta(self):
        """cos_theta(self) -> PyObject"""
        return _shtns.sht_cos_theta(self)

    def idx(self, *args):
        """idx(self, int l, int m) -> int"""
        return _shtns.sht_idx(self, *args)

    def spat_to_SH(self, *args):
        """spat_to_SH(self, PyObject Vr, PyObject Qlm)"""
        return _shtns.sht_spat_to_SH(self, *args)

    def SH_to_spat(self, *args):
        """SH_to_spat(self, PyObject Qlm, PyObject Vr)"""
        return _shtns.sht_SH_to_spat(self, *args)

    def spat_to_SHsphtor(self, *args):
        """spat_to_SHsphtor(self, PyObject Vt, PyObject Vp, PyObject Slm, PyObject Tlm)"""
        return _shtns.sht_spat_to_SHsphtor(self, *args)

    def SHsphtor_to_spat(self, *args):
        """SHsphtor_to_spat(self, PyObject Slm, PyObject Tlm, PyObject Vt, PyObject Vp)"""
        return _shtns.sht_SHsphtor_to_spat(self, *args)

    def SHsph_to_spat(self, *args):
        """SHsph_to_spat(self, PyObject Slm, PyObject Vt, PyObject Vp)"""
        return _shtns.sht_SHsph_to_spat(self, *args)

    def SHtor_to_spat(self, *args):
        """SHtor_to_spat(self, PyObject Tlm, PyObject Vt, PyObject Vp)"""
        return _shtns.sht_SHtor_to_spat(self, *args)

    def spat_to_SHqst(self, *args):
        """
        spat_to_SHqst(self, PyObject Vr, PyObject Vt, PyObject Vp, PyObject Qlm, 
            PyObject Slm, PyObject Tlm)
        """
        return _shtns.sht_spat_to_SHqst(self, *args)

    def SHqst_to_spat(self, *args):
        """
        SHqst_to_spat(self, PyObject Qlm, PyObject Slm, PyObject Tlm, PyObject Vr, 
            PyObject Vt, PyObject Vp)
        """
        return _shtns.sht_SHqst_to_spat(self, *args)

    def SH_to_point(self, *args):
        """SH_to_point(self, PyObject Qlm, double cost, double phi) -> double"""
        return _shtns.sht_SH_to_point(self, *args)

    def SHqst_to_point(self, *args):
        """
        SHqst_to_point(self, PyObject Qlm, PyObject Slm, PyObject Tlm, double cost, 
            double phi)
        """
        return _shtns.sht_SHqst_to_point(self, *args)

    def SH_Zrotate(self, *args):
        """SH_Zrotate(self, PyObject Qlm, double alpha, PyObject Rlm)"""
        return _shtns.sht_SH_Zrotate(self, *args)

    def SH_Yrotate(self, *args):
        """SH_Yrotate(self, PyObject Qlm, double alpha, PyObject Rlm)"""
        return _shtns.sht_SH_Yrotate(self, *args)

    def SH_Yrotate90(self, *args):
        """SH_Yrotate90(self, PyObject Qlm, PyObject Rlm)"""
        return _shtns.sht_SH_Yrotate90(self, *args)

    def SH_Xrotate90(self, *args):
        """SH_Xrotate90(self, PyObject Qlm, PyObject Rlm)"""
        return _shtns.sht_SH_Xrotate90(self, *args)

sht_swigregister = _shtns.sht_swigregister
sht_swigregister(sht)


def nlm_calc(*args):
  """nlm_calc(long lmax, long mmax, long mres) -> long"""
  return _shtns.nlm_calc(*args)
# This file is compatible with both classic and new-style classes.

