# Python setup

from distutils.core import setup, Extension
from numpy import get_include

numpy_inc = get_include()		#  NumPy include path.
shtns_o = ['SHT.o', 'sht_std.o', 'sht_ltr.o', 'sht_m0.o', 'sht_m0ltr.o']

shtns_module = Extension('_shtns', sources=['shtns_numpy_wrap.c'],
	extra_objects=shtns_o, depends=shtns_o,
	include_dirs=[numpy_inc],
	libraries=['fftw3', 'fftw3_omp'])

setup(name='SHTns',
	version='2.2',
	description='High performance Spherical Harmonic Transform',
	license='CeCILL',
	author='Nathanael Schaeffer',
	author_email='nschaeff@ujf-grenoble.fr',
	url='https://bitbucket.org/nschaeff/shtns',
	ext_modules=[shtns_module],
	py_modules=["shtns"],
	)
