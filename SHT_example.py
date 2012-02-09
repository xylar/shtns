#!/usr/bin/python2
#
# SHTns Python interface example
#

import numpy		# numpy for arrays
import shtns		# shtns interface (must be build with "make python" and copied in the current directory)

lmax = 15			# maximum degree of spherical harmonic representation.
mmax = 3			# maximum order of spherical harmonic representation.

sh = shtns.sht(lmax, mmax)		# create sht object with given lmax and mmax.
# mres = 2						# use 2-fold symmetry in phi
# norm = sht_schmidt			# use schmidt semi-normalized harmonics
# sh = shtns.sht(lmax, mmax, mres, norm)	# advanced creation of sht object.

sh.set_grid_auto()				# build grid with default size
print(sh.nlat, sh.nphi)			# displays the latitudinal and longitudinal grid sizes.

cost = sh.cos_theta()			# latitudinal coordinates of the grid as cos(theta)
el = sh.l()						# array of size sh.nlm giving the spherical harmonic degree l for any sh coefficient
l2 = el*(el+1)					# array l(l+1) that is useful for computing laplacian

# nlat = lmax*2
# nphi = mmax*3
# sh.set_grid(nlat, nphi, sht_gauss_fly|SHT_PHI_CONTIGUOUS, 1.0e-10)		# use advanced options to create a gauss grid with optimized on-the-fly transforms.

ylm = numpy.zeros(sh.nlm, dtype=complex)		# a spherical harmonic spectral array

vr = numpy.zeros((sh.nphi, sh.nlat))			# a spatial array matching the grid (with SHT_THETA_CONTIGUOUS by default)
# vr = numpy.zeros((sh.nlat, sh.nphi))			# a spatial array matching a grid build with SHT_PHI_CONTIGUOUS


ylm[sh.idx(1,0)] = 1.0		# set sh coefficient l=1, m=0 to value 1

ylm = ylm * l2				# multiply by l(l+1)

sh.SH_to_spat(ylm, vr)		# transform sh description ylm into spatial representation vr (vr is overwritten)

print(vr)		# display spatial field


zlm = numpy.zeros(sh.nlm, dtype=complex)		# a spherical harmonic spectral array
sh.spat_to_SH(vr, zlm)							# transform the spatial field back to spectral

print(zlm)
