#!/usr/bin/python
#plot a slice
from pylab import *     # matplotlib
import sys              # acces a la ligne de commande.

ir0=1
rg=0.35

argn = len(sys.argv)

for i in range(1, argn):
	print 'loading',sys.argv[i]
	a=load(sys.argv[i])
	s = a.shape

	ct = a[0,1:s[1]]
	st = sqrt(1-ct*ct)
	r = transpose(matrix(a[ir0:s[0],0]))

	a = a[ir0:s[0],1:s[1]]
	x = r*matrix(st)
	y = r*matrix(ct)

	#convert Up to angular velocity
	#a=a/array(x)
	#a[:,s[1]-2] = a[:,s[1]-3]	#remove nan
	#a[:,0] = a[:,1]			#remove nan

	m=amax(abs(a))
	print 'max value=',m
	figure()
	#pcolor(array(x),array(y),a,shading='interp')
	#colormaps : cm.PuOr, cm.RdBu, cm.RdGy
	contourf(array(x),array(y),a,20,cmap=cm.RdBu)
	theta = linspace(-pi/2,pi/2,100)
	plot(rg*cos(theta),rg*sin(theta),color='gray')
	plot(cos(theta),sin(theta),color='gray')

	axis('equal')
	axis('off')
#	subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.98, wspace=0.1, hspace=0.1)
	colorbar()
	clim(-m,m)

show()
