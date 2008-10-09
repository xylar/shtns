#!/usr/bin/python
#plot a slice
from pylab import *     # matplotlib
import sys              # acces a la ligne de commande.

ir0=1

argn = len(sys.argv)
print argn

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

	print a.shape
	
	#convert Up to angular velocity
	#a=a/array(x)
	
	figure()
	#pcolor(array(x),array(y),a,shading='interp')
	contourf(array(x),array(y),a,15,cmap=cm.PuOr)
	axis('equal')
	colorbar()
	show()
