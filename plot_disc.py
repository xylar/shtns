#!/usr/bin/python
#plot a slice
from pylab import *     # matplotlib
import sys              # acces a la ligne de commande.

ir0=0
rg=0.35

argn = len(sys.argv)
if argn < 2:
	print "usage: " + sys.argv[0] + " <xspp_disc_file> [0|1|2] [1|2|0] ..."
	print "  optional parameters 0,1 or 2 chooses component (x,y or z)"
	exit()

print 'loading',sys.argv[1]
a=load(sys.argv[1],comments='%')
s = a.shape

np = (s[1]-1)/3+1
ip = mgrid[0:np]
ip[len(ip)-1] = 0	# ip loops around
phi = ip*2*pi/(np-1)
r = transpose(matrix(a[ir0:,0]))
x = r*matrix(cos(phi))
y = r*matrix(sin(phi))

title = array(['r','phi','z'])
comp = range(0,3)

if argn > 2:
	for i in range(2,argn):
		comp[i-2] = int(sys.argv[i])
	comp = comp[0:argn-2]

print title[comp]
for i in comp:
	figure()
	b = a[:, 1+i +ip*3]
	m=amax(abs(b))
	print 'max value=',m
	
	contourf(array(x),array(y),b,20,cmap=cm.PuOr)
	axis('equal')
	axis('off')
	clim(-m,m)
	subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.98, wspace=0.1, hspace=0.1)
	colorbar()
	figtext(0.05, 0.9, title[i])
	theta = linspace(-pi/2,pi/2,100)
	plot(rg*cos(theta),rg*sin(theta),color='gray')
	plot(cos(theta),sin(theta),color='gray')
	plot(-rg*cos(theta),rg*sin(theta),color='gray')
	plot(-cos(theta),sin(theta),color='gray')

show()
