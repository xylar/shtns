#!/usr/bin/python
#plot a slice
from pylab import *     # matplotlib
import sys              # acces a la ligne de commande.

argn = len(sys.argv)
if argn < 2:
	print "usage: " + sys.argv[0] + " <xspp_surf_file> [0|1|2] [1|2|0] ..."
	print "  optional parameters 0,1 or 2 chooses component (x,y or z)"
	exit()

print 'loading',sys.argv[1]
a=load(sys.argv[1],comments='%')
s = a.shape

t = array(a[0,1:s[1]:3])
t = t * 180./pi
p = array(a[1:s[0],0])
p = p * 180./pi

title = array(['r','phi','z'])
comp = range(0,3)

if argn > 2:
	for i in range(2,argn):
		comp[i-2] = int(sys.argv[i])
	comp = comp[0:argn-2]

print title[comp]
for i in comp:
	b = array(transpose(matrix(a[1:s[0],(1+i)::3])))
	m=amax(abs(b))
	print 'max value=',m
	figure()
	contourf(array(p),array(t),b,20,cmap=cm.PuOr)

#	axis('equal')
#	axis('off')
#	subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.98, wspace=0.1, hspace=0.1)
	colorbar()
	clim(-m,m)
	xlim(0,360)
	ylim(180,0)
	figtext(0.05, 0.9, title[i])

show()