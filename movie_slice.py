#!/usr/bin/python
#plot a slice
from pylab import *     # matplotlib
from subprocess import *
import sys              # acces a la ligne de commande.
import glob
import time

ir0=1
rg=0.35

i=1;	# snapshot numbering

base=sys.argv[1]
job=sys.argv[2]
files=sorted(glob.glob( base + '_*.' + job))

if len(sys.argv) >3:
	i = int(sys.argv[3])

for f in files[i-1:] :
	retcode = call("./xspp " + f + " axi", shell=True)
	#retcode = os.system("./xspp " + f + " axi")
	if retcode != 0:
		print 'error from xspp again, abort.'
		exit()

	a=load('o_Vp',comments='%')
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
	contourf(array(x),array(y),a,15,cmap=cm.RdBu)
	theta = linspace(-pi/2,pi/2,100)
	plot(rg*cos(theta),rg*sin(theta),color='gray')
	plot(cos(theta),sin(theta),color='gray')

	axis('equal')
	axis('off')
	colorbar()
	clim(-m,m)
	subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.98, wspace=0.1, hspace=0.1)
	figtext(0.05, 0.9, str(i))

#	savefig( "slide_" + str(i) + ".png" )
	savefig( "slide_%04d.png" % i )
	print '>>> slide', i, "saved.\n"
	i+=1
	close()

#retcode = os.system( 'mencoder "mf://slide_*.png" -mf fps=10 -o movie.avi -ovc lavc -lavcopts vcodec=mjpeg' )
