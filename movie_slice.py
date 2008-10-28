#!/usr/bin/python
#plot a slice
from pylab import *     # matplotlib
import os
import sys              # acces a la ligne de commande.

ir0=1
rg=0.35

i=1;	# snapshot numbering

for f in sys.argv[1:] :
	retcode = os.system("./xspp " + f + " axi |grep [load]")
	if retcode != 0:
		print 'error from xspp'
		exit()
	
	a=load('o_Vp')
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

#	savefig( "slide_" + str(i) + ".png" )
	savefig( "slide_%04d.png" % i )
	print '>>> slide', i, "saved.\n"
	i+=1
	close()

#retcode = os.system( 'mencoder "mf://slide_*.png" -mf fps=10 -o movie.avi -ovc lavc -lavcopts vcodec=mjpeg' )
