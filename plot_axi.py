#!/usr/bin/python
#plot a slice
from pylab import *     # matplotlib
from subprocess import *	# lance une ligne de commande
import sys              # acces a la ligne de commande.

ir0=1

retcode = call("./xspp " + sys.argv[1] + " axi", shell=True)
if retcode != 0:
	print 'error from xspp'
	exit()

#Up
#print 'loading',sys.argv[1]
#a=load(sys.argv[1])
a=load('o_Vp')
s = a.shape

ct = a[0,1:s[1]]
st = sqrt(1-ct*ct)
r = transpose(matrix(a[ir0:s[0],0]))

a = a[ir0:s[0],1:s[1]]
x = r*matrix(st)
y = r*matrix(ct)

#Poloidal scalar
#print 'loading',sys.argv[2]
#b = load(sys.argv[2])
b = load('o_Vpol')
b = b[ir0:s[0],1:s[1]]

#convert Up to angular velocity
c=a/array(x)
c[:,s[1]-2] = c[:,s[1]-3]	#remove nan
c[:,0] = c[:,1]			#remove nan

m=amax(abs(c))
print 'max angular velocity (absolute value) =',m

#colormap for phi component
#pcolor(array(x),array(y),c,shading='interp')
contourf(array(x),array(y),c,15,cmap=cm.RdBu)
colorbar()
clim(-m,m)

#contour for poloidal scalar
m=amax(abs(b))
contour(array(x),array(y),b,arange(m/6,m,m/3),colors='k')
contour(array(x),array(y),b,arange(-m/6,-m,-m/3),colors='k')
#axvline(x=0, ymin=-1, ymax=1, color='k')

axis('equal')
axis('off')
xlim(0,1)
ylim(0,1)
savefig('axisym.png')
#savefig('axisym.pdf')
show()
