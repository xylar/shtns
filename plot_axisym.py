#!/usr/bin/python
#plot a slice
from pylab import *     # matplotlib
import sys              # acces a la ligne de commande.

ir0=1

#Up
print 'loading',sys.argv[1]
a=load(sys.argv[1])
s = a.shape
print s

ct = a[0,1:s[1]]
st = sqrt(1-ct*ct)
r = transpose(matrix(a[ir0:s[0],0]))

a = a[ir0:s[0],1:s[1]]
x = r*matrix(st)
y = r*matrix(ct)

#Poloidal scalar
print 'loading',sys.argv[2]
b = load(sys.argv[2])
b = b[ir0:s[0],1:s[1]]

#convert Up to angular velocity
c=a/array(x)
c[:,s[1]-2] = c[:,s[1]-3]	#remove nan
c[:,0] = c[:,1]			#remove nan

#colormap for phi component
#pcolor(array(x),array(y),c,shading='interp')
contourf(array(x),array(y),c,15,cmap=cm.Spectral)
colorbar()

#contour for poloidal scalar
contour(array(x),array(y),b,15,colors='k')

axis('equal')
axis('off')
xlim(0,1)
ylim(0,1)
savefig('axisym.png')
#savefig('axisym.pdf')
show()
