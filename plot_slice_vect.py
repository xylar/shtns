#!/usr/bin/python
#plot a slice
from pylab import *     # matplotlib
import sys              # acces a la ligne de commande.

ir0=1

#Ur
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

#Ut
print 'loading',sys.argv[2]
b = load(sys.argv[2])
b = b[ir0:s[0],1:s[1]]

#Up
print 'loading',sys.argv[3]
c = load(sys.argv[3])
c = c[ir0:s[0],1:s[1]]

print c.shape

#convert Up to angular velocity
c=c/array(x)

#colormap for phi component
#pcolor(array(x),array(y),c,shading='interp')
contourf(array(x),array(y),c,15,cmap=cm.Spectral)
colorbar()

#quiver for r,theta
l = s[0]
m = s[1]
s = l/11
t = (m+15)/16
print s,t

#cartesian coordinates
ux = a*st + b*ct
uy = a*ct - b*st

quiver(array(x[s:l:s,0:m:t]),array(y[s:l:s,0:m:t]),ux[s:l:s,0:m:t],uy[s:l:s,0:m:t])
#quiver(array(x[l/3:l:s,:]),array(y[l/3:l:s,:]),ux[l/3:l:s,:],uy[l/3:l:s,:])

axis('equal')
savefig('slice.png')
#savefig('slice.pdf')
show()
