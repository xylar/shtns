#!/usr/bin/python
#plot a slice
from pylab import *     # matplotlib
import sys              # acces a la ligne de commande.

#Ur
print 'loading',sys.argv[1]
a=load(sys.argv[1])
s = a.shape

ct = a[0,1:s[1]]
st = sqrt(1-ct*ct)
r = transpose(matrix(a[1:s[0],0]))

a = a[1:s[0],1:s[1]]
x = r*matrix(st)
y = r*matrix(ct)

#Ut
print 'loading',sys.argv[2]
b = load(sys.argv[2])
b = b[1:s[0],1:s[1]]

#Up
print 'loading',sys.argv[3]
c = load(sys.argv[3])
c = c[1:s[0],1:s[1]]

print c.shape

#colormap for phi component
pcolor(array(x),array(y),c,shading='interp')
colorbar()

#quiver for r,theta
l = s[0]
m = s[1]
s = l/11
t = m/16
print s,t

#cartesian coordinates
ux = a*st + b*ct
uy = a*ct - b*st

quiver(array(x[s:l:s,0:m:t]),array(y[s:l:s,0:m:t]),ux[s:l:s,0:m:t],uy[s:l:s,0:m:t])
#quiver(array(x[l/3:l:s,:]),array(y[l/3:l:s,:]),ux[l/3:l:s,:],uy[l/3:l:s,:])

axis('equal')
show()
