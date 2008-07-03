#!/usr/bin/python
#plot a slice
from pylab import *     # matplotlib
import sys              # acces a la ligne de commande.

print 'loading',sys.argv[1]
a=load(sys.argv[1])
s = a.shape

ct = a[0,1:s[1]]
st = sqrt(1-ct*ct)
r = transpose(matrix(a[1:s[0],0]))

a = a[1:s[0],1:s[1]]
x = r*matrix(st)
y = r*matrix(ct)

print a.shape

pcolor(array(x),array(y),a,shading='interp')
axis('equal')
colorbar()
show()
