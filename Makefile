# le compilateur :
## generic gcc
cmd = gcc -O3
## recent gcc with native support
#cmd = gcc -O3 -march=native -mfpmath=sse
## gcc k8 (lgitl3)
#cmd = gcc -O3 -march=k8
## gcc core2 (calcul1&2)
#cmd = gcc -O3 -march=core2 -mfpmath=sse
## icare 64bits opteron
#cmd = cc -fast -xarch=amd64 -I/users/nschaeff/include -L/users/nschaeff/lib
#cmdp = cc -fast -xarch=amd64 -xopenmp=parallel -I/users/nschaeff/include -L/users/nschaeff/lib


default: xshells

xshells : xshells.c SHT.c SHT.h grid.c Makefile
	$(cmd) xshells.c -lfftw3 -lgsl -lgslcblas -lm -o xshells
sphere : sphere.c SHT.c SHT.h Makefile
	$(cmd) sphere.c -lfftw3 -lgsl -lgslcblas -lm -o sphere
time_SHT : time_SHT.c SHT.c SHT.h Makefile
	$(cmd) time_SHT.c -lfftw3 -lgsl -lgslcblas -lm -o time_SHT
sphere2 : sphere2.c SHTfast.c SHT.h Makefile
	$(cmd) sphere2.c -lfftw3 -lgsl -lgslcblas -lm -o sphere2
sphshell : sphshell.c SHT.c SHT.h grid.c Makefile
	$(cmd) sphshell.c -lfftw3 -lgsl -lgslcblas -lm -o sphshell
dyncin : dyncin.c SHT.c SHT.h grid.c Makefile
	$(cmd) dyncin.c -lfftw3 -lgsl -lgslcblas -lm -o dyncin

