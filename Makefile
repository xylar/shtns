# definit le chemin contenant les programmes
#p = ~/Recherche/Progs/VxCart

# le compilateur :
cmd = gcc -O3

sphere : sphere.c SHT.c SHT.h Makefile
	$(cmd) sphere.c -lfftw3 -lgsl -lgslcblas -lm -o sphere
time_SHT : time_SHT.c SHT.c SHT.h Makefile
	$(cmd) time_SHT.c -lfftw3 -lgsl -lgslcblas -lm -o time_SHT
sphere2 : sphere2.c SHTfast.c SHT.h Makefile
	$(cmd) sphere2.c -lfftw3 -lgsl -lgslcblas -lm -o sphere2
sphshell : sphshell.c SHT.c SHT.h grid.c Makefile
	$(cmd) sphshell.c -lfftw3 -lgsl -lgslcblas -lm -o sphshell

