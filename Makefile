# definit le chemin contenant les programmes
#p = ~/Recherche/Progs/VxCart

# le compilateur :
cmd = gcc -O3

sphere : sphere.c SHT.c SHT.h Makefile
	$(cmd) sphere.c -lfftw3 -lgsl -lgslcblas -lm -o sphere
