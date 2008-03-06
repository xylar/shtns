# definit le chemin contenant les programmes
#p = ~/Recherche/Progs/VxCart

# le compilateur :
cmd = gcc -O3

sphere : sphere.c Makefile
	$(cmd) sphere.c -lgsl -lgslcblas -lm -o sphere
