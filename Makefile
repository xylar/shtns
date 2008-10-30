# le compilateur :
## generic gcc
cmd = gcc -O3
## profiling
cmd = gcc -O3 -p -fno-inline
## recent gcc with native support
#cmd = gcc -O3 -march=native -mfpmath=sse
## gcc k8 (lgitl3)
#cmd = gcc -O3 -march=k8
## gcc core2 (calcul1&2)
#cmd = gcc -O3 -march=core2 -mfpmath=sse
## icare 64bits opteron
#cmd = cc -fast -xarch=amd64 -I/users/nschaeff/include -L/users/nschaeff/lib
#cmdp = cc -fast -xarch=amd64 -xopenmp=parallel -I/users/nschaeff/include -L/users/nschaeff/lib

shtfiles = SHT.c SHT/SH_to_spat.c SHT/spat_to_SH.c SHT/SH_to_spat.gen.c SHT/Makefile
ini = inc_B0ini.c inc_U0ini.c

default : xshells

SHT/SH_to_spat.c : SHT/SH_to_spat.gen.c SHT/Makefile
	$(MAKE) -C SHT
SHT/spat_to_SH.c : SHT/spat_to_SH.gen.c SHT/Makefile
	$(MAKE) spat_to_SH.c -C SHT

xshells : xshells.c SHT.h grid.c xshells_fields.c xshells_io.c Makefile $(shtfiles) $(ini)
	$(cmd) xshells.c -lfftw3 -lgsl -lgslcblas -lm -o xshells
xshells_imp : xshells.c SHT.h grid.c xshells_fields.c xshells_io.c Makefile $(shtfiles) $(ini)
	$(cmd) xshells.c -D_IMPULSE_ -lfftw3 -lgsl -lgslcblas -lm -o xshells_imp
xspp : xspp.c grid.c xshells_fields.c xshells_io.c Makefile $(shtfiles)
	$(cmd) xspp.c -lfftw3 -lgsl -lgslcblas -lm -o xspp

sphere : sphere.c SHT.c SHT.h Makefile
	$(cmd) sphere.c -lfftw3 -lgsl -lgslcblas -lm -o sphere
time_SHT : time_SHT.c SHT.h Makefile $(shtfiles)
	$(cmd) time_SHT.c -lfftw3 -lgsl -lgslcblas -lm -o time_SHT
sphere2 : sphere2.c SHTfast.c SHT.h Makefile
	$(cmd) sphere2.c -lfftw3 -lgsl -lgslcblas -lm -o sphere2
sphshell : sphshell.c SHT.c SHT.h grid.c Makefile
	$(cmd) sphshell.c -lfftw3 -lgsl -lgslcblas -lm -o sphshell
dyncin : dyncin.c SHT.c SHT.h grid.c Makefile
	$(cmd) dyncin.c -lfftw3 -lgsl -lgslcblas -lm -o dyncin

