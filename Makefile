# le compilateur :
## generic gcc
cmd = gcc -O3 -ffast-math
cmdp = gcc -O3 -fopenmp -ffast-math
## profiling
#cmd = gcc -O3 -p -fno-inline
## recent gcc with native support
#cmd = gcc -O3 -march=native -mfpmath=sse
## gcc k8 (lgitl3)
#cmd = gcc -O3 -march=k8
## gcc core2 (calcul1&2)
#cmd = gcc -O3 -march=core2 -mfpmath=sse
## icare 64bits opteron
#cmd = cc -fast -xarch=amd64 -I/users/nschaeff/include -L/users/nschaeff/lib
#cmdp = cc -fast -xarch=amd64 -xopenmp=parallel -I/users/nschaeff/include -L/users/nschaeff/lib
## r2d2
#cmd = gcc -march=core2 -O3 -ffast-math -m64 -I/home/ciment/nschaeff/include -L/home/ciment/nschaeff/lib
#cmdp = gcc -march=core2 -O3 -fopenmp -ffast-math -m64 -I/home/ciment/nschaeff/include -L/home/ciment/nschaeff/lib


NTH=8 #number of threads for parallel version (can be overwritten with command line)
shtfiles = SHT.c SHT/SH_to_spat.c SHT/spat_to_SH.c SHT/dct_SH_to_spat.c SHT/SH_to_spat.gen.c SHT/spat_to_SH.gen.c SHT/dct_SH_to_spat.gen.c SHT/Makefile
ini = xshells.h

default : xshells

SHT/SH_to_spat.c : SHT/SH_to_spat.gen.c SHT/Makefile
	$(MAKE) SH_to_spat.c -C SHT
SHT/spat_to_SH.c : SHT/spat_to_SH.gen.c SHT/Makefile
	$(MAKE) spat_to_SH.c -C SHT
SHT/dct_SH_to_spat.c : SHT/dct_SH_to_spat.gen.c SHT/Makefile
	$(MAKE) dct_SH_to_spat.c -C SHT

xshells : xshells.c SHT.h grid.c xshells_fields.c xshells_io.c Makefile $(shtfiles) $(ini)
	$(cmd) xshells.c -lfftw3 -lgsl -lgslcblas -lm -o xshells
pxshells : xshells.c SHT.h grid.c xshells_fields.c xshells_io.c Makefile $(shtfiles) $(ini)
	$(cmdp) xshells.c -D_NTH_=$(NTH) -lfftw3 -lgsl -lgslcblas -lm -o pxshells

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


#fftw compiling options :
#-O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math -fno-schedule-insns -fno-web -fno-loop-optimize --param inline-unit-growth=1000 --param large-function-growth=1000
