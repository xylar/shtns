## path prefix for make install (installs to $(PREFIX)/lib and $(PREFIX)/include)
PREFIX=$(HOME)

## comment to disable parallel spherical harmonic transforms using openmp
#OMP=-fopenmp -lfftw3_omp

# compile for a given architecture : native is the best (gcc >= 4.0), generic if you don't know.
#march = -march=core2
march = -march=native

## global options for gcc
## there should be -ffast-math or at least -fcx-limited-range to produce fast code.
go= $(march) -pipe -fomit-frame-pointer -std=gnu99 -D_GNU_SOURCE -fpic $(OMP)

# intel compiler may be used for codelets
#shtcc = icc -axT -xT -O3 -prec-div -complex-limited-range
# gcc + vector intrinsic leads to faster code (with _GCC_VEC_ set to 1 in sht_config.h)
# gcc compiler command with options for the sht codelets
shtcc = gcc $(go) -O3 -ffast-math
# gcc compiler command with options for other source (initialization, ...)
cc = gcc $(go) -O2 -ffast-math

## "version" identification string
HGID=`hg id -ti`

hfiles = sht_private.h sht_config.h shtns.h

default : libshtns.a

libshtns.a : Makefile SHT.o sht_std.o sht_ltr.o sht_m0.o sht_m0ltr.o
	ar rcs libshtns.a SHT.o sht_std.o sht_ltr.o sht_m0.o sht_m0ltr.o
	@echo " "
	@cat COPYRIGHT

install :
	@mkdir -p $(PREFIX)/lib/
	@mkdir -p $(PREFIX)/include/
	cp libshtns.a $(PREFIX)/lib/
	cp shtns.h $(PREFIX)/include/
	cp shtns.f $(PREFIX)/include/
	@echo " "
	@cat COPYRIGHT

codelets1 = SHT/SH_to_spat.c SHT/SH_to_spat_fly.c SHT/spat_to_SH.c SHT/spat_to_SH_fly.c SHT/sht_generic.c
ifeq ($(strip $(OMP)),)
  codelets = $(codelets1)
else
  codelets = $(codelets1) SHT/SH_to_spat_omp.c SHT/spat_to_SH_omp.c
endif

# codelets :
SHT/SH_to_spat.c : SHT/hyb_SH_to_spat.gen.c
	$(MAKE) SH_to_spat.c -C SHT 
SHT/spat_to_SH.c : SHT/hyb_spat_to_SH.gen.c
	$(MAKE) spat_to_SH.c -C SHT 
SHT/SH_to_spat_fly.c : SHT/fly_SH_to_spat.gen.c 
	$(MAKE) SH_to_spat_fly.c -C SHT SFX=fly
SHT/spat_to_SH_fly.c : SHT/fly_spat_to_SH.gen.c 
	$(MAKE) spat_to_SH_fly.c -C SHT SFX=fly
SHT/SH_to_spat_omp.c : SHT/omp_SH_to_spat.gen.c
	$(MAKE) SH_to_spat_omp.c -C SHT SFX=omp
SHT/spat_to_SH_omp.c : SHT/omp_spat_to_SH.gen.c
	$(MAKE) spat_to_SH_omp.c -C SHT SFX=omp

# objects :
SHT.o : SHT.c Makefile sht_legendre.c $(hfiles) cycle.h
	$(cc) -D_HGID_="\"$(HGID)\"" -c SHT.c -o SHT.o
	@echo "DONE SHT.o"

sht_std.o : sht_std.c Makefile $(hfiles) $(codelets)
	$(shtcc) -c sht_std.c -o sht_std.o
	@echo "DONE sht_std"
sht_ltr.o : sht_ltr.c Makefile $(hfiles) $(codelets)
	$(shtcc) -c sht_ltr.c -o sht_ltr.o
	@echo "DONE sht_ltr"
sht_m0.o : sht_m0.c Makefile $(hfiles) $(codelets)
	$(shtcc) -c sht_m0.c -o sht_m0.o
	@echo "DONE sht_m0"
sht_m0ltr.o : sht_m0ltr.c Makefile $(hfiles) $(codelets)
	$(shtcc) -c sht_m0ltr.c -o sht_m0ltr.o
	@echo "DONE sht_m0ltr"

# programs :
time_SHT : shtns.h time_SHT.c libshtns.a Makefile
	$(cc) time_SHT.c -I$(PREFIX)/include -L$(PREFIX)/lib ./libshtns.a -lfftw3 -lm -lrt -o time_SHT
test_rot : shtns.h test_rot.c libshtns.a Makefile
	$(cc) test_rot.c -I$(PREFIX)/include -L$(PREFIX)/lib ./libshtns.a -lfftw3 -lm -o test_rot

SHT_example : SHT_example.c libshtns.a Makefile shtns.h
	$(cc) -I$(PREFIX)/include -L$(PREFIX)/lib SHT_example.c ./libshtns.a -lfftw3 -lm -o SHT_example

SHT_fort_ex : SHT_example.f libshtns.a Makefile shtns.f
	gfortran -fdefault-real-8 -I$(PREFIX)/include -L$(PREFIX)/lib SHT_example.f -lshtns -lfftw3 -lm -lc -o SHT_fort_ex

#documentation :
docs :
	doxygen doxygen.conf

clean :
	$(MAKE) clean -C SHT
	rm -f *.o
	rm -rf doc/html
	rm -rf doc/latex

# "make python" builds a python interface using SWIG.
# use it with "from shtns import *" in a python program/shell
python : shtns_numpy_wrap.c Makefile $(hfiles) SHT.o sht_std.o sht_ltr.o sht_m0.o sht_m0ltr.o
	gcc $(march) -O2 -fpic -I/usr/include/python2.7 -c shtns_numpy_wrap.c
	gcc $(march) -O2 -fpic -shared $(OMP) -lfftw3 SHT.o sht_std.o sht_ltr.o sht_m0.o sht_m0ltr.o shtns_numpy_wrap.o -o _shtns.so

# generate python and c glue code with SWIG.
shtns_numpy_wrap.c : shtns_numpy.i sht_private.h shtns.h
	swig -python shtns_numpy.i

# update the copyright notice
updatecpy : COPYRIGHT
	./update-copyright.sh shtns.h
	./update-copyright.sh SHT.c
	./update-copyright.sh sht_legendre.c
	./update-copyright.sh sht_config.h
	./update-copyright.sh sht_private.h
	./update-copyright.sh sht_std.c
	./update-copyright.sh sht_ltr.c
	./update-copyright.sh sht_m0.c
	./update-copyright.sh sht_m0ltr.c
	./update-copyright.sh SHT/sht_generic.c
	./update-copyright.sh SHT/hyb_SH_to_spat.gen.c
	./update-copyright.sh SHT/hyb_spat_to_SH.gen.c
	./update-copyright.sh SHT/fly_SH_to_spat.gen.c
	./update-copyright.sh SHT/fly_spat_to_SH.gen.c
	./update-copyright.sh SHT/omp_SH_to_spat.gen.c
	./update-copyright.sh SHT/omp_spat_to_SH.gen.c
	./update-copyright.sh time_SHT.c
	./update-copyright.sh SHT_example.c
	./update-copyright.sh -fortran SHT_example.f
	./update-copyright.sh -fortran shtns.f
	./update-copyright.sh shtns_numpy.i
	./update-copyright.sh -python SHT_example.py

.PHONY : updatecpy install clean docs python


#fftw compiling options :
#-O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math -fno-schedule-insns -fno-web -fno-loop-optimize --param inline-unit-growth=1000 --param large-function-growth=1000

