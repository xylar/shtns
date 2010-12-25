## "version" identification string
HGID=`hg id -ti`
## path prefix for make install (installs to $(PREFIX)/lib and $(PREFIX)/include)
PREFIX=$(HOME)

## global options for gcc
## there should be -ffast-math or at least -fcx-limited-range to produce fast code.
go= -O3 -fomit-frame-pointer -std=gnu99 -ffast-math -D_GNU_SOURCE

## compiler :
## generic gcc
cmd = gcc $(go)
## profiling
#cmd = gcc $(go) -p -fno-inline
## recent gcc with native support
cmd = gcc $(go) -march=native
## gcc k8 (lgitl3)
#cmd = gcc $(go) -march=k8
## gcc core2 (calcul1&2)
#cmd = gcc $(go) -march=core2 -I$(PREFIX)/include -L$(PREFIX)/lib
## icare 64bits opteron
#cmd = cc -fast -xarch=amd64 -I$(PREFIX)/include -L$(PREFIX)/lib
## r2d2
#cmd = gcc $(go) -march=core2 -m64 -I$(PREFIX)/include -L$(PREFIX)/lib

# intel compiler may be used for codelets
#shtcc = icc -axT -xT -O3 -prec-div -complex-limited-range -D_HGID_="\"$(HGID)\""
# gcc + vector intrinsic leads to faster code (with _GCC_VEC_ set to 1 in sht_config.h)
shtcc = $(cmd)

shtfiles = SHT/SH_to_spat_fly.c SHT/fly_SH_to_spat.gen.c SHT/SH_to_spat.c SHT/spat_to_SH.c SHT/SHeo_to_spat.c SHT/spat_to_SHeo.c SHT/hyb_SH_to_spat.gen.c SHT/hyb_spat_to_SH.gen.c SHT/sparse_spat_to_SH.gen.c SHT/sparse_SH_to_spat.gen.c SHT/Makefile sht_legendre.c

hfiles = sht_private.h sht_config.h shtns.h

default : libshtns.a

libshtns.a : Makefile SHT.o sht_std.o sht_ltr.o sht_m0.o sht_eo.o sht_m0ltr.o
	ar rcs libshtns.a SHT.o sht_std.o sht_ltr.o sht_m0.o sht_eo.o sht_m0ltr.o
	@echo " "
	@cat COPYRIGHT

install :
	cp libshtns.a $(PREFIX)/lib/
	cp shtns.h $(PREFIX)/include/
	cp shtns.f $(PREFIX)/include/
	@echo " "
	@cat COPYRIGHT

# codelets :
SHT/SH_to_spat_fly.c : SHT/fly_SH_to_spat.gen.c SHT/Makefile
	$(MAKE) SH_to_spat_fly.c -C SHT
SHT/SH_to_spat.c : SHT/hyb_SH_to_spat.gen.c SHT/Makefile
	$(MAKE) SH_to_spat.c -C SHT
SHT/spat_to_SH.c : SHT/hyb_spat_to_SH.gen.c SHT/Makefile
	$(MAKE) spat_to_SH.c -C SHT
SHT/SHeo_to_spat.c : SHT/sparse_SH_to_spat.gen.c SHT/Makefile
	$(MAKE) SHeo_to_spat.c -C SHT
SHT/spat_to_SHeo.c : SHT/sparse_spat_to_SH.gen.c SHT/Makefile
	$(MAKE) spat_to_SHeo.c -C SHT

# objects :
SHT.o : SHT.c Makefile sht_legendre.c $(hfiles) cycle.h
	$(cmd) -D_HGID_="\"$(HGID)\"" -c SHT.c -o SHT.o

sht_std.o : sht_std.c Makefile $(hfiles) SHT/sht_generic.c SHT/SH_to_spat.c SHT/spat_to_SH.c SHT/SH_to_spat_fly.c
	$(shtcc) -c sht_std.c -o sht_std.o
sht_ltr.o : sht_ltr.c Makefile $(hfiles) SHT/sht_generic.c SHT/SH_to_spat.c SHT/spat_to_SH.c SHT/SH_to_spat_fly.c
	$(shtcc) -c sht_ltr.c -o sht_ltr.o
sht_m0.o : sht_m0.c Makefile $(hfiles) SHT/sht_generic.c SHT/SH_to_spat.c SHT/spat_to_SH.c SHT/SH_to_spat_fly.c
	$(shtcc) -c sht_m0.c -o sht_m0.o
sht_m0ltr.o : sht_m0ltr.c Makefile $(hfiles) SHT/sht_generic.c SHT/SH_to_spat.c SHT/spat_to_SH.c SHT/SH_to_spat_fly.c
	$(shtcc) -c sht_m0ltr.c -o sht_m0ltr.o
sht_eo.o : sht_eo.c Makefile $(hfiles) SHT/SHeo_to_spat.c SHT/spat_to_SHeo.c
	$(shtcc) -c sht_eo.c -o sht_eo.o

# programs :
time_SHT : shtns.h time_SHT.c libshtns.a Makefile
	$(cmd) time_SHT.c libshtns.a -lfftw3 -lm -o time_SHT

SHT_example : SHT_example.c libshtns.a Makefile shtns.h
	$(cmd) -I$(PREFIX)/include -L$(PREFIX)/lib SHT_example.c -lshtns -lfftw3 -lm -o SHT_example

SHT_fort_ex : SHT_example.f libshtns.a Makefile shtns.f
	gfortran -fdefault-real-8 -I$(PREFIX)/include -L$(PREFIX)/lib SHT_example.f -lshtns -lfftw3 -lm -lc -o SHT_fort_ex

#documentation :
docs :
	doxygen doxygen.conf

clean :
	$(MAKE) clean -C SHT
	rm -f *.o
	rm -rf doc/

# build a python interface using SWIG.
# use it with "from shtns import *" in a python program/shell
python : shtns.h shtns.i
	swig -python shtns.i
	gcc -fpic -I/usr/include/python2.6 -c shtns_wrap.c 
	gcc -shared /usr/lib/libfftw3.so SHT.o sht_*.o shtns_wrap.o -o _shtns.so

# update the copyright notice
updatecpy : copyright
	./update-copyright.sh shtns.h
	./update-copyright.sh SHT.c
	./update-copyright.sh sht_legendre.c
	./update-copyright.sh sht_config.h
	./update-copyright.sh sht_private.h
	./update-copyright.sh sht_std.c
	./update-copyright.sh sht_ltr.c
	./update-copyright.sh sht_m0.c
	./update-copyright.sh sht_m0ltr.c
	./update-copyright.sh sht_eo.c
	./update-copyright.sh SHT/sht_generic.c
	./update-copyright.sh SHT/hyb_SH_to_spat.gen.c
	./update-copyright.sh SHT/hyb_spat_to_SH.gen.c
	./update-copyright.sh SHT/sparse_SH_to_spat.gen.c
	./update-copyright.sh SHT/sparse_spat_to_SH.gen.c
	./update-copyright.sh time_SHT.c
	./update-copyright.sh SHT_example.c
	./update-copyright.sh -fortran SHT_example.f
	./update-copyright.sh -fortran shtns.f
	./update-copyright.sh shtns.i

#fftw compiling options :
#-O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math -fno-schedule-insns -fno-web -fno-loop-optimize --param inline-unit-growth=1000 --param large-function-growth=1000
