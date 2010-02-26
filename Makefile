## "version" identification string
HGID=`hg id -ti`
## path prefix for make install (installs to $(PREFIX)/lib and $(PREFIX)/include)
PREFIX=$(HOME)

## global options for gcc
## there should be -ffast-math or at least -fcx-limited-range to produce fast code.
go= -O3 -std=gnu99 -ffast-math -D_GNU_SOURCE

## compiler :
## generic gcc
cmd = gcc $(go)
## profiling
#cmd = gcc $(go) -p -fno-inline
## recent gcc with native support
#cmd = gcc $(go) -march=native -mfpmath=sse
## gcc k8 (lgitl3)
#cmd = gcc $(go) -march=k8
## gcc core2 (calcul1&2)
#cmd = gcc $(go) -march=core2 -mfpmath=sse
## icare 64bits opteron
#cmd = cc -fast -xarch=amd64 -I/users/nschaeff/include -L/users/nschaeff/lib
## r2d2
#cmd = gcc $(go) -march=core2 -m64 -I/home/ciment/nschaeff/include -L/home/ciment/nschaeff/lib

# intel compiler is used for codelets (leads to faster code)
shtcc = icc -axT -xT -O3 -prec-div -complex-limited-range -D_HGID_="\"$(HGID)\""
# if icc is not available, use gcc instead by uncommenting the following line.
shtcc = $(cmd)

shtfiles = SHT/SH_to_spat.c SHT/spat_to_SH.c SHT/SHeo_to_spat.c SHT/spat_to_SHeo.c SHT/hyb_SH_to_spat.gen.c SHT/hyb_spat_to_SH.gen.c SHT/sparse_spat_to_SH.gen.c SHT/sparse_SH_to_spat.gen.c SHT/Makefile sht_legendre.c

hfiles = sht_private.h sht_config.h shtns.h

default : libshtns.a

libshtns.a : Makefile SHT.o sht_std.o sht_ltr.o sht_m0.o sht_eo.o sht_m0ltr.o
	ar rcs libshtns.a SHT.o sht_std.o sht_ltr.o sht_m0.o sht_eo.o sht_m0ltr.o

install :
	cp libshtns.a $(PREFIX)/lib
	cp shtns.h $(PREFIX)/include
	cp shtns.f $(PREFIX)/include

# codelets :
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

sht_std.o : sht_std.c Makefile $(hfiles) SHT/sht_generic.c SHT/SH_to_spat.c SHT/spat_to_SH.c
	$(shtcc) -c sht_std.c -o sht_std.o
sht_ltr.o : sht_ltr.c Makefile $(hfiles) SHT/sht_generic.c SHT/SH_to_spat.c SHT/spat_to_SH.c
	$(shtcc) -c sht_ltr.c -o sht_ltr.o
sht_m0.o : sht_m0.c Makefile $(hfiles) SHT/sht_generic.c SHT/SH_to_spat.c SHT/spat_to_SH.c
	$(shtcc) -c sht_m0.c -o sht_m0.o
sht_m0ltr.o : sht_m0ltr.c Makefile $(hfiles) SHT/sht_generic.c SHT/SH_to_spat.c SHT/spat_to_SH.c
	$(shtcc) -c sht_m0ltr.c -o sht_m0ltr.o
sht_eo.o : sht_eo.c Makefile $(hfiles) SHT/SHeo_to_spat.c SHT/spat_to_SHeo.c
	$(shtcc) -c sht_eo.c -o sht_eo.o

# programs :
time_SHT : shtns.h time_SHT.c libshtns.a Makefile
	$(cmd) time_SHT.c libshtns.a -lfftw3 -lm -o time_SHT

SHT_example : SHT_example.c libshtns.a Makefile shtns.h
	$(cmd) SHT_example.c libshtns.a -lfftw3 -lm -o SHT_example

SHT_fort_ex : SHT_example.f libshtns.a Makefile shtns.f
	gfortran -fdefault-real-8 SHT_example.f libshtns.a -lfftw3 -lm -lc -o SHT_fort_ex

#documentation :
docs :
	doxygen doxygen.conf

clean :
	$(MAKE) clean -C SHT
	rm -f *.o
	rm -rf doc/

#fftw compiling options :
#-O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math -fno-schedule-insns -fno-web -fno-loop-optimize --param inline-unit-growth=1000 --param large-function-growth=1000
