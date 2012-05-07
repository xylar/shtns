#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([SHTns],[2.2],[],[shtns],[https://bitbucket.org/nschaeff/shtns])
AC_CONFIG_SRCDIR([SHT.c])
AC_CONFIG_HEADERS([sht_config.h])

# optional features with --enable-XXX
AC_ARG_ENABLE([verbose],
	AS_HELP_STRING([--enable-verbose=0-3], [define verbosity level for the library: 0=silent, 1=default, 2=debug, 3=full]),
	[],[enable_verbose=1])
AC_ARG_ENABLE([openmp],
	AS_HELP_STRING([--enable-openmp], [Enable multi-threading with OpenMP]))
AC_ARG_ENABLE([long-double],
	AS_HELP_STRING([--enable-long-double], [Use long double during initialization for (maybe) slightly better precision]))
AC_ARG_ENABLE([f77],
	AS_HELP_STRING([--disable-f77], [Do not include F77 wrapper to call SHTns library from Fortran]))
AC_ARG_ENABLE([simd],
	AS_HELP_STRING([--disable-simd], [Do not use SSE2 or AVX extensions]))

# Checks for programs.
if test "x$CFLAGS" = "x"; then
	CFLAGS="-O2"
fi
AC_PROG_CC
AC_PROG_SED
if test "$SED" = :; then
	AC_MSG_ERROR([sed program required.])
fi

# define macrao AX_CHECK_COMPILE_FLAG
AC_DEFUN([AX_CHECK_COMPILE_FLAG],
[AS_VAR_PUSHDEF([CACHEVAR],[ax_cv_check_[]_AC_LANG_ABBREV[]flags_$4_$1])dnl
AC_CACHE_CHECK([whether _AC_LANG compiler accepts $1], CACHEVAR, [
  ax_check_save_flags=$[]_AC_LANG_PREFIX[]FLAGS
  _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $4 $1"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
    [AS_VAR_SET(CACHEVAR,[yes])],
    [AS_VAR_SET(CACHEVAR,[no])])
  _AC_LANG_PREFIX[]FLAGS=$ax_check_save_flags])
AS_IF([test x"AS_VAR_GET(CACHEVAR)" = xyes],
  [m4_default([$2], :)],
  [m4_default([$3], :)])
AS_VAR_POPDEF([CACHEVAR])dnl
])dnl AX_CHECK_COMPILE_FLAGS

# add gcc compile options if supported.
AX_CHECK_COMPILE_FLAG([-march=native],[CC="$CC -march=native"])
AX_CHECK_COMPILE_FLAG([-ffast-math],[CFLAGS="$CFLAGS -ffast-math"])
AX_CHECK_COMPILE_FLAG([-fomit-frame-pointer],[CFLAGS="$CFLAGS -fomit-frame-pointer"])
AX_CHECK_COMPILE_FLAG([-std=gnu99],[CFLAGS="$CFLAGS -std=gnu99"])

# Checks for libraries.
AC_CHECK_LIB([m],[cos],,AC_MSG_ERROR([math library not found.]))
AC_CHECK_LIB([fftw3],[fftw_plan_many_dft],,AC_MSG_ERROR([FFTW3 library required.]))

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h stdio.h string.h math.h complex.h fftw3.h])

# Checks related to OpenMP
AS_IF([test "x$enable_openmp" = "xyes"], [
	AC_OPENMP
	AC_CHECK_LIB([fftw3_omp], [fftw_init_threads])
	dnl AC_CHECK_HEADERS([omp.h])
])

# Checks related to long double
AS_IF([test "x$enable_long_double" = "xyes"], [
	AC_TYPE_LONG_DOUBLE_WIDER])

# Disable Fortran interface ?
AS_IF([test "x$enable_f77" != "xno"], [
	AC_DEFINE([SHT_F77_API],[1],[Compile the Fortran API])
])

# Disable SIMD ?
AS_IF([test "x$enable_simd" != "xno"], [
	AC_DEFINE([_GCC_VEC_],[1],[I compile with GCC 4 or later, and I would like fast vectorized code (if SSE2 or AVX is supported) !])
])


# Verbosity setting
AS_IF([test "x$enable_verbose" == "xno"], [enable_verbose=0])
AS_IF([test "x$enable_verbose" == "xyes"], [enable_verbose=1])
AC_DEFINE_UNQUOTED([SHT_VERBOSE],$enable_verbose,[0:no output, 1:output info to stdout, 2:more output (debug info), 3:also print fftw plans.])
echo "verbose level = $enable_verbose"

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
#AC_FUNC_MALLOC
#AC_FUNC_REALLOC
AC_CHECK_FUNCS([fftw_cost])

### for cycle.h ###
   AC_C_INLINE
   AC_HEADER_TIME
   AC_CHECK_HEADERS([sys/time.h c_asm.h intrinsics.h mach/mach_time.h])

   AC_CHECK_TYPE([hrtime_t],[AC_DEFINE(HAVE_HRTIME_T, 1, [Define to 1 if hrtime_t is defined in <sys/time.h>])],,[#if HAVE_SYS_TIME_H
#include <sys/time.h>
#endif])

   AC_CHECK_FUNCS([gethrtime read_real_time time_base_to_time clock_gettime mach_absolute_time])

   dnl Cray UNICOS _rtc() (real-time clock) intrinsic
   AC_MSG_CHECKING([for _rtc intrinsic])
   rtc_ok=yes
   AC_TRY_LINK([#ifdef HAVE_INTRINSICS_H
#include <intrinsics.h>
#endif], [_rtc()], [AC_DEFINE(HAVE__RTC,1,[Define if you have the UNICOS _rtc() intrinsic.])], [rtc_ok=no])
   AC_MSG_RESULT($rtc_ok)
### end cycle.h ###

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
