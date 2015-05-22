dnl Process this file with autoconf to produce a configure script.
AC_INIT([denise], [0.1.0])
AM_INIT_AUTOMAKE([denise], [0.1.0])

dnl Check for options.
AC_ARG_ENABLE(fft, AC_HELP_STRING([--enable-fft],
   [Compile FFT support (default=no)]), , enable_fft=no)
AM_CONDITIONAL(ENABLE_FFT, test x$enable_fft = xyes)

AC_ARG_ENABLE(nwp, AC_HELP_STRING([--enable-nwp],
   [Compile NWP support (default=yes)]), , enable_nwp=yes)
AM_CONDITIONAL(ENABLE_NWP, test x$enable_nwp = xyes)

AC_ARG_ENABLE(cairomm, AC_HELP_STRING([--enable-cairomm],
   [Compile CAIROMM support (default=yes)]), , enable_cairomm=yes)
AM_CONDITIONAL(ENABLE_CAIROMM, test x$enable_cairomm = xyes)

AC_ARG_ENABLE(gtkmm, AC_HELP_STRING([--enable-gtkmm],
   [Compile GTKMM support (default=yes)]), , enable_gtkmm=yes)
AM_CONDITIONAL(ENABLE_GTKMM, test x$enable_gtkmm = xyes)

dnl Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LIBTOOL
AC_C_BIGENDIAN
AC_LANG(C++)
AC_CHECK_LIBM

dnl Checks for libraries.
PKG_CHECK_MODULES([ZLIB], [zlib])
PKG_CHECK_MODULES([GSL], [gsl])

DENISE_REQUIRES="$DENISE_REQUIRES zlib gsl"
DENISE_CFLAGS="$ZLIB_CFLAGS $GSL_CFLAGS $DENISE_CFLAGS"
DENISE_LIBS="$ZLIB_LIBS $GSL_LIBS $DENISE_LIBS"

if test x"$enable_cairomm" = xyes
then
   PKG_CHECK_MODULES([CAIROMM_1_0], [cairomm-1.0])
   DENISE_REQUIRES="$DENISE_REQUIRES cairomm-1.0"
   DENISE_CFLAGS="$CAIROMM_1_0_CFLAGS $DENISE_CFLAGS"
   DENISE_LIBS="$CAIROMM_1_0_LIBS $DENISE_LIBS"
fi

if test x"$enable_gtkmm" = xyes
then
   PKG_CHECK_MODULES([GTKMM_3_0], [gtkmm-3.0])
   DENISE_REQUIRES="$DENISE_REQUIRES gtkmm-3.0"
   DENISE_CFLAGS="$GTKMM_3_0_CFLAGS $DENISE_CFLAGS"
   DENISE_LIBS="$GTKMM_3_0_LIBS $DENISE_LIBS"
fi

if test x"$enable_fft" = xyes
then
   PKG_CHECK_MODULES([FFTW3], [fftw3])
   DENISE_REQUIRES="$DENISE_REQUIRES fftw3"
   DENISE_CFLAGS="$FFTW3_CFLAGS $DENISE_CFLAGS"
   DENISE_LIBS="$FFTW3_LIBS $DENISE_LIBS"
fi

if test x"$enable_nwp" = xyes
then
   AC_CHECK_HEADERS(jasper/jasper.h)
   AC_CHECK_LIB(jasper, main, [DENISE_LIBS="-ljasper $DENISE_LIBS"],
      AC_MSG_ERROR(jasper not found), )
   PKG_CHECK_MODULES([NETCDF], [netcdf])
   DENISE_REQUIRES="$DENISE_REQUIRES netcdf"
   DENISE_CFLAGS="$DENISE_CFLAGS $NETCDF_CFLAGS"
   DENISE_LIBS="$DENISE_LIBS $NETCDF_LIBS"
fi

dnl Checks for header files.
AC_CHECK_HEADERS(ieeefp.h)
AC_CHECK_HEADERS(cstdint)

dnl Checks for typedefs, structures, and compiler characteristics.

dnl Checks for library functions.
AC_CHECK_FUNCS(exp10, [DENISE_CFLAGS="-DHAVE_EXP10=1 $DENISE_CFLAGS"])
AC_CHECK_FUNCS(timezone)

CPPFLAGS="$DENISE_CFLAGS $CPPFLAGS"
LIBS="$DENISE_LIBS $LIBS"

AC_SUBST(DENISE_REQUIRES)
AC_SUBST(DENISE_CFLAGS)
AC_SUBST(DENISE_LIBS)

AC_CONFIG_FILES([denise.pc Makefile denise/Makefile])
AC_OUTPUT
