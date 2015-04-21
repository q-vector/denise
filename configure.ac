dnl Process this file with autoconf to produce a configure script.
AC_INIT([denise], [808])
AM_INIT_AUTOMAKE([denise], [808])

dnl Check for options.
dnl AC_ARG_ENABLE(cairo_pdf, AC_HELP_STRING([--enable-cairo-pdf],
dnl    [Compile CAIRO_PDF support (default=no)]), , enable_cairo_pdf=no)
dnl AM_CONDITIONAL(ENABLE_CAIRO_PDF, test x$enable_cairo_pdf = xyes)

dnl AC_ARG_ENABLE(cairo_ps, AC_HELP_STRING([--enable-cairo-ps],
dnl    [Compile CAIRO_PS support (default=no)]), , enable_cairo_ps=no)
dnl AM_CONDITIONAL(ENABLE_CAIRO_PS, test x$enable_cairo_ps = xyes)

dnl AC_ARG_ENABLE(cairo_svg, AC_HELP_STRING([--enable-cairo-svg],
 dnl   [Compile CAIRO_SVG support (default=no)]), , enable_cairo_svg=no)
dnl AM_CONDITIONAL(ENABLE_CAIRO_SVG, test x$enable_cairo_svg = xyes)

dnl AC_ARG_ENABLE(cairo_xlib, AC_HELP_STRING([--enable-cairo-xlib],
dnl    [Compile CAIRO_XLIB support (default=no)]), , enable_cairo_xlib=no)
dnl AM_CONDITIONAL(ENABLE_CAIRO_XLIB, test x$enable_cairo_xlib = xyes)

AC_ARG_ENABLE(gtkmm, AC_HELP_STRING([--enable-gtkmm],
   [Compile GTKMM support (default=yes)]), , enable_gtkmm=yes)
AM_CONDITIONAL(ENABLE_GTKMM, test x$enable_gtkmm = xyes)

AC_ARG_ENABLE(fftw3, AC_HELP_STRING([--enable-fftw3],
   [Compile FFTW3 support (default=yes)]), , enable_fftw3=yes)
AM_CONDITIONAL(ENABLE_FFTW3, test x$enable_fftw3 = xyes)

AC_ARG_ENABLE(sdl, AC_HELP_STRING([--enable-sdl],
   [Compile SDL support (default=no)]), , enable_sdl=no)
AM_CONDITIONAL(ENABLE_SDL, test x$enable_sdl = xyes)

dnl Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LIBTOOL
AC_C_BIGENDIAN
AC_LANG(C++)

AC_CHECK_LIBM

dnl Checks for libraries.
AC_CHECK_LIB(jasper, main, [DENISE_LIBS="-ljasper $DENISE_LIBS"],
   AC_MSG_ERROR(jasper not found), )
PKG_CHECK_MODULES([ZLIB], [zlib])
PKG_CHECK_MODULES([GSL], [gsl])
PKG_CHECK_MODULES([CAIRO], [cairo])
PKG_CHECK_MODULES([CAIRO_FT], [cairo-ft])
PKG_CHECK_MODULES([CAIRO_PNG], [cairo-png])
PKG_CHECK_MODULES([CAIROMM_1_0], [cairomm-1.0])

DENISE_REQUIRES="$DENISE_REQUIRES zlib gsl cairo cairo-ft cairo-png cairomm-1.0"
DENISE_CFLAGS="$CAIRO_CFLAGS $ZLIB_CFLAGS $CAIRO_FT_CFLAGS $CAIRO_PNG_CFLAGS $CAIROMM_1_0_CFLAGS $GSL_CFLAGS $DENISE_CFLAGS"
DENISE_LIBS="$CAIRO_LIBS $ZLIB_LIBS $CAIRO_FT_LIBS $CAIRO_PNG_LIBS $CAIROMM_1_0_LIBS $GSL_LIBS $DENISE_LIBS"

dnl if test x"$enable_cairo_pdf" = xyes
dnl then
dnl    PKG_CHECK_MODULES([CAIRO_PDF], [cairo-pdf])
dnl    DENISE_REQUIRES="$DENISE_REQUIRES cairo-pdf"
dnl    DENISE_CFLAGS="$CAIRO_PDF_CFLAGS $DENISE_CFLAGS"
dnl    DENISE_LIBS="$CAIRO_PDF_LIBS $DENISE_LIBS"
dnl fi

dnl if test x"$enable_cairo_ps" = xyes
dnl then
dnl    PKG_CHECK_MODULES([CAIRO_PS], [cairo-ps])
dnl    DENISE_REQUIRES="$DENISE_REQUIRES cairo-ps"
dnl    DENISE_CFLAGS="$CAIRO_PS_CFLAGS $DENISE_CFLAGS"
dnl    DENISE_LIBS="$CAIRO_PS_LIBS $DENISE_LIBS"
dnl fi

dnl if test x"$enable_cairo_svg" = xyes
dnl then
dnl    PKG_CHECK_MODULES([CAIRO_SVG], [cairo-svg])
dnl    DENISE_REQUIRES="$DENISE_REQUIRES cairo-svg"
dnl    DENISE_CFLAGS="$CAIRO_SVG_CFLAGS $DENISE_CFLAGS"
dnl    DENISE_LIBS="$CAIRO_SVG_LIBS $DENISE_LIBS"
dnl fi

dnl if test x"$enable_cairo_xlib" = xyes
dnl then
dnl    PKG_CHECK_MODULES([CAIRO_XLIB], [cairo-xlib])
dnl    DENISE_REQUIRES="$DENISE_REQUIRES cairo-xlib"
dnl    DENISE_CFLAGS="$CAIRO_XLIB_CFLAGS $DENISE_CFLAGS"
dnl    DENISE_LIBS="$CAIRO_XLIB_LIBS $DENISE_LIBS"
dnl fi

if test x"$enable_gtkmm" = xyes
then
   PKG_CHECK_MODULES([GTKMM_3_0], [gtkmm-3.0])
   DENISE_REQUIRES="$DENISE_REQUIRES gtkmm-3.0"
   DENISE_CFLAGS="$GTKMM_3_0_CFLAGS $DENISE_CFLAGS"
   DENISE_LIBS="$GTKMM_3_0_LIBS $DENISE_LIBS"
fi

if test x"$enable_fftw3" = xyes
then
   PKG_CHECK_MODULES([FFTW3], [fftw3])
   DENISE_REQUIRES="$DENISE_REQUIRES fftw3"
   DENISE_CFLAGS="$FFTW3_CFLAGS $DENISE_CFLAGS"
   DENISE_LIBS="$FFTW3_LIBS $DENISE_LIBS"
fi

dnl if test x"$enable_sdl" = xyes
dnl then
dnl    PKG_CHECK_MODULES([SDL], [sdl])
dnl    DENISE_REQUIRES="$DENISE_REQUIRES sdl"
dnl    DENISE_CFLAGS="$SDL_CFLAGS $DENISE_CFLAGS"
dnl    DENISE_LIBS="$SDL_LIBS $DENISE_LIBS"
dnl fi

dnl Checks for header files.
AC_CHECK_HEADERS(ieeefp.h)
AC_CHECK_HEADERS(cstdint)
AC_CHECK_HEADERS(jasper/jasper.h)

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
