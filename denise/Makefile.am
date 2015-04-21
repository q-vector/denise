INCLUDES = -I$(top_srcdir) -I$(top_builddir)
lib_LTLIBRARIES	= libdenise.la

#pkginclude_HEADERS = analysis.h basics.h aviation.h cairo.h cluster.h \
#   ds.h dstring.h dtime.h exception.h geodesy.h geometry.h kalman.h \
#   met.h stat.h streamline.h thermo.h util.h visualize.h
pkginclude_HEADERS = analysis.h astronomy.h basics.h aviation.h \
   cairoable.h cluster.h dataset.h ds.h dstring.h dtime.h exception.h \
   fractal.h geodesy.h geometry.h gis.h graphics.h grib.h gtkmm.h \
   gts.h histogram.h linalg.h kalman.h marine.h met.h metcode.h \
   stat.h streamline.h symbolic.h thermo.h transform.h util.h \
   visualize.h

libdenise_la_LDFLAGS = -release $(VERSION)
#libdenise_la_SOURCES = analysis.cc basics.cc aviation.cc cairo.cc cluster.cc \
#   ds.cc dstring.cc dtime.cc exception.cc geodesy.cc geometry.cc kalman.cc \
#   met.cc stat.cc streamline.cc thermo.cc util.cc visualize.cc
libdenise_la_SOURCES = analysis.cc astronomy.cc basics.cc aviation.cc \
   cairoable.cc cluster.cc dataset.cc ds.cc dstring.cc dtime.cc exception.cc \
   fractal.cc geodesy.cc geometry.cc gis.cc graphics.cc grib.cc gtkmm.cc \
   gts.cc histogram.cc linalg.cc kalman.cc marine.cc met.cc metcode.cc \
   stat.cc streamline.cc symbolic.cc thermo.cc transform.cc util.cc \
   visualize.cc

if ENABLE_NWP
pkginclude_HEADERS += nwp.h
libdenise_la_SOURCES += nwp.cc
endif

#if ENABLE_CAIRO_PDF
#pkginclude_HEADERS += cairo-pdf.h
#libdenise_la_SOURCES += cairo-pdf.cc
#endif

#if ENABLE_CAIRO_PS
#pkginclude_HEADERS += cairo-ps.h
#libdenise_la_SOURCES += cairo-ps.cc
#endif

#if ENABLE_CAIRO_SVG
#pkginclude_HEADERS += cairo-svg.h
#libdenise_la_SOURCES += cairo-svg.cc
#endif

#if ENABLE_CAIRO_XLIB
#pkginclude_HEADERS += cairo-xlib.h
#libdenise_la_SOURCES += cairo-xlib.cc
#endif

if ENABLE_FFT
pkginclude_HEADERS += fft.h
libdenise_la_SOURCES += fft.cc
endif

#if ENABLE_SDL
#pkginclude_HEADERS += sdl.h
#libdenise_la_SOURCES += sdl.cc
#endif

