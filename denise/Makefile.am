AM_CPPFLAGS = -I$(top_srcdir) -I$(top_builddir)
lib_LTLIBRARIES	= libdenise.la

pkginclude_HEADERS = analysis.h astronomy.h basics.h aviation.h \
   cluster.h dataset.h ds.h dstring.h dtime.h exception.h fractal.h \
   geodesy.h geometry.h gis.h gts.h gzstream.h histogram.h linalg.h \
   kalman.h marine.h met.h metcode.h stat.h streamline.h \
   symbolic.h tc.h thermo.h transform.h util.h visualize.h

libdenise_la_LDFLAGS = -release $(VERSION)
libdenise_la_SOURCES = analysis.cc astronomy.cc basics.cc aviation.cc \
   cluster.cc dataset.cc ds.cc dstring.cc dtime.cc exception.cc fractal.cc \
   geodesy.cc geometry.cc gis.cc gts.cc gzstream.cc histogram.cc linalg.cc \
   kalman.cc marine.cc met.cc metcode.cc stat.cc streamline.cc \
   symbolic.cc tc.cc thermo.cc transform.cc util.cc visualize.cc

if ENABLE_CAIROMM
pkginclude_HEADERS += cairoable.h graphics.h
libdenise_la_SOURCES += cairoable.cc graphics.cc
endif

if ENABLE_GTKMM
pkginclude_HEADERS += gtkmm.h
libdenise_la_SOURCES += gtkmm.cc
endif

if ENABLE_NWP
pkginclude_HEADERS += nwp.h
libdenise_la_SOURCES += nwp.cc
endif

if ENABLE_FFT
pkginclude_HEADERS += fft.h
libdenise_la_SOURCES += fft.cc
endif

