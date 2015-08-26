INCLUDES	= -I$(top_builddir) -I$(top_srcdir) \
		  -I$(top_builddir)/denise -I$(top_srcdir)/denise

noinst_HEADERS	= andrea.h package.h \
		  geodesy.h geodetic_mesh.h geodetic_transform.h \
		  gshhs.h journey.h sounding.h surface.h tc.h
bin_PROGRAMS	= andrea

andrea_LDADD	= $(top_builddir)/denise/libdenise.la
andrea_SOURCES	= main.cc andrea.cc package.cc \
		  geodesy.cc geodetic_mesh.cc geodetic_transform.cc \
		  gshhs.cc journey.cc sounding.cc surface.cc tc.cc

