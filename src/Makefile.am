INCLUDES	= -I$(top_builddir) -I$(top_srcdir) \
		  -I$(top_builddir)/denise -I$(top_srcdir)/denise

noinst_HEADERS	= andrea.h
bin_PROGRAMS	= andrea

andrea_LDADD	= $(top_builddir)/denise/libdenise.la
andrea_SOURCES	= main.cc andrea.cc

