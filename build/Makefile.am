AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libtofreco.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(OFFLINE_MAIN)/include/eigen3 \
  -I$(ROOTSYS)/include

pkginclude_HEADERS = \
  ../src/TOFReco.h

if ! MAKEROOT6
  ROOT5_DICTS = \
    TOFReco_Dict.cc
endif

libtofreco_la_SOURCES = \
  $(ROOT5_DICTS) \
  ../src/TOFReco.cc 

libtofreco_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lfun4all \
  -lphool \
  -lphgeom \
  -lphfield \
  -lPHGenFit \
  -lg4trackfastsim

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
