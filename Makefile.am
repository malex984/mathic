ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -I${top_srcdir}/
# -I${top_builddir} 

lib_LTLIBRARIES = libmathic.la

#libmathic_la_CFLAGS   = 
libmathic_la_LDFLAGS = -lmemtailor

libmathic_la_SOURCES = \
mathic/Timer.cpp mathic/ColumnPrinter.cpp mathic/DivMask.cpp \
mathic/Action.cpp mathic/BoolParameter.cpp mathic/CliParameter.cpp \
mathic/CliParser.cpp mathic/ElementDeleter.cpp mathic/NameFactory.cpp \
mathic/error.cpp mathic/HelpAction.cpp mathic/IntegerParameter.cpp \
mathic/StringParameter.cpp mathic/display.cpp

nobase_include_HEADERS = mathic.h mathic/stdinc.h mathic/Timer.h \
mathic/ColumnPrinter.h mathic/ElementDeleter.h mathic/error.h \
mathic/DivList.h mathic/KDTree.h mathic/TourTree.h mathic/StlSet.h \
mathic/Heap.h mathic/Geobucket.h mathic/Action.h mathic/DivMask.h \
mathic/BoolParameter.h mathic/CliParameter.h mathic/CliParser.h \
mathic/HelpAction.h mathic/IntegerParameter.h mathic/StringParameter.h \
mathic/display.h mathic/NameFactory.h
 
EXTRA_DIST = autogen.sh

# SUBDIRS=test
# TESTDIRS=test
