ACLOCAL_AMFLAGS = -I ../m4

AM_CPPFLAGS = -I${top_srcdir}/
AM_LDFLAGS=${top_builddir}/libmathic.la 

# TESTS=divsim
# check_PROGRAMS=$(TESTS)

EXTRA_PROGRAMS=divsim

divsim_SOURCES=divMain.cpp Simulation.cpp \
DivListModel.h divMain.h KDTreeModel.h \
Monomial.h Simulation.h stdinc.h
