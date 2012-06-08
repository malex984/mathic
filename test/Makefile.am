ACLOCAL_AMFLAGS = -I ../m4

AM_CPPFLAGS = -I${top_srcdir}/

TESTS=test
check_PROGRAMS=$(TESTS)

test_SOURCES=DivFinder.cpp
test_LDFLAGS=-lgtest -lgtest_main -lmemtailor ${top_builddir}/libmathic.la 
