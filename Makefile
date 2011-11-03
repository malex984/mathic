# ***** Variables

rawSources := Timer.cpp ColumnPrinter.cpp DivMask.cpp libs/memtailor.cpp
rawDivSources := divsim/Simulation.cpp divsim/divMain.cpp
rawPqSources := pqsim/Item.cpp pqsim/Model.cpp pqsim/Simulator.cpp	\
  pqsim/pqMain.cpp

rawTestSources = libs/gtest.cpp test/DivFinder.cpp

GTEST_DIR = libs/gtest/
GTEST_VERSION = 1.6.0

ifndef ldflags
  ldflags = $(LDFLAGS)
endif

ifndef CXX
  CXX      = "g++"
endif

ifndef BIN_INSTALL_DIR
  BIN_INSTALL_DIR = "/usr/local/bin/"
endif

cflags = $(CFLAGS) $(CPPFLAGS) -Wall -Isrc/ -Wno-uninitialized	\
  -Wno-unused-parameter -Ilibs/memtailor/include -isystem $(GTEST_DIR)	\
  -isystem $(GTEST_DIR)include
pqProgram = pq
divProgram = div
testProgram = matest

ifndef MODE
 MODE=release
endif

MATCH=false
ifeq ($(MODE), release)
  outdir = bin/release/
  cflags += -O2
  MATCH=true
endif
ifeq ($(MODE), debug)
  rawSources := $(rawSources)
  outdir = bin/debug/
  cflags += -g -D DEBUG -fno-inline -Werror -Wextra -Wno-uninitialized \
            -Wno-unused-parameter
  MATCH=true
endif
ifeq ($(MODE), shared)
  outdir = bin/shared/
  cflags += -O2 -fPIC
  MATCH=true
endif
ifeq ($(MODE), profile)
  outdir = bin/profile/
  cflags += -g -pg -O2 -D PROFILE
  ldflags += -pg
  MATCH=true
  benchArgs = _profile $(FROBBYARGS)
endif
ifeq ($(MODE), analysis)
  rawSources := $(rawSources)
  outdir = bin/analysis/
  cflags += -Wextra -fsyntax-only -O1 -Wfloat-equal -Wundef				\
  -Wno-endif-labels -Wshadow -Wlarger-than-1000 -Wpointer-arith			\
  -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wsign-compare	\
  -Waggregate-return -Wmissing-noreturn -Wmissing-format-attribute		\
  -Wno-multichar -Wno-deprecated-declarations -Wpacked					\
  -Wno-redundant-decls -Wunreachable-code -Winline						\
  -Wno-invalid-offsetof -Winvalid-pch -Wlong-long						\
  -Wdisabled-optimization -D DEBUG -Werror
  MATCH=true
endif

ifeq ($(MATCH), false)
  $(error Unknown value of MODE: "$(MODE)")
endif

sources    = $(patsubst %.cpp, src/%.cpp, $(rawSources))
pqSources  = $(patsubst %.cpp, src/%.cpp, $(rawSources) $(rawPqSources))
divSources = $(patsubst %.cpp, src/%.cpp, $(rawSources) $(rawDivSources))

objs = $(patsubst %.cpp, $(outdir)%.o, $(rawSources))
allObjs = $(patsubst %.cpp, $(outdir)%.o, \
  $(rawSources) $(rawPqSources) $(rawDivSources) $(rawTestSources))
pqObjs = $(patsubst %.cpp, $(outdir)%.o, $(rawSources) $(rawPqSources))
divObjs = $(patsubst %.cpp, $(outdir)%.o, $(rawSources) $(rawDivSources))
testObjs = $(patsubst %.cpp, $(outdir)%.o, $(rawSources) $(rawTestSources))

# ***** Compilation

.PHONY: all depend clean bin/$(divProgram) bin/$(pqProgram) test distribution clear fixspace bin/$(testProgram)

all: bin/$(pqProgram) $(outdir)$(pqProgram) bin/$(divProgram) $(outdir)$(divProgram) bin/$(testProgram) test bin/$(testProgram)

# Make symbolic link to program from bin/
bin/$(divProgram): $(outdir)$(divProgram)
ifneq ($(MODE), analysis)
	mkdir -p $(dir $@); rm -f $@; ln -s ../$< $@
endif

# Link object files into executable
$(outdir)$(divProgram): $(divObjs) | $(outdir)
	@mkdir -p $(dir $@)
ifeq ($(MODE), analysis)
	echo > $@
endif
ifneq ($(MODE), analysis)
	$(CXX) $(divObjs) $(ldflags) -o $@
	if [ -f $@.exe ]; then \
      mv -f $@.exe $@; \
	fi
endif
ifeq ($(MODE), release)
	strip $@
endif

# Make symbolic link to program from bin/
bin/$(pqProgram): $(outdir)$(pqProgram)
ifneq ($(MODE), analysis)
	mkdir -p $(dir $@); rm -f $@; ln -s ../$< $@
endif

# Link object files into executable
$(outdir)$(pqProgram): $(pqObjs) | $(outdir)
	@mkdir -p $(dir $@)
ifeq ($(MODE), analysis)
	echo > $@
endif
ifneq ($(MODE), analysis)
	$(CXX) $(pqObjs) $(ldflags) -o $@
	if [ -f $@.exe ]; then \
      mv -f $@.exe $@; \
	fi
endif
ifeq ($(MODE), release)
	strip $@
endif

test: $(outdir)$(testProgram)
	@$(outdir)$(testProgram)


# Make symbolic link to program from bin/
bin/$(testProgram): $(outdir)$(testProgram)
ifneq ($(MODE), analysis)
	mkdir -p $(dir $@); rm -f $@; ln -s ../$< $@
endif

# Link object files into executable
$(outdir)$(testProgram): $(testObjs) | $(outdir)
	@mkdir -p $(dir $@)
ifeq ($(MODE), analysis)
	echo > $@
endif
ifneq ($(MODE), analysis)
	$(CXX) $(testObjs) $(ldflags) -o $@
	if [ -f $@.exe ]; then \
      mv -f $@.exe $@; \
	fi
endif
ifeq ($(MODE), release)
	strip $@
endif


# Compile and output object files.
# In analysis mode no file is created, so create one
# to allow dependency analysis to work.
$(outdir)%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) ${cflags} -c $< -o $@
	$(CXX) $(cflags) -MM -c $< > $(@:.o=.d).tmp
# using /usr/bin/env echo to get the non-built-in echo on OS X, since
# the built-in one does not understand the parameter -n.
	@/usr/bin/env echo -n "$(dir $@)" > $(@:.o=.d)
	@cat $(@:.o=.d).tmp >> $(@:.o=.d)
	@sed -e 's/.*://' -e 's/\\$$//' < $(@:.o=.d).tmp | fmt -1 | \
	  sed -e 's/^ *//' -e 's/$$/:/' >> $(@:.o=.d)
	@rm -f $(@:.o=.d).tmp
ifeq ($(MODE), analysis)
	  echo > $@
endif

-include $(allObjs:.o=.d)

gtest: libs/gtest
libs/gtest:
	rm -rf bin/tmp/ libs/gtest-$(GTEST_VERSION) libs/gtest
	@mkdir -p bin/tmp/
	@mkdir -p libs/
	(cd bin/tmp; wget http://googletest.googlecode.com/files/gtest-$(GTEST_VERSION).zip);
	cd bin/tmp; unzip gtest-$(GTEST_VERSION).zip;
	rm -rf bin/tmp/gtest-$(GTEST_VERSION).zip;
	mv bin/tmp/gtest-$(GTEST_VERSION) libs/
	cd libs; ln -s gtest-$(GTEST_VERSION) gtest

clean: tidy
	rm -rf bin

# ***** Miscellaneous

tidy:
	find .|grep -x ".*~\|.*/\#.*\#|.*\.stackdump\|gmon\.out\|.*\.orig\|.*/core\|core"|xargs rm -f

# Fixes various white space related issues.
fixspace:
	find src/ -type f|xargs ./fixspace;

commit: test
	echo
	hg commit -m "$(MSG)"

# ***** Distribution

remoteUrl = ssh://daimi/projs/momodel
pull:
	hg pull $(remoteUrl)
push:
	hg push $(remoteUrl)

distribution:
ifndef VER
	echo "Please specify version of Frobby distribution using VER=x.y.z";
	exit 1;
endif
	rm -fr frobby_v$(VER).tar.gz frobby_v$(VER)
	mkdir frobby_v$(VER)
	cp -r frobgrob COPYING Makefile src test doc frobby_v$(VER)
	mkdir frobby_v$(VER)/4ti2
	tar --create --gzip --file=frobby_v$(VER).tar.gz frobby_v$(VER)/
	rm -fr frobby_v$(VER)	
	ls -l frobby_v$(VER).tar.gz

spkg: tidy depend
ifndef VER
	echo "Please specify version of Frobby spkg using VER=x.y.z";
	exit 1;
endif
	if [ "$$SAGE_LOCAL" = "" ]; then \
	  echo "SAGE_LOCAL undefined ... exiting"; \
	  echo "Maybe run 'sage -sh?'" \
	  exit 1; \
	fi

	if [ ! -d sage/ ]; then echo "sage/ directory not found."; exit 1; fi
# Ensure that previous builds have been cleaned up
	rm -rf bin/sagetmp bin/frobby-$(VER) bin/frobby-$(VER).spkg

	hg clone sage bin/sagetmp

	mkdir bin/sagetmp/src
	cp -r COPYING Makefile src test bin/sagetmp/src

	mv bin/sagetmp bin/frobby-$(VER)
	cd bin/; $(SAGE_ROOT)/sage -pkg `pwd`/frobby-$(VER)
	rm -rf bin/frobby-$(VER)
