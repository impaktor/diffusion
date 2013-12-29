EXECUTABLE = prog
CC = g++

CFLAGS = -c -pedantic -Wall -Wextra -std=c++0x
LDFLAGS =
#LDFLAGS = -lboost_thread
# -L/usr/local/lib

ifeq ($(mode),debug)
   $(info COMPILING DEBUG MODE)
   FLAGS = $(CFLAGS) -O0 -g3 -Wno-unknown-pragmas
else ifeq ($(mode),profile)
   $(info COMPILING PROFILING MODE (use gprof))
   FLAGS = $(CFLAGS) -O2 -D NDEBUG -Wno-unknown-pragmas
   LDFLAGS += -pg -fprofile-arcs -ftest-coverage
else ifeq ($(mode),serial)
   $(info COMPILING SERIAL MODE)
   FLAGS = $(CFLAGS) -O2 -D NDEBUG -Wno-unknown-pragmas
else ifeq ($(mode),)
   $(info COMPILING STANDARD RELEASE MODE (no debug, optimized, parallel))
   mode = release
   FLAGS = $(CFLAGS) -O3 -D NDEBUG -fopenmp
   LDFLAGS += -lgomp
else
   $(info ERROR: WRONG mode SELECTED)
   LDFLAGS = -BLABLABLA
   FLAGS = -YADDAYADDA
endif
#	@echo "Building on "$(mode)" mode"

#  #flags:
#  -O0        no optimization
#  -O2        optimize quite a bit
#  -g3        compile debug info
#  -std       use new c++ standard
#  -D NDEBUG  NoDebug, disables all asserts()
#  -pg        link with  for gprof profiling.

########################### END #######################################

MY_HEADERS = classes.h auxiliary.h lattice.h save.h superInteraction.h

SOURCES_BASIC  = lattice.cpp classes.cpp auxiliary.cpp save.cpp
SOURCES_SUPER  = superInteraction.cpp main_super.cpp
SOURCES_MAIN   = main.cpp
SOURCES_EXTRA  = read_data.cpp

OBJECTS_BASIC  = $(SOURCES_BASIC:.cpp=.o)
OBJECTS_SUPER  = $(SOURCES_SUPER:.cpp=.o)
OBJECTS_MAIN   = $(SOURCES_MAIN:.cpp=.o)
OBJECTS_EXTRA  = read_data.o


ALL_FILES = $(SOURCES_SUPER) $(SOURCES_MAIN) $(SOURCES_BASIC) $(MY_HEADERS)\
	nr/* simpleini/* Makedep.rule Makefile README.txt input.dat *.py

all: $(EXECUTABLE)

.cpp.o:
	$(CC) $(FLAGS) -o $@ $<

super: $(OBJECTS_BASIC) $(OBJECTS_SUPER)
	$(CC) $(LDFLAGS) -o $@ $+

depend:
	$(CC) $(FLAGS) -MM $(SOURCES_BASIC) $(SOURCES_SUPER) $(SOURCES_MAIN) $(SOURCES_EXTRA) >Makedep.rule

Makedep.rule: Makefile
	@echo "The makefile has changed since the last 'make depend'."

# Dependencies for executable
$(EXECUTABLE): $(OBJECTS_BASIC) $(OBJECTS_MAIN)
	$(CC) $(LDFLAGS) -o $@ $+

clean:
	touch -d20020101 Makedep.rule
	make depend
	\rm -v $(OBJECTS_MAIN) $(OBJECTS_BASIC) $(OBJECTS_SUPER) $(EXECUTABLE) $(OBJECTS_EXTRA) super read_data

#use gz:
tar:
	tar -vzcf prog.tar.gz $(ALL_FILES)

zip:
	zip prog.zip $(ALL_FILES)

#For Emacs-lovers
tags:
	etags *.h *.cpp

lines:
	wc -l *.h *.cpp

read_data: read_data.o save.o auxiliary.o classes.o
	$(CC) $(LDFLAGS) read_data.o save.o classes.o auxiliary.o  -o $@

# Dependencies for object files
include Makedep.rule
