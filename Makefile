EXECUTABLE = prog
CC   = g++

override CFLAGS := -c -pedantic -Wall -Wno-switch $(CFLAGS)

LIBS =
LDFLAGS =
#LDFLAGS = $(LIBS) -pg -fprofile-arcs -ftest-coverage
# -L/usr/local/lib

###SLOW: -O0 = no optimization, -g3 compile debug info, std = use new c++ standard
#======================================================================
FLAGS = $(CFLAGS) -O0 -g3 -std=c++0x

###FAST: -O2 = optimize, use new standard, NDEBUG = NoDebug, disables asserts()
#======================================================================
#FLAGS = $(CFLAGS) -O2 -std=c++0x -D NDEBUG

########################### END #######################################

MY_HEADERS = classes.h auxiliary.h lattice.h save.h superInteraction.h global.h

SOURCES_BASIC  = lattice.cpp classes.cpp auxiliary.cpp save.cpp
SOURCES_SUPER  = superInteraction.cpp main_super.cpp
SOURCES_THREAD = main_thread.cpp
SOURCES_MAIN   = main.cpp

OBJECTS_BASIC  = $(SOURCES_BASIC:.cpp=.o)
OBJECTS_SUPER  = $(SOURCES_SUPER:.cpp=.o)
OBJECTS_THREAD = $(SOURCES_THREAD:.cpp=.o)
OBJECTS_MAIN   = $(SOURCES_MAIN:.cpp=.o)

ALL_FILES = $(SOURCES_SUPER) $(SOURCES_THREAD) $(SOURCES_MAIN) $(SOURCES_BASIC) $(MY_HEADERS) nr/* Makefile README.txt input.dat

all: $(EXECUTABLE)

.cpp.o:
	$(CC) $(FLAGS) -o $@ $<

thread: $(OBJECTS_BASIC) $(OBJECTS_THREAD)
	$(CC) -lboost_thread -o $@ $+

super: $(OBJECTS_BASIC) $(OBJECTS_SUPER)
	$(CC) -o $@ $+

depend:
	$(CC) $(FLAGS) -MM $(SOURCES_BASIC) $(SOURCES_SUPER) $(SOURCES_THREAD) $(SOURCES_MAIN) >Makedep.rule

Makedep.rule: Makefile
	@echo "The makefile has changed since the last 'make depend'."

# Dependencies for executable
$(EXECUTABLE): $(OBJECTS_BASIC) $(OBJECTS_MAIN)
	$(CC) $(LDFLAGS) -o $@ $+

# Dependencies for object files
include Makedep.rule

clean:
	\rm -v *.o *~ a.out $(EXECUTABLE) super thread
	touch -d20020101 Makedep.rule
	make depend
	make tags

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
