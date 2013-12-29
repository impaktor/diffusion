CC = g++

#The -c flag makes .o files, that will speed up future
#compilations in multi-source projects
CFLAGS = -c -O3

LDFLAGS = -Wall -pedantic
DEBUG = -g
SOURCES = prog.cpp lattice.cpp classes.cpp auxiliary.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = prog

ALL_FILES = prog.cpp lattice.cpp lattice.h classes.h classes.cpp \
auxiliary.cpp auxiliary.h Makefile nr/* README.txt TESTinteraction.cpp $(SUPER)
MY_HEADERS = classes.h auxiliary.h lattice.h

#These files are only included for easy tar.gz/zip-ing:
SUPER = superInteraction.cpp superInteraction.h superProg.cpp make.super

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

all: $(SOURCES) $(EXECUTABLE)

clean:
	\rm -v *.o *~

debug:
	$(CC) $(DEBUG) $(SOURCES) -o $(EXECUTABLE)

tar:
	tar -vzcf prog.tar.gz $(ALL_FILES)

zip:
	prog.zip $(ALL_FILES)

lines:
	wc -l $(SOURCES) $(MY_HEADERS)

simple:
	g++ -O2 prog.cpp classes.cpp auxiliary.cpp -o prog

tags:
	etags $(SOURCES)

super:
	make -f make.super


# ETAGS = etags
# rm TAGS
# find . -name '*.cpp' -o -name '*.h' -o -name '*.c' -print0 \
# | xargs $(ETAGS) --extra=+q --fields=+fksaiS --c++-kinds=+px --append





# Automatic variables are the variable which are set by Make, once a rule is matched. Some of the important automatic variables are:

#     * $@     this is the name of the target
#     * $^     the name of all the prequisites (with duplicate file names removed)
#     * $?     names of all the prequisites which are newer than the target
#     * $<     name of the first prequisite
#     * $+     same as $^, but duplicates are not removed
#     * $*     the stem of the target name

#         Refer to the Makefile in example ex01. The make file has been modified to use automatic variables:

# app: main.o operation.o
#         gcc $^ -o myapp.out

# main.o: main.c operation.h
#         gcc -c $< -o $@

# operation.o: operation.c
#         gcc -c $< -o $@

# clean:
#         rm -f *.o *.out

#         Listing.2 Makefile of example ex01 