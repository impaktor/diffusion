CC = g++

#The -c flag makes .o files, that will speed up future 
#compilations in multi-source projects
CFLAGS = -c -O3

LDFLAGS = -Wall -pedantic
DEBUG = -g 
SOURCES = prog.cpp classes.cpp auxilary.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = prog

ALL_FILES = prog.cpp classes.h classes.cpp \
auxilary.cpp auxilary.h Makefile ran.h nr3.h ludcmp.h README.txt
MY_HEADERS = classes.h auxilary.h 

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@


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
	g++ -O2 prog.cpp classes.cpp auxilary.cpp -o prog

tags:
	etags $(SOURCES)


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