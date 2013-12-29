CC = g++
CFLAGS = -c $(OPTFLAGS)
LFLAGS = -Wall -pedantic $(OPTFLAGS)
DEBUG = -g
OPTFLAGS = -O2
EXECUTABLE = prog

BASIC_SOURCES = lattice.cpp classes.cpp auxiliary.cpp save.cpp
SUPER = superInteraction.cpp superInteraction.h superProg.cpp
SOURCES = $(BASIC_SOURCES) main.cpp
OBJECTS = $(SOURCES:.cpp=.o)
SOURCES2 = $(BASIC_SOURCES) main_thread.cpp
OBJECTS2 = $(SOURCES2:.cpp=.o)
MY_HEADERS = classes.h auxiliary.h lattice.h save.h
ALL_FILES = $(SOURCES) $(MY_HEADERS) $(SUPER) nr/* Makefile README.txt TESTinteraction.cpp main_*


all: $(EXECUTABLE)

thread: $(OBJECTS2)
	$(CC) $(LFLAGS) -lboost_thread $(OBJECTS2) -o prog_thread

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o $@

save.o: save.h save.cpp auxiliary.o
	$(CC) $(CFLAGS) save.cpp

lattice.o: lattice.cpp lattice.h classes.o auxiliary.o
	$(CC) $(CFLAGS) lattice.cpp

classes.o: classes.cpp classes.h
	$(CC) $(CFLAGS) classes.cpp

auxiliary.o: auxiliary.cpp auxiliary.h
	$(CC) $(CFLAGS) auxiliary.cpp

super: $(BASIC_SOURCES:.cpp=.o) $(SUPER:.cpp=.o)
	$(CC) $(LFLAGS) $(BASIC_SOURCES:.cpp=.o) $(SUPER:.cpp=.o) -o $@

superProg.o: superProg.cpp save.o auxiliary.o classes.o superInteraction.o
	$(CC) $(CFLAGS) superProg.cpp

superInteraction.o: superInteraction.cpp superInteraction.h lattice.o classes.o auxiliary.o
	$(CC) $(CFLAGS) superInteraction.cpp

tags:
	etags $(SOURCES) $(SUPER)

clean:
	\rm -v *.o *~

tar:
	tar -vzcf prog.tar.gz $(ALL_FILES)

zip:
	zip prog.zip $(ALL_FILES)

lines:
	wc -lc $(SOURCES) $(MY_HEADERS) $(SUPER)
