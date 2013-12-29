CC = g++

#The -c flag makes .o files, that will speed up future 
#compilations in multi-source projects
CFLAGS = -c -O3

LDFLAGS = -Wall -pedantic
DEBUG = -g -Wall                 #+ remove -O3 manualy!
SOURCES = prog.cpp classes.cpp auxilary.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = prog
ALL_FILES = prog.cpp classes.h classes.cpp \
auxilary.cpp auxilary.h Makefile ran.h nr3.h README.txt


all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@


clean:     
	\rm -v *.o *~

tar:
	tar -vzcf prog.tar.gz $(ALL_FILES)

zip:
	prog.zip $(ALL_FILES)