CC=g++
#The -c flag makes .o files, that will speed up future 
#compilations in multi-source projects
CFLAGS= -c -O3

#For debugging
#CFLAGS=-g -Wall
LDFLAGS=
SOURCES=prog.cpp classes.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=prog

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
