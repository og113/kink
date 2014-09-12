# makefile for false vacuum program
CC = g++
CFLAGS = -static -ansi -g -Wall -pedantic -std=c++0x -ggdb #-ggdb is for debugging; static is because gsl-shared wasn't working
LFLAGS = 
LIBS = -lm -lgsl -lgslcblas -DHAVE_INLINE

HEADERS =
INCLUDES = -I/home/og113/Documents/kink/ -I/home/og113/Documents/c++/eigen_build/eigen/ -I/home/usr/local/lib/ -I/home/usr/local/include/ #includes these locations
OBJS = $(x:.cc=.o)
MAIN = $(x:.cc=)

.PHONY: depend clean

all:	$(MAIN)
	@echo Simple compiler named $(MAIN) has been compiled

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

.cc.o:  
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	$(RM) *.o~ $(MAIN)

depend: $(x) $(HEADERS)
	makedepend $(INCLUDES) $^
# DO NOT DELETE
