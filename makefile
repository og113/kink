# makefile for false vacuum program
CC = g++
CFLAGS = -std=c++0x -static -g -Wall
#CFLAGS EXPLAINED:
#-std=c++0x : added so that to_string and similar functions can be used
#-g : for using gdb to debug
#-ggdb : for using gdb to debug - not sure how this differs from above
#-Wall : turns on most compiler warnings
#-O : does some optimizations for running (compiling takes longer though)
#-Os : optimises as long as code size isn't increased
#-O2 : does some more optimizations that -O
#-O3 : does all the optimizations for running
#-static : linking to libraries statically
#-ansi -pedantic : turns off some extensions of g++ which are incompatible with the ANSI language

LFLAGS = 
LIBS = -lm -lgsl -lgslcblas -DHAVE_INLINE

HEADERS =
INCLUDES = -I/home/og113/Documents/kink/ -I/home/og113/Documents/c++/eigen_build/eigen -I/home/usr/local/lib/ -I/home/usr/local/include/ -I/home/og113/Documents/c++/gnuplot/gnuplot-cpp #includes these locations
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
