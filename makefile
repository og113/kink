# makefile for false vacuum program
CC = g++
CFLAGS = -static -std=c++0x -g -Wall# -ansi -pedantic
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
