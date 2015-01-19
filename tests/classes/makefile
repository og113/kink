CC 			= g++
OPTIM 		= 
CFLAGS 		= -Wall -g
LFLAGS 		= 
INCLUDES	= 
LIBS 		= 
SRC			= simple.cc error.cc bag.cc test.cc
HEADERS 	= simple.h error.h bag.h
OBJS 		= $(SRC:.cc=.o)
MAIN 		= main

all:	$(MAIN)
	@echo Simple compiler named $(MAIN) has been compiled

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) $^ -o $@ $(LIBS)
    
%.o: %.cc
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	$(RM) *.o~ $(MAIN)
