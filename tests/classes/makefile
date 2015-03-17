SDIR		= src
HDIR		= include
ODIR		= objects
CC 			= g++
OPTIM 		= 
CFLAGS 		= -Wall -g
LFLAGS 		= 
INCLUDES	= -I$(HDIR)
LIBS 		= 
_SRC		= simple.cc error.cc bag.cc test.cc
_HEADERS 	= simple.h error.h bag.h
_OBJS 		= $(_SRC:.cc=.o)
SRC			= $(patsubst %,$(SDIR)/%,$(_SRC))
HEADERS		= $(patsubst %,$(HDIR)/%,$(_HEADERS))
OBJS		= $(patsubst %,$(ODIR)/%,$(_OBJS))
MAIN 		= main

all:	$(MAIN)
	@echo Simple compiler named $(MAIN) has been compiled

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) $^ -o $@ $(LIBS)
    
$(ODIR)/%.o: $(SDIR)/%.cc
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	$(RM) *.o~ $(MAIN)
