CFLAGS = -Wall -pg -g3
LDFLAGS = 
CC = g++
LD = g++
EXEC = tausort.x
LIBS =
OBJS = tausort.o 

######### Suffix rules ########################################
.SUFFIXES :    .o .cc

.cc.o:
	$(CC) $(CFLAGS) -c $<

all: 
	$(MAKE) $(OBJS)
	$(LD) -o $(EXEC) $(OBJS) $(LIBS) $(LDFLAGS)

clean:
	rm -f $(EXEC) $(OBJS)
