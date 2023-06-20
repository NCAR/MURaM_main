CFLAGS = -Wall -O3
LDFLAGS = 
CC = icpc
LD = icpc
EXEC = eos.x
LIBS =
OBJS = atmdatb.o \
	   molecb.o \
	   partf.o \
	   saha.o \
	   io.o \
	   Wittmann.o

######### Suffix rules ########################################
.SUFFIXES :    .o .cc

.cc.o:
	$(CC) $(CFLAGS) -c $<

all: 
	$(MAKE) $(OBJS)
	$(LD) -o $(EXEC) $(OBJS) $(LIBS) $(LDFLAGS)

clean:
	rm -f $(EXEC) $(OBJS)
