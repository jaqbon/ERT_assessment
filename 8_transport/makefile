#--------------------------------------------
CC = gcc
CFLAGS =

#--------------------------------------------
PROGRAM = transport
OBJECTS = transport.o flux_corrected_transport.o

#--------------------------------------------
all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $(OBJECTS)

$(OBJECTS): *.c
	$(CC) $(CFLAGS) -c $*.c

clean:
	rm -f $(PROGRAM) $(OBJECTS)
