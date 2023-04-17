all: call_iri

# local dir
MKFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
MKFILE_DIR  := $(dir $(MKFILE_PATH))

#--------------------------------------------
CC = gcc
FC = gfortran
CFLAGS =
FFLAGS = -std=legacy -fpic
LBFLAGS = -shared

LIBLOC = global

ifeq ($(LIBLOC), local)
  LIBLOCAL = -Wl,-rpath,$(MKFILE_DIR)
else
  LIBLOCAL =
endif

#--------------------------------------------
PROGRAM = call_iri
IRILIB  = iri_lib.so
COBJECTS = call_iri.o
FOBJECTS = cira.o igrf.o iridreg.o iriflip.o irifun.o irisub.o iritec.o

#--------------------------------------------
all: $(PROGRAM)

$(PROGRAM): $(COBJECTS) $(IRILIB)
	$(FC) -o $@ $(COBJECTS) -L. -l:$(IRILIB) $(LIBLOCAL)

$(IRILIB): $(FOBJECTS)
	$(FC) $(LBFLAGS) -o $@ $+

$(COBJECTS): *.c
	$(CC) $(CFLAGS) -c $*.c

$(FOBJECTS): *.for
	$(FC) $(FFLAGS) -c $*.for

clean:
	rm -f $(PROGRAM) $(IRILIB) $(COBJECTS) $(FOBJECTS)
