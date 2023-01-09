############################ -*- Mode: Makefile -*- ###########################
## Makefile --- Cours MPI : TP3 : Communications collectives
##
## Auteur          : Dimitri Lecas (CNRS/IDRIS) <Dimitri.Lecas@idris.fr>
##
###############################################################################

# Compilateurs, options de compilation et d'�dition de liens
include ../arch/make_inc

OBJS  = pi.o

# R�gle implicite de compilation
.SUFFIXES: .o .f90
.f90.o:
	$(CF95) -c $(FFLAGS_TP3) $<

all: pi

pi: $(OBJS)
	$(CF95) -o $@ $(LDFLAGS_TP3) $(OBJS)
	$(MPIEXEC_TP3) ./pi

clean:
	rm -f $(OBJS) pi core
