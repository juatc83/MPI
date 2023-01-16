############################ -*- Mode: Makefile -*- ###########################
## Makefile --- Cours MPI : TP4 : Type de données derivées (type_transpose)
##
## Auteur          : Isabelle Dupays (CNRS/IDRIS) <Isabelle.Dupays@idris.fr>
## 
###############################################################################

# Compilateurs, options de compilation et d'édition de liens
include ../arch/make_inc

OBJS = transpose.o

# Règle implicite de compilation
.SUFFIXES: .o .f90
.f90.o:
	$(CF95) -c $(FFLAGS_TP4) $<

all : transpose

transpose: $(OBJS)
	$(CF95) -o $@ $(LDFLAGS_TP4) $(OBJS)
	$(MPIEXEC_TP4) ./transpose

clean :
	rm -f $(OBJS) transpose core
