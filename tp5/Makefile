############################ -*- Mode: Makefile -*- ###########################
## Makefile --- Cours MPI : TP5 : Communications collectives
##
## Auteur          : Jalel Chergui (CNRS/IDRIS)
## 
###############################################################################

# Compilateurs, options de compilation et d'édition de liens
include ../arch/make_inc

OBJS = produit_matrices.o

# Règle implicite de compilation
.SUFFIXES: .o .f90
.f90.o:
	$(CF95) -c $(FFLAGS_TP5) $<

all : produit_matrices

produit_matrices: $(OBJS)
	$(CF95) -o $@ $(LDFLAGS_TP5) $(OBJS)
	$(MPIEXEC_TP5) ./produit_matrices

clean:
	rm -f $(OBJS) produit_matrices core
