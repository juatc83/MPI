!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! transpose.f90  --- Utilisation d'un type derive (type_transpose)
!!                    pour transposer une matrice.
!!
!!
!! Auteur          : Isabelle DUPAYS (CNRS/IDRIS - France)
!!                   <Isabelle.Dupays@idris.fr>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PROGRAM transpose
  USE MPI
  IMPLICIT NONE
  INTEGER, PARAMETER                     :: nb_lignes=5,nb_colonnes=4,&
                                            etiquette=1000
  INTEGER                                :: code,rang,type_transpose,&
                                            taille_reel,i,j
  REAL, DIMENSION(nb_lignes,nb_colonnes) :: A
  REAL, DIMENSION(nb_colonnes,nb_lignes) :: AT
  INTEGER(kind=MPI_ADDRESS_KIND)         :: pas
  INTEGER, DIMENSION(MPI_STATUS_SIZE)    :: statut


  !Initialisation de MPI
  CALL MPI_INIT(code)

  !-- Savoir qui je suis
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)

  !-- Initialisation de la matrice AT
  AT(:,:) = 0.

  ! Définition du type type_transpose
  call MPI_TYPE_VECTOR(nb_colonnes,1,nb_lignes, MPI_REAL ,type_transpose,code)
  ! Validation du type type_transpose
  call MPI_TYPE_COMMIT(type_transpose,code)

  !Construction du type derive type_transpose pour transposer la
  !matrice A composee de nb_lignes et de nb_colonnes
  call MPI_TYPE_VECTOR(nb_colonnes,1,nb_lignes, MPI_REAL ,type_transpose,code)
  !Validation du type cree type_transpose
  call MPI_TYPE_COMMIT(type_transpose,code)

  IF (rang == 0) THEN
    !Initialisation de la matrice A sur le processus 0
    A(:,:) = RESHAPE( (/ (i,i=1,nb_lignes*nb_colonnes) /), &
                      (/ nb_lignes,nb_colonnes /) )

    PRINT *,'Matrice A'
    DO i=1,nb_lignes
      PRINT *,A(i,:)
    END DO

    do i=1,nb_lignes
      !Envoi de la matrice A au processus 1 avec le type type_transpose
      call MPI_SEND (a(i,1),1,type_transpose,1,etiquette, MPI_COMM_WORLD,code)
    enddo
  ELSE
    do i=1,nb_lignes
      !Reception pour le processus 1 dans la matrice AT
      call MPI_RECV (at(1,i),nb_lignes, MPI_REAL ,0,etiquette,&
    MPI_COMM_WORLD,statut,code)
    enddo

    PRINT *,'Matrice transposee AT'
    DO i=1,nb_colonnes
      PRINT *,AT(i,:)
    END DO

  END IF

  call mpi_type_free(type_transpose,code)

  !Sortie de MPI
  CALL MPI_FINALIZE(code)

END PROGRAM transpose
