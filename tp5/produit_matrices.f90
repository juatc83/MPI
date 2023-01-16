!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! produit_matrices.f90 --- TP5 : produit de matrices
!
! Auteur          : Jalel Chergui (CNRS/IDRIS - France) <Jalel.Chergui@idris.fr>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Remarques :
! ---------
!
!   * On veut r�aliser le produit de matrices C = A * B en parall�le.
!
!   * On suppose que ces matrices sont carr�es et que leur ordre N
!     est divisible par le nombre Nprocs de processus.
!
!   * Le processus 0 initialise les matrices A et B qu'il distribue
!     ensuite aux autres processus.
!
!   * La distribution de A se fait par bandes horizontales.
!     La distribution de B se fait par bandes verticales.
!
!   * Chaque processus poss�de une bande des matrices A et B.
!
!   * Chaque processus calcule ainsi un bloc de la diagonale principale
!     de C avec les �l�ments qu'il poss�de. Le calcul des blocs
!     extra-diagonaux n�cessite des communications avec les autres
!     processus.
!
!   * En fait, l'op�ration se reduit ici � un produit de matrices par bloc.

program produit_matrices

  use mpi
  implicit none

  integer, parameter                  :: etiquette=1000, Nprocs_max=16
  integer                             :: rang, Nprocs, N, NL, code, k, type_temp, type_tranche
  integer                             :: rang_suivant, rang_precedent, taille_type_reel
  integer, dimension(MPI_STATUS_SIZE) :: statut
  real                                :: Emax
  real, allocatable, dimension(:,:)   :: A, B, C, CC, E
  real, allocatable, dimension(:,:)   :: AL, BL, CL, TEMP
  integer(kind=MPI_ADDRESS_KIND)      :: borne_inferieure=0, taille_deplacement_type_tranche

  ! Initialisation de MPI
  call MPI_INIT(code)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, Nprocs, code)

  if (rang == 0) then
    print *, 'Entrez l''ordre N global des matrices :'
    read *,N

    ! Il faut que N soit divisible par Nprocs
    if ( mod(N, Nprocs) == 0 ) then
      NL = N / Nprocs

      ! Le processus 0 diffuse N � tous les autres processus
  
    else
      print *, 'N n''est pas divisible par Nprocs'
      ! On arrete l'execution
      
    end if
  end if

  ! Le processus 0 initialise les matrices A et B
  if (rang == 0) then
    ! Allocation dynamique de m�moire, entre autres, des matrices A, B et C
    allocate( A(N,N), B(N,N), C(N,N), CC(N,N) )

    ! Initialisation de A et B
    call RANDOM_NUMBER(A)
    call RANDOM_NUMBER(B)

    ! Calcul monoprocesseur du produit matriciel A*B
    CC(:,:) = matmul(A(:,:), B(:,:))

  end if

  ! Allocation dynamique de m�moire des divers tableaux locaux
  allocate( AL(NL,N), BL(N,NL), CL(N,NL), TEMP(NL,N) )

  call MPI_TYPE_SIZE(MPI_REAL, taille_type_reel, code)

  ! Construction du type qui correspond a 1 bloc de NL lignes et N colonnes


  ! Le processus 0 distribue dans AL les tranches horizontales de la matrice A


  ! Le processus 0 distribue dans BL les tranches verticales de la matrice B


  ! Calcul des blocs diagonaux de la matrice resultante.
  CL(rang*NL+1:(rang+1)*NL,:) = matmul( AL(:,:), BL(:,:) )

  ! Premier algorithme (deux fois plus co�teux que le second)
  do k = 0, Nprocs-1
    ! Chaque processus ENVOIE sa tranche AL au processus k
    ! et RE�OIT dans TEMP la tranche AL du processus k
    if (rang /= k) then


      ! Chaque processus calcule les blocs situ�s au-dessus
      ! et en dessous du bloc de la diagonale principale
      CL(k*NL+1:(k+1)*NL,:)=matmul(TEMP(:,:),BL(:,:))
    end if
  end do

  ! Second algorithme
!  rang_precedent = mod(Nprocs+rang-1,Nprocs)
!  rang_suivant   = mod(rang+1,Nprocs)
!  do k = 1, Nprocs-1
!    ! Chaque processus ENVOIE sa tranche AL au processus pr�c�dent
!    ! et RE�OIT la tranche AL du processus suivant (mais les contenus changent)
! 
!
!    ! Chaque processus calcule les blocs situ�s au-dessus
!    ! et en dessous du bloc de la diagonale principale
!    CL(mod(rang+k,Nprocs)*NL+1:(mod(rang+k,Nprocs)+1)*NL,:)=matmul(AL(:,:),BL(:,:))
!  end do

  ! Le processus 0 collecte les tranches CL de tous les processus
  ! pour former la matrice r�sultante C
  

  ! Les tableaux locaux sont d�sormais inutiles
  deallocate( AL, BL, CL, TEMP )

  ! V�rification des r�sultats
  if (rang == 0) then
    allocate( E(N,N) )
    E(:,:) = abs(C(:,:) - CC(:,:))
    Emax   = maxval( E(:,:) ) / N**2
    deallocate( A, B, C, CC, E )

    if ( Emax <= epsilon(1.0) ) then
      print'(/,40X,"Bravo !",/,  &
            & 20X,"Le produit matriciel A*B calcule en parallele",/, &
            & 20X,"est bien egal a celui calcule en monoprocesseur")'
    else
      print'(/,33X,"Resultat incorrect !",/, &
            & 20X,"Le produit matriciel A*B calcule en parallele",/, &
            & 20X,"est different de celui calcule en monoprocesseur")'
    end if

  end if

  call MPI_FINALIZE(code)

end program produit_matrices
