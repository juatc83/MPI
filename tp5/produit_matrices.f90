!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! produit_matrices.f90 --- TP5 : produit de matrices
!
! Auteur          : Jalel Chergui (CNRS/IDRIS - France) <Jalel.Chergui@idris.fr>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Remarques :
! ---------
!
!   * On veut realiser le produit de matrices C = A * B en parallele.
!
!   * On suppose que ces matrices sont carrees et que leur ordre N
!     est divisible par le nombre Nprocs de processus.
!
!   * Le processus 0 initialise les matrices A et B qu'il distribue
!     ensuite aux autres processus.
!
!   * La distribution de A se fait par bandes horizontales.
!     La distribution de B se fait par bandes verticales.
!
!   * Chaque processus possede une bande des matrices A et B.
!
!   * Chaque processus calcule ainsi un bloc de la diagonale principale
!     de C avec les elements qu'il possede. Le calcul des blocs
!     extra-diagonaux necessite des communications avec les autres
!     processus.
!
!   * En fait, l'operation se reduit ici e un produit de matrices par bloc.

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
    print *, 'Entrez l''ordre N global des matrices : (non : ca sera 64)'
    N = 32
    
    ! Il faut que N soit divisible par Nprocs
    if ( mod(N, Nprocs) == 0 ) then
      NL = N / Nprocs
      
      ! Le processus 0 diffuse N et NL a tous les autres processus
      call MPI_BCAST(N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, code)
      call MPI_BCAST(NL, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, code)
    else
      print *, 'N n''est pas divisible par Nprocs'
      
      ! On arrete l'execution
      call MPI_ABORT(MPI_COMM_WORLD, 1, code)
    
    end if
  end if

  ! Le processus 0 initialise les matrices A et B
  if (rang == 0) then
    ! Allocation dynamique de memoire, entre autres, des matrices A, B et C
    allocate( A(N,N), B(N,N), C(N,N), CC(N,N) )
    
    ! Initialisation de A et B
    call RANDOM_NUMBER(A)
    call RANDOM_NUMBER(B)

    ! Calcul monoprocesseur du produit matriciel A*B
    CC(:,:) = matmul(A(:,:), B(:,:))

      ! Allocation dynamique de memoire des divers tableaux locaux
    allocate( AL(NL,N), BL(N,NL), CL(N,NL), TEMP(NL,N) )
    print*,'ok_alloc'
  end if



  call MPI_TYPE_SIZE(MPI_REAL, taille_type_reel, code)

  ! Construction du type qui correspond a 1 bloc de NL lignes et N colonnes
  call MPI_TYPE_VECTOR(N,NL,N,MPI_REAL,type_tranche,code)
  ! Validation du type type_tranche
  call MPI_TYPE_COMMIT(type_tranche,code)

  if(rang==0)then
    ! Le processus 0 distribue dans AL les tranches horizontales de la matrice A
    do k = 0, Nprocs-1
      call MPI_SEND(A(k*NL+1,:), 1, type_tranche, k, etiquette, MPI_COMM_WORLD, code)
    end do
    print*,'send_ok'
    ! Le processus 0 distribue dans BL les tranches verticales de la matrice B
    do k = 0, Nprocs-1
      call MPI_SEND(B(:,k*NL+1), 1, type_tranche, k, etiquette, MPI_COMM_WORLD, code)
    end do

  ! Les autres processus recoivent les informations
  else
    call MPI_RECV(AL, 1, type_tranche, 0, etiquette, MPI_COMM_WORLD, statut, code)
    call MPI_RECV(BL, 1, type_tranche, 0, etiquette, MPI_COMM_WORLD, statut, code)
  end if


  ! Calcul des blocs diagonaux de la matrice resultante.
  CL(rang*NL+1:(rang+1)*NL,:) = matmul( AL(:,:), BL(:,:) )

  ! Premier algorithme (deux fois plus coeteux que le second)
  do k = 0, Nprocs-1
    ! Chaque processus ENVOIE sa tranche AL au processus k
    ! et RECOIT dans TEMP la tranche AL du processus k
    if (rang /= k) then
      call MPI_SEND(AL, 1, type_tranche, k, etiquette, MPI_COMM_WORLD, code)
      call MPI_RECV(TEMP, 1, type_tranche, k, etiquette, MPI_COMM_WORLD, statut, code)

      ! Chaque processus calcule les blocs situes au-dessus
      ! et en dessous du bloc de la diagonale principale
      CL(k*NL+1:(k+1)*NL,:)=matmul(TEMP(:,:),BL(:,:))
    end if
  end do

  ! Second algorithme
!  rang_precedent = mod(Nprocs+rang-1,Nprocs)
!  rang_suivant   = mod(rang+1,Nprocs)
!  do k = 1, Nprocs-1
!    ! Chaque processus ENVOIE sa tranche AL au processus precedent
!    ! et REeOIT la tranche AL du processus suivant (mais les contenus changent)
! 
!
!    ! Chaque processus calcule les blocs situes au-dessus
!    ! et en dessous du bloc de la diagonale principale
!    CL(mod(rang+k,Nprocs)*NL+1:(mod(rang+k,Nprocs)+1)*NL,:)=matmul(AL(:,:),BL(:,:))
!  end do

  ! Le processus 0 collecte les tranches CL de tous les processus
  ! pour former la matrice resultante C
  

  ! Les tableaux locaux sont desormais inutiles
  deallocate( AL, BL, CL, TEMP )

  ! Verification des resultats
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
