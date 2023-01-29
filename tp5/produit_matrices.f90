!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! produit_matrices.f90 --- TP5 : produit de matrices
!
! Auteur          : Jalel Chergui (CNRS/IDRIS - France) <Jalel.Chergui@idris.fr>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Remarques :
! ---------
!
!   * On veut réaliser le produit de matrices C = A * B en parallèle.
!
!   * On suppose que ces matrices sont carrées et que leur ordre N
!     est divisible par le nombre Nprocs de processus.
!
!   * Le processus 0 initialise les matrices A et B qu'il distribue
!     ensuite aux autres processus.
!
!   * La distribution de A se fait par bandes horizontales.
!     La distribution de B se fait par bandes verticales.
!
!   * Chaque processus possède une bande des matrices A et B.
!
!   * Chaque processus calcule ainsi un bloc de la diagonale principale
!     de C avec les éléments qu'il possède. Le calcul des blocs
!     extra-diagonaux nécessite des communications avec les autres
!     processus.
!
!   * En fait, l'opération se reduit ici à un produit de matrices par bloc.

! Début du programme
program produit_matrices
  ! On utilise MPI
  use mpi

  ! On ne veut pas de déclaration implicite
  implicit none

  ! Déclaration des variables
  integer, parameter                  :: etiquette=1000, Nprocs_max=16, nb_tests=10,dp = kind(1.d0)
  integer                             :: rang, Nprocs, N, NL, code, k, type_temp, type_tranche
  integer                             :: rang_suivant, rang_precedent, taille_real, i
  integer, dimension(nb_tests)        :: tableau_N
  integer, dimension(MPI_STATUS_SIZE) :: statut
  real                                :: Emax
  real(kind=dp)                       :: temps_debut,temps_fin
  real, allocatable, dimension(:,:)   :: A, B, C, CC, E
  real, allocatable, dimension(:,:)   :: AL, BL, CL, TEMP
  integer(kind=MPI_ADDRESS_KIND)      :: borne_inferieure=0, taille_tranche
  character(len=256)                  :: arg

  ! Initialisation de MPI
  call MPI_INIT(code)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, Nprocs, code)

  ! Lecture de l'ordre des matrices
  if (rang.eq.0) then
    ! print *, 'Entrez l''ordre N global des matrices :'
    ! read(*,'(i4)')N
    print'(/,20X,"On essaye plusieurs valeurs de N à la suite : ",/, &
    & 20X,"N = 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096",/, &
    & 20X,"avec un total de huit processeurs.")'
    print*,' '
    tableau_N = (/  8, 16, 32, 64, 128, 256, &
                     512, 1024, 2048, 4096 /)
  end if

  do i=1,nb_tests
    if (rang.eq.0) then
      N = tableau_N(i)
    end if
    ! Le processus 0 diffuse N à tous les autres processus
    call MPI_BCAST(N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, code)

    ! Il faut que N soit divisible par Nprocs
    if ( mod(N, Nprocs) == 0) then
      ! Si c'est bon : NL = N / Nprocs
      NL = N / Nprocs
    else
      ! Sinon : le processus 0 affiche un message d'erreur et on arrete l'execution
      print *, 'N=',N, ' n''est pas divisible par Nprocs=', Nprocs
      ! On arrete l'execution
      call MPI_ABORT(MPI_COMM_WORLD, 1, code)
    end if

    ! Le processus 0 initialise les matrices A et B
    if (rang == 0) then
      ! Allocation dynamique de mémoire, entre autres, des matrices A, B et C
      allocate( A(N,N), B(N,N), C(N,N), CC(N,N) )

      ! Initialisation de A et B
      call RANDOM_NUMBER(A)
      call RANDOM_NUMBER(B)

      ! Calcul monoprocesseur du produit matriciel A*B
      CC(:,:) = matmul(A(:,:), B(:,:))

    end if

    ! Allocation dynamique de mémoire des divers tableaux locaux
    allocate( AL(NL,N), BL(N,NL), CL(N,NL), TEMP(NL,N) )

    call MPI_TYPE_SIZE(MPI_REAL, taille_real, code)

    ! Construction du type qui correspond a 1 bloc de NL lignes et N colonnes
    call MPI_TYPE_VECTOR(N, NL, N, MPI_REAL, type_temp, code)
    taille_tranche = taille_real*NL

    ! On resize le type
    call MPI_TYPE_CREATE_RESIZED(type_temp, borne_inferieure, taille_tranche,type_tranche, code)

    ! On doit appeler MPI_TYPE_COMMIT pour que le type soit effectivement créé
    call MPI_TYPE_COMMIT(type_tranche, code)

    ! Le processus 0 distribue dans AL les tranches horizontales de la matrice A
    call MPI_SCATTER(A, 1, type_tranche, AL, N*NL, MPI_REAL, 0, MPI_COMM_WORLD, code)

    ! Le processus 0 distribue dans BL les tranches verticales de la matrice B

    ! Remarque :
    ! ---------
    !
    !   * Puisque le programme lit les matrices colonnes par colonnes, et non pas ligne par ligne,
    !     on peut se passer de la création d'un nouveau type MPI.
    !
    call MPI_SCATTER(B, N*NL, MPI_REAL, BL, N*NL, MPI_REAL, 0, MPI_COMM_WORLD, code)

    ! Calcul des blocs diagonaux de la matrice resultante.
    CL(rang*NL+1:(rang+1)*NL,:) = matmul( AL(:,:), BL(:,:) )

    ! Calcul des blocs extra-diagonaux

    temps_debut=MPI_WTIME()

  ! Premier algorithme (deux fois plus coûteux que le second)
   if (rang.eq.0 .and. i.eq.1)  then
   print*, "PREMIER ALGORITHME"
   end if
   do k = 0, Nprocs-1
     ! Chaque processus ENVOIE sa tranche AL au processus k
     ! et REÇOIT dans TEMP la tranche AL du processus k
     if (rang /= k) then
      ! On peut utiliser MPI_SENDRECV car les messages sont de même taille
      ! et que les processus envoient et reçoivent en même temps
       call MPI_SENDRECV(AL,   NL*N, MPI_REAL, k, etiquette, &
                         TEMP, NL*N, MPI_REAL, k, etiquette, MPI_COMM_WORLD, statut, code)

       ! Chaque processus calcule les blocs situés au-dessus
       ! et en dessous du bloc de la diagonale principale
       CL(k*NL+1:(k+1)*NL,:)=matmul(TEMP(:,:),BL(:,:))
     end if
   end do

  !   ! Second algorithme
  !  if (rang.eq.0 .and. i.eq.1) then
  !  print*, "SECOND ALGORITHME"
  !  end if
  !   rang_precedent = mod(Nprocs+rang-1,Nprocs)
  !   rang_suivant   = mod(rang+1,Nprocs)
  !   do k = 1, Nprocs-1
  !     ! Chaque processus ENVOIE sa tranche AL au processus précédent
  !     ! et REÇOIT la tranche AL du processus suivant (mais les contenus changent)

  !     ! On peut utiliser MPI_SENDRECV_REPLACE car les messages sont de même taille
  !     ! et que les processus envoient et reçoivent en même temps (mais les contenus changent)
  !     call MPI_SENDRECV_REPLACE(AL, NL*N, MPI_REAL, rang_precedent, etiquette, &
  !                               rang_suivant, etiquette, MPI_COMM_WORLD, statut, code)

  !     ! Chaque processus calcule les blocs situés au-dessus
  !     ! et en dessous du bloc de la diagonale principale
  !     CL(mod(rang+k,Nprocs)*NL+1:(mod(rang+k,Nprocs)+1)*NL,:)=matmul(AL(:,:),BL(:,:))
  !   end do

    temps_fin=MPI_WTIME()

    ! Le processus 0 collecte les tranches CL de tous les processus
    ! pour former la matrice résultante C
    call MPI_GATHER(CL, NL*N, MPI_REAL, C, NL*N, MPI_REAL, 0, MPI_COMM_WORLD, code)

    ! Les tableaux locaux sont désormais inutiles
    deallocate( AL, BL, CL, TEMP )

    ! De même, les types MPI ne sont plus nécessaires
    call MPI_TYPE_FREE(type_temp,code)
    call MPI_TYPE_FREE(type_tranche,code)

    ! --------------------------
    ! Vérification des résultats
    ! --------------------------
    if (rang == 0) then
      allocate( E(N,N) )
      E(:,:) = abs(C(:,:) - CC(:,:))
      Emax   = maxval( E(:,:) ) / N**2
      deallocate( A, B, C, CC, E )

      if ( Emax <= epsilon(1.0) ) then
        ! print'(/,40X,"Bravo !",/,  &
        !       & 20X,"Le produit matriciel A*B calcule en parallele",/, &
        !       & 20X,"est bien egal a celui calcule en monoprocesseur")'
        print ('("Pour la valeur N = ",i4, " le calcul est juste" &
        & " et effectue en ",f8.6," secondes")'),N,temps_fin-temps_debut
      else
        ! print'(/,33X,"Resultat incorrect !",/, &
        !       & 20X,"Le produit matriciel A*B calcule en parallele",/, &
        !       & 20X,"est different de celui calcule en monoprocesseur")'
        print ('("Pour la valeur N = ",i4, " le calcul est faux")'),N
      end if
    end if
  enddo
  ! Nettoyage MPI 
  call MPI_FINALIZE(code)

  ! Fin du programme
end program produit_matrices