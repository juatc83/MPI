program reduce
  use mpi
  implicit none
  integer :: nb_procs,rang,code
  integer, parameter :: dp = kind(1.d0)
  integer, parameter :: li = selected_int_kind(15)
  integer(kind=li):: nbbloc,i
  real(kind=dp) :: largeur,somme,valeur,x,longueur_tranche,pi

  call MPI_INIT (code)
  call MPI_COMM_SIZE ( MPI_COMM_WORLD ,nb_procs,code)
  call MPI_COMM_RANK ( MPI_COMM_WORLD ,rang,code)

  ! Nombre d'intervalles
  nbbloc = 3*1000*1000_li*100
  ! largeur des intervalles
  largeur = 1._dp / real(nbbloc,dp)

  longueur_tranche = nbbloc/real(nbbloc,dp)

  somme = 0._dp

  do i = longueur_tranche*(rang)+1, longueur_tranche*(rang+1)
    ! Point au milieu de l'intervalle
    x = largeur*(i-0.5_dp)
    ! Calcul de l'aire
    somme = somme + largeur*(4._dp / (1._dp + x*x))
  end do


  call MPI_REDUCE (somme,pi,1, MPI_DOUBLE_PRECISION , MPI_SUM ,0, MPI_COMM_WORLD,code)

  if(rang==0)then
    print *, "Pi =", pi  
  endif

  call MPI_FINALIZE (code)
  end program reduce