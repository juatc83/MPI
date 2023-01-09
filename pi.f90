program pi
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer, parameter :: li = selected_int_kind(15)
  integer(kind=li):: nbbloc,i
  real(kind=dp) :: largeur,somme,x

  ! Nombre d'intervalles
  nbbloc = 3*1000*1000_li*100
  ! largeur des intervalles
  largeur = 1._dp / real(nbbloc,dp)

  somme = 0._dp

  do i=1, nbbloc
    ! Point au milieu de l'intervalle
    x = largeur*(i-0.5_dp)
    ! Calcul de l'aire
    somme = somme + largeur*(4._dp / (1._dp + x*x))
  end do

  print *, "Pi =", somme
end program