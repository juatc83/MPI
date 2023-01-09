!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ping_pong_1.f90 --- TP2 : Communications point � point :
!!                           envoi d'un message du processus 0 au processus 1
!! 
!! Auteur          : Denis GIROU (CNRS/IDRIS - France) <Denis.Girou@idris.fr>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program ping_pong_1
  USE MPI
  implicit none

  integer, dimension(MPI_STATUS_SIZE)  :: statut
  integer, parameter                   :: nb_valeurs=1000,etiquette=99
  integer                              :: rang,code
  integer, parameter                   :: dp = kind(1.d0)
  real(kind=dp), dimension(nb_valeurs) :: valeurs

  call MPI_INIT (code)
  call MPI_COMM_RANK ( MPI_COMM_WORLD ,rang,code)

  if (rang == 0) then
    call random_number(valeurs)
    call MPI_SEND (valeurs,nb_valeurs, MPI_DOUBLE_PRECISION ,1,etiquette, MPI_COMM_WORLD,code)
  elseif (rang == 1) then
    call MPI_RECV (valeurs,nb_valeurs, MPI_DOUBLE_PRECISION ,0,etiquette, MPI_COMM_WORLD,statut,code)
    print ('("Moi, processus 1, j''ai recu ",i4," valeurs (derniere = ", &
    & f4.2,") du processus 0.")'), nb_valeurs,valeurs(nb_valeurs)
  end if
  
  call MPI_FINALIZE (code)

end program ping_pong_1


