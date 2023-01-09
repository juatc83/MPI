!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ping_pong_3.f90 --- TP2 : Communications point � point :
!!                           ping-pong pour des tailles variables de messages
!! 
!! Auteur          : Denis GIROU (CNRS/IDRIS - France) <Denis.Girou@idris.fr>
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program ping_pong_3
  USE MPI
  implicit none

  integer, dimension(MPI_STATUS_SIZE)          :: statut
  integer, parameter                           :: nb_valeurs_max=7000000, &
       nb_tests=9,etiquette=99
  integer, dimension(nb_tests)                 :: nb_valeurs
  integer                                      :: rang,code,i
  integer, parameter                           :: dp = kind(1.d0)
  real(kind=dp), dimension(0:nb_valeurs_max-1) :: valeurs
  real(kind=dp)                                :: temps_debut,temps_fin

  call MPI_INIT (code)
  call MPI_COMM_RANK ( MPI_COMM_WORLD ,rang,code)

  nb_valeurs = (/ 0,1,10,100,1000,10000,100000,1000000,7000000 /)

  do i=1, nb_tests
    if (rang == 0) then
      call random_number(valeurs)
      temps_debut=MPI_WTIME()
      call MPI_SEND (valeurs,nb_valeurs(i), MPI_DOUBLE_PRECISION ,1,etiquette, MPI_COMM_WORLD,code)
      call MPI_RECV (valeurs,nb_valeurs(i), MPI_DOUBLE_PRECISION ,1,etiquette, MPI_COMM_WORLD,statut,code)  
      temps_fin=MPI_WTIME()
      call affichage
    elseif (rang == 1) then
      call MPI_RECV (valeurs,nb_valeurs(i), MPI_DOUBLE_PRECISION ,0,etiquette, MPI_COMM_WORLD,statut,code)
      call MPI_SEND (valeurs,nb_valeurs(i), MPI_DOUBLE_PRECISION ,0,etiquette, MPI_COMM_WORLD,code)
    end if
  enddo
  
  call MPI_FINALIZE (code) 

contains
  subroutine affichage

    if (nb_valeurs(i)/=0) then
      print ('("Moi, processus 0, j''ai envoye et recu ",i8, &
            & " valeurs (derniere = ",f4.2,") du processus 1", &
            & " en ",f8.6," secondes, soit avec un debit de ",f7.2, &
            & " Mo/s.")'), &
            nb_valeurs(i),valeurs(nb_valeurs(i)-1),temps_fin-temps_debut, &
            real(2*nb_valeurs(i)*8)/1000000./(temps_fin-temps_debut)
    else
      print ('("Moi, processus 0, j''ai envoye et recu ",i8, &
            & " valeurs en ",f8.6," secondes, soit avec un debit de ",f7.2, &
            & " Mo/s.")'), &
            nb_valeurs(i),temps_fin-temps_debut, &
            real(2*nb_valeurs(i)*8)/1000000./(temps_fin-temps_debut)
    end if
  end subroutine affichage

end program ping_pong_3
