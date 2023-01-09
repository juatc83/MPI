
call MPI_INIT (code)
call MPI_COMM_RANK ( MPI_COMM_WORLD ,rang,code)
if (rang == 2) then
valeur=1000
call MPI_SEND (valeur,1, MPI_INTEGER ,5,etiquette, MPI_COMM_WORLD,code)
elseif (rang == 5) then
call MPI_RECV (valeur,1, MPI_INTEGER ,2,etiquette, MPI_COMM_WORLD,statut,code)
print *,’Moi, processus 5, ai reçu ’,valeur,’ du processus 2.’
end if
call MPI_FINALIZE (code)
end program point_a_point