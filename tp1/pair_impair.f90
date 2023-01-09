
program qui_je_suis

    use mpi

    implicit none
    integer :: nb_procs,rang,code


    call MPI_INIT (code)
    call MPI_COMM_SIZE ( MPI_COMM_WORLD ,nb_procs,code)
    call MPI_COMM_RANK ( MPI_COMM_WORLD ,rang,code)

    if (mod(rang,2)==0)then
        print *,'Je suis le processus pair de rang ',rang
    else
        print *,'Je suis le processus impair de rang ',rang
    endif
    
    call MPI_FINALIZE (code)

end program qui_je_suis
