module MOD_INIT
    use mpi
    implicit none

contains
    ! INITIALIZE PARTICLES IN A SIMPLE CUBIC CONFIGURATION AND DISTRIBUTE PARTICLES !
    subroutine do_SCC(Nat, L, r_xyz, atoms_list, nprocs, rank, fileout_scc,pos_to_transfer,start_atom,end_atom,displs)
        ! Nat: number of atoms
        ! L: size of one side of the simulation box
        ! r_xyz: OUTPUT, position of the particles
        ! atoms_list: OUTPUT, list of the indexes of the atoms corresponding to the processor
        ! nprocs: number of processors
        ! rank: number associated at the current processor
        ! fileout_scc: name of output file where to print positions
        ! pos_to_transfer: array with number of data to transfer/receive from each processors in ALLGATHER
        ! start_atom: first atom index to consider for this processor
        ! end_atom: last atom index to consider for this processor
        ! displs: displacement matrix for ALLGATHER
        implicit none
        integer, intent(in) :: Nat, nprocs, rank
        integer, intent(out), allocatable :: atoms_list(:)
        real(8), intent(in) :: L
        real(8), intent(out) :: r_xyz(Nat,3)
        real(8):: pos(Nat,3)
        integer :: pos_to_transfer(nprocs),displs(nprocs)
        character(len=*), intent(in) :: fileout_scc

        integer :: M,i, j, k, indx, atoms_per_proc, start_atom, end_atom, ierror, n_atoms_remaining
        real(8) :: a

        ! Calculate lattice parameter
        M=int(Nat**(1.d0/3.d0))
        a=L/real(M,8)
        if (rank==0) print *, L, M, a

        ! Calculate positions for this proc's atoms
            indx = 0
            do i = 1, nint(Nat ** (1.0d0 / 3.0d0))
                do j = 1, nint(Nat ** (1.0d0 / 3.0d0))
                    do k = 1, nint(Nat ** (1.0d0 / 3.0d0))
                        indx = indx + 1
                        if (indx >= start_atom .and. indx <= end_atom) then ! save only the atoms associated to processor
                            r_xyz(indx, :) = (/(i - 1)*a ,(j - 1)*a ,(k - 1)*a /)
                        end if
                    end do
                end do
            end do
        
        ! Synchronize before writing to file
        call MPI_ALLGATHERV(r_xyz(start_atom:end_atom, 1), pos_to_transfer(rank+1), MPI_REAL8,  r_xyz(:,1), pos_to_transfer, displs, MPI_REAL8, MPI_COMM_WORLD, ierror)
        call MPI_ALLGATHERV(r_xyz(start_atom:end_atom, 2), pos_to_transfer(rank+1), MPI_REAL8,  r_xyz(:,2), pos_to_transfer, displs, MPI_REAL8, MPI_COMM_WORLD, ierror)
        call MPI_ALLGATHERV(r_xyz(start_atom:end_atom, 3), pos_to_transfer(rank+1), MPI_REAL8,  r_xyz(:,3), pos_to_transfer, displs, MPI_REAL8, MPI_COMM_WORLD, ierror)

        if (rank == 0) then
            open(unit=25, file=fileout_scc, status='replace')
            write(25, *) Nat
            write(25, *) 'Initial SCC configuration'
            do i = 1, Nat
                write(25, "(a4,3f12.6)") 'Atom', r_xyz(i, 1), r_xyz(i, 2), r_xyz(i, 3) !"(a4,3f12.6)"
            end do
            close(25)
        end if
        return

    end subroutine do_SCC


    ! APPLY PERIODIC BOUNDARY CONDITIONS ! 
    subroutine PBC(Nat,L,r_xyz)
        ! Nat: number of atoms
        ! L: size of one side of the simulation box
        ! r_xyz: INPUT & OUTPUT, position of the particles
        implicit none
        integer, intent(in) :: Nat
        real(8), intent(in) :: L
        real(8), dimension(Nat,3), intent(inout) :: r_xyz

        integer :: i, j

        do i=1,Nat
            do j=1,3
                if (r_xyz(i,j)>L/2.d0) then
                        r_xyz(i,j) = r_xyz(i,j)-dble(L)

                else if (r_xyz(i,j)<-L/2.d0) then
                        r_xyz(i,j) = r_xyz(i,j)+dble(L)
                endif
            enddo
        enddo
   end subroutine PBC


end module MOD_INIT

