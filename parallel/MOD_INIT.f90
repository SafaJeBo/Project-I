module MOD_INIT
    use mpi
    implicit none

contains
    ! INITIALIZE PARTICLES IN A SIMPLE CUBIC CONFIGURATION AND DISTRIBUTE PARTICLES !
    subroutine do_SCC(Nat, L, r_xyz, atoms_list, nprocs, rank, fileout_scc)
        ! Nat: number of atoms
        ! L: size of one side of the simulation box
        ! r_xyz: OUTPUT, position of the particles
        ! atoms_list: OUTPUT, list of the indexes of the atoms corresponding to the processor
        ! nprocs: number of processors
        ! rank: number associated at the current processor
        ! fileout_scc: name of output file where to print positions
        implicit none
        integer, intent(in) :: Nat, nprocs, rank
        integer, intent(out), allocatable :: atoms_list(:)
        real(8), intent(in) :: L
        real(8), intent(out) :: r_xyz(Nat,3)
        real(8):: pos(Nat,3)
        integer :: pos_to_transfer(nprocs),displs(nprocs)
        character(len=*), intent(in) :: fileout_scc

        integer :: M,i, j, k, indx, atoms_per_proc, start_atom, end_atom, ierror
        real(8) :: a

        ! Calculate lattice parameter
        M=int(Nat**(1.d0/3.d0))
        a=L/real(M,8)

        if (rank==0) print *, L, M, a

        ! Determine number of atoms per proc and assign start and end indices
        atoms_per_proc = Nat / nprocs
        start_atom = rank * atoms_per_proc + 1
        end_atom = start_atom + atoms_per_proc - 1

        if (rank == nprocs - 1) then
            atoms_per_proc = atoms_per_proc + Nat-end_atom
            end_atom = Nat  ! Ensure last proc gets any extra atoms    
        end if

        ! Save indexes in list
        allocate(atoms_list(atoms_per_proc))
        indx = 1
        do i = start_atom, end_atom
            atoms_list(indx) = i
            indx = indx + 1
        end do

        ! Generate an array with all the number of positions that will be sent later
        call MPI_ALLGATHER(atoms_per_proc,1,MPI_INT,pos_to_transfer,1,MPI_INT,MPI_COMM_WORLD, ierror)

        ! Calculate displs
        displs(1) = 0
        do i = 2, nprocs
            displs(i) = displs(i-1)+pos_to_transfer(i)
        end do
        print *, "I", rank, "have",pos_to_transfer, displs

        r_xyz = 0
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
        print *, "I", rank, "finished with my particles", atoms_per_proc
        
        ! Synchronize before writing to file

        call MPI_ALLGATHERV(r_xyz(start_atom:end_atom, 1), pos_to_transfer(rank+1), MPI_REAL8,  &
        r_xyz(:,1), pos_to_transfer, displs, MPI_REAL8, MPI_COMM_WORLD, ierror)
        call MPI_ALLGATHERV(r_xyz(start_atom:end_atom, 2), pos_to_transfer(rank+1), MPI_REAL8,  &
        r_xyz(:,2), pos_to_transfer, displs, MPI_REAL8, MPI_COMM_WORLD, ierror)
        call MPI_ALLGATHERV(r_xyz(start_atom:end_atom, 3), pos_to_transfer(rank+1), MPI_REAL8,  &
        r_xyz(:,3), pos_to_transfer, displs, MPI_REAL8, MPI_COMM_WORLD, ierror)


        if (rank == 0) then
            open(unit=21, file=fileout_scc, status='replace')
            write(21, *) Nat
            write(21, *) 'Initial SCC configuration'
            do i = 1, Nat
                write(21, "(a4,3f12.6)") 'Atom', r_xyz(i, 1), r_xyz(i, 2), r_xyz(i, 3) !"(a4,3f12.6)"
            end do
            close(21)
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
        real(16), intent(in) :: L
        real(16), dimension(Nat,3), intent(inout) :: r_xyz

        integer :: i, j

        do i=1,Nat
            do j=1,3
                if (r_xyz(i,j)>L/2.d0) then
                        r_xyz(i,j) = r_xyz(i,j)-dble(L)

                else if (r_xyz(i,j)<L/2.d0) then
                        r_xyz(i,j) = r_xyz(i,j)+dble(L)
                endif
            enddo
        enddo
   end subroutine PBC


end module MOD_INIT

