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
        integer:: my_arr(3,3), arr(3,3), transf(2)
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
        print *, "I", rank, "have",pos_to_transfer
        r_xyz = 0
        pos = 0
        ! Calculate positions for this proc's atoms
            indx = 0
            do i = 1, nint(Nat ** (1.0d0 / 3.0d0))
                do j = 1, nint(Nat ** (1.0d0 / 3.0d0))
                    do k = 1, nint(Nat ** (1.0d0 / 3.0d0))
                        indx = indx + 1
                        if (indx >= start_atom .and. indx <= end_atom) then ! save only the atoms associated to processor
                            r_xyz(indx, :) = (/(i - 1) ,(j - 1) ,(k - 1) /)
                            ! r_xyz(indx, 2) = (j - 1) * a
                            ! r_xyz(indx, 3) = (k - 1) * a
                            print*, rank, indx,"---", r_xyz(indx,1),r_xyz(indx,2),r_xyz(indx,3)
                        end if
                    end do
                end do
            end do
        print *, "I", rank, "finished with my particles", atoms_per_proc
        ! Synchronize before writing to file
        

        displs = 0
        if (rank==1) then
            print *, pos_to_transfer
            print *, r_xyz
            print *, "----"
            print *, r_xyz(125,:)
            ! do i=1,Nat
            !     print *, r_xyz(i,:)
            ! end do
        end if
        print*, displs


        my_arr = rank+2
        arr = 5
        print*, "IIIIIII", rank, my_arr(rank+1,:)
        transf=(/3,6/)
        start_atom = rank+1
        end_atom = rank*2+1
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)
        ! call MPI_ALLGATHERV(r_xyz(start_atom:end_atom, :), atoms_per_proc, MPI_REAL8, pos(start_atom:end_atom,:), pos_to_transfer, (/0,186/) ,MPI_REAL8, MPI_COMM_WORLD, ierror)
        call MPI_ALLGATHERV(my_arr(start_atom:end_atom,:), transf(rank+1), MPI_INT, arr, transf,(/0,3/) ,MPI_INT, MPI_COMM_WORLD, ierror)
        ! print *, "I", rank, "have",arr
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)


        if (rank == 0) then
            print *, "iiiiiiiiiiiiiii"
            print*, arr
            do i = 1,3
                do j = 1,3
                    print*,arr(i,j)
                enddo
            enddo
            print *, "iiiiiiiiiiiiiii"



            open(unit=21, file=fileout_scc, status='replace')
            write(21, *) Nat
            write(21, *) 'Initial SCC configuration'
            do i = 1, Nat
                write(21, "(a4,3f12.6)") 'Atom', pos(i, 1), pos(i, 2), pos(i, 3) !"(a4,3f12.6)"
            end do
            close(21)
        end if

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

