module forces
    !use mpi
    use MOD_INIT, only: PBC
    implicit none
    contains    
    subroutine find_force_LJ(pos,N,d,L,force,cutoff,Upot)
        implicit none
        include 'mpif.h'
        integer :: N,d,i,j,rank,size,ierror,start,end,k,n_atoms_per_proc,n_atoms_this_proc,n_atoms_remaining
        integer, allocatable :: displs(:),counts(:)
        real(8) :: pos(N,d),force(N,d),dx,dy,dz,dr2,cutoff,L,cf2,Upot,dr6,dr12,potcut,fij,local_Upot,dr(3)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
        local_Upot=0d0
        Upot=0d0
        force=0d0
        cf2=cutoff*cutoff
        potcut=4.d0*(1.d0/cf2**6-1.d0/cf2**3)
        n_atoms_per_proc = N / size
        start = rank * n_atoms_per_proc + 1
        end = min((rank + 1) * n_atoms_per_proc, N)
        n_atoms_remaining = mod(N, size)

        if (rank < n_atoms_remaining) then
            n_atoms_this_proc = n_atoms_per_proc + 1
        else
            n_atoms_this_proc = n_atoms_per_proc
        end if
        do i=start,end
            do j=1,N
                if (i.ne.j) then
                    dx=pos(i,1)-pos(j,1)
                    dy=pos(i,2)-pos(j,2)
                    dz=pos(i,3)-pos(j,3)
                    dr(1)=dx; dr(2)=dy; dr(3)=dz
                    call pbc(N,L,dr)
                    dr2=dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
                    if (dr2.lt.cf2) then
                        dr6=dr2*dr2*dr2
                        dr12=dr6*dr6
                        fij=48.d0/(dr6*dr6*dr2)-24.d0/(dr2*dr6)
                        force(i,1)=force(i,1)+fij*dr(1)
                        force(i,2)=force(i,2)+fij*dr(2)
                        force(i,3)=force(i,3)+fij*dr(3)
                        local_Upot=local_Upot+4.d0*(1.d0/dr12-1.d0/dr6)-potcut
                        
                    endif
                endif
            enddo
        enddo
        local_Upot=local_Upot/2d0
        allocate(counts(size), displs(size))

        do i = 0, size - 1
            if (i < n_atoms_remaining) then
                counts(i + 1) = n_atoms_per_proc + 1
            else
                counts(i + 1) = n_atoms_per_proc
            end if
        end do

        displs(1) = 0
        do i = 2, size
            displs(i) = displs(i - 1) + counts(i - 1)
        end do
        call MPI_ALLGATHERV(force(start:end,1), counts(rank+1), MPI_DOUBLE_PRECISION, force(:,1), counts, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_ALLGATHERV(force(start:end,2), counts(rank+1), MPI_DOUBLE_PRECISION, force(:,2), counts, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_ALLGATHERV(force(start:end,3), counts(rank+1), MPI_DOUBLE_PRECISION, force(:,3), counts, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allreduce(local_Upot,Upot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
        if (rank.eq.1) then
            do k=1,N
                write(20,*)rank,k,force(k,1),force(k,2),force(k,3)
            enddo
        endif
        deallocate(displs,counts)
    return
    end
end module forces