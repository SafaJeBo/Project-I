module forces
    use mpi
    use MOD_INIT, only: PBC
    implicit none
    contains

    ! CALCULATES FORCES BETWEEN ALL PAIRS FOR THE PARTICLES ASSOCIATED TO THE PROCESSOR ! 
    subroutine find_force_LJ(nprocs, pos,N,d,L,force,cutoff,Upot,pos_to_transfer,start_atom,end_atom,displs)
        ! nprocs: number of total processors
        ! pos: positions array for all particles
        ! N: number of particles
        ! d: number of dimensions in system
        ! L: size of one side of the simulation box
        ! force: OUTPUT, total force from all processors
        ! cutoff: maximum particle distance where there is interaction
        ! Upot: OUTPUT, total potential energy from all processors
        ! pos_to_transfer: array with number of data to transfer/receive from each processors in ALLGATHER
        ! start_atom: first atom index to consider for this processor
        ! end_atom: last atom index to consider for this processor
        ! displs: displacement matrix for ALLGATHER
        implicit none
        integer, intent(in) :: nprocs
        integer :: N, d,i,j,rank,ierror,start_atom,end_atom,k
        integer, dimension(nprocs) :: displs, pos_to_transfer
        real(8) :: pos(N,d),force(N,d),dx,dy,dz,dr2,cutoff,L,cf2,Upot,dr6,dr12,potcut,fij,local_Upot,dr(3)
        
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
        local_Upot=0d0
        Upot=0d0
        force=0d0
        cf2=cutoff*cutoff
        potcut=4.d0*(1.d0/cf2**6-1.d0/cf2**3)

        ! Calculate forces and potential energy for assigned atoms in this process
        do i=start_atom,end_atom
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
        
        ! Distribute forces to all processors and sum all the potential energies on MASTER processor
        local_Upot=local_Upot/2d0
        call MPI_ALLGATHERV(force(start_atom:end_atom,1), pos_to_transfer(rank+1), MPI_REAL8, force(:,1), pos_to_transfer, displs, MPI_REAL8, MPI_COMM_WORLD, ierror)
        call MPI_ALLGATHERV(force(start_atom:end_atom,2), pos_to_transfer(rank+1), MPI_REAL8, force(:,2), pos_to_transfer, displs, MPI_REAL8, MPI_COMM_WORLD, ierror)
        call MPI_ALLGATHERV(force(start_atom:end_atom,3), pos_to_transfer(rank+1), MPI_REAL8, force(:,3), pos_to_transfer, displs, MPI_REAL8, MPI_COMM_WORLD, ierror)
        call MPI_Allreduce(local_Upot,Upot,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierror)
        
    return
    end
end module forces
