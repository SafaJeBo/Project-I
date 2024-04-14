module integrate
    !use mpi
    use forces, only: find_force_LJ
    use MOD_INIT, only: PBC,pbc2
    implicit none
    contains
    subroutine time_step_vVerlet(pos,N,d,L,vel,dt,cutoff,nu,sigma,Upot)
        
        implicit none
        include 'mpif.h'
        integer :: N,d,i,rank,size,start,end,k,ierror,n_atoms_per_proc,n_atoms_remaining,n_atoms_this_proc
        integer, allocatable :: displs(:),counts(:)
        !integer :: displs
        real(8) :: pos(N,d),force(N,d),dt,L,cutoff,vel(N,d),nu,sigma,Upot,dr(3)
        real(8),allocatable :: local_velx(:),local_vely(:),local_velz(:),local_posx(:),local_posy(:),local_posz(:)

        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)

        n_atoms_per_proc = N / size
        start = rank * n_atoms_per_proc + 1
        end = min((rank + 1) * n_atoms_per_proc, N)
        n_atoms_remaining = mod(N, size)

        ! Calculate the number of atoms this process will handle
        if (rank < n_atoms_remaining) then
            n_atoms_this_proc = n_atoms_per_proc + 1
        else
            n_atoms_this_proc = n_atoms_per_proc
        end if

        call find_force_LJ(pos,N,d,L,force,cutoff,Upot)
        

        do i=start,end
            pos(i,1)=pos(i,1)+vel(i,1)*dt+0.5d0*force(i,1)*dt*dt
            pos(i,2)=pos(i,2)+vel(i,2)*dt+0.5d0*force(i,2)*dt*dt
            pos(i,3)=pos(i,3)+vel(i,3)*dt+0.5d0*force(i,3)*dt*dt
            vel(i,1)=vel(i,1)+force(i,1)*0.5d0*dt
            vel(i,2)=vel(i,2)+force(i,2)*0.5d0*dt
            vel(i,3)=vel(i,3)+force(i,3)*0.5d0*dt
        enddo

        !print*,'asd1'
        !allocate(displs(size), counts(size))
        !displs=rank*N/size+1
        allocate(counts(size), displs(size))

        ! Set the counts array
        do i = 0, size - 1
            if (i < n_atoms_remaining) then
                counts(i + 1) = n_atoms_per_proc + 1
            else
                counts(i + 1) = n_atoms_per_proc
            end if
        end do

        ! Set the displacements array
        displs(1) = 0
        do i = 2, size
            displs(i) = displs(i - 1) + counts(i - 1)
        end do

        call MPI_Allgatherv(vel(start:end,1), counts(rank+1), MPI_DOUBLE_PRECISION, vel(:,1), counts, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(vel(start:end,2), counts(rank+1), MPI_DOUBLE_PRECISION, vel(:,2), counts, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(vel(start:end,3), counts(rank+1), MPI_DOUBLE_PRECISION, vel(:,3), counts, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)

        !if (rank.eq.0) then

        do i=start,end
            dr(1)=pos(i,1); dr(2)=pos(i,2);dr(3)=pos(i,3)
            call pbc(N,L,dr)
            pos(i,1)=dr(1);pos(i,2)=dr(2);pos(i,3)=dr(3)
        enddo

        call MPI_Allgatherv(pos(start:end,1), n_atoms_this_proc, MPI_DOUBLE_PRECISION, pos(:,1), n_atoms_this_proc, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(pos(start:end,2), n_atoms_this_proc, MPI_DOUBLE_PRECISION, pos(:,2), n_atoms_this_proc, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(pos(start:end,3), n_atoms_this_proc, MPI_DOUBLE_PRECISION, pos(:,3), n_atoms_this_proc, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)

        if (rank.eq.0) then
            do k=1,N
                write(21,*)rank,k,pos(k,1),pos(k,2),pos(k,3)
                write(21,*)rank,k,vel(k,1),vel(k,2),vel(k,3)
            enddo
        endif


        call find_force_LJ(pos,N,d,L,force,cutoff,Upot)
        do i=start,end
            vel(i,1)=vel(i,1)+force(i,1)*0.5d0*dt
            vel(i,2)=vel(i,2)+force(i,2)*0.5d0*dt
            vel(i,3)=vel(i,3)+force(i,3)*0.5d0*dt
        enddo
        
        call MPI_Allgatherv(vel(start:end,1), n_atoms_this_proc, MPI_DOUBLE_PRECISION, vel(:,1), n_atoms_this_proc, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(vel(start:end,2), n_atoms_this_proc, MPI_DOUBLE_PRECISION, vel(:,2), n_atoms_this_proc, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(vel(start:end,3), n_atoms_this_proc, MPI_DOUBLE_PRECISION, vel(:,3), n_atoms_this_proc, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        

!        call andersen_thermo(vel,N,d,nu,sigma)
    return
    end



!    subroutine andersen_thermo(vel,N,d,nu,sigma)
!        implicit none
!        integer :: N,d,i
!        real(8) :: sigma,vel(N,d),nu
    
!        do i=1,N
!            if (rand().lt.nu) then
!                call gauss(0d0,sigma,vel(i,1))
!                call gauss(0d0,sigma,vel(i,2))
!                call gauss(0d0,sigma,vel(i,3))
!            endif
!        enddo
!    return
!    end
    
!    subroutine gauss(mean,var,val)
!        implicit none
!        real(8) :: chi1,chi2,pi,var,mean,val
!        parameter(pi=4.d0*atan(1.d0))
!        chi1=rand()
!        chi2=rand()
!        val=var*dsqrt(-2.d0*dlog(1.d0-chi1))*dcos(2.d0*pi*chi2)+mean
!    return
!    end

end module integrate