module integrate

    use mpi
    use forces, only: find_force_LJ
    use MOD_INIT, only: PBC
    implicit none

    contains

     subroutine time_step_vVerlet(nprocs,pos,N,d,L,vel,dt,cutoff,nu,sigma,Upot,force,pos_to_transfer,start_atom,end_atom,displs)
        implicit none

        !include 'mpif.h'
        integer, intent(in) :: nprocs
        integer :: Nat, N,d,i,rank,start_atom,end_atom,k,ierror
        integer, dimension(nprocs) :: displs, pos_to_transfer
        real(8) :: pos(N,d),force(N,d),dt,L,cutoff,vel(N,d),nu,sigma,Upot,dr(3)
        real(8),allocatable :: local_velx(:),local_vely(:),local_velz(:),local_posx(:),local_posy(:),local_posz(:)
       
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
        !call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)

!       n_atoms_per_proc = N / size
!       n_atoms_remaining = mod(N, size)
!
!       ! Calculate the number of atoms this process will handle
!   if (rank < n_atoms_remaining) then
!       n_atoms_this_proc = n_atoms_per_proc + 1
!       start = rank * n_atoms_this_proc + 1
!   else
!       n_atoms_this_proc = n_atoms_per_proc
!       start = n_atoms_remaining * (n_atoms_per_proc + 1) + (rank - n_atoms_remaining) * n_atoms_per_proc + 1
!   end if
!
!   end = start + n_atoms_this_proc - 1
!
!       ! Calculate counts
!   !print*,'start-end',rank,start,end
!   allocate(counts(size))
!   do i = 1, size
!   if (i <= n_atoms_remaining) then
!       counts(i) = (n_atoms_per_proc + 1) 
!   else
!       counts(i) = n_atoms_per_proc 
!   end if
!   end do
!
!   ! Calculate displs
!   allocate(displs(size))
!   displs(1) = 0
!   do i = 2, size
!   if (i <= n_atoms_remaining) then
!       displs(i) = displs(i-1) + (n_atoms_per_proc + 1) 
!   else
!       displs(i) = displs(i-1) + n_atoms_per_proc 
!   end if
!   end do
        
!       allocate(pos_to_transfer(nprocs),displs(nprocs)) 
!       ! Ensure last proc gets any extra atoms 
!       n_atoms_remaining = mod(N, nprocs)
!         
!       !print*,'check velocity',N,nprocs, n_atoms_remaining
!       if (rank < n_atoms_remaining) then
!           atoms_per_proc = N / nprocs + 1
!           start_atom = rank * atoms_per_proc + 1
!           end_atom = start_atom + atoms_per_proc - 1 
!       else
!           atoms_per_proc = N / nprocs
!           start_atom = n_atoms_remaining * (atoms_per_proc + 1) + (rank - n_atoms_remaining) * atoms_per_proc + 1
!           end_atom = start_atom + atoms_per_proc - 1
!       end if
!
!       !print *, "Rank-start-app-end", start_atom,atoms_per_proc,end_atom
!
!       ! Save indexes in list
!       allocate(atoms_list(atoms_per_proc))
!       indx = 1
!       do i = start_atom, end_atom
!           atoms_list(indx) = i
!           indx = indx + 1
!       end do
!       !print*,'beforeallgather'
!
!       ! Generate an array with all the number of positions that will be sent later
!       call MPI_ALLGATHER(atoms_per_proc,1,MPI_INT,pos_to_transfer,1,MPI_INT,MPI_COMM_WORLD, ierror)
!       !print*,'ididalgather'
!       ! Calculate displs
!
!       displs(1) = 0
!       do i = 2, nprocs
!           displs(i) = displs(i-1)+pos_to_transfer(i-1)
!       end do
!       !print*,'thisisdisps'
!       !print *, "I", rank, "have",pos_to_transfer, displs
!


        call find_force_LJ(nprocs,pos,N,d,L,force,cutoff,Upot,pos_to_transfer,start_atom,end_atom,displs)
        

        do i=start_atom,end_atom
            pos(i,1)=pos(i,1)+vel(i,1)*dt+0.5d0*force(i,1)*dt*dt
            pos(i,2)=pos(i,2)+vel(i,2)*dt+0.5d0*force(i,2)*dt*dt
            pos(i,3)=pos(i,3)+vel(i,3)*dt+0.5d0*force(i,3)*dt*dt
            vel(i,1)=vel(i,1)+force(i,1)*0.5d0*dt
            vel(i,2)=vel(i,2)+force(i,2)*0.5d0*dt
            vel(i,3)=vel(i,3)+force(i,3)*0.5d0*dt
        enddo

        !allocate(displs(size), counts(size))
        !displs=rank*N/size+1
        !allocate(counts(size), displs(size))


        call MPI_Allgatherv(vel(start_atom:end_atom,1), pos_to_transfer(rank+1), MPI_REAL8, vel(:,1), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(vel(start_atom:end_atom,2), pos_to_transfer(rank+1), MPI_REAL8, vel(:,2), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(vel(start_atom:end_atom,3), pos_to_transfer(rank+1), MPI_REAL8, vel(:,3), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)

        !if (rank.eq.0) then

        do i=start_atom,end_atom
            dr(1)=pos(i,1); dr(2)=pos(i,2);dr(3)=pos(i,3)
            call pbc(N,L,dr)
            pos(i,1)=dr(1);pos(i,2)=dr(2);pos(i,3)=dr(3)
        enddo

        call MPI_Allgatherv(pos(start_atom:end_atom,1), pos_to_transfer(rank+1), MPI_REAL8, pos(:,1), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(pos(start_atom:end_atom,2), pos_to_transfer(rank+1), MPI_REAL8, pos(:,2), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(pos(start_atom:end_atom,3), pos_to_transfer(rank+1), MPI_REAL8, pos(:,3), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)

        !call MPI_BARRIER(ierror)
        !if (rank.eq.0) then
        !    do k=1,N
        !        write(21,*)rank,k,pos(k,1),pos(k,2),pos(k,3)
        !        write(21,*)rank,k,vel(k,1),vel(k,2),vel(k,3)
        !    enddo
        !endif


        call find_force_LJ(nprocs,pos,N,d,L,force,cutoff,Upot,pos_to_transfer,start_atom,end_atom,displs)
        do i=start_atom,end_atom
            vel(i,1)=vel(i,1)+force(i,1)*0.5d0*dt
            vel(i,2)=vel(i,2)+force(i,2)*0.5d0*dt
            vel(i,3)=vel(i,3)+force(i,3)*0.5d0*dt
        enddo
        
        call MPI_Allgatherv(vel(start_atom:end_atom,1), pos_to_transfer(rank+1), MPI_REAL8, vel(:,1), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(vel(start_atom:end_atom,2), pos_to_transfer(rank+1), MPI_REAL8, vel(:,2), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(vel(start_atom:end_atom,3), pos_to_transfer(rank+1), MPI_REAL8, vel(:,3), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)

        call MPI_Allgatherv(pos(start_atom:end_atom,1), pos_to_transfer(rank+1), MPI_REAL8, pos(:,1), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(pos(start_atom:end_atom,2), pos_to_transfer(rank+1), MPI_REAL8, pos(:,2), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(pos(start_atom:end_atom,3), pos_to_transfer(rank+1), MPI_REAL8, pos(:,3), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        
        !if (mod(iter,100).eq.0) then
        !    do i=1,N    
        !        write(22,*)iter,i,force(i,:)
        !    enddo
        !endif

        call andersen_thermo(vel,N,d,nu,sigma)
       ! deallocate(pos_to_transfer,displs)
    return
    end

    subroutine time_step_vVerlet2(nprocs,pos,N,d,L,vel,dt,cutoff,nu,sigma,Upot,force,pos_to_transfer,start_atom,end_atom,displs)
        implicit none

        !include 'mpif.h'
        integer, intent(in) :: nprocs
        integer :: Nat, N,d,i,rank,start_atom,end_atom,k,ierror
        integer, dimension(nprocs) :: displs, pos_to_transfer
        real(8) :: pos(N,d),force(N,d),dt,L,cutoff,vel(N,d),nu,sigma,Upot,dr(3)
        real(8),allocatable :: local_velx(:),local_vely(:),local_velz(:),local_posx(:),local_posy(:),local_posz(:)
       
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
        !call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)

!       n_atoms_per_proc = N / size
!       n_atoms_remaining = mod(N, size)
!
!       ! Calculate the number of atoms this process will handle
!   if (rank < n_atoms_remaining) then
!       n_atoms_this_proc = n_atoms_per_proc + 1
!       start = rank * n_atoms_this_proc + 1
!   else
!       n_atoms_this_proc = n_atoms_per_proc
!       start = n_atoms_remaining * (n_atoms_per_proc + 1) + (rank - n_atoms_remaining) * n_atoms_per_proc + 1
!   end if
!
!   end = start + n_atoms_this_proc - 1
!
!       ! Calculate counts
!   !print*,'start-end',rank,start,end
!   allocate(counts(size))
!   do i = 1, size
!   if (i <= n_atoms_remaining) then
!       counts(i) = (n_atoms_per_proc + 1) 
!   else
!       counts(i) = n_atoms_per_proc 
!   end if
!   end do
!
!   ! Calculate displs
!   allocate(displs(size))
!   displs(1) = 0
!   do i = 2, size
!   if (i <= n_atoms_remaining) then
!       displs(i) = displs(i-1) + (n_atoms_per_proc + 1) 
!   else
!       displs(i) = displs(i-1) + n_atoms_per_proc 
!   end if
!   end do
        
!       allocate(pos_to_transfer(nprocs),displs(nprocs)) 
!       ! Ensure last proc gets any extra atoms 
!       n_atoms_remaining = mod(N, nprocs)
!         
!       !print*,'check velocity',N,nprocs, n_atoms_remaining
!       if (rank < n_atoms_remaining) then
!           atoms_per_proc = N / nprocs + 1
!           start_atom = rank * atoms_per_proc + 1
!           end_atom = start_atom + atoms_per_proc - 1 
!       else
!           atoms_per_proc = N / nprocs
!           start_atom = n_atoms_remaining * (atoms_per_proc + 1) + (rank - n_atoms_remaining) * atoms_per_proc + 1
!           end_atom = start_atom + atoms_per_proc - 1
!       end if
!
!       !print *, "Rank-start-app-end", start_atom,atoms_per_proc,end_atom
!
!       ! Save indexes in list
!       allocate(atoms_list(atoms_per_proc))
!       indx = 1
!       do i = start_atom, end_atom
!           atoms_list(indx) = i
!           indx = indx + 1
!       end do
!       !print*,'beforeallgather'
!
!       ! Generate an array with all the number of positions that will be sent later
!       call MPI_ALLGATHER(atoms_per_proc,1,MPI_INT,pos_to_transfer,1,MPI_INT,MPI_COMM_WORLD, ierror)
!       !print*,'ididalgather'
!       ! Calculate displs
!
!       displs(1) = 0
!       do i = 2, nprocs
!           displs(i) = displs(i-1)+pos_to_transfer(i-1)
!       end do
!       !print*,'thisisdisps'
!       !print *, "I", rank, "have",pos_to_transfer, displs
!


        call find_force_LJ(nprocs,pos,N,d,L,force,cutoff,Upot,pos_to_transfer,start_atom,end_atom,displs)
        

        do i=start_atom,end_atom
            pos(i,1)=pos(i,1)+vel(i,1)*dt+0.5d0*force(i,1)*dt*dt
            pos(i,2)=pos(i,2)+vel(i,2)*dt+0.5d0*force(i,2)*dt*dt
            pos(i,3)=pos(i,3)+vel(i,3)*dt+0.5d0*force(i,3)*dt*dt
            vel(i,1)=vel(i,1)+force(i,1)*0.5d0*dt
            vel(i,2)=vel(i,2)+force(i,2)*0.5d0*dt
            vel(i,3)=vel(i,3)+force(i,3)*0.5d0*dt
        enddo

        !allocate(displs(size), counts(size))
        !displs=rank*N/size+1
        !allocate(counts(size), displs(size))


        call MPI_Allgatherv(vel(start_atom:end_atom,1), pos_to_transfer(rank+1), MPI_REAL8, vel(:,1), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(vel(start_atom:end_atom,2), pos_to_transfer(rank+1), MPI_REAL8, vel(:,2), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(vel(start_atom:end_atom,3), pos_to_transfer(rank+1), MPI_REAL8, vel(:,3), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)

        !if (rank.eq.0) then

        do i=start_atom,end_atom
            dr(1)=pos(i,1); dr(2)=pos(i,2);dr(3)=pos(i,3)
            call pbc(N,L,dr)
            pos(i,1)=dr(1);pos(i,2)=dr(2);pos(i,3)=dr(3)
        enddo

        call MPI_Allgatherv(pos(start_atom:end_atom,1), pos_to_transfer(rank+1), MPI_REAL8, pos(:,1), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(pos(start_atom:end_atom,2), pos_to_transfer(rank+1), MPI_REAL8, pos(:,2), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(pos(start_atom:end_atom,3), pos_to_transfer(rank+1), MPI_REAL8, pos(:,3), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)

        !call MPI_BARRIER(ierror)
        !if (rank.eq.0) then
        !    do k=1,N
        !        write(21,*)rank,k,pos(k,1),pos(k,2),pos(k,3)
        !        write(21,*)rank,k,vel(k,1),vel(k,2),vel(k,3)
        !    enddo
        !endif


        call find_force_LJ(nprocs,pos,N,d,L,force,cutoff,Upot,pos_to_transfer,start_atom,end_atom,displs)
        do i=start_atom,end_atom
            vel(i,1)=vel(i,1)+force(i,1)*0.5d0*dt
            vel(i,2)=vel(i,2)+force(i,2)*0.5d0*dt
            vel(i,3)=vel(i,3)+force(i,3)*0.5d0*dt
        enddo
        
        call MPI_Allgatherv(vel(start_atom:end_atom,1), pos_to_transfer(rank+1), MPI_REAL8, vel(:,1), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(vel(start_atom:end_atom,2), pos_to_transfer(rank+1), MPI_REAL8, vel(:,2), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(vel(start_atom:end_atom,3), pos_to_transfer(rank+1), MPI_REAL8, vel(:,3), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)

        call MPI_Allgatherv(pos(start_atom:end_atom,1), pos_to_transfer(rank+1), MPI_REAL8, pos(:,1), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(pos(start_atom:end_atom,2), pos_to_transfer(rank+1), MPI_REAL8, pos(:,2), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(pos(start_atom:end_atom,3), pos_to_transfer(rank+1), MPI_REAL8, pos(:,3), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        
        !if (mod(iter,100).eq.0) then
        !    do i=1,N    
        !        write(22,*)iter,i,force(i,:)
        !    enddo
        !endif

       ! call andersen_thermo(vel,N,d,nu,sigma)
       ! deallocate(pos_to_transfer,displs)
    return
    end




    subroutine andersen_thermo(vel,N,d,nu,sigma)
        implicit none
        integer :: N,d,i
        real(8) :: sigma,vel(N,d),nu, rand
    
        do i=1,N
            call RANDOM_NUMBER(rand)
            if (rand.lt.nu) then
                call gauss(0d0,sigma,vel(i,1))
                call gauss(0d0,sigma,vel(i,2))
                call gauss(0d0,sigma,vel(i,3))
            endif
        enddo
    return
    end
    
    subroutine gauss(mean,var,val)
        implicit none
        real(8) :: chi1,chi2,pi,var,mean,val
        parameter(pi=4.d0*atan(1.d0))
        call RANDOM_NUMBER(chi1)
        call RANDOM_NUMBER(chi2)
        val=var*dsqrt(-2.d0*dlog(1.d0-chi1))*dcos(2.d0*pi*chi2)+mean
    return
    end

end module integrate
