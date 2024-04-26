module integrate

    use mpi
    use forces, only: find_force_LJ,force_Verlet
    use MOD_INIT, only: PBC
    implicit none

    contains
     ! VELOCITY VERLET INTEGRATOR FOR 1 TIMESTEP !
     subroutine time_step_vVerlet(nprocs,pos,N,d,L,vel,dt,cutoff,nu,sigma,Upot,force,pos_to_transfer,start_atom,end_atom,displs,list,nlist)
        ! nprocs: number of processors
        ! pos: INPUT & OUTPUT, positions array 
        ! N: total number of atoms
        ! d: number of spatial dimensions
        ! L: size of one side of the simulation box
        ! vel:  INPUT & OUTPUT, velocities array
        ! dt: size of timestep
        ! cutoff: maximum particle distance where there is interaction
        ! nu: Andersen thermostat's associated probability
        ! sigma: deviation of the temperature
        ! Upot: OUTPUT, potential energy
        ! force: OUTPUT, array of forces
        ! pos_to_transfer: array with number of data to transfer/receive from each processors in ALLGATHER
        ! start_atom: first atom index to consider for this processor
        ! end_atom: last atom index to consider for this processor
        ! displs: displacement matrix for ALLGATHER
        implicit none
        integer, intent(in) :: nprocs
        integer :: Nat, N,d,i,rank,start_atom,end_atom,k,ierror, list(N),nlist(N,N)
        integer, dimension(nprocs) :: displs, pos_to_transfer
        real(8) :: pos(N,d),force(N,d),dt,L,cutoff,vel(N,d),nu,sigma,Upot,dr(3)
        real(8),allocatable :: local_velx(:),local_vely(:),local_velz(:),local_posx(:),local_posy(:),local_posz(:)
       
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

        ! Calculate forces between particles
        !call find_force_LJ(nprocs,pos,N,d,L,force,cutoff,Upot,pos_to_transfer,start_atom,end_atom,displs)
        call force_Verlet(nprocs,N,d,L,pos,force,cutoff,Upot,list,nlist,pos_to_transfer,start_atom,end_atom,displs)
        
        ! Calculate positions and velocities for each of the assigned atoms
        do i=start_atom,end_atom
            pos(i,1)=pos(i,1)+vel(i,1)*dt+0.5d0*force(i,1)*dt*dt
            pos(i,2)=pos(i,2)+vel(i,2)*dt+0.5d0*force(i,2)*dt*dt
            pos(i,3)=pos(i,3)+vel(i,3)*dt+0.5d0*force(i,3)*dt*dt
            vel(i,1)=vel(i,1)+force(i,1)*0.5d0*dt
            vel(i,2)=vel(i,2)+force(i,2)*0.5d0*dt
            vel(i,3)=vel(i,3)+force(i,3)*0.5d0*dt
        enddo

        ! Distribute velocities between processors
        call MPI_Allgatherv(vel(start_atom:end_atom,1), pos_to_transfer(rank+1), MPI_REAL8, vel(:,1), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(vel(start_atom:end_atom,2), pos_to_transfer(rank+1), MPI_REAL8, vel(:,2), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(vel(start_atom:end_atom,3), pos_to_transfer(rank+1), MPI_REAL8, vel(:,3), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)

        ! Apply PBC to atomic positions before distributing
        do i=start_atom,end_atom
            dr(1)=pos(i,1); dr(2)=pos(i,2);dr(3)=pos(i,3)
            call pbc(N,L,dr)
            pos(i,1)=dr(1);pos(i,2)=dr(2);pos(i,3)=dr(3)
        enddo

        ! Distribute positions between processors
        call MPI_Allgatherv(pos(start_atom:end_atom,1), pos_to_transfer(rank+1), MPI_REAL8, pos(:,1), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(pos(start_atom:end_atom,2), pos_to_transfer(rank+1), MPI_REAL8, pos(:,2), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(pos(start_atom:end_atom,3), pos_to_transfer(rank+1), MPI_REAL8, pos(:,3), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)


        ! Calculate forces and second part of velocities
        !call find_force_LJ(nprocs,pos,N,d,L,force,cutoff,Upot,pos_to_transfer,start_atom,end_atom,displs)
        call force_Verlet(nprocs,N,d,L,pos,force,cutoff,Upot,list,nlist,pos_to_transfer,start_atom,end_atom,displs)
        do i=start_atom,end_atom
            vel(i,1)=vel(i,1)+force(i,1)*0.5d0*dt
            vel(i,2)=vel(i,2)+force(i,2)*0.5d0*dt
            vel(i,3)=vel(i,3)+force(i,3)*0.5d0*dt
        enddo
        
        ! Distribute final velocities between processors
        call MPI_Allgatherv(vel(start_atom:end_atom,1), pos_to_transfer(rank+1), MPI_REAL8, vel(:,1), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(vel(start_atom:end_atom,2), pos_to_transfer(rank+1), MPI_REAL8, vel(:,2), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
        call MPI_Allgatherv(vel(start_atom:end_atom,3), pos_to_transfer(rank+1), MPI_REAL8, vel(:,3), pos_to_transfer, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)

        ! Apply Andersen thermostate
        call andersen_thermo(vel,N,d,nu,sigma)
    return
    end

    ! APPLY ANDERSEN THERMOSTATE !
    subroutine andersen_thermo(vel,N,d,nu,sigma)
        ! vel: INPUT & OUTPUT, velocity of particles
        ! N: number of particles
        ! d: number of spatial dimentions
        ! nu: probability of changing velocity
        ! sigma: deviation of temperature
        implicit none
        integer :: N,d,i
        real(8) :: sigma,vel(N,d),nu, rand
        
        ! With a probability nu, change the velocity of that particle 
        ! so it fits a normal distribution with variance equal to the temperature
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

    ! GENERATE A NUMBER WITHIN A GAUSSIAN DISTRIBUTION !
    subroutine gauss(mean,var,val)
        ! mean: mean value of the distribution
        ! var: standard deviation of the distribution
        ! val: OUTPUT, sampled value from the gaussian distribution
        implicit none
        real(8) :: chi1,chi2,pi,var,mean,val
        parameter(pi=4.d0*atan(1.d0))
        call RANDOM_NUMBER(chi1)
        call RANDOM_NUMBER(chi2)
                
        ! Sample random value from gaussian distribution using Box-Muller method
        val=var*dsqrt(-2.d0*dlog(1.d0-chi1))*dcos(2.d0*pi*chi2)+mean
    return
    end

end module integrate
