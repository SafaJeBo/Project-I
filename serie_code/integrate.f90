! This module includes the Velocity Verlet and Andersen thermostat subroutines.
! it requires force and pbc modules for full functionality.
module integrate
    use forces, only: find_force_LJ, force_Verlet,new_vlist
    use initialize, only: pbc
    contains

    ! VELOCITY VERLET INTEGRATOR FOR 1 TIMESTEP !
    subroutine time_step_vVerlet(pos,N,d,L,vel,dt,cutoff,list,nlist,nu,sigma,Upot)
        ! pos: INPUT & OUTPUT, position array
        ! N: number of particles
        ! d: number of spatial dimensions
        ! L: size of one side of the simulation box
        ! vel: INPUT & OUTPUT, velocity array
        ! dt: size of timestep
        ! cutoff: maximum particle distance where there is interaction
        ! list: Verlet list. Array with the indexes of the particles within cutoff radius of another 
        ! nlist: array containing the number of particles within cutoff radius of another
        ! nu: Andersen thermostat's associated probability
        ! sigma: deviation of the temperature
        ! Upot: OUTPUT, potential energy
        implicit none
        integer :: N,d,i,list(N,N),nlist(N)
        real(8) :: pos(N,d),force(N,d),dt,L,cutoff,vel(N,d),nu,sigma,Upot
    
        ! Calculate force between particles using Verlet lists
        !call find_force_LJ(pos,N,d,L,force,cutoff,Upot)
        call force_Verlet(N,d,L,pos,force,cutoff,Upot,list,nlist)
        
        ! Calculate postions and one part of velocities of the following timestep
        pos=pos+vel*dt+0.5d0*force*dt*dt
        vel=vel+force*0.5d0*dt
    
        ! Apply periodic boundary conditions
        do i=1,N
            call pbc(pos(i,1),L,pos(i,1))
            call pbc(pos(i,2),L,pos(i,2))
            call pbc(pos(i,3),L,pos(i,3))
        enddo
        ! Calculate the other part of the velocities
        !call find_force_LJ(pos,N,d,L,force,cutoff,Upot)
        call force_Verlet(N,d,L,pos,force,cutoff,Upot,list,nlist)
        vel=vel+force*0.5d0*dt

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
        real(8) :: sigma,vel(N,d),nu,rand
        
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
