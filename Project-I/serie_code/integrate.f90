! this module includes the Velocity Verlet subroutine.
! it requires force and pbc modules for full functionality.
module integrate
    use forces, only: find_force_LJ
    use initialize, only: pbc
    contains
    subroutine time_step_vVerlet(pos,N,d,L,vel,dt,cutoff,nu,sigma,Upot)
        implicit none
        integer :: N,d,i
        real(8) :: pos(N,d),force(N,d),dt,L,cutoff,vel(N,d),nu,sigma,Upot
    
        call find_force_LJ(pos,N,d,L,force,cutoff,Upot)
    
        pos=pos+vel*dt+0.5d0*force*dt*dt
        vel=vel+force*0.5d0*dt
    
        do i=1,N
            call pbc(pos(i,1),L,pos(i,1))
            call pbc(pos(i,2),L,pos(i,2))
            call pbc(pos(i,3),L,pos(i,3))
        enddo
        call find_force_LJ(pos,N,d,L,force,cutoff,Upot)
        vel=vel+force*0.5d0*dt
        call andersen_thermo(vel,N,d,nu,sigma)
    return
    end
    
    subroutine andersen_thermo(vel,N,d,nu,sigma)
        implicit none
        integer :: N,d,i
        real(8) :: sigma,vel(N,d),nu
    
        do i=1,N
            if (rand().lt.nu) then
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
        chi1=rand()
        chi2=rand()
        val=var*dsqrt(-2.d0*dlog(1.d0-chi1))*dcos(2.d0*pi*chi2)+mean
    return
    end
end module integrate
