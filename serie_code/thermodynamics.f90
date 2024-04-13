! this module includes different observables: Instantaneous temperature, kinetic energy, pressure, msd and g(r).

module thermodynamics
    use initialize, only: pbc
    contains
    real(8) function temp_inst(ke,N)
        implicit none
        integer :: N
        real(8) :: ke,kb
        kb=1d0
        temp_inst=2d0*ke/((3d0*N-3d0)*kb)
    return
    end

    subroutine kin_energy(vel,N,d,ene)
        implicit none
        integer :: N,d,i
        real(8) :: vel(N,d),ene,vel_mod
        ene=0.d0
        do i=1,N
            vel_mod=vel(i,1)**2+vel(i,2)**2+vel(i,3)**2
            ene=ene+0.5d0*vel_mod
        enddo
    return
    end

    subroutine pression(pos,N,d,L,cutoff,press)
        implicit none
        integer :: N,d,i,j
        real(8) :: pos(N,d),L,press,cutoff,dx,dy,dz,cf2,dr2,dr6,fij
        press=0d0
        cf2=cutoff*cutoff
        do i=1,N
            do j=i+1,N
                dx=pos(i,1)-pos(j,1)
                call pbc(dx,L,dx)
                dy=pos(i,2)-pos(j,2)
                call pbc(dy,L,dy)
                dz=pos(i,3)-pos(j,3)
                call pbc(dz,L,dz)
                dr2=dx*dx+dy*dy+dz*dz
                if (dr2.lt.cf2) then
                    dr6=dr2*dr2*dr2
                    fij=48.d0/(dr6*dr6*dr2)-24.d0/(dr2*dr6)
                    press=press+fij*dsqrt(dr2)
                endif
            enddo
        enddo
        press=press/(3d0*L**3)
    return
    end

    subroutine msd(pos,N,d,pos0,L,val)
        implicit none
        integer :: N,d,i
        real(8) :: pos(N,d),pos0(N,d),val,dx,dy,dz,L
        val=0d0
        do i=1,N
            dx=pos(i,1)-pos0(i,1)
            call pbc(dx,L,dx)
            dy=pos(i,2)-pos0(i,2)
            call pbc(dy,L,dy)
            dz=pos(i,3)-pos0(i,3)
            call pbc(dz,L,dz)
            val=val+dx*dx+dy*dy+dz*dz
        enddo
        val=val/N
    return
    end

    ! ######## CALCULATE g(r) ######## !
    subroutine gr(pos,N,d,numdr,L,rdf)       
        implicit none
        integer :: N,d,i,j,numdr,x,y,z
        real(8) :: pos(N,d),deltar,L,rdf(numdr),dx,dy,dz,dr,limit
        
        ! Set the limits of the g(r) to be half the simulation box and calculate deltar 
        limit=0.5d0*L
        deltar=limit/numdr

        ! Calculate distance between each pair of particles
        do i=1,N
            do j=i+1,N
                dx=pos(i,1)-pos(j,1)
                call pbc(dx,L,dx)
                dy=pos(i,2)-pos(j,2)
                call pbc(dy,L,dy)
                dz=pos(i,3)-pos(j,3)
                call pbc(dz,L,dz)
                
                dr=dsqrt(dx**2+dy**2+dz**2)
                
                ! Classify distance by how many deltar it comprises
                if (dr.le.limit) then
                    rdf(int(dr/deltar)+1)=rdf(int(dr/deltar)+1)+1d0
                endif
             enddo
        enddo
    return
    end
end module thermodynamics
