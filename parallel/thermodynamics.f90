! This module includes different observables: Instantaneous temperature, kinetic energy, pressure, msd and g(r).
module thermodynamics
    use mpi
    use MOD_INIT, only: pbc
    contains

    ! CALCULATE INSTANTANEOUS TEMPERATURE !
    real(8) function temp_inst(ke,N)
        ! ke: kinetic energy of the timestep
        ! N: number of particles
        implicit none
        integer :: N
        real(8) :: ke,kb
        kb=1d0 ! Boltzmann constant in reduced units
        temp_inst=2d0*ke/((3d0*N-3d0)*kb)
    return
    end

    ! CALCULATE KINETIC ENERGY !
    subroutine kin_energy(vel,N,d,start_atom,end_atom,ene)
        ! vel: array of velocities for each particle and dimension
        ! N: number of particles
        ! d: nummber of spatial dimensions
        ! start_atom: first atom index to consider for this processor
        ! end_atom: last atom index to consider for this processor
        ! ene: OUTPUT, kinetic energy
        implicit none
        integer :: N,d,i,ierror,start_atom,end_atom
        real(8) :: vel(N,d),ene,vel_mod,local_ene
        ene=0d0
        local_ene=0d0
        ! Calculate Kinetic energy for 
        do i=start_atom,end_atom
            vel_mod=vel(i,1)**2+vel(i,2)**2+vel(i,3)**2
            local_ene=local_ene+0.5d0*vel_mod
        enddo
        call MPI_Allreduce(local_ene,ene,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierror)
    return
    end

    ! CALCULATE PRESSURE (ONLY SUMMATION TERM) !
    subroutine pression(pos,N,d,L,cutoff,start_atom,end_atom,press)
        ! pos: position array
        ! N: number of particles
        ! d: number of spatial dimmensions
        ! L: size of one side of the simulation box
        ! cutoff: maximum particle distance where there is interaction
        ! start_atom: first atom index to consider for this processor
        ! end_atom: last atom index to consider for this processor
        ! press: OUTPUT, pressure of the system
        implicit none
        integer :: N,d,i,j,ierror,start_atom,end_atom
        real(8) :: pos(N,d),L,press,cutoff,dx,dy,dz,cf2,dr2,dr6,fij,local_press,dr(3)
        press=0d0
        local_press=0d0
        cf2=cutoff*cutoff
        ! Calculate distances between pairs of particles for current processor
        do i=start_atom,end_atom
            do j=1,N
                if (i.ne.j) then
                    dx=pos(i,1)-pos(j,1)
                    dy=pos(i,2)-pos(j,2)
                    dz=pos(i,3)-pos(j,3)
                    dr(1)=dx; dr(2)=dy; dr(3)=dz
                    call pbc(N,L,dr)
                    dr2=dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)

                    ! Calculate pressure for particles that interact
                    if (dr2.lt.cf2) then
                        dr6=dr2*dr2*dr2
                        fij=48.d0/(dr6*dr6*dr2)-24.d0/(dr2*dr6)
                        local_press=local_press+fij*dsqrt(dr2)
                    endif
                endif
            enddo
        enddo
        local_press=local_press/2d0
        call MPI_Allreduce(local_press,press,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierror)
        press=press/(3d0*L**3)
    return
    end

    ! CALCULATE MEAN SQUARED DISPLACEMENT ! 
    subroutine msd(pos,N,d,pos0,L,start_atom,end_atom,val)
        ! pos: postion array at time t_{i}
        ! N: number of particles
        ! pos0: position array at time t_{i-1}
        ! L: size of one side of the simulation box
        ! start_atom: first atom index to consider for this processor
        ! end_atom: last atom index to consider for this processor
        ! val: OUTPUT, mean squared displacement
        implicit none
        integer :: N,d,i,start_atom,end_atom,ierror
        real(8) :: pos(N,d),pos0(N,d),val,dx,dy,dz,L,dr(3),local_val
        local_val=0d0
        val=0d0
        
        ! Calculate squared displacement between 2 consecutive timesteps for current processor's particles
        do i=start_atom,end_atom
            dx=pos(i,1)-pos0(i,1)
            dy=pos(i,2)-pos0(i,2)
            dz=pos(i,3)-pos0(i,3)
            dr(1)=dx; dr(2)=dy; dr(3)=dz
            call pbc(N,L,dr)
            local_val=local_val+dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
        enddo
        call MPI_Allreduce(local_val,val,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierror)
        ! Average the displacement
        val=val/N
    return
    end

    ! CALCULATE PAIRWISE RADIAL DISTRIBUTION FUNCTION (RDF) !
    subroutine gr(pos,N,d,numdr,L,start_atom,end_atom,local_rdf)       
        ! pos: postion array
        ! N: number of particles
        ! d: number of spatial dimensions
        ! numdr: number of bins for rdf
        ! L: size of one side of the simulation box
        ! start_atom: first atom index to consider for this processor
        ! end_atom: last atom index to consider for this processor
        ! rdf: OUTPUT, array containing number of particles in each distance bin 
        implicit none
        integer :: N,d,i,j,numdr,x,y,z,start_atom,end_atom,ierror
        real(8) :: pos(N,d),deltar,L,local_rdf(numdr),dx,dy,dz,dr,limit,drnew(3)
        
        ! Set the limits of the g(r) to be half the simulation box and calculate deltar 
        limit=0.5d0*L
        deltar=limit/numdr

        ! Calculate distance between each pair of particles of current processor
        do i=start_atom,end_atom
            do j=1,N
                if (i.ne.j) then
                    dx=pos(i,1)-pos(j,1)
                    dy=pos(i,2)-pos(j,2)
                    dz=pos(i,3)-pos(j,3)
                    drnew(1)=dx; drnew(2)=dy; drnew(3)=dz
                    call pbc(N,L,drnew)
                    dr=dsqrt(drnew(1)*drnew(1)+drnew(2)*drnew(2)+drnew(3)*drnew(3))
                    ! Classify distance by how many deltar it comprises
                    if (dr.le.limit) then
                        local_rdf(int(dr/deltar)+1)=local_rdf(int(dr/deltar)+1)+1d0
                    endif
                endif
            enddo
        enddo
    return
    end
end module thermodynamics
