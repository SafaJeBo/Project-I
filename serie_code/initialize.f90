! This module includes a initialize and periodic boundary conditions subroutine.

module initialize
    contains

    ! INITIALIZE PARTICLES POSITIONS !
    subroutine ini_pos_sc(N,a,M,d,pos)
        ! N: number of particles
        ! a: separation of particles 
        ! M: number of unit cells in each dimension
        ! d: number of spatial dimensions
        ! pos: OUTPUT, array of particle positions
        implicit none
        integer :: N,M,i,j,k,d,num_part
        real(8) :: pos(N,d),a

        ! Put N particles, M in each dimension, spaced a distance "a" in each dimension
        pos=0
        num_part=1
        do i=0,M-1
            do j=0,M-1
                do k=0,M-1
                    pos(num_part,1)=i*a
                    pos(num_part,2)=j*a
                    pos(num_part,3)=k*a
                    num_part=num_part+1
                enddo
            enddo
        enddo
    return
    end

    ! APPLY PERIODIC BOUNDARY CONDITIONS !
    subroutine pbc(xold,L,x)
        ! xold: initial position
        ! L: size of one side of the simulation box, which is centered at 0
        ! x: OUTPUT, position after applying PBC
        implicit none
        real(8) :: x,L,xold
        if (xold.gt.(L/2.d0)) then
            x=xold-L
        else if (xold.lt.(-L/2.d0)) then
            x=xold+L
        else
            x=xold
        endif
    return
    end
end module initialize
