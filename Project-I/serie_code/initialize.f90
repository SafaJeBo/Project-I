! this module includes a initialize and periodic boundary conditions subroutine.

module initialize
    contains
    subroutine ini_pos_sc(N,a,M,d,pos)
        implicit none
        integer :: N,M,i,j,k,d,num_part
        real(8) :: pos(N,d),a
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
    subroutine pbc(xold,L,x)
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
