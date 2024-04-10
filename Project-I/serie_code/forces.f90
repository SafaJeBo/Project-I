! module used to calculate the forces acting on a particle.
! requires module pbc for full functionality.

module forces
    use initialize, only: pbc
    contains
    subroutine find_force_LJ(pos,N,d,L,force,cutoff,Upot)
        implicit none
        integer :: N,d,i,j
        real(8) :: pos(N,d),force(N,d),dx,dy,dz,dr2,cutoff,L,cf2,Upot,dr6,dr12,potcut,fij
        force=0d0
        Upot=0d0
        cf2=cutoff*cutoff
        potcut=4.d0*(1.d0/cf2**6-1.d0/cf2**3)
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
                    dr12=dr6*dr6
                    fij=48.d0/(dr6*dr6*dr2)-24.d0/(dr2*dr6)
                    force(i,1)=force(i,1)+fij*dx
                    force(i,2)=force(i,2)+fij*dy
                    force(i,3)=force(i,3)+fij*dz
                    force(j,1)=force(j,1)-fij*dx
                    force(j,2)=force(j,2)-fij*dy
                    force(j,3)=force(j,3)-fij*dz
                    Upot=Upot+4.d0*(1.d0/dr12-1.d0/dr6)-potcut
                endif
            enddo
        enddo
    return
    end
end module forces
