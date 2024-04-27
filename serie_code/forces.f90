! Module used to calculate the forces acting on a particle.
! requires module pbc for full functionality.

module forces
    use initialize, only: pbc
    contains

    ! CALCULATE THE FORCE BETWEEN ALL PAIRS OF PARTICLES FOLLOWING A LENNARD-JONES POTENTIAL !
    subroutine find_force_LJ(pos,N,d,L,force,cutoff,Upot)
        ! pos: position array
        ! N: number of particles
        ! d: number of spatial dimensions
        ! L: size of one side of the simulation box
        ! force: OUTPUT, array of forces that applies that particle
        ! cutoff: maximum particle distance where there is interaction
        ! Upot: OUTPUT, potential energy of the system
        implicit none
        integer :: N,d,i,j
        real(8) :: pos(N,d),force(N,d),dx,dy,dz,dr2,cutoff,L,cf2,Upot,dr6,dr12,potcut,fij
        force=0d0
        Upot=0d0
        cf2=cutoff*cutoff
        potcut=4.d0*(1.d0/cf2**6-1.d0/cf2**3) 

        do i=1,N
            do j=i+1,N
                ! Calculate distance
                dx=pos(i,1)-pos(j,1)
                call pbc(dx,L,dx)
                dy=pos(i,2)-pos(j,2)
                call pbc(dy,L,dy)
                dz=pos(i,3)-pos(j,3)
                call pbc(dz,L,dz)
                dr2=dx*dx+dy*dy+dz*dz

                ! Determine if particles are within interactive radius
                if (dr2.lt.cf2) then
                    dr6=dr2*dr2*dr2
                    dr12=dr6*dr6
                    ! Calculate force between pair of particles
                    fij=48.d0/(dr6*dr6*dr2)-24.d0/(dr2*dr6)

                    ! Add that force to each spatial component
                    force(i,1)=force(i,1)+fij*dx
                    force(i,2)=force(i,2)+fij*dy
                    force(i,3)=force(i,3)+fij*dz
                    force(j,1)=force(j,1)-fij*dx
                    force(j,2)=force(j,2)-fij*dy
                    force(j,3)=force(j,3)-fij*dz

                    ! Calculate potential energy between particles
                    Upot=Upot+4.d0*(1.d0/dr12-1.d0/dr6)-potcut
                endif
            enddo
        enddo
    return
    end
    
    ! CALCULATE VERLET LIST !
    subroutine new_vlist(N,d,L,pos,list,nlist,cutoff)
        ! N: number of particles
        ! d: number of spatial dimensions
        ! L: size of one side of the simulation box
        ! pos: position array
        ! list: OUTPUT, Verlet list. Array with the indexes of the particles within cutoff radius of another
        ! nlist: OUTPUT, array containing the number of particles within cutoff radius of another
        ! cutoff: maximum particle distance where there is interaction
        implicit none
        integer :: i,j,N,d,nlist(N),list(N,N)
        real(8) :: pos(N,d),L,dx,dy,dz,cutoff,cf2
        real(8) :: dr2
        nlist=0
        cf2=cutoff**2

        ! Calculate distances between particles
        do i=1,N
            do j=i+1,N
                if (i.ne.j) then
                dx=pos(i,1)-pos(j,1)
                call pbc(dx,L,dx)
                dy=pos(i,2)-pos(j,2)
                call pbc(dy,L,dy)
                dz=pos(i,3)-pos(j,3)
                call pbc(dz,L,dz)
                dr2=dx*dx+dy*dy+dz*dz
                ! Store the number of particles within cutoff radius of "i"
                ! and store also the index of those particles
                if (dr2.le.cf2) then
                    nlist(i)=nlist(i)+1
                    nlist(j)=nlist(j)+1
                    list(i,nlist(i))=j
                    list(j,nlist(j))=i
                endif
            endif
            enddo
        enddo
        return
        end
        
        ! CALCULATE THE FORCE BETWEEN PAIRS OF CLOSE (USING VERLET LISTS) PARTICLES FOLLOWING A LENNARD-JONES POTENTIAL !
        subroutine force_Verlet(N,d,L,pos,force,cutoff,Upot,list,nlist)
            ! N: number of particles
            ! d: number of spatial dimensions
            ! L: size of one side of the simulation box
            ! pos: position array
            ! force: OUTPUT, array of forces that applies that particle
            ! cutoff: maximum particle distance where there is interaction
            ! Upot: OUTPUT, potential energy of the system
            ! list: Verlet list. Array with the indexes of the particles within cutoff radius of another
            ! nlist: array containing the number of particles within cutoff radius of another
            implicit none
            integer :: N,d,i,nlist(N),jj,j,list(N,N)
            real(8) :: dy,fij,L,dz,cutoff,cf2,dr2,dr6,dr12,Upot
            real(8) :: force(N,d),dx,pos(N,d),potcut
            
            force=0d0
            Upot=0d0
            cf2=cutoff*cutoff
            potcut=4.d0*(1.d0/cf2**6-1.d0/cf2**3) 

            ! Calculate force for particles in Verlet list
            do i=1,N
                do jj=1,nlist(i)
                    j=list(i,jj) ! Get particle index
                    if (j.gt.i) then ! Avoid repetition of pairs
                        ! Calculate distance for pair of particles within 
                        dx=pos(i,1)-pos(j,1)
                        call pbc(dx,L,dx)
                        dy=pos(i,2)-pos(j,2)
                        call pbc(dy,L,dy)
                        dz=pos(i,3)-pos(j,3)
                        call pbc(dz,L,dz)
                        dr2=dx*dx+dy*dy+dz*dz
        
                        ! Check if particles are still within interactive radius
                        if (dr2.lt.cf2) then
                            dr6=dr2*dr2*dr2
                            dr12=dr6*dr6
                            ! Calculate force between pair of particles
                            fij=48.d0/(dr6*dr6*dr2)-24.d0/(dr2*dr6)
        
                            ! Add that force to each spatial component
                            force(i,1)=force(i,1)+fij*dx
                            force(i,2)=force(i,2)+fij*dy
                            force(i,3)=force(i,3)+fij*dz
                            force(j,1)=force(j,1)-fij*dx
                            force(j,2)=force(j,2)-fij*dy
                            force(j,3)=force(j,3)-fij*dz
        
                            ! Calculate potential energy between particles
                            Upot=Upot+4.d0*(1.d0/dr12-1.d0/dr6)-potcut
                        endif
                    endif
                enddo
            enddo
            return
            end
end module forces
