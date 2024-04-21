implicit none
integer,parameter :: d=3,N=125,nsim_temp=5000,nsim_tot=50000,numdr=1000
integer :: M,i
real(8) :: density,L,a,pos(N,d),vel(N,d),dt,cutoff,temp,nu,sigma,pot
real(8) :: timeini,timefin,rdf(numdr),r,deltar,volumdr,pi,ke,pos0(N,d)
real(8) :: msdval,press,temperatura,temp_inst
pi=4d0*datan(1d0)

open(15,file='thermodynamics_test.dat')
open(16,file='resultsrdflong_def_test.dat')

call cpu_time(timeini)
dt=0.0005d0
cutoff=2.5d0
density=0.8d0

L=(dble(N)/density)**(1.d0/3.d0)
M=int(N**(1.d0/3.d0))+1
a=L/dble(M)

temp=300d0
sigma=sqrt(temp)
nu=0.1d0

call ini_pos_sc(N,a,M,d,pos)

do i=1,nsim_temp
    call time_step_vVerlet(pos,N,d,L,vel,dt,cutoff,nu,sigma,pot)
enddo

temp=2d0
sigma=sqrt(temp)
nu=0.1d0
vel=0d0
rdf=0d0
pos0=pos

do i=1,nsim_tot
    call time_step_vVerlet(pos,N,d,L,vel,dt,cutoff,nu,sigma,pot)
    if (mod(i,100).eq.0) then
        !Mesures
        call kin_energy(vel,N,d,ke)
        call msd(pos,N,d,pos0,L,msdval)
        call pression(pos,N,d,L,cutoff,press)
        temperatura=temp_inst(ke,N)
        write(15,*)i*dt,ke,pot,pot+ke,temperatura,msdval,press
    endif
    if (mod(i,10000).eq.0) then
        print*,i
    endif
    if (i.gt.1e4) then
        call gr(pos,N,d,numdr,L,rdf)
    endif
enddo

r=0d0
deltar=1.5d0*L/numdr
do i=1,numdr
    r=(2*i-1)*deltar
    volumdr=4d0*pi*((r+deltar/2d0)**3-(r-deltar/2d0)**3)/3d0
    write(16,*)r,rdf(i)/(sum(rdf)*volumdr*density)
enddo

call cpu_time(timefin)
print*,'FINAL cputime = ',timefin-timeini

end


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

subroutine gr(pos,N,d,numdr,L,rdf)
    implicit none
    integer :: N,d,i,j,numdr,x,y,z
    real(8) :: pos(N,d),deltar,L,rdf(numdr),dx,dy,dz,limit,dxnew,dynew,dznew,drnew
    limit=1.5d0*L
    deltar=limit/numdr
    do i=1,N
        do j=i+1,N
            dx=pos(i,1)-pos(j,1)
            call pbc(dx,L,dx)
            dy=pos(i,2)-pos(j,2)
            call pbc(dy,L,dy)
            dz=pos(i,3)-pos(j,3)
            call pbc(dz,L,dz)
            do x=-1,1,1
                do y=-1,1,1
                    do z=-1,1,1
                        dxnew=dx+x*L; dynew=dy+y*L; dznew=dz+z*L
                        drnew=dsqrt(dxnew**2+dynew**2+dznew**2)
                        if (drnew.le.limit) then
                            rdf(int(drnew/deltar)+1)=rdf(int(drnew/deltar)+1)+1d0
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo
return
end

