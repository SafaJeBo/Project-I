! main program
program main
    ! using the provided modules
    use initialize
    use forces
    use integrate
    use thermodynamics
    use binning_gestor

    ! variable declaration
    implicit none
    integer,parameter :: d=3
    integer :: N,nsim_temp,nsim_tot,numdr,verlet_step
    integer :: M,i,j
    real(8) :: density,L,a,dt,cutoff,temp1,temp2,nu,sigma,temperatura,ke,pot
    real(8) :: timeini,timefin,msdval,r,deltar,volumdr,pi,press
    real(8), allocatable :: pos(:,:), vel(:,:), pos0(:,:), rdf(:)
    integer,allocatable :: list(:,:),nlist(:)
    integer :: size_seed, seed
    integer, allocatable :: seed2(:)
    character(15) :: dummy
    double precision , dimension(:), allocatable :: ekin_arr, epot_arr, etot_arr, temp_arr, msd_arr, press_arr

    
    ! Read parameters
    read(*,*) dummy, N
    read(*,*) dummy, nsim_temp
    read(*,*) dummy, nsim_tot
    read(*,*) dummy, verlet_step
    read(*,*) dummy, numdr
    read(*,*) dummy, dt
    read(*,*) dummy, cutoff
    read(*,*) dummy, density
    read(*,*) dummy, temp1
    read(*,*) dummy, temp2
    read(*,*) dummy, nu

    pi=4d0*datan(1d0)

    ! Allocate variables
    allocate(pos(N,d))
    allocate(vel(N,d))
    allocate(pos0(N,d))
    allocate(rdf(numdr))
    allocate(list(N,N),nlist(N))
    allocate(ekin_arr((nsim_tot-1000)/100), epot_arr((nsim_tot-1000)/100), etot_arr((nsim_tot-1000)/100))
    allocate(temp_arr((nsim_tot-1000)/100), msd_arr((nsim_tot-1000)/100), press_arr((nsim_tot-1000)/100))
    
    ! Opening files to save results
    open(14,file='trajectory.xyz')
    open(15,file='thermodynamics.dat')
    open(16,file='resultsrdflong_def.dat')


    ! Getting the initial time to account for total simulation time
    call cpu_time(timeini)
    write(*,*)timeini
    
    !call srand(1)

    !  Initialize random number generator according to system clock (different results each time)  !
    call random_seed(size=size_seed)
    allocate (seed2(size_seed))
    call system_clock(count=seed)
    seed2 = seed
    call random_seed(put=seed2)
    
    ! Initialize system
    L=(dble(N)/density)**(1.d0/3.d0)
    M=int(N**(1.d0/3.d0))+1
    a=L/dble(M)

    print *, L, cutoff, M, a
    
    call ini_pos_sc(N,a,M,d,pos)
    
    !Thermalization
    sigma = sqrt(temp1)
    call new_vlist(N,d,L,pos,list,nlist,cutoff)
    vel=0d0
    do i=1,nsim_temp
        if (mod(i,verlet_step).eq.0) then
            call new_vlist(N,d,L,pos,list,nlist,cutoff)
        endif
        call time_step_vVerlet(pos,N,d,L,vel,dt,cutoff,list,nlist,nu,sigma,pot)
    enddo
    print *, "Finished Thermalization"    


    !Starting production run
    sigma=sqrt(temp2)
    pos0=pos
    rdf=0d0
    vel=0d0

    call new_vlist(N,d,L,pos,list,nlist,cutoff)

    do i=1,nsim_tot
        if (mod(i,verlet_step).eq.0) then
            call new_vlist(N,d,L,pos,list,nlist,cutoff)
        endif
        call time_step_vVerlet(pos,N,d,L,vel,dt,cutoff,list,nlist,nu,sigma,pot)
        if (mod(i,100).eq.0) then
            ! Mesure energy, pressure and msd
            call kin_energy(vel,N,d,ke)
            call msd(pos,N,d,pos0,L,msdval)
            call pression(pos,N,d,L,cutoff,press)
            temperatura=temp_inst(ke,N)
            write(15,*)i*dt,ke,pot,pot+ke,temperatura,msdval,press+temperatura*density
            print*,i,ke,pot,pot+ke
            if (mod(i,50000).eq.0) then ! Control state of simulation
                print*,i
            endif

            ! Write trajectory in file "trajectory.xyz"
            write(14,'(I5)')N
            write(14, *) 
            do j = 1, 125
                write(14, '(A, 3F12.6)')'A', pos(j, 1), pos(j, 2), pos(j, 3)
            end do
            
            ! Mesure RDF after a certain timestep
            if (i.gt.1e3) then
                call gr(pos,N,d,numdr,L,rdf)
                ! Save results in arrays 
                ekin_arr((i-1000)/100) = dble(ke)
                epot_arr((i-1000)/100) = dble(pot)
                etot_arr((i-1000)/100) = dble(pot+ke)
                temp_arr((i-1000)/100) = dble(temperatura)
                msd_arr((i-1000)/100) = dble(msdval)
                press_arr((i-1000)/100) = dble(press+temperatura*density)
            endif
        endif
    enddo
    
    print *, "Gonna bin"
    call binning(ekin_arr, (nsim_tot-1000)/100, "Ekin_mean.dat")
    call binning(epot_arr, (nsim_tot-1000)/100, "Epot_mean.dat")
    call binning(etot_arr, (nsim_tot-1000)/100, "Etot_mean.dat")
    call binning(temp_arr, (nsim_tot-1000)/100, "Temp_mean.dat")
    call binning(msd_arr, (nsim_tot-1000)/100, "MSD_mean.dat")
    call binning(press_arr, (nsim_tot-1000)/100, "Press_mean.dat")

    ! Normalization of RDF
    r=0d0
    deltar=0.5d0*L/numdr
    do i=1,numdr
        r=(i-1)*deltar
        volumdr=4d0*pi*((r+deltar/2d0)**3-(r-deltar/2d0)**3)/3d0
        write(16,*)r,rdf(i)/(sum(rdf)*volumdr*density)
    enddo
    
    call cpu_time(timefin)
    print*,'FINAL time = ',timefin-timeini

    deallocate(seed2)
    end program main
    
