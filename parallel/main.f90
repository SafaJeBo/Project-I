! main program
program main
    ! Using the provided modules
    use mpi

    use MOD_INIT
    ! use [integration]
    ! use [force]
    ! use [thermodynamics]

    ! MPI variables
    integer :: ierror, rank, nprocs

    ! variable declaration
    integer,parameter :: d=3
    integer :: N,nsim_temp,nsim_tot,numdr
    integer :: M,i,j
    integer, allocatable :: atoms_list(:)
    real(8) :: density,L,a,dt,cutoff,temp1,temp2,nu,sigma,temperatura,ke,pot
    real(8) :: timeini,timefin,msdval,r,deltar,volumdr,pi,press
    real(8), allocatable :: pos(:,:), vel(:,:), pos0(:,:), rdf(:)
    character(15) :: dummy

    pi=4d0*datan(1d0)
    
    ! ------------------------------------------------------------------ !
    ! Initialize MPI
    call MPI_INIT(ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)

    print *, "Hi! I am", rank
    if (rank.eq.0) then
        ! Read parameters (only master)
        read(*,*) dummy, N
        read(*,*) dummy, nsim_temp
        read(*,*) dummy, nsim_tot
        read(*,*) dummy, numdr
        read(*,*) dummy, dt
        read(*,*) dummy, cutoff
        read(*,*) dummy, density
        read(*,*) dummy, temp1
        read(*,*) dummy, temp2
        read(*,*) dummy, nu
        print *, "I am ", rank, "I finished reading inputs. Going to broadcast"
    end if

    ! Broadcast values
    call MPI_BCAST(N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(nsim_temp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(nsim_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(numdr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(dt, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(cutoff, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(density, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(temp1, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(temp2, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(nu, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierror)

    ! Allocate variables
    allocate(pos(N,d))
    allocate(vel(N,d))
    allocate(pos0(N,d))
    allocate(rdf(numdr))

    ! ! Opening files to save results
    ! open(14,file='trajectory.xyz')
    ! open(15,file='thermodynamics.dat')
    ! open(16,file='resultsrdflong_def.dat')

    call MPI_BARRIER(MPI_COMM_WORLD,ierror) ! Barrier to start program at the same time

    ! Getting the initial time to account for total simulation time
    if (rank.eq.0) then 
        call cpu_time(timeini)
        write(*,*)timeini
    end if
    
    call srand(1)

    ! Initialize system
    L=(real(N,8)/density)**(1.d0/3.d0)
    call do_SCC(N, L, pos, atoms_list ,nprocs, rank, "SCCconf_init.xyz")

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    if (rank.eq.0) then 
        print *, rank
        print *, atoms_list
    end if 

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    if (rank.eq.1) then 
        print *, rank
        print *, atoms_list
    end if 

    ! HERE GOES SEPARATION OF PARTICLES BY PROC
    
    ! !Thermalization
    ! sigma = sqrt(temp1)
    ! do i=1,nsim_temp
    !     call time_step_vVerlet(pos,N,d,L,vel,dt,cutoff,nu,sigma,pot)
    !     write(14,'(A)') trim(adjustl(''))  ! Blank line
    !     write(14,'(I0)') N  
    !     write(14,'(A)') trim(adjustl(''))  ! Blank line
 	!     write(14,*) pos
    ! enddo
    ! print *, "Finished Thermalization"    


    ! !Starting production run
    ! sigma=sqrt(temp2)
    ! pos0=pos
    ! rdf=0d0
    
    ! do i=1,nsim_tot
    !     call time_step_vVerlet(pos,N,d,L,vel,dt,cutoff,nu,sigma,pot)
    !     if (mod(i,100).eq.0) then
    !         ! Mesure energy, pressure and msd
    !         call kin_energy(vel,N,d,ke)
    !         call msd(pos,N,d,pos0,L,msdval)
    !         call pression(pos,N,d,L,cutoff,press)
    !         temperatura=temp_inst(ke,N)
    !         write(15,*)i*dt,ke,pot,pot+ke,temperatura,msdval,press+temperatura*density

    !         if (mod(i,50000).eq.0) then ! Control state of simulation
    !             print*,i
    !         endif

    !         ! Mesure RDF after a certain timestep
    !         if (i.gt.1e3) then
    !             call gr(pos,N,d,numdr,L,rdf)
    !         endif
    !     endif
    !     AQUI PRBLY TOCA POSAR UNA BARRERA
    ! enddo
    
    ! ! Normalization of RDF
    ! r=0d0
    ! deltar=0.5d0*L/numdr
    ! do i=1,numdr
    !     r=(i-1)*deltar
    !     volumdr=4d0*pi*((r+deltar/2d0)**3-(r-deltar/2d0)**3)/3d0
    !     write(16,*)r,rdf(i)/(sum(rdf)*volumdr*density)
    ! enddo

    call MPI_BARRIER(MPI_COMM_WORLD,ierror) ! Final barrier to get time

    if (rank.eq.0) then
        call cpu_time(timefin)
        print*,'FINAL time = ',timefin-timeini
    end if
    ! Finalize MPI
    call MPI_FINALIZE(ierror)
    
    end program main
    
