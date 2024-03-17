program binning_gestor
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer(8) :: dat_num, ii, mcs, num_lines, iunit
    integer :: file_status
    double precision, dimension(:), allocatable :: times, ekin, epot, etot, temp, msdval, press
    double precision :: time
    !!real(dp) :: ek, ep, et, t, mm, p, total_energy, total_magnetization
    character(len=100) :: filename

    ! Ask user for the file name 
    print*, "Enter the name of the file:"
    read(*,*) filename

    ! Open the file
    open(unit=10, file=filename, status='old', action='read', iostat=iunit)
    !Check the file opened successfully
    if (iunit /= 0) then
        print*, "Error opening file"
        stop
    endif

    ! Read the data and count the number of lines 
    num_lines = 0
    do
        read(10, *, iostat=iunit) time
        if (iunit /= 0) exit
        num_lines = num_lines + 1
    enddo

    ! Close the file
    close(10)
    print*, "Number of lines:", num_lines

    ! read data
    allocate(times(num_lines))
    allocate(ekin(num_lines))
    allocate(epot(num_lines))
    allocate(etot(num_lines))
    allocate(temp(num_lines))
    allocate(msdval(num_lines))
    allocate(press(num_lines))

    ! Call binning subroutine
    call binning(ekin, num_lines, "ekin_mean.dat", file_status)
    call binning(epot, num_lines, "epot_mean.dat", file_status)
    call binning(etot, num_lines, "etot_mean.dat", file_status)
    call binning(temp, num_lines, "temp_mean.dat", file_status)
    call binning(msdval, num_lines, "msdval_mean.dat", file_status)
    call binning(press, num_lines, "press_mean.dat", file_status)

contains

    subroutine binning(data_arr, num, file_name, file_status)
        implicit none
        character(len=*), intent(in) :: file_name
        integer(8), intent(in) :: num
        integer, intent(out) :: file_status
        double precision, dimension(:), intent(in) :: data_arr

        ! Binning code
        integer(8) :: ii, mm, max_m, block_length, block_num
        double precision :: mean, sigma, variance
        double precision, dimension(:,:), allocatable :: binning_mat
        double precision, dimension(:), allocatable :: means_arr
        double precision :: med_mean, med_std

        ! Calculate max binning length
        max_m = int(log(dble(num)) / log(2.0d0)) ! Max block length
        allocate(binning_mat(max_m+1, 3))

        do mm = 0, max_m ! Iterate over block length
             block_length = 2**mm
             block_num = 0
             allocate(means_arr(ceiling(num/dble(block_length))))
             ! For each block length
             ! Store mean of each sub block
             do ii = 1, num, block_length
                 block_num = block_num + 1
                 if ((ii+block_length-1).gt.num) then
                     means_arr(block_num) = sum(data_arr(ii:num))/size(data_arr(ii:num))
                 else
                     means_arr(block_num) = sum(data_arr(ii:ii+block_length-1))/dble(block_length)
                 end if
             end do
             ! Calculate mean of means
             mean = sum(means_arr)/dble(block_num)

             ! Calculate standard deviation and correct it
             variance = 0.
             do ii = 1, block_num
                 variance = variance + (means_arr(ii)-mean)**2
             end do
             sigma = sqrt(variance)/sqrt(dble(block_num)*dble(block_num-1))

             ! Store in binning_mat
             binning_mat(mm+1,1) = block_length
             binning_mat(mm+1,2) = mean
             binning_mat(mm+1,3) = sigma
             deallocate(means_arr)
        end do

        ! Print middle value for mean and std (apaño)
        med_mean = binning_mat(block_num/2,2)
        med_std = binning_mat(block_num/2,3)

        ! Write a file with results
        open(1, file=file_name, status='replace', iostat=file_status)
        if (file_status /= 0) then
             print*, "Error opening file for writing"
             return
        endif

        write(1,'(A, E20.10, A, E20.10)') "La mitjana estadistica es:", med_mean, &
                                " i la desviació estandar es:", med_std

        close(1)
    end subroutine binning

end program binning_gestor


