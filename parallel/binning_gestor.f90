module binning_gestor
    implicit none
        contains
        subroutine binning(data_arr, num, output_file)
             implicit none 
             character(12), intent(in) :: output_file
             integer, intent(in) :: num ! number of elements in data_arr
             double precision, dimension(:), intent(in) :: data_arr
             integer(8) :: ii, mm, max_m, block_length, block_num
             double precision :: mean, sigma, variance, mean_median, std_median
             double precision, dimension(:), allocatable :: means_arr, mean_results, std_results, my_data
             !double precision :: med_mean, med_std
             integer :: file_status, ierror, comm, rank, numproc
             integer, parameter :: MASTER = 0
             include 'mpif.h'
             double precision :: med_mean, med_std
             integer :: status(MPI_STATUS_SIZE)

             ! MPI initialization
             call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
             call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, ierror)

             ! Calculate max binning length
             allocate(mean_results(numproc))
             allocate(std_results(numproc))
             
             block_length = INT(num/(numproc))
             allocate(my_data(block_length))

             call MPI_Scatter(data_arr, block_length, MPI_DOUBLE_PRECISION, my_data,&
                block_length, MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD)

             if (rank /= MASTER) then
             !do i = numproc, 1 ! Iterate over block length
                  block_length = INT(num/(rank))
                  block_num = 0
                  allocate(means_arr(ceiling(num/dble(block_length))))
                  ! For each block length
                  ! Store mean of each sub block
                  do ii = 1, num, block_length
                      block_num = block_num + 1
                      if ((ii+block_length-1).gt.num) then
                          means_arr(block_num) = sum(my_data(ii:num))/size(my_data(ii:num))
                      else
                          means_arr(block_num) = sum(my_data(ii:ii+block_length-1))/dble(block_length)
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

                  deallocate(means_arr)
             end if

             print *, 'Mean:', med_mean
             print *, 'Standard deviation:', med_std

             call MPI_Gather(mean, 1, MPI_DOUBLE_PRECISION, mean_results, 1, MPI_DOUBLE_PRECISION, &
                MASTER, MPI_COMM_WORLD, ierror)

        call MPI_Gather(sigma, 1, MPI_DOUBLE_PRECISION, std_results, 1, MPI_DOUBLE_PRECISION, &
                MASTER, MPI_COMM_WORLD, ierror)
        
        ! Calculate median
        mean_median = calculate_median(mean_results, numproc)
        std_median = calculate_median(std_results, numproc)

        ! Write file with results
        if (rank == MASTER) then
             open(2, file=output_file, status='replace', action='write', iostat=file_status)
                 if (file_status /= 0) then
                      print*, "Error opening file for writing"
                      call MPI_Finalize(ierror)
                      stop
                 endif
        
                 write(2,*) "The mean value is:",mean_median , "and standard  deviation is:",std_median
             close(2)
        endif
    return
    end 

    real function calculate_median(arr, n)
        integer(4), intent(in) :: n
        double precision, dimension(:), intent(in) :: arr
        real(kind=kind(0.0d0)) :: median
        integer :: mid_index
        
        mid_index = n / 2
        
        if (mod(n, 2) == 0) then
            median = real(arr(mid_index) + arr(mid_index + 1), kind=kind(0.0d0)) / 2.0
        else
            median = real(arr(mid_index + 1), kind=kind(0.0d0))
        end if
        
        calculate_median = median
        
    return    
    end 

end module binning_gestor
