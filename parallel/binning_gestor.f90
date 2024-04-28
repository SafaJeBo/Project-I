module binning_gestor
    implicit none
    contains
    ! DO BINNING OF DATA !
    subroutine binning(data_arr, num, output_file, output_block)
        ! data_arr: array containing the data to bin
        ! num: size of data_arr
        ! file_name: name of the output file
        implicit none 
        character(len=*), intent(in) :: output_file, output_block
        integer, intent(in) :: num ! number of elements in data_arr
        double precision, dimension(:), intent(in) :: data_arr
        integer(8) :: ii,i, mm, max_m, block_length, block_num
        double precision :: mean, sigma, variance, mean_median, std_median
        double precision, dimension(:), allocatable :: means_arr, mean_results, std_results, block_results
        double precision, dimension(:,:), allocatable :: binning_mat
        integer, dimension(:), allocatable :: pos_to_transfer, displs
        integer :: file_status, ierror, comm, rank, nprocs, n_blocks_remaining, blocks_per_proc, end_block, start_block
        integer, parameter :: MASTER = 0
        include 'mpif.h'
        double precision :: med_mean, med_std
        integer :: status(MPI_STATUS_SIZE)

        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
        allocate(pos_to_transfer(nprocs))
        allocate(displs(nprocs))
            
        ! Calculate max binning length
        max_m = int(log(dble(num)) / log(2.0d0)) ! Max block length exponent
        allocate(binning_mat(max_m+1,3))
        allocate(mean_results(max_m+1))
        allocate(std_results(max_m+1))
        allocate(block_results(max_m+1))
            
        ! Distribute blocks between processors
        n_blocks_remaining = mod(max_m+1, nprocs)

        if (rank < n_blocks_remaining) then
            blocks_per_proc = (max_m+1) / nprocs + 1
            start_block = rank * blocks_per_proc
            end_block = start_block + blocks_per_proc - 1 
        else
            blocks_per_proc = (max_m+1) / nprocs
            start_block = n_blocks_remaining * (blocks_per_proc + 1) + (rank - n_blocks_remaining) * blocks_per_proc
            end_block = start_block + blocks_per_proc - 1
        end if

        ! Generate an array with all the number of positions that will be sent later
        call MPI_ALLGATHER(blocks_per_proc,1,MPI_INT,pos_to_transfer,1,MPI_INT,MPI_COMM_WORLD, ierror)

        displs(1) = 0
        do i = 2, nprocs
            displs(i) = displs(i-1)+pos_to_transfer(i-1)
        end do
        
        ! Start doing binning
        do mm = start_block, end_block ! Iterate over block length
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
            binning_mat(mm+1,1) = mean
            binning_mat(mm+1,2) = sigma
            binning_mat(mm+1,3) = block_length
            deallocate(means_arr)
        end do
        close(32)
        ! Share results to master
        call MPI_Gatherv(binning_mat(start_block+1:end_block+1,1), pos_to_transfer(rank+1), MPI_DOUBLE_PRECISION, &
        mean_results, pos_to_transfer,displs, MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD, ierror)
        call MPI_Gatherv(binning_mat(start_block+1:end_block+1,2), pos_to_transfer(rank+1), MPI_DOUBLE_PRECISION, &
        std_results, pos_to_transfer,displs, MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD, ierror)
        call MPI_Gatherv(binning_mat(start_block+1:end_block+1,3), pos_to_transfer(rank+1), MPI_DOUBLE_PRECISION, &
        block_results, pos_to_transfer,displs, MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD, ierror) 
        
        ! Calculate median
        mean_median = calculate_median(mean_results, max_m+1)
        std_median = calculate_median(std_results, max_m+1)

        ! Write file with results
        if (rank == MASTER) then
            open(2, file=output_file, status='replace', action='write', iostat=file_status)
                if (file_status /= 0) then
                    print*, "Error opening file for writing"
                    stop
                endif
                write(2,*) "The mean value is:",mean_median , "and standard  deviation is:",std_median
            close(2)    
            open(3, file=output_block, iostat=file_status)  
                if (file_status /= 0) then
                     print*, "Error opening file for writing"
                     stop  
                endif 
                write(3,*) "Block length                 Mean                  Average"
                do i = 1, max_m+1
                     write(3,*) block_results(i), mean_results(i), std_results(i)
                enddo
            close(3)
        endif
    return
    end 
    
    ! CALCULATE MEDIAN OF DATA !
    real function calculate_median(arr, n)
        ! arr: array with data to calculate median 
        ! n: length of array
        integer(8), intent(in) :: n
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
