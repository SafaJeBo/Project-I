module binning_gestor
implicit none
    contains
    ! DO BINNING OF DATA !
    subroutine binning(data_arr, num, file_name)
        ! data_arr: array containing the data to bin
        ! num: size of data_arr
        ! file_name: name of the output file
        implicit none
        character(len=*), intent(in) :: file_name
        integer, intent(in) :: num
        integer :: file_status
        double precision, dimension(:), intent(in) :: data_arr
        integer :: ii, mm, max_m, block_length, block_num
        double precision :: mean, sigma, variance
        double precision, dimension(:,:), allocatable :: binning_mat
        double precision, dimension(:), allocatable :: means_arr
        double precision :: med_mean, med_std

        ! Calculate max binning length
        max_m = int(log(dble(num)) / log(2.0d0)) ! Max block length
        allocate(binning_mat(max_m+1, 3))
        
        ! Start doing binning
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

        ! Calculate median value for mean and std
        med_mean = calculate_median(binning_mat(:,2),max_m+1)
        med_std = calculate_median(binning_mat(:,3),max_m+1)

        ! Write a file with results
        open(1, file=file_name, status='replace', iostat=file_status)
        if (file_status /= 0) then
            print*, "Error opening file for writing"
            return
        endif

        write(1,*) "Mean value is:", med_mean, &
                                " and standard deviation is:", med_std

        close(1)
    end subroutine binning

    ! CALCULATE MEDIAN OF DATA !
    real function calculate_median(arr, n)
        ! arr: array with data to calculate median 
        ! n: length of array
        integer, intent(in) :: n
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

! end program binning_gestor


