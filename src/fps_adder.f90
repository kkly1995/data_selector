program fps_adder

    use tables
    use misc
    implicit none

    real(dp), allocatable       :: initial(:,:), candidate(:,:), total(:,:), &
        scaled(:,:), initial_scaled(:,:), candidate_scaled(:,:), &
        self_distance(:,:), cross_distance(:,:), nearest_neighbor(:)
    integer, allocatable        :: selected(:)
    integer                     :: placeholder(1) ! for maxloc
    character(len=80)           :: fname, dummy
    integer                     :: L, M, N, number_to_select, i, j
    real(dp)                    :: time_start, time_finish

    call get_command_argument(1, fname)
    call read_data(trim(fname), initial)
    call get_command_argument(2, fname)
    call read_data(trim(fname), candidate)
    call get_command_argument(3, dummy)
    read(dummy, *) number_to_select

    call cpu_time(time_start)

    ! can now allocate and calculate
    L = size(candidate, dim=2)
    M = size(candidate, dim=1)
    N = size(initial, dim=2)

    ! combine data sets for scaling
    allocate(total(M,N+L))
    total(:,:N) = initial
    total(:,N+1:) = candidate
    allocate(scaled(M,N+L))
    call standard_scale(total, scaled)
    ! separate again
    allocate(initial_scaled(M,N))
    allocate(candidate_scaled(M,L))
    initial_scaled = scaled(:,:N)
    candidate_scaled = scaled(:,N+1:)
    print *, 'standard scaling of data complete!'

    allocate(self_distance(L,L))
    allocate(cross_distance(N,L))
    call self_table(candidate_scaled, self_distance)
    call cross_table(initial_scaled, candidate_scaled, cross_distance)
    print *, 'calculation of distances complete!'

    allocate(selected(number_to_select))
    print *, 'going to select ', number_to_select, ' points'
    ! initialize nearest neighbor array
    allocate(nearest_neighbor(L))
    nearest_neighbor = minval(cross_distance, dim=1)

    call cpu_time(time_finish)
    print *, 'initialization time: ', time_finish - time_start, 's'

    call cpu_time(time_start)

    do i = 1, number_to_select
        placeholder = maxloc(nearest_neighbor) ! since this returns array
        selected(i) = placeholder(1)
        ! update nearest neighbors
        do j = 1, L
            if (self_distance(selected(i),j) < nearest_neighbor(j)) then
                nearest_neighbor(j) = self_distance(selected(i), j)
            end if
        end do
    end do

    call cpu_time(time_finish)
    print *, 'search time: ', time_finish - time_start, 's'
    
    ! print out selections
    print *, 'selected indices (starting from 0, for python):'
    do i = 1, number_to_select
        print *, selected(i) - 1
    end do
    print *, 'selected data:'
    do i = 1, number_to_select
        print *, candidate(:,selected(i))
    end do

end program fps_adder
