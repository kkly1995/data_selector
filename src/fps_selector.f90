program fps_selector

    use tables
    use misc
    implicit none

    real(dp), allocatable       :: candidate(:,:), scaled(:,:), distance(:,:), &
        nearest_neighbor(:)
    integer, allocatable        :: selected(:)
    character(len=80)           :: fname, dummy
    integer                     :: M, N, number_to_select, i, j, initial, &
        placeholder(1) ! for maxloc
    real(dp)                    :: time_start, time_finish

    call get_command_argument(1, fname)
    call read_data(trim(fname), candidate)
    call get_command_argument(2, dummy)
    read(dummy, *) number_to_select
    call get_command_argument(3, dummy)
    read(dummy, *) initial

    call cpu_time(time_start)

    ! can now allocate and calculate
    M = size(candidate, dim=1)
    N = size(candidate, dim=2)

    allocate(scaled(M,N))
    call standard_scale(candidate, scaled)
    print *, 'standard scaling of data complete!'

    allocate(distance(N,N))
    call self_table(scaled, distance)
    print *, 'calculation of distances complete!'

    allocate(selected(number_to_select))
    allocate(nearest_neighbor(N))
    print *, 'going to select ', number_to_select, ' points'
    print *, 'initially picking point ', initial
    selected(1) = initial
    nearest_neighbor = distance(:,initial) ! only one neighbor at the start

    call cpu_time(time_finish)
    print *, 'initialization time: ', time_finish - time_start, 's'

    call cpu_time(time_start)

    do i = 2, number_to_select
        placeholder = maxloc(nearest_neighbor)
        selected(i) = placeholder(1)
        ! update nearest neighbors
        do j = 1, N
            if (distance(selected(i), j) < nearest_neighbor(j)) then
                nearest_neighbor(j) = distance(selected(i), j)
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

end program fps_selector
