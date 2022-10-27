program adder

    use tables
    use misc
    implicit none

    real(dp), allocatable       :: initial(:,:), candidate(:,:), total(:,:), &
        scaled(:,:), initial_scaled(:,:), candidate_scaled(:,:), &
        self_distance(:,:), cross_distance(:,:)
    logical, allocatable        :: selected_mask(:)
    integer, allocatable        :: selected(:), not_selected(:)
    character(len=80)           :: fname, dummy
    integer                     :: L, M, N, number_to_select, i, j, &
        i_remove, i_add, placeholder
    real(dp)                    :: temperature, energy, acceptance_rate, &
        dE, probability, dT
    real(dp)                    :: time_start, time_finish

    call get_command_argument(1, fname)
    call read_data(trim(fname), initial)
    call get_command_argument(2, fname)
    call read_data(trim(fname), candidate)
    call get_command_argument(3, dummy)
    read(dummy, *) number_to_select
    call get_command_argument(4, dummy)
    read(dummy, *) temperature

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
    call self_energy(candidate_scaled, self_distance)
    call cross_energy(initial_scaled, candidate_scaled, cross_distance)
    print *, 'calculation of distances complete!'

    allocate(selected_mask(L))
    allocate(selected(number_to_select))
    allocate(not_selected(L - number_to_select))
    print *, 'going to select ', number_to_select, ' points'
    print *, 'initially picking the first ', number_to_select, ' points'
    ! initialize selection
    do i = 1, number_to_select
        selected(i) = i
    end do
    do i = 1, L - number_to_select
        not_selected(i) = i + number_to_select
    end do
    selected_mask = .false.
    selected_mask(:number_to_select) = .true.
    energy = total_energy(self_distance, cross_distance, selected_mask)
    print *, 'initial energy: ', energy

    call cpu_time(time_finish)
    print *, 'initialization time: ', time_finish - time_start, 's'

    ! fixed: cooling in 100 steps, 10*L samples at each temp
    dT = temperature / 100.0_dp

    call cpu_time(time_start)

    print *, 'begin annealing at temperature', temperature
    cooling: do i = 1, 100
        acceptance_rate = 0
        sampling: do j = 1, 10*L
            i_remove = randrange(number_to_select)
            i_add = randrange(L - number_to_select)
            dE = change_in_energy(self_distance, cross_distance, &
                selected(i_remove), not_selected(i_add), selected_mask)
            call random_number(probability)
            if (probability < exp(-dE / temperature)) then
                ! accept and perform the swap
                placeholder = selected(i_remove)
                selected(i_remove) = not_selected(i_add)
                not_selected(i_add) = placeholder
                ! update mask
                selected_mask(selected(i_remove)) = .true.
                selected_mask(not_selected(i_add)) = .false.
                energy = energy + dE
                acceptance_rate = acceptance_rate + 1
            end if
        end do sampling
        print *, temperature, energy, acceptance_rate/(10*L)
        temperature = temperature - dT
    end do cooling
    print *, 'done annealing!'

    call cpu_time(time_finish)
    print *, 'annealing time: ', time_finish - time_start, 's'

    ! final sanity check
    print *, 'quick sanity check...'
    print *, 'final energy: ', energy
    print *, 'expected: ', &
        total_energy(self_distance, cross_distance, selected_mask)
    
    ! print out selections
    print *, 'selected indices (starting from 0, for python):'
    do i = 1, number_to_select
        print *, selected(i) - 1
    end do
    print *, 'selected data:'
    do i = 1, number_to_select
        print *, candidate(:,selected(i))
    end do

end program adder
