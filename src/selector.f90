program selector

    use tables
    use misc
    implicit none

    real, allocatable           :: candidate(:,:), scaled(:,:), distance(:,:)
    logical, allocatable        :: selected_mask(:)
    integer, allocatable        :: selected(:), not_selected(:)
    character(len=80)           :: fname, dummy
    integer                     :: M, N, number_to_select, i, j, &
        i_remove, i_add, placeholder
    real                        :: temperature, energy, acceptance_rate, &
        dE, probability, dT
    real                        :: time_start, time_finish

    call get_command_argument(1, fname)
    call read_data(trim(fname), candidate)
    call get_command_argument(2, dummy)
    read(dummy, *) number_to_select
    call get_command_argument(3, dummy)
    read(dummy, *) temperature

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

    allocate(selected_mask(N))
    allocate(selected(number_to_select))
    allocate(not_selected(N - number_to_select))
    print *, 'going to select ', number_to_select, ' points'
    print *, 'initially picking the first ', number_to_select, ' points'
    ! initialize selection
    do i = 1, number_to_select
        selected(i) = i
    end do
    do i = 1, N - number_to_select
        not_selected(i) = i + number_to_select
    end do
    selected_mask = .false.
    selected_mask(:number_to_select) = .true.
    energy = total_self_energy(distance, selected_mask)
    print *, 'initial energy: ', energy

    call cpu_time(time_finish)
    print *, 'initialization time: ', time_finish - time_start, 's'

    ! fixed: cooling in 1000 steps, 10*N samples at each temp
    dT = temperature / 1000

    call cpu_time(time_start)

    print *, 'begin annealing at temperature', temperature
    cooling: do i = 1, 1000
        acceptance_rate = 0
        sampling: do j = 1, 10*N
            i_remove = randrange(number_to_select)
            i_add = randrange(N - number_to_select)
            dE = change_in_self_energy(distance, selected(i_remove), &
                not_selected(i_add), selected_mask)
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
        print *, temperature, energy, acceptance_rate/(10*N)
        temperature = temperature - dT
    end do cooling
    print *, 'done annealing!'

    call cpu_time(time_finish)
    print *, 'annealing time: ', time_finish - time_start, 's'

    ! final sanity check
    print *, 'quick sanity check...'
    print *, 'final energy: ', energy
    print *, 'expected: ', total_self_energy(distance, selected_mask)
    
    ! print out selections
    print *, 'selected indices (starting from 0, for python):'
    do i = 1, number_to_select
        print *, selected(i) - 1
    end do
    print *, 'selected data:'
    do i = 1, number_to_select
        print *, candidate(:,selected(i))
    end do

end program selector
