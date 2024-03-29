module tables

    use misc, only: dp
    implicit none

    real(dp), parameter :: large_number = 10000.0_dp, small_number = 0.00001_dp

    contains

        subroutine self_table(x, y)
            ! x is an array of size (M, N)
            ! representing N data points in M-dimensional space
            ! y is an array of size (N, N)
            ! where y(i,j) represents euclidean distance
            ! between x(i) and x(j)
            real(dp), intent(in)        :: x(:,:)
            real(dp), intent(in out)    :: y(:,:)
            integer                     :: N, i, j

            N = size(x, dim=2)
            y = 0
            do j = 1, N
                do i = 1, j - 1
                    y(i,j) = norm2(x(:,i) - x(:,j))
                    y(j,i) = y(i,j)
                end do
            end do
        end subroutine self_table

        subroutine self_energy(x, y)
            ! same as self_table but calculates an energy instead
            ! i.e. y(i,j) represents the RECIPROCAL of the euclidean distance
            ! between x(i) and x(j)
            !
            ! note: if x(i) and x(j) are sufficiently close,
            ! their entry will be set to large_number
            real(dp), intent(in)        :: x(:,:)
            real(dp), intent(in out)    :: y(:,:)
            integer                     :: N, i, j
            real(dp)                    :: r

            N = size(x, dim=2)
            y = 0
            !$omp parallel do shared(y) private(i, r)
            do j = 1, N
                do i = 1, j - 1
                    r = norm2(x(:,i) - x(:,j))
                    if (r > small_number) then
                        y(i,j) = 1.0_dp / r
                    else
                        y(i,j) = large_number
                    end if
                    y(j,i) = y(i,j)
                end do
            end do
            !$omp end parallel do
        end subroutine self_energy

        subroutine cross_table(x, y, z)
            ! computes distances between two sets of data points x and y
            ! x has size (M, N)
            ! y has size (M, L), i.e. data has same dimension
            ! z has size (N, L) and will store the final table
            ! z(i,j) represents the distance between x(:,i) and y(:,j)
            real(dp), intent(in)        :: x(:,:), y(:,:)
            real(dp), intent(in out)    :: z(:,:)
            
            integer                     :: N, L, i, j

            N = size(x, dim=2)
            L = size(y, dim=2)
            do i = 1, N
                do j = 1, L
                    z(i,j) = norm2(x(:,i) - y(:,j))
                end do
            end do
        end subroutine cross_table

        subroutine cross_energy(x, y, z)
            ! same as cross_table but computes energy instead
            ! e.g. z(i,j) represents the reciprocal distance between
            ! x(:,i) and y(:,j)
            !
            ! again, if these two points are suffciently close, 
            ! their entry will be set to large_number
            real(dp), intent(in)        :: x(:,:), y(:,:)
            real(dp), intent(in out)    :: z(:,:)
            
            integer                     :: N, L, i, j
            real(dp)                    :: r

            N = size(x, dim=2)
            L = size(y, dim=2)
            !$omp parallel do shared(z) private(j, r)
            do i = 1, N
                do j = 1, L
                    r = norm2(x(:,i) - y(:,j))
                    if (r > small_number) then
                        z(i,j) = 1.0_dp / r
                    else
                        z(i,j) = large_number
                    end if
                end do
            end do
            !$omp end parallel do
        end subroutine cross_energy

        function total_self_energy(x, indices) result(E)
            ! calculates total energy for a selected set of data points
            ! x is the table of distances, having size (N, N)
            ! indices is a mask indicating which points are selected,
            ! having size N
            real(dp), intent(in):: x(:,:)
            logical, intent(in) :: indices(:)
            real(dp)            :: E

            integer             :: i, N
            
            E = 0
            N = size(x, dim=1)
            !$omp parallel do reduction(+:E)
            do i = 1, N
                if (indices(i)) then
                    E = E + sum(x(:,i), mask=indices)
                end if
            end do
            !$omp end parallel do
            E = 0.5_dp * E
        end function total_self_energy

        function total_energy(x, y, indices) result(E)
            ! calculates total energy
            ! for a fixed data set + candidate data set
            ! i.e. x is a table of self distances for the candidate set
            ! and y is a table of distances between the fixed and candidate
            ! indices is a mask indicating which candidates are selected
            !
            ! NOTE: if the candidate set has L points and the fixed set has N,
            ! x will have size (L,L)
            ! and y will have size (N,L)
            ! i.e. cross_table should be called s.t. dim 1 of the cross table
            ! corresponds to the fixed set
            real(dp), intent(in):: x(:,:), y(:,:)
            logical, intent(in) :: indices(:)
            real(dp)            :: E

            integer             :: N, L, i

            L = size(y, dim=2)
            N = size(y, dim=1)
            
            E = 0
            !$omp parallel do reduction (+:E)
            do i = 1, L
                if (indices(i)) then
                    E = E + sum(x(:,i), mask=indices)
                end if
            end do
            !$omp end parallel do
            E = 0.5_dp * E
            !$omp parallel do reduction (+:E)
            do i = 1, N
                E = E + sum(y(i,:), mask=indices)
                ! getting y(i,:) instead of y(:,i) is obviously slower
                ! but this routine is only called a few times
                ! whereas change_in_energy, which will grab y(:,i)
                ! will be called many times
            end do
            !$omp end parallel do
        end function total_energy

        function change_in_self_energy(x, i, j, indices) result(dE)
            ! calculates the change in energy if index i is removed
            ! and j is added
            ! indices is the current mask
            ! i.e. indices(i) must be true and indices(j) must be false
            real(dp), intent(in):: x(:,:)
            integer, intent(in) :: i, j
            logical, intent(in) :: indices(:)
            real(dp)            :: dE

            dE = -sum(x(:,i), mask=indices)
            dE = dE + sum(x(:,j), mask=indices)
!            N = size(x, dim=1)
!            dE = 0
!            !$omp simd reduction (+:dE)
!            do k = 1, N
!                if (indices(k)) dE = dE - x(k,i) + x(k,j)
!            end do
!            !$omp end simd
            dE = dE - x(i,j)
        end function change_in_self_energy

        pure function change_in_energy(x, y, i, j, indices) result(dE)
            ! calculates the change in total energy
            ! if data point i is removed from the candidate set
            ! and data point j is added
            ! x is the table of self distances for the candidate set
            ! y is the table of cross distances between the fixed set
            ! and the candidate set
            !
            ! again, the convention is that if there are L candidate points
            ! x should have size (L,L)
            ! and y should have size (N,L)
            real(dp), intent(in):: x(:,:), y(:,:)
            integer, intent(in) :: i, j
            logical, intent(in) :: indices(:)
            real(dp)            :: dE

            dE = -sum(x(:,i), mask=indices) - sum(y(:,i))
            dE = dE + sum(x(:,j), mask=indices) + sum(y(:,j))
            dE = dE - x(i,j)
        end function change_in_energy

end module tables
