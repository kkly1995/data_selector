module misc
    ! collection of things not related to distance tables

    implicit none

    integer, parameter  :: dp = kind(0.d0) ! double precision

    contains

        subroutine read_data(fname, x)
            ! read data from fname into unallocated x
            ! i.e. the size of x will be determined directly from file
            ! file is a table of number, where the first line is
            ! '# N M', indicating a N by M table of numbers
            character(*), intent(in)                    :: fname
            real(dp), allocatable, intent(in out)       :: x(:,:)
            integer                                     :: fileunit, N, M, total
            character(len=1)                            :: hashtag
            integer                                     :: exit_status

            open(newunit=fileunit, file=fname)

            read(fileunit, *) hashtag, N, M
            if (hashtag .ne. '#') then
                print *, 'warning: expected # as first character, ', &
                    'instead got ', hashtag
                call exit(exit_status)
            end if
            print *, 'reading', N, ' data points of dimension', M, &
                ' from ', fname

            ! can now allocate and read
            allocate(x(M,N))
            read(fileunit, *) x

            close(fileunit)
        end subroutine read_data

        function randrange(N)
            ! name and functionality borrowed from python's random.randrange()
            ! except instead of returning integer between [0, N),
            ! returns integer between [1, N] (notice brackets)
            integer, intent(in) :: N
            integer             :: randrange

            real                :: x
            
            call random_number(x)
            x = N*x
            randrange = floor(x)
            randrange = randrange + 1
        end function randrange

        subroutine standard_scale(x, y)
            ! rescales data x so that it has mean 0 and variance 1
            ! the rescaled data is placed into y
            real(dp), intent(in)        :: x(:,:)
            real(dp), intent(in out)    :: y(:,:)

            real(dp)                    :: mean, std
            real(dp), allocatable       :: dx(:)
            integer                     :: M, N, i

            M = size(x, dim=1)
            N = size(x, dim=2)
            allocate(dx(N))

            ! do each dimension one by one
            !$OMP PARALLEL DO SHARED(y) PRIVATE(mean, dx, std)
            do i = 1, M
                mean = sum(x(i,:)) / N
                dx = x(i,:) - mean
                std = sum(dx**2) / N
                std = sqrt(std)
                y(i,:) = (x(i,:) - mean) / std
            end do
            !$OMP END PARALLEL DO

        end subroutine standard_scale
    
end module misc
