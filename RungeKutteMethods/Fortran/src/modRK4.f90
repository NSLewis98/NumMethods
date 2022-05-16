! Simple Module that solves a N-D Linear ODE using the Runge-Kutta 4 method
!.. Note that the function pointer past to the solver 

module RK4
implicit none
    interface
        function func(t,y)
            real(kind(1d0)), optional, intent(in)  :: t
            real(kind(1d0)), intent(in)            :: y(:)
            real(kind(1d0)), dimension(size(y))    :: func
        end function func
    end interface 
    public func
    
    contains 
    ! Calculate a single Runge-Kutta 4 step
    function rk4_step(fPtr, t, y, h)
        procedure(func), pointer              :: fPtr
        real(kind(1d0)), optional, intent(in) :: t
        real(kind(1d0))          , intent(in) :: y(:)
        real(kind(1d0))          , intent(in) :: h
        
        real(kind(1d0)), dimension(size(y)) :: rk4_step
        real(kind(1d0)), dimension(size(y)) :: k1, k2, k3, k4

        k1 = fPtr(t, y)
        k2 = fPtr(t+ 0.5d0 *h, y + h*0.5d0*k1)
        k3 = fPtr(t+ 0.5d0 *h, y + h*0.5d0*k2)
        k4 = fPtr(t+ h, y+ h*k3)
        rk4_step = y + h*1d0/6d0 * (k1 + 2d0*k2 + 2d0*k3 + k4)
    end function rk4_step

    ! Solve a system of linear ODEs using a Runge-Kutta 4 method
    function rk4_solve(fPtr, t0, tf, y0, numpts)
        procedure(func), pointer    :: fPtr
        real(kind(1d0)), intent(in) :: t0, tf
        real(kind(1d0)), intent(in) :: y0(:)
        integer        , intent(in) :: numpts

        real(kind(1d0)), dimension(size(y0), numpts) :: rk4_solve
        real(kind(1d0)) :: h, tn, y_prev(size(y0))
        integer :: n
        
        rk4_solve = 0d0
        rk4_solve(:, 1) = y0
        h = (tf - t0) / dble(numpts)
        print*, "h = ", h
        do n = 1, numpts - 1
            tn = h*dble(n) + t0
            y_prev = rk4_solve(:, n)
            rk4_solve(:, n+1) = rk4_step(fPtr, tn, y_prev, h)
        enddo
    end function rk4_solve
end module RK4
