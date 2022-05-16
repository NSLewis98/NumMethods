module ODESolvers
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
    function euler_step(fPtr, t, y, h)
        procedure(func), pointer              :: fPtr
        real(kind(1d0)), optional, intent(in) :: t
        real(kind(1d0))          , intent(in) :: y(:)
        real(kind(1d0))          , intent(in) :: h
        real(kind(1d0)), dimension(size(y)) :: euler_step
        euler_step = y + h*fPtr(t, y)
    end function euler_step

    function euler_solver(fPtr, t0, tf, y0, numpts)
        procedure(func), pointer    :: fPtr
        real(kind(1d0)), intent(in) :: t0, tf
        real(kind(1d0)), intent(in) :: y0(:)
        integer        , intent(in) :: numpts

        real(kind(1d0)), dimension(size(y0), numpts) :: euler_solver

        real(kind(1d0)) :: h, tn, y_prev(size(y0))
        integer :: n

        euler_solver = 0d0
        euler_solver(:, 1) = y0
        h = (tf - t0) / dble(numpts)
        do n = 1, numpts - 1
            tn = h*dble(n) + t0
            y_prev = euler_solver(:, n)
            euler_solver(:, n+1) = euler_step(fPtr, tn, y_prev, h)
        enddo
    end function euler_solver

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
    function rk4_solver(fPtr, t0, tf, y0, numpts)
        procedure(func), pointer    :: fPtr
        real(kind(1d0)), intent(in) :: t0, tf
        real(kind(1d0)), intent(in) :: y0(:)
        integer        , intent(in) :: numpts

        real(kind(1d0)), dimension(size(y0), numpts) :: rk4_solver
        real(kind(1d0)) :: h, tn, y_prev(size(y0))
        integer :: n
        
        rk4_solver = 0d0
        rk4_solver(:, 1) = y0
        h = (tf - t0) / dble(numpts)
        do n = 1, numpts - 1
            tn = h*dble(n) + t0
            y_prev = rk4_solver(:, n)
            rk4_solver(:, n+1) = rk4_step(fPtr, tn, y_prev, h)
        enddo
    end function rk4_solver
end module ODESolvers

program main
    use ODESolvers
    implicit none
    real(kind(1d0)), parameter :: PI = 3.14159265358979323846d0
    integer        , parameter :: N = 1000

    procedure(func), pointer :: fPtr
    real(kind(1d0)) :: y0(1), t0, tf, h
    real(kind(1d0)), dimension(size(y0), N) :: y, w

    real(kind(1d0)) :: xi
    integer         :: i

    fPtr => ODE

    y0 = 0.5d0
    t0 = 0d0
    tf = 2d0*PI
    
    y =  euler_solver(fPtr, t0, tf, y0, N)
    w = rk4_solver(fPtr, t0, tf, y0, N)
    h = (tf - t0) / DBLE(N)
    print*, h
    do i = 1, N 
        xi = h * DBLE(i) + t0
        write(10, *) xi, y(:, i), w(:, i)
    enddo

    contains

    function ODE(t,y)
        real(kind(1d0)), optional, intent(in)  :: t
        real(kind(1d0)), intent(in)            :: y(:)
        real(kind(1d0)), dimension(size(y))    :: ODE
        ODE = -sin(y)
    end function ODE

    
end program main