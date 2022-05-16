module RK45
    implicit none

    real(kind(1d0)), parameter :: AK1 = 0d0
    real(kind(1d0)), parameter :: AK2 = 2d0/9d0
    real(kind(1d0)), parameter :: AK3 = 1d0 / 3d0
    real(kind(1d0)), parameter :: AK4 = 0.75d0
    real(kind(1d0)), parameter :: AK5 = 1d0
    real(kind(1d0)), parameter :: AK6 = 5d0 /6d0

    real(kind(1d0)), parameter :: B21 = 2d0 /9d0
    real(kind(1d0)), parameter :: B31 = 1d0 / 12d0
    real(kind(1d0)), parameter :: B41 = 69d0 / 128d0
    real(kind(1d0)), parameter :: B51 = -17d0 / 12d0
    real(kind(1d0)), parameter :: B61 = 65d0 / 432d0

    real(kind(1d0)), parameter :: B32 = 0.25d0 
    real(kind(1d0)), parameter :: B42 = -243d0 / 128d0 
    real(kind(1d0)), parameter :: B52 = 27d0 / 4d0  
    real(kind(1d0)), parameter :: B62 = -5d0 /16d0 

    real(kind(1d0)), parameter :: B43 = 135d0 / 64d0 
    real(kind(1d0)), parameter :: B53 = -27d0 / 5d0 
    real(kind(1d0)), parameter :: B63 = 13d0 / 16d0

    real(kind(1d0)), parameter :: B54 = 16d0 / 15d0 
    real(kind(1d0)), parameter :: B64 = 4d0/ 27d0

    real(kind(1d0)), parameter :: B65 =  5d0 / 144d0

    real(kind(1d0)), parameter :: CK1 = 1d0/9d0 
    real(kind(1d0)), parameter :: CK2 = 0d0 
    real(kind(1d0)), parameter :: CK3 = 9d0 /20d0
    real(kind(1d0)), parameter :: CK4 = 16d0 / 45d0
    real(kind(1d0)), parameter :: CK5 = 1d0 / 12d0

    real(kind(1d0)), dimension(6), parameter :: CTl = &
                    [-1d0/150d0 ,0d0, 3d0/100d0, -16d0 /75d0, -1d0/20d0, 6d0/25d0] 
    real(kind(1d0)), dimension(6), parameter :: CHm = &
                    [47d0/450d0 , 0d0, 3d0 /100d0, -16d0 /75d0, -1d0/20d0, 6d0/25d0]
    real(kind(1d0)), dimension(5), parameter :: Cki = &
                    [1d0/9d0 ,  0d0 , 9d0/20d0, 16d0/45d0, 1d0/12d0]

    interface
        function func(t,y)
            real(kind(1d0)), optional, intent(in)  :: t
            real(kind(1d0)), intent(in)            :: y(:)
            real(kind(1d0)), dimension(size(y))    :: func
        end function func
    end interface 
    public func
    contains 
    subroutine RK45Step(fPtr, t, y0, yn, h, hBounds, tolerance)
        procedure(func), pointer       :: fPtr
        real(kind(1d0)), intent(inout) :: t, h
        real(kind(1d0)), intent(in)    :: hBounds(2), tolerance, y0(:)
        real(kind(1d0)), dimension(size(y0)), intent(out) :: yn
        

        real(kind(1d0)), dimension(size(y0), 6) :: k
        real(kind(1d0)) :: truncError, h0
        integer         :: m

        do m = 1, 100
            k(:, 1) = h*fPtr(t, y0)
            k(:, 2) = h*fPtr( t+AK2*h, y0+B21*k(:,1))
            k(:, 3) = h*fPtr( t+AK3*h, y0+B31*k(:,1) + B32*k(:,2))
            k(:, 4) = h*fPtr( t+AK4*h, y0+B41*k(:,1)+B42*k(:,2)+B43*k(:,3))
            k(:, 5) = h*fPtr( t+AK5*h, y0+B51*k(:,1)+B52*k(:,2)+B53*k(:,3)+B54*k(:,4))
            k(:, 6) = h*fPtr( t+AK6*h, y0+B61*k(:,1)+B62*k(:,2)+B63*k(:,3)+B64*k(:,4)+B65*k(:,5))
            truncError = EstTruncError(k)
            h0 = h
            h = UpdateStepSize(h0, hBounds, truncError, tolerance)
            if(h0.eq.h) exit ! Means h wasnt changed
        enddo

        yn = y0
        do m = 1, 6
            yn = yn + k(:, m)*CHm(m)
        enddo
        
        t = t + h

    end subroutine RK45Step
    
    real(kind(1d0)) function EstTruncError(k) result(TE)
        real(kind(1d0)), intent(in) :: k(:,:)
        real(kind(1d0)), dimension(size(k(:,:), 1)) :: kSum
        integer :: i

        TE = 0d0
        kSum = 0d0

        do i =1, 6
            kSum = kSum + k(:, i)*CTl(i)
        enddo
        TE = SQRT( dot_product(kSum, kSum) )
        write(*,*) "TE = ", TE
    end function 

    real(kind(1d0)) function UpdateStepSize(h, hBounds, error, tolerance) result(hNew)
        real(kind(1d0)), intent(in) :: h, hBounds(2), error, tolerance
        hNew = h
        if(error.le.tolerance) return
        hNew = 0.9d0 * h * (tolerance / error)**(1d0/5d0)
        if(hNew.lt.hBounds(1)) then
            hNew = hBounds(1)
            write(*,*) "hNew = ", hNew
            return
        elseif(hNew.gt.hBounds(2)) then
            hNew = hBounds(2)
            write(*,*) "hNew = ", hNew
            return
        endif

    end function UpdateStepSize
end module

program main
    use RK45
    implicit none
    real(kind(1d0)), parameter :: TOLERANCE = 1d-10 
    procedure(func), pointer :: fPtr
    real(kind(1d0)) :: y0(1), y(1), t0, hBounds(2), h0

    fPtr => ODE
    y0(1) = 0d0
    t0 = 0d0
    
    hBounds(1) = TOLERANCE - 1d-6
    hBounds(2) = 0.1d0
    h0 = hBounds(2)
    print*, t0, y0, h0
    call RK45Step(fPtr, t0, y0, y, h0, hBounds, TOLERANCE)
    print*, t0, y, h0
    contains
    function ODE(t, y)
        real(kind(1d0)), optional, intent(in)  :: t
        real(kind(1d0)), intent(in)            :: y(:)
        real(kind(1d0)), dimension(size(y))    :: ODE
        ODE = cos(t)
    end function ODE
end program main