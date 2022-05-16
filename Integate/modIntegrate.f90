! Module containing numerical integration methods.
! The LebdevIntegral integral routine is the most efficient 3D integration method.
! It uses a Lebedev pt Quadrrature method for integration over \mathcal{S}^2 and 
! trapezoiad rule for the radial integration

module modIntegral
    implicit none
    interface
    !.. Expects a function f: \mathcal{R}^3 -> \mathcal{C}
    !.. The variables are in sphercial coordinates: (r, \theta, \phi)
        complex(kind(1d0)) function func(r,theta, phi, parvec) result(res)
            real(kind(1d0)), optional, intent(in)  :: r, theta, phi, parvec(*)        
        end function func
    end interface 
    public func

    public  :: OneDimIntegral, ThreeDimIntegral, LebdevIntegral
    private :: getLebPoints, weights
    contains
    
    ! Numericaly Integrate a 1D complex function, f,  using trapezoidal rule
    !.. f: R -> C
    complex(kind(1d0)) function OneDimIntegral(fPtr, upperBound, lowerBound, numPoints, parvec) result(Integral)
        procedure(func) , pointer               :: fPtr
        real(kind(1d0)) , intent(in)            :: upperBound, lowerBound
        integer         , intent(in)            :: numPoints
        real(kind(1d0)) , optional, intent(in)  :: parvec(*)

        real(kind(1d0)) :: xi
        integer :: i

        Integral = (0d0, 0d0)
        xi = (0d0, 0d0)
        Integral = 0.5d0 *( fPtr(r=upperBound, parvec=parvec) + fPtr(r=lowerBound, parvec=parvec) )

        do i = 1, numPoints - 1
            xi = (upperBound - lowerBound) / dble(numPoints) * dble(i) + lowerBound
            Integral = Integral + fPtr(r=xi, parvec=parvec)
        enddo
        Integral = Integral * (upperBound - lowerBound) / (numPoints)
    end function OneDimIntegral

    ! Numericaly Integrate a 3D complex function, f,  using trapezoidal rule
    !.. f: R^3 -> C; f(r, theta, phi)
    complex(kind(1d0)) function ThreeDimIntegral(fPtr, r_max, numPoints, parvec) result(Integral)
        procedure(func),           pointer    :: fPtr
        real(kind(1d0)),           intent(in) :: r_max
        integer        ,           intent(in) :: numPoints(3)
        real(kind(1d0)), optional, intent(in) :: parvec(*)
        DoublePrecision, parameter :: PI = 3.14159265358979323846d0

        real(kind(1d0)) :: rNumPoints(3), phik, thetaj, wTheta, ri, wR
        integer         :: i, j, k

        rNumPoints = DBLE(numPoints)
        Integral = (0d0, 0d0)

        do k =1, numPoints(3)
            phik = 2d0*PI/rNumPoints(3) *DBLE(k -1)
            do j =1, numPoints(2)
                thetaj = PI / (rNumPoints(2) - 1d0) *DBLE(j-1)
                wTheta = weights(j, 1, numPoints(2))
                do i = 1, numPoints(1)
                    ri = r_max / (rNumPoints(1) - 1d0) * DBLE(i-1)
                    wR = weights(i, 1, numPoints(1))
                    Integral = Integral + ri**2 * sin(thetaj)*fPtr(r=ri,theta=thetaj,phi=phik, parvec=parvec)&
                                *wTheta*wR
                enddo
            enddo
        enddo
        Integral = Integral * 2d0 * PI**2 * r_max &
        /((rNumPoints(1) - 1)*(rNumPoints(2)-1)*rNumPoints(1) )
        return
    end function ThreeDimIntegral

    real(kind(1d0)) function weights(n, N_MIN, N_MAX) result(res)
        integer, intent(in) :: n, N_MIN, N_MAX
        res = 1d0
        if( (n.eq.N_MIN).or.(n.eq.N_MAX)) then
            res = 0.5d0
            return
        endif
    end function weights

    ! *************** Lebedev Quadrature ***************

    !.. Make sure works
    ! Numericaly Integrate a 3D complex function, f,
    ! using Lebdev rule for integral over \mathcal{S}^2 and trapezoidal for the radial integral
    !.. f: R^3 -> C; f(r, theta, phi) 
    complex(kind(1d0)) function LebdevIntegral(fPtr, max_r, numRpts, numLebPts) result(res)
        procedure(func)                   , pointer    :: fPtr
        real(kind(1d0))                   , intent(in) :: max_r
        integer                           , intent(in) :: numLebPts, numRpts

        DoublePrecision, parameter :: PI      = 3.14159265358979323846d0
        real(kind(1d0)), dimension(numLebPts) :: L_wts, theta, phi
        real(kind(1d0)) :: ri, rWt
        integer :: i, j

        call getLebPoints(numLebPts, L_wts, theta, phi)

        res = (0d0, 0d0)
        do i = 1, numRpts
            ri = max_r / (numRpts - 1) * dble(i - 1)
            rWt = weights(i, 1, numRpts)
            do j = 1, numLebPts
                res = res + L_wts(j)*rWt*fPtr(ri, theta(j), phi(j))*ri**2
            enddo
        enddo
        res =4d0*PI*max_r*res/(numRpts - 1)
    end function LebdevIntegral
    
    !.. Make sure works
    ! Opens a txt file given the numPts passed and fills the lebdev weights, and solid angle vals
    subroutine getLebPoints(numPts, weights, theta, phi)
        integer                           , intent(in)  :: numPts
        real(kind(1d0)), dimension(numPts), intent(out) :: weights, theta, phi

        DoublePrecision, parameter :: PI      = 3.14159265358979323846d0

        integer          :: i
        character(len=5) :: x1
        character(len=8) :: fmt ! format descriptor

        fmt = '(I4.4)' ! an integer of width 5 with zeros at the left

        write (x1,fmt) numPts ! converting integer to string using a 'internal file'

        weights = 0d0 
        theta   = 0d0 
        phi     = 0d0 
        ! write(*,*) "Reading File: ", 'Leb'//trim(x1)//'pt.txt'
        open( 10, file = '../src/LebPts/Leb'//trim(x1)//'pt.txt', status='old')
        do i = 1, numPts
            read( 10, *) phi(i), theta(i), weights(i)
        end do
        theta = PI / 180d0 * theta 
        phi   = PI / 180d0 * phi + PI/2d0 
        close(10)
        return
    end subroutine getLebPoints

end module modIntegral

program main
    use modIntegral
    use ModuleAngularMomentum
    implicit none
    DoublePrecision, parameter :: PI      = 3.14159265358979323846d0
    real(kind(1d0)), parameter :: MAX_R   = 1d0
    real(kind(1d0)), parameter :: MIN_R   = 0d0
    integer        , parameter :: NUM_PTS = 300
    integer        , parameter :: LEB_PTS = 5810
    integer        , parameter :: R_PTS   = 100

    procedure(func), pointer :: fPtr
    complex(kind(1d0))       :: Integral
    integer, dimension(3)    :: numPts
    
    numPts = NUM_PTS
    fPtr => f
    print*, "1D Integral of f"
    Integral = OneDimIntegral(fPtr, MAX_R, MIN_R, NUM_PTS)
    print*, Integral
    
    fPtr => g
    print*, "3D Integral of g using trap"
    print*, "Num of pts: ", NUM_PTS**3
    Integral = ThreeDimIntegral(fPtr, MAX_R, numPts)
    print*, Integral
    
    Integral = LebdevIntegral(fPtr, MAX_R, R_PTS, LEB_PTS)
    print*, "3D Integral of g using leb and trap"
    print*, "Num of pts: ", R_PTS*LEB_PTS
    print*, Integral

    print*, "The trap method uses ", dble(NUM_PTS**3)  / dble(R_PTS*LEB_PTS) , "x more pts"
    contains
    
    complex(kind(1d0)) function f(r,theta, phi, parvec) result(res)
        real(kind(1d0)), optional, intent(in)  :: r, theta, phi, parvec(*)
        res = (1d0 ,0d0)*cos(r)
    end function f
    
    complex(kind(1d0)) function g(r,theta, phi, parvec) result(res)
        real(kind(1d0)), optional, intent(in)  :: r, theta, phi, parvec(*)
        !Ylm(theta,phi,l,m)
        res = (1d0 ,0d0) * exp(r*conjg(Ylm(PI/4d0, PI/4d0, 2, 1)*Ylm(theta, phi, 2, 1)))
    end function g
end program main