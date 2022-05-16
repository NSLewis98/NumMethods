!.. The Lebdev method seems to be working for the test functions.
!.. The test functions used have been linear combinations of Spherical Harmonics
!.. By integrating over \mathcad{S}^2 a function of a sum of spherical harmonics will always equal 0
!.. unless the Y_{00} harmonic is included

program LebedevQuadrature
    use ModuleAngularMomentum
    use SphericalIntegral
    implicit none

    real(kind(1d0)), parameter :: PI = 3.14159265358979323844d0
    integer        , parameter :: NUM_PTS = 5810
    integer        , parameter :: N = INT(SQRT(DBLE(NUM_PTS)))


    real(kind(1d0)), dimension(NUM_PTS) :: weights, theta, phi
    real(kind(1d0)), dimension(N)       :: alpha, beta
    complex(kind(1d0)) :: integral

    call getLebPoints(NUM_PTS, weights, theta, phi)

    integral = LebdevIntegral(NUM_PTS, weights, theta, phi)
    print*, "Lebdev Integral: "
    print*, integral

    integral = (0d0 ,0d0)

    write(*,*)

    integral = TrapRuleIntegral(N)
    print*, "Trapezoidal Integral: "
    print*, integral

    contains

    complex(kind(1d0)) function LebdevIntegral(numPts, weights, theta, phi) result(res)
        integer                           , intent(in) :: numPts
        real(kind(1d0)), dimension(numPts), intent(in) :: weights, theta, phi

        integer :: i, j
        res = (0d0, 0d0)
        do i = 1, numPts
            res = res + weights(i)*funcYlm(theta(i), phi(i))
        enddo
        res =4d0 * PI*res
    end function LebdevIntegral

    complex(kind(1d0)) function TrapRuleIntegral(numPts) result(res)
        integer                           , intent(in) :: numPts
        real(kind(1d0)), dimension(numPts) :: theta, phi
        real(kind(1d0)) :: wTheta, wPhi
        integer :: i, j
        res = (0d0, 0d0)
        call fillPoints(numPts, theta, phi)
        do i = 1, numPts
            wTheta = 0.5d0
            if((i.eq.1).or.(i.eq.numPts)) wTheta = 1
            do j = 1, numPts
                res = res + wTheta  * funcYlm(theta(i), phi(i) )
            enddo
        enddo
       res = 2d0 * PI**2 / (DBLE (numPts)*DBLE(numPts -1) ) * res
    end function TrapRuleIntegral

    subroutine fillPoints(numPts, theta, phi)
        integer                           , intent(in)  :: numPts
        real(kind(1d0)), dimension(numPts), intent(out) :: theta, phi

        integer :: i

        theta = 0d0
        phi   = 0d0

        do i = 1, numPts
            theta(i) = PI / dble(numPts - 1) * dble(i-1)
            phi(i)   = 2d0*PI / DBLE(numPts) * dble(i - 1)
        enddo 
        return

    end subroutine fillPoints

    subroutine getLebPoints(numPts, weights, theta, phi)
        integer                           , intent(in)  :: numPts
        real(kind(1d0)), dimension(numPts), intent(out) :: weights, theta, phi

        integer          :: error, i
        character(len=5) :: x1
        character(len=8) :: fmt ! format descriptor

        fmt = '(I4.4)' ! an integer of width 5 with zeros at the left

        write (x1,fmt) numPts ! converting integer to string using a 'internal file'

        weights = 0d0 
        theta   = 0d0 
        phi     = 0d0 
        write(*,*) "Reading File: ", 'Leb'//trim(x1)//'pt.txt'
        open( 10, file = 'Leb'//trim(x1)//'pt.txt', status='old')
        do i = 1, numPts
            read( 10, *) phi(i), theta(i), weights(i)
        end do
        theta = PI / 180d0 * theta 
        phi   = PI / 180d0 * phi + PI/2d0 
        close(10)
        return
    end subroutine getLebPoints
    ! f = \sum_{\ell, m} Y_{\ell, m}
    ! The integral of this function should be zero for all orders except Y_{00}
    complex(kind(1d0)) function funcYlm(theta, phi) result(res)
        real(kind(1d0)), intent(in) :: theta, phi
        integer :: l, m, i
        res = (0d0, 0d0)
        ! do l = 1,10
        !     do m = -l, l
        !         res = res + Ylm(theta, phi, l, m)
        !     enddo
        !     ! res = res + Ylm(theta, phi, l, 0)
        ! enddo

        res = COS(theta)**2d0 / (4d0*PI)
        
    end function funcYlm
end program LebedevQuadrature