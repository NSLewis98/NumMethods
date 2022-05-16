module SphericalIntegral
    
implicit none
    DoublePrecision, private, parameter :: PI    = 3.14159265358979323846d0
    interface
        real(kind(1d0)) function func(r,theta,phi, parvec) result(res)
            real(kind(1d0)), optional, intent(in)  :: r, theta, phi
            real(kind(1d0)), optional, intent(in)  :: parvec(*)    
        end function func
    end interface 
    public func

    interface 
        complex(kind(1d0)) function imgFunc(r, theta, phi, parvec) result(res)
            real(kind(1d0)), optional, intent(in) :: r, theta, phi
            real(kind(1d0)), optional, intent(in) :: parvec(*)
        end function imgFunc
    end interface

    public :: ThreeDimSphIntegral, rThreeDimSphIntegral, imgrThreeDimSphIntegral
    private :: weights

    contains

    ! Computes the 3D Integral in Spherical Coordinates using the trapezoidal rule
    real(kind(1d0)) function ThreeDimSphIntegral(fPtrRThPh, upperBound, lowerBound, numPoints, parvec ) result(Integral)
        procedure(func),           pointer    :: fPtrRThPh
        real(kind(1d0)),           intent(in) :: upperBound(3), lowerBound(3)
        integer        ,           intent(in) :: numPoints(3)
        real(kind(1d0)), optional, intent(in) :: parvec(*)

        real(kind(1d0)) :: rNumPoints(3), intrvl(3), ri, thetaj, phik, wR, wTheta, wPhi, step_size(3)
        integer         :: i, j, k

        rNumPoints = DBLE(numPoints)
        intrvl = ( upperBound - lowerBound )
        step_size = intrvl / rNumPoints
        Integral = 0d0

        do k =0, numPoints(3)
            phik = step_size(3) * k + lowerBound(3)
            wPhi = weights(k, 0, numPoints(3))
            do j =0, numPoints(2)
                thetaj = step_size(2) * j + lowerBound(2)
                wTheta = weights(k, 0, numPoints(2))
                do i = 0, numPoints(1)
                    ri = step_size(1) * i + lowerBound(1)
                    wR = weights(k, 0, numPoints(1))
                    Integral = Integral + ri**2 * sin(thetaj)*fPtrRThPh(r=ri,theta=thetaj,phi=phik)&
                                *wPhi*wTheta*wR
                enddo
            enddo
        enddo
        Integral = Integral * intrvl(1)*intrvl(2)*intrvl(3)&
                /( rNumPoints(3)*rNumPoints(2)*rNumPoints(1))
    end function ThreeDimSphIntegral


    ! Integrates over unit sphere and radially from 0 to r_max.
    ! The accepted function is of f | R^3 -> R
    real(kind(1d0)) function rThreeDimSphIntegral(funcPtr, r_max, numPoints, parvec) result(Integral)
        procedure(func),           pointer    :: funcPtr
        real(kind(1d0)),           intent(in) :: r_max
        integer        ,           intent(in) :: numPoints(3)
        real(kind(1d0)), optional, intent(in) :: parvec(*)

        real(kind(1d0)) :: rNumPoints(3), rStep, phik, wPhi,thetaj, wTheta, ri, wR
        integer         :: i, j, k

        rNumPoints = DBLE(numPoints)
        rStep = r_max / (rNumPoints(1) - 1d0)
        do k =1, numPoints(3)
            phik = 2d0*PI/rNumPoints(3) *DBLE(k -1)
            wPhi = weights(k, 1, numPoints(3)) ! Dont want to count this boundary..
            do j =1, numPoints(2)
                thetaj = PI / (rNumPoints(2) - 1d0) *DBLE(j-1)
                wTheta = weights(k, 1, numPoints(2))
                do i = 1, numPoints(1)
                    ri = r_max / (rNumPoints(1) - 1d0) * DBLE(i-1)
                    wR = weights(k, 1, numPoints(1))
                    Integral = Integral + ri**2 * sin(thetaj)*funcPtr(r=ri,theta=thetaj,phi=phik, parvec=parvec)&
                                *wPhi*wTheta*wR
                enddo
            enddo
        enddo
        Integral = Integral * 2d0 * PI**2 * r_max &
        /((rNumPoints(1) - 1)*(rNumPoints(2)-1)*rNumPoints(1) )
    end function rThreeDimSphIntegral

    ! Integrates over unit sphere and radially from 0 to r_max.
    ! The accepted function is of f :| R^3 -> C
    complex(kind(1d0)) function imgrThreeDimSphIntegral(funcPtr, r_max, numPoints, parvec) result(Integral)
        procedure(imgFunc),        pointer    :: funcPtr
        real(kind(1d0)),           intent(in) :: r_max
        integer        ,           intent(in) :: numPoints(3)
        real(kind(1d0)), optional, intent(in) :: parvec(*)

        real(kind(1d0)) :: rNumPoints(3), rStep, phik, wPhi,thetaj, wTheta, ri, wR
        integer         :: i, j, k

        rNumPoints = DBLE(numPoints)
        rStep = r_max / (rNumPoints(1) - 1d0)
        Integral = (0d0, 0d0)

        do k =1, numPoints(3)
            phik = 2d0*PI/rNumPoints(3) *DBLE(k -1)
            ! wPhi = 1d0 ! Get rid of this weight
            do j =1, numPoints(2)
                thetaj = PI / (rNumPoints(2) - 1d0) *DBLE(j-1)
                wTheta = weights(j, 1, numPoints(2))
                do i = 1, numPoints(1)
                    ri = r_max / (rNumPoints(1) - 1d0) * DBLE(i-1)
                    wR = weights(i, 1, numPoints(1))
                    Integral = Integral + ri**2 * sin(thetaj)*funcPtr(r=ri,theta=thetaj,phi=phik, parvec=parvec)&
                                *wTheta*wR
                    ! print*, "Integral: ", REAL(Integral) , ri, thetaj, phik
                enddo
            enddo
        enddo
        Integral = Integral * 2d0 * PI**2 * r_max &
        /((rNumPoints(1) - 1)*(rNumPoints(2)-1)*rNumPoints(1) )
        return
    end function imgrThreeDimSphIntegral

    real(kind(1d0)) function weights(n, N_MIN, N_MAX) result(res)
        integer, intent(in) :: n, N_MIN, N_MAX
        res = 1d0
        if( (n.eq.N_MIN).or.(n.eq.N_MAX)) then
            res = 0.5d0
            return
        endif
    end function weights

end module SphericalIntegral