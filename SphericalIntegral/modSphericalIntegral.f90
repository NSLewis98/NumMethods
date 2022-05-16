module SphericalIntegral
implicit none
    interface
        real(kind(1d0)) function func(r,theta,phi, parvec) result(res)
            real(kind(1d0)), optional, intent(in)  :: r, theta, phi
            real(kind(1d0)), optional, intent(in)  :: parvec(*)        
        end function func
    end interface 
    public func

    interface
        real(kind(1d0)) function funcXY(x,y, parvec) result(res)
            real(kind(1d0)), optional, intent(in)  :: x, y
            real(kind(1d0)), optional, intent(in)  :: parvec(*)        
        end function funcXY
    end interface 
    public funcXY
    contains
    real(kind(1d0)) function OneDimIntegral(fPtr, upperBound, lowerBound, numPoints, parvec) result(Integral)
        procedure(func) , pointer               :: fPtr
        real(kind(1d0)) , intent(in)            :: upperBound, lowerBound
        integer         , intent(in)            :: numPoints
        real(kind(1d0)) , optional, intent(in)  :: parvec(*)

        real(kind(1d0)) :: xi
        integer :: i

        Integral = 0.5d0 *( fPtr(r=upperBound, parvec=parvec) + fPtr(r=lowerBound, parvec=parvec) )

        do i = 1, numPoints - 1
            xi = (upperBound - lowerBound) / dble(numPoints) * dble(i) + lowerBound
            Integral = Integral + fPtr(r=xi, parvec=parvec)
        enddo
        Integral = Integral * (upperBound - lowerBound) / (numPoints)
    end function OneDimIntegral

    real(kind(1d0)) function TwoDimIntegral(fPtrXY, upperBound, lowerBound, numPoints, parvec) result(Integral)
        procedure(funcXY) , pointer             :: fPtrXY
        real(kind(1d0)) , intent(in)            :: upperBound(2), lowerBound(2)
        integer         , intent(in)            :: numPoints(2)
        real(kind(1d0)) , optional, intent(in)  :: parvec(*)

        DoublePrecision, parameter :: PI = 3.14159265358979323846d0

        real(kind(1d0)) :: xi, yj, xInt, yInt, rNumPoints(2), wy, wx
        integer :: i, j

        rNumPoints = DBLE(numPoints)

        xInt = upperBound(1) - lowerBound(1)
        yInt = upperBound(2) - lowerBound(2)

        print*, xInt*yInt / (rNumPoints(1)*rNumPoints(2))

        do j = 0, numPoints(2)
            yj = yInt / dble(rNumPoints(2)) * dble(j) + lowerBound(2)
            wy = weights(j, 0, numPoints(2) )
            do i = 0, numPoints(1)
                xi = xInt/ dble(rNumPoints(1)) * dble(i) + lowerBound(1)
                wx = weights(i, 0, numPoints(1))
                Integral = Integral + fPtrXY(x=xi, y=yj, parvec=parvec)*wy*wx
            enddo
        enddo
        Integral = Integral * xInt*yInt / (rNumPoints(1)*rNumPoints(2))
    end function TwoDimIntegral

    real(kind(1d0)) function TwoDimPolarInt(fPtrRPhi, upperBound, lowerBound, numPoints, parvec) result(Integral)
        procedure(func) , pointer             :: fPtrRPhi
        real(kind(1d0)) , intent(in)            :: upperBound(2), lowerBound(2)
        integer         , intent(in)            :: numPoints(2)
        real(kind(1d0)) , optional, intent(in)  :: parvec(*)

        DoublePrecision, parameter :: PI = 3.14159265358979323846d0

        real(kind(1d0)) :: ri, phij, rInt, phiInt, rNumPoints(2), wPhi, wR
        integer :: i, j

        rNumPoints = DBLE(numPoints)

        rInt = upperBound(1) - lowerBound(1)
        phiInt = upperBound(2) - lowerBound(2)
        Integral = 0d0

        do j = 0, numPoints(2)
            phij = phiInt / dble(rNumPoints(2)) * dble(j) + lowerBound(2)
            wPhi = weights(j, 0, numPoints(2) )
            do i = 0, numPoints(1)
                ri = rInt/ dble(rNumPoints(1)) * dble(i) + lowerBound(1)
                wR = weights(i, 0, numPoints(1))
                Integral = Integral + ri*fPtrRPhi(r=ri, phi=phij, parvec=parvec)*wPhi*wR
            enddo
        enddo
        Integral = Integral * rInt*phiInt / (rNumPoints(1)*rNumPoints(2))
    end function TwoDimPolarInt

    real(kind(1d0)) function ThreeDimCartIntegral(fPtrXYZ, upperBound, lowerBound, numPoints, parvec) result(Integral)
        procedure(func),           pointer    :: fPtrXYZ
        real(kind(1d0)),           intent(in) :: upperBound(3), lowerBound(3)
        integer        ,           intent(in) :: numPoints(3)
        real(kind(1d0)), optional, intent(in) :: parvec(*)

        real(kind(1d0)) :: rNumPoints(3), intrvl(3), xi, yj, zk, wx, wy, wz, step_size(3)
        real(kind(1d0)) :: xInt, yInt, zInt
        integer         :: i, j, k

        rNumPoints = DBLE(numPoints)
        intrvl = ( upperBound - lowerBound )
        step_size = intrvl / rNumPoints

        Integral = 0d0
        
        do k =0,  numPoints(3)
            zk = step_size(3) * k + lowerBound(3)
            wz = weights(k, 0, numPoints(3))
            do j =0, numPoints(2)
                yj = step_size(2) * j + lowerBound(2)
                wy = weights(k, 0, numPoints(2))
                do i = 0, numPoints(1)
                    xi = step_size(1) * i + lowerBound(1)
                    wx = weights(k, 0, numPoints(1))
                    Integral = Integral + fPtrXYZ(r=xi, theta=yj, phi=zk)*wz*wy*wx
                enddo
            enddo
        enddo
        Integral = Integral * intrvl(1)*intrvl(2)*intrvl(3)&
                /( rNumPoints(3)*rNumPoints(2)*rNumPoints(1))

    end function ThreeDimCartIntegral

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
            wPhi = weights(k, 0, numPoints(3)) ! Dont want to count this boundary..
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

    real(kind(1d0)) function weights(n, N_MIN, N_MAX) result(res)
        integer, intent(in) :: n, N_MIN, N_MAX
        res = 1d0
        if( (n.eq.N_MIN).or.(n.eq.N_MAX)) then
            res = 0.5d0
            return
        endif
    end function weights
end module SphericalIntegral

program main
    use SphericalIntegral
    implicit none

    DoublePrecision, parameter :: PI = 3.14159265358979323846d0

    real(kind(1d0)), parameter :: UPBOUND = 0d0
    real(kind(1d0)), parameter :: LWBOUND = -PI
    integer        , parameter :: NUMPTS  = 200
    procedure(func), pointer :: funcPtr
    procedure(funcXY), pointer :: fXYPtr
    

    real(kind(1d0)) :: IntegralOne , lowerbounds(3), upperbounds(3), step_size
    integer :: nmpts(3), i

    funcPtr => SphFunc

    lowerbounds = 0d0
    upperBounds(1) = 10d0
    upperBounds(2) = PI
    upperBounds(3) = 2d0 * PI
    nmpts = 2000
    step_size = (upperBounds(1) - lowerBounds(1) ) / nmpts(1)
    print*, "Step Size: ", step_size

    IntegralOne = ThreeDimSphIntegral(funcPtr, upperBounds, lowerBounds, nmpts)

    write(*,*) "Approx: ", IntegralOne , nmpts(:)

    contains
    real(kind(1d0)) function SphFunc(r, theta, phi, parvec) result(res)
        real(kind(1d0)), optional, intent(in)  :: r, theta, phi
        real(kind(1d0)), optional, intent(in)  :: parvec(*)
        res = PI**(-3d0/2d0) * EXP(-r**2)
    end function SphFunc
end program main