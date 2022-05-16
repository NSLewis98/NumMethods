! This module uses the recourrence relation to compute the Spherical Bessel Functions of the first
! kind for complex inputs. ie z = ix,  \forall x \in \mathbb{R}
! NOTE: Since the functions are found recursively, there are current numerically instabilites for
!       k > 9. Currently looking into Miller's Algorithm or other algorithms to reduce instabilties
module SphericalBesselOnImgAxis
    implicit none
    
  contains
    
  ! Performs the following Recurrence Relations for Bessel Funcs:
  !... j_{n+1}(z) = -j_{n-1}(z) + \frac{2n+1}{z} j_n (z)
  !... Where z = ix for x \in R
  !.. j_0(ix) = \frac{sinh(x)}{x} \qquad j_1(ix) = \frac{sinh}{ix^2} - \frac{cosh(x)}{ix}
  pure recursive function j_k(k, z) result(res)
    integer,            intent(in) :: k
    complex(kind(1d0)), intent(in) :: z
    complex(kind(1d0)) :: res
    complex(kind(1d0)), parameter :: Zi = (0, 1d0)

    if(k == 0) then 
        res = SINH(abs(z))/ abs(z) 
    else if(k == 1) then
        res =  (-SINH(abs(z)) /abs(z)**2 + COSH(abs(z)) / abs(z) ) / Zi
    else if(k < 0) then
        res = 0
    else
        res = - j_k(k -2, z) + (2*k-1) /z * j_k(k-1, z)
    endif
  end function j_k


  ! Performs the following Recurrence Relations for Bessel Funcs:
  !... j_{n+1}(z) = -j_{n-1}(z) + \frac{2n+1}{z} j_n (z)
  !... Where z = ix for x \in R
  !.. j_0(ix) = \frac{sinh(x)}{x} \qquad j_1(ix) = \frac{sinh}{ix^2} - \frac{cosh(x)}{ix}
  pure function Besselj(k,z) result(res)
    integer,            intent(in) :: k
    complex(kind(1d0)), intent(in) :: z
    complex(kind(1d0)) :: res, zB0, zB1
    integer :: ik
    real(kind(1d0)), parameter :: THRESHOLD = 1.d-12
    complex(kind(1d0)), parameter :: Zi = (0, 1d0)

    res=(0.d0,0.d0)
    if(k<0)return
    
    if(abs(z)<THRESHOLD)then
      res = Besselj_z0limit(k,z)
      return
    endif
    
    zB0 = SINH(AIMAG(z))/ AIMAG(z) 
    if(k==0)then
      res = zB0
      return
    endif
    
    zB1 = Zi * ( COSH( AIMAG(z) )/AIMAG(z) - SINH(AIMAG(z))/ AIMAG(z)**2 )
    if(k==1)then
      res= zB1
      return
    endif
    res = (0d0, 0d0)
    do ik=2,k
      res = -zB0 + (2.d0*ik-1.d0)*zB1/z 
      zB0 = zB1
      zB1 = res
    enddo
  end function Besselj
  

  !.. Computes spherical Bessel functions in the 
  !   neighborhood of the origin
  !..
  pure function Besselj_z0limit(k, z) result(res)
    integer,            intent(in) :: k
    complex(kind(1d0)), intent(in) :: z
    complex(kind(1d0)) :: res
    integer            :: ik
    res=1.d0
    do ik=1,k
        res = res * z / dble(2*ik+1)
    enddo
  end function Besselj_z0limit
end module SphericalBesselOnImgAxis
