subroutine l_polynomial ( m, n, x, v )
!*****************************************************************************80
!
!! L_POLYNOMIAL evaluates the Laguerre polynomial L(n,x).
!
!  First terms:
!
!      1
!     -X     +  1
!   (  X^2 -  4 X      +  2 ) / 2
!   ( -X^3 +  9 X^2 -  18 X    +    6 ) / 6
!   (  X^4 - 16 X^3 +  72 X^2 -   96 X +      24 ) / 24
!   ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -   600 X    +  120 ) / 120
!   (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 +  5400 X^2 -  4320 X     + 720 ) 
!     / 720
!   ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3 + 52920 X^2 - 35280 X 
!     + 5040 ) / 5040
!
!  Recursion:
!
!    L(0,X) = 1
!    L(1,X) = 1 - X
!    L(N,X) = (2*N-1-X)/N * L(N-1,X) - (N-1)/N * L(N-2,X)
!
!  Orthogonality:
!
!    Integral ( 0 <= X < oo ) exp ( - X ) * L(N,X) * L(M,X) dX = delta ( M, N )
!
!  Relations:
!
!    L(N,X) = (-1)^N / N! * exp ( x ) * (d/dx)^n ( exp ( - x ) * x^n )  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) V(M,0:N), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ) v(m,0:n)
  real ( kind = 8 ) x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  v(1:m,1) = 1.0D+00 - x(1:m)
 
  do j = 2, n

    v(1:m,j) = ( ( real ( 2 * j - 1, kind = 8 ) - x(1:m) ) * v(1:m,j-1)   &
                 - real (     j - 1, kind = 8 )            * v(1:m,j-2) ) &
                 / real (     j,     kind = 8 )

  end do

  return
end