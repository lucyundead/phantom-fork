program test
implicit none

real*8 x,cc
real*8 p(3)

x = 4.3

call eif(x,p)

write(*,*) p


end program

subroutine eif(x,ei)

!*****************************************************************************80
!
!! EIX computes the exponential integral Ei(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.
!    However, 
!    they give permission to incorporate this routine into a user
!    progra`m 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    10 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) EI, the function value.
!
  implicit none

  real ( kind = 8 ) ei(3)
  real ( kind = 8 ) ga
  integer ( kind = 4 ) k
  real ( kind = 8 ) r
  real ( kind = 8 ), intent(in) :: x

  if ( x == 0.0D+00 ) then

    ei = -1.0D+300

  else if ( x <= 40.0D+00 ) then

    ei = 1.0D+00
    r = 1.0D+00
    do k = 1, 100
      r = r * k * x / ( k + 1.0D+00 )**2
      ei = ei + r
      if ( abs ( r / ei(1) ) <= 1.0D-15 ) then
        exit
      end if
    end do

    ga = 0.5772156649015328D+00
    ei = ga + log ( x ) + x * ei

  else

    ei = 1.0D+00
    r = 1.0D+00
    do k = 1, 20
      r = r * k / x
      ei = ei + r
    end do
    ei(1) = exp ( x ) / x * ei(1)
    ei(2) = 0.
    ei(3) = 1.
  end if

  return 

end
