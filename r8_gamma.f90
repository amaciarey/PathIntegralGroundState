function r8_gamma ( x )

!*****************************************************************************80
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This function was originally named DGAMMA.
!
!    However, a number of Fortran compilers now include a library
!    function of this name.  To avoid conflicts, this function was
!    renamed R8_GAMMA.
!
!    This routine calculates the GAMMA function for a real argument X.
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the GAMMA
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Modified:
!
!    18 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_GAMMA, the value of the function.
!
  implicit none
!
!  Coefficients for minimax approximation over (12, INF).
!
  real    ( kind = 8 ), dimension ( 7 ) :: c = (/ &
   -1.910444077728D-03, &
    8.4171387781295D-04, &
   -5.952379913043012D-04, &
    7.93650793500350248D-04, &
   -2.777777777777681622553D-03, &
    8.333333333333333331554247D-02, &
    5.7083835261D-03 /)
  real    ( kind = 8 ) eps
  real    ( kind = 8 ) fact
  real    ( kind = 8 ) half
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real    ( kind = 8 ) one
  real    ( kind = 8 ) p(8)
  logical parity
  real    ( kind = 8 ) pi
  real    ( kind = 8 ) q(8)
  real    ( kind = 8 ) r8_gamma
  real    ( kind = 8 ) res
  real    ( kind = 8 ) sqrtpi
  real    ( kind = 8 ) sum
  real    ( kind = 8 ) twelve
  real    ( kind = 8 ) two
  real    ( kind = 8 ) x
  real    ( kind = 8 ) xbig
  real    ( kind = 8 ) xden
  real    ( kind = 8 ) xinf
  real    ( kind = 8 ) xminin
  real    ( kind = 8 ) xnum
  real    ( kind = 8 ) y
  real    ( kind = 8 ) y1
  real    ( kind = 8 ) ysq
  real    ( kind = 8 ) z
  real    ( kind = 8 ) zero
!
!  Mathematical constants
!
  data one /1.0D+00 /
  data half /0.5D+00/
  data twelve /12.0D+00/
  data two /2.0D+00 /
  data zero /0.0D+00/
  data sqrtpi /0.9189385332046727417803297D+00/
  data pi /3.1415926535897932384626434D+00/
!
!  Machine dependent parameters
!
  data xbig / 171.624D+00 /
  data xminin / 2.23D-308 /
  data eps / 2.22D-16 /
  data xinf /1.79D+308/
!
!  Numerator and denominator coefficients for rational minimax
!  approximation over (1,2).
!
  data p / -1.71618513886549492533811D+0,   2.47656508055759199108314D+01, &
           -3.79804256470945635097577D+02,  6.29331155312818442661052D+02, &
            8.66966202790413211295064D+02, -3.14512729688483675254357D+04, &
           -3.61444134186911729807069D+04,  6.64561438202405440627855D+04 /

  data q / -3.08402300119738975254353D+01,  3.15350626979604161529144D+02, &
           -1.01515636749021914166146D+03, -3.10777167157231109440444D+03, &
            2.25381184209801510330112D+04,  4.75584627752788110767815D+03, &
           -1.34659959864969306392456D+05, -1.15132259675553483497211D+05 /

  parity = .false.
  fact = one
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= zero ) then

    y = - x
    y1 = aint ( y )
    res = y - y1

    if ( res /= zero ) then

      if ( y1 /= aint ( y1 * half ) * two ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * res )
      y = y + one

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Argument is positive.
!
  if ( y < eps ) then
!
!  Argument < EPS.
!
    if ( xminin <= y ) then
      res = one / y
    else
      res = xinf
      r8_gamma = res
      return
    end if

  else if ( y < twelve ) then

    y1 = y
!
!  0.0 < argument < 1.0.
!
    if ( y < one ) then

      z = y
      y = y + one
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
    else

      n = int ( y ) - 1
      y = y - real ( n, kind = 8 )
      z = y - one

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = zero
    xden = one
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    res = xnum / xden + one
!
!  Adjust result for case  0.0 < argument < 1.0.
!
    if ( y1 < y ) then

      res = res / y1
!
!  Adjust result for case 2.0 < argument < 12.0.
!
    else if ( y < y1 ) then

      do i = 1, n
        res = res * y
        y = y + one
      end do

    end if

  else
!
!  Evaluate for 12.0 <= argument.
!
    if ( y <= xbig ) then

      ysq = y * y
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum / y - y + sqrtpi
      sum = sum + ( y - half ) * log ( y )
      res = exp ( sum )

    else

      res = huge ( res )
      r8_gamma = res
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    res = - res
  end if

  if ( fact /= one ) then
    res = fact / res
  end if

  r8_gamma = res

  return

end function r8_gamma

