!DIR$ EXPORT pbdv
module fortfuncs
    use, intrinsic :: iso_c_binding
    implicit none
    private
    public :: gamma, dvsa, dvla, vvla, pbdv, chgm

contains

    subroutine gamma ( x, ga ) bind(C, name = "gamma_")

    !*****************************************************************************80
    !
    !! GAMMA evaluates the Gamma function.
    !
    !  Licensing:
    !
    !    The original FORTRAN77 version of this routine is copyrighted by
    !    Shanjie Zhang and Jianming Jin.  However, they give permission to
    !    incorporate this routine into a user program that the copyright
    !    is acknowledged.
    !
    !  Modified:
    !
    !    08 September 2007
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Shanjie Zhang, Jianming Jin,
    !    Computation of Special Functions,
    !    Wiley, 1996,
    !    ISBN: 0-471-11963-6,
    !    LC: QA351.C45
    !
    !  Parameters:
    !
    !    Input, real ( kind = rk ) X, the argument.
    !    X must not be 0, or any negative integer.
    !
    !    Output, real ( kind = rk ) GA, the value of the Gamma function.
    !
      integer, parameter :: rk = kind ( 1.0D+00 )

      real ( kind = rk ), dimension ( 26 ) :: g = (/ &
        1.0D+00, &
        0.5772156649015329D+00, &
       -0.6558780715202538D+00, &
       -0.420026350340952D-01, &
        0.1665386113822915D+00, &
       -0.421977345555443D-01, &
       -0.96219715278770D-02, &
        0.72189432466630D-02, &
       -0.11651675918591D-02, &
       -0.2152416741149D-03, &
        0.1280502823882D-03, &
       -0.201348547807D-04, &
       -0.12504934821D-05, &
        0.11330272320D-05, &
       -0.2056338417D-06, &
        0.61160950D-08, &
        0.50020075D-08, &
       -0.11812746D-08, &
        0.1043427D-09, &
        0.77823D-11, &
       -0.36968D-11, &
        0.51D-12, &
       -0.206D-13, &
       -0.54D-14, &
        0.14D-14, &
        0.1D-15 /)
      real(c_double) :: ga
      real ( kind = rk ) gr
      integer k
      integer m
      integer m1
      real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
      real ( kind = rk ) r
      real(c_double) :: x
      real ( kind = rk ) z

      if ( x == aint ( x ) ) then

        if ( 0.0D+00 < x ) then
          ga = 1.0D+00
          m1 = int ( x ) - 1
          do k = 2, m1
            ga = ga * k
          end do
        else
          ga = 1.0D+300
        end if

      else

        if ( 1.0D+00 < abs ( x ) ) then
          z = abs ( x )
          m = int ( z )
          r = 1.0D+00
          do k = 1, m
            r = r * ( z - real ( k, kind = rk ) )
          end do
          z = z - real ( m, kind = rk )
        else
          z = x
        end if

        gr = g(26)
        do k = 25, 1, -1
          gr = gr * z + g(k)
        end do

        ga = 1.0D+00 / ( gr * z )

        if ( 1.0D+00 < abs ( x ) ) then
          ga = ga * r
          if ( x < 0.0D+00 ) then
            ga = - pi / ( x* ga * sin ( pi * x ) )
          end if
        end if

      end if
    end subroutine gamma

    subroutine dvsa ( va, x, pd ) bind(C, name = "dvsa_")

    !*****************************************************************************80
    !
    !! DVSA computes parabolic cylinder functions Dv(x) for small argument.
    !
    !  Licensing:
    !
    !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
    !    they give permission to incorporate this routine into a user program
    !    provided that the copyright is acknowledged.
    !
    !  Modified:
    !
    !    07 July 2012
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
    !    Input, real ( kind = rk ) VA, the order.
    !
    !    Input, real ( kind = rk ) X, the argument.
    !
    !    Output, real ( kind = rk ) PD, the function value.
    !
      integer, parameter :: rk = kind ( 1.0D+00 )

      real ( kind = rk ) a0
      real ( kind = rk ) ep
      real ( kind = rk ) eps
      real ( kind = rk ) g0
      real ( kind = rk ) g1
      real ( kind = rk ) ga0
      real ( kind = rk ) gm
      integer m
      real(c_double) :: pd
      real ( kind = rk ) pi
      real ( kind = rk ) r
      real ( kind = rk ) r1
      real ( kind = rk ) sq2
      real(c_double) :: va
      real ( kind = rk ) va0
      real ( kind = rk ) vm
      real ( kind = rk ) vt
      real(c_double) :: x

      eps = 1.0D-15
      pi = 3.141592653589793D+00
      sq2 = sqrt ( 2.0D+00 )
      ep = exp ( -0.25D+00 * x * x )
      va0 = 0.5D+00 * ( 1.0D+00 - va )

      if ( va == 0.0D+00 ) then

        pd = ep

      else

        if ( x == 0.0D+00 ) then
          if ( va0 <= 0.0D+00 .and. va0 == int ( va0 ) ) then
            pd = 0.0D+00
          else
            call gamma ( va0, ga0 )
            pd = sqrt ( pi ) / ( 2.0D+00 ** ( -0.5D+00 * va ) * ga0 )
          end if

        else

          call gamma ( -va, g1 )
          a0 = 2.0D+00 ** ( -0.5D+00 * va - 1.0D+00 ) * ep / g1
          vt = -0.5D+00 * va
          call gamma ( vt, g0 )
          pd = g0
          r = 1.0D+00
          do m = 1, 250
            vm = 0.5D+00 * ( m - va )
            call gamma ( vm, gm )
            r = -r * sq2 * x / m
            r1 = gm * r
            pd = pd + r1
            if ( abs ( r1 ) < abs ( pd ) * eps ) then
              exit
            end if
          end do

          pd = a0 * pd

        end if

      end if
    end subroutine dvsa

    subroutine dvla ( va, x, pd ) bind(C, name = "dvla_")

    !*****************************************************************************80
    !
    !! DVLA computes parabolic cylinder functions Dv(x) for large argument.
    !
    !  Licensing:
    !
    !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
    !    they give permission to incorporate this routine into a user program
    !    provided that the copyright is acknowledged.
    !
    !  Modified:
    !
    !    06 July 2012
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
    !    Input, real ( kind = rk ) X, the argument.
    !
    !    Input, real ( kind = rk ) VA, the order.
    !
    !    Output, real ( kind = rk ) PD, the function value.
    !
      integer, parameter :: rk = kind ( 1.0D+00 )

      real ( kind = rk ) a0
      real ( kind = rk ) ep
      real ( kind = rk ) eps
      real ( kind = rk ) gl
      integer k
      real(c_double) :: pd
      real ( kind = rk ) pi
      real ( kind = rk ) r
      real(c_double) :: va
      real ( kind = rk ) vl
      real(c_double) :: x
      real ( kind = rk ) x1

      pi = 3.141592653589793D+00
      eps = 1.0D-12
      ep = exp ( -0.25D+00 * x * x )
      a0 = abs ( x ) ** va * ep
      r = 1.0D+00
      pd = 1.0D+00
      do k = 1, 16
        r = -0.5D+00 * r * ( 2.0D+00 * k - va - 1.0D+00 ) &
          * ( 2.0D+00 * k - va - 2.0D+00 ) / ( k * x * x )
        pd = pd + r
        if ( abs ( r / pd ) < eps ) then
          exit
        end if
      end do

      pd = a0 * pd

      if ( x < 0.0D+00 ) then
        x1 = - x
        call vvla ( va, x1, vl )
        call gamma ( -va, gl )
        pd = pi * vl / gl + cos ( pi * va ) * pd
      end if
    end subroutine dvla

    subroutine vvla ( va, x, pv ) bind(C, name = "vvla_")

    !*****************************************************************************80
    !
    !! VVLA computes parabolic cylinder function Vv(x) for large arguments.
    !
    !  Licensing:
    !
    !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
    !    they give permission to incorporate this routine into a user program
    !    provided that the copyright is acknowledged.
    !
    !  Modified:
    !
    !    04 July 2012
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
    !    Input, real ( kind = rk ) X, the argument.
    !
    !    Input, real ( kind = rk ) VA, the order nu.
    !
    !    Output, real ( kind = rk ) PV, the value of V(nu,x).
    !
      integer, parameter :: rk = kind ( 1.0D+00 )

      real ( kind = rk ) a0
      real ( kind = rk ) dsl
      real ( kind = rk ) eps
      real ( kind = rk ) gl
      integer k
      real ( kind = rk ) pdl
      real ( kind = rk ) pi
      real(c_double) :: pv
      real ( kind = rk ) qe
      real ( kind = rk ) r
      real(c_double) :: va
      real(c_double) :: x
      real ( kind = rk ) x1

      pi = 3.141592653589793D+00
      eps = 1.0D-12
      qe = exp ( 0.25D+00 * x * x )
      a0 = abs ( x ) ** ( -va - 1.0D+00 ) * sqrt ( 2.0D+00 / pi ) * qe

      r = 1.0D+00
      pv = 1.0D+00
      do k = 1, 18
        r = 0.5D+00 * r * ( 2.0D+00 * k + va - 1.0D+00 ) &
          * ( 2.0D+00 * k + va ) / ( k * x * x )
        pv = pv + r
        if ( abs ( r / pv ) < eps ) then
          exit
        end if
      end do

      pv = a0 * pv

      if ( x < 0.0D+00 ) then
        x1 = -x
        call dvla ( va, x1, pdl )
        call gamma ( -va, gl )
        dsl = sin ( pi * va ) * sin ( pi * va )
        pv = dsl * gl / pi * pdl - cos ( pi * va ) * pv
      end if
    end subroutine vvla

    subroutine pbdv ( v, x, dv, dp, pdf, pdd ) bind(C, name = "pbdv_")

    !*****************************************************************************80
    !
    !! PBDV computes parabolic cylinder functions Dv(x) and derivatives.
    !
    !  Licensing:
    !
    !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
    !    they give permission to incorporate this routine into a user program
    !    provided that the copyright is acknowledged.
    !
    !  Modified:
    !
    !    29 July 2012
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
    !    Input, real ( kind = rk ) V, the order.
    !
    !    Input, real ( kind = rk ) X, the argument.
    !
    !    Output, real ( kind = rk ) DV(0:*), DP(0:*), the values of
    !    Dn+v0(x), Dn+v0'(x).
    !
    !    Output, real ( kind = rk ) PDF, PDD, the values of Dv(x) and Dv'(x).
    !
      integer, parameter :: rk = kind ( 1.0D+00 )

      real(c_double) :: dp(0:*)
      real(c_double) :: dv(0:*)
      real ( kind = rk ) ep
      real ( kind = rk ) f
      real ( kind = rk ) f0
      real ( kind = rk ) f1
      integer ja
      integer k
      integer l
      integer m
      integer na
      integer nk
      integer nv
      real ( kind = rk ) pd
      real ( kind = rk ) pd0
      real ( kind = rk ) pd1
      real(c_double) :: pdd
      real(c_double) :: pdf
      real ( kind = rk ) s0
      real(c_double) :: v
      real ( kind = rk ) v0
      real ( kind = rk ) v1
      real ( kind = rk ) v2
      real ( kind = rk ) vh
      real(c_double) :: x
      real ( kind = rk ) xa

      xa = abs ( x )
      vh = v
      v = v + sign ( 1.0D+00, v )
      nv = int ( v )
      v0 = v - nv
      na = abs ( nv )
      ep = exp ( -0.25D+00 * x * x )

      if ( 1 <= na ) then
        ja = 1
      end if

      if ( 0.0D+00 <= v ) then
        if ( v0 == 0.0D+00 ) then
          pd0 = ep
          pd1 = x * ep
        else
          do l = 0, ja
            v1 = v0 + l
            if ( xa <= 5.8D+00 ) then
              call dvsa ( v1, x, pd1 )
            else
              call dvla ( v1, x, pd1 )
            end if
            if ( l == 0 ) then
              pd0 = pd1
            end if
          end do
        end if

        dv(0) = pd0
        dv(1) = pd1
        do k = 2, na
          pdf = x * pd1 - ( k + v0 - 1.0D+00 ) * pd0
          dv(k) = pdf
          pd0 = pd1
          pd1 = pdf
        end do

      else

        if ( x <= 0.0D+00 ) then

          if ( xa <= 5.8D+00 )  then
            call dvsa ( v0, x, pd0 )
            v1 = v0 - 1.0D+00
            call dvsa ( v1, x, pd1 )
          else
            call dvla ( v0, x, pd0 )
            v1 = v0 - 1.0D+00
            call dvla ( v1, x, pd1 )
          end if

          dv(0) = pd0
          dv(1) = pd1
          do k = 2, na
            pd = ( - x * pd1 + pd0 ) / ( k - 1.0D+00 - v0 )
            dv(k) = pd
            pd0 = pd1
            pd1 = pd
          end do

        else if ( x <= 2.0D+00 ) then

          v2 = nv + v0
          if ( nv == 0 ) then
            v2 = v2 - 1.0D+00
          end if

          nk = int ( - v2 )
          call dvsa ( v2, x, f1 )
          v1 = v2 + 1.0D+00
          call dvsa ( v1, x, f0 )
          dv(nk) = f1
          dv(nk-1) = f0
          do k = nk - 2, 0, -1
            f = x * f0 + ( k - v0 + 1.0D+00 ) * f1
            dv(k) = f
            f1 = f0
            f0 = f
          end do

        else

          if ( xa <= 5.8D+00 ) then
            call dvsa ( v0, x, pd0 )
          else
            call dvla ( v0, x, pd0 )
          end if

          dv(0) = pd0
          m = 100 + na
          f1 = 0.0D+00
          f0 = 1.0D-30
          do k = m, 0, -1
            f = x * f0 + ( k - v0 + 1.0D+00 ) * f1
            if ( k <= na ) then
              dv(k) = f
            end if
            f1 = f0
            f0 = f
          end do
          s0 = pd0 / f
          do k = 0, na
            dv(k) = s0 * dv(k)
          end do

        end if

      end if

      do k = 0, na - 1
        v1 = abs ( v0 ) + k
        if ( 0.0D+00 <= v ) then
          dp(k) = 0.5D+00 * x * dv(k) - dv(k+1)
        else
          dp(k) = -0.5D+00 * x * dv(k) - v1 * dv(k+1)
        end if
      end do

      pdf = dv(na-1)
      pdd = dp(na-1)
      v = vh
    end subroutine pbdv

    subroutine chgm ( a, b, x, hg ) bind(C, name = "chgm_")

    !*****************************************************************************80
    !
    !! CHGM computes the confluent hypergeometric function M(a,b,x).
    !
    !  Licensing:
    !
    !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
    !    they give permission to incorporate this routine into a user program
    !    provided that the copyright is acknowledged.
    !
    !  Modified:
    !
    !    27 July 2012
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
    !    Input, real ( kind = rk ) A, B, parameters.
    !
    !    Input, real ( kind = rk ) X, the argument.
    !
    !    Output, real ( kind = rk ) HG, the value of M(a,b,x).
    !
      integer, parameter :: rk = kind ( 1.0D+00 )


      real(c_double) :: a
      real ( kind = rk ) a0
      real ( kind = rk ) a1
      real(c_double) :: b
      real(c_double) :: hg
      real ( kind = rk ) hg1
      real ( kind = rk ) hg2
      integer i
      integer j
      integer k
      integer la
      integer m
      integer n
      integer nl
      real ( kind = rk ) pi
      real ( kind = rk ) r
      real ( kind = rk ) r1
      real ( kind = rk ) r2
      real ( kind = rk ) rg
      real ( kind = rk ) sum1
      real ( kind = rk ) sum2
      real ( kind = rk ) ta
      real ( kind = rk ) tb
      real ( kind = rk ) tba
      real(c_double) :: x
      real ( kind = rk ) x0
      real ( kind = rk ) xg
      real ( kind = rk ) y0
      real ( kind = rk ) y1

      pi = 3.141592653589793D+00
      a0 = a
      a1 = a
      x0 = x
      hg = 0.0D+00

      if ( b == 0.0D+00 .or. b == - abs ( int ( b ) ) ) then
        hg = 1.0D+300
      else if ( a == 0.0D+00 .or. x == 0.0D+00 ) then
        hg = 1.0D+00
      else if ( a == -1.0D+00 ) then
        hg = 1.0D+00 - x / b
      else if ( a == b ) then
        hg = exp ( x )
      else if ( a - b == 1.0D+00 ) then
        hg = ( 1.0D+00 + x / b ) * exp ( x )
      else if ( a == 1.0D+00 .and. b == 2.0D+00 ) then
        hg = ( exp ( x ) - 1.0D+00 ) / x
      else if ( a == int ( a ) .and. a < 0.0D+00 ) then
        m = int ( - a )
        r = 1.0D+00
        hg = 1.0D+00
        do k = 1, m
          r = r * ( a + k - 1.0D+00 ) / k / ( b + k - 1.0D+00 ) * x
          hg = hg + r
        end do
      end if

      if ( hg /= 0.0D+00 ) then
        return
      end if

      if ( x < 0.0D+00 ) then
        a = b - a
        a0 = a
        x = abs ( x )
      end if

      if ( a < 2.0D+00 ) then
        nl = 0
      end if

      if ( 2.0D+00 <= a ) then
        nl = 1
        la = int ( a )
        a = a - la - 1.0D+00
      end if

      do n = 0, nl

        if ( 2.0D+00 <= a0 ) then
          a = a + 1.0D+00
        end if

        if ( x <= 30.0D+00 + abs ( b ) .or. a < 0.0D+00 ) then

          hg = 1.0D+00
          rg = 1.0D+00
          do j = 1, 500
            rg = rg * ( a + j - 1.0D+00 ) &
              / ( j * ( b + j - 1.0D+00 ) ) * x
            hg = hg + rg
            if ( abs ( rg / hg ) < 1.0D-15 ) then
              exit
            end if
          end do

        else

          call gamma ( a, ta )
          call gamma ( b, tb )
          xg = b - a
          call gamma ( xg, tba )
          sum1 = 1.0D+00
          sum2 = 1.0D+00
          r1 = 1.0D+00
          r2 = 1.0D+00
          do i = 1, 8
            r1 = - r1 * ( a + i - 1.0D+00 ) * ( a - b + i ) / ( x * i )
            r2 = - r2 * ( b - a + i - 1.0D+00 ) * ( a - i ) / ( x * i )
            sum1 = sum1 + r1
            sum2 = sum2 + r2
          end do
          hg1 = tb / tba * x ** ( - a ) * cos ( pi * a ) * sum1
          hg2 = tb / ta * exp ( x ) * x ** ( a - b ) * sum2
          hg = hg1 + hg2

        end if

        if ( n == 0 ) then
          y0 = hg
        else if ( n == 1 ) then
          y1 = hg
        end if

      end do

      if ( 2.0D+00 <= a0 ) then
        do i = 1, la - 1
          hg = ( ( 2.0D+00 * a - b + x ) * y1 + ( b - a ) * y0 ) / a
          y0 = y1
          y1 = hg
          a = a + 1.0D+00
        end do
      end if

      if ( x0 < 0.0D+00 ) then
        hg = hg * exp ( x0 )
      end if

      a = a1
      x = x0

      return
    end subroutine chgm

end module fortfuncs
