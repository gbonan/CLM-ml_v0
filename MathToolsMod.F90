module MathToolsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Math tools
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl, only : iulog
  use abortutils, only : endrun
  use CanopyFluxesMultilayerType, only : mlcanopy_type
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: hybrid             ! Solve for the root of a function using secant and Brent's methods
  public :: zbrent             ! Use Brent's method to find the root of a function
  public :: quadratic          ! Solve a quadratic equation for its two roots
  public :: tridiag            ! Solve a tridiagonal system of equations
  public :: beta_function      ! Evaluate the beta function at p and q: B(p,q)
  public :: log_gamma_function ! Evaluate the log natural of the gamma function at x: ln(G(x))

  interface
    function xfunc (p, iv, il, mlcanopy_inst, x) result(f)
    use shr_kind_mod, only : r8 => shr_kind_r8
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    integer :: p, iv, il
    real(r8) :: x, f
    type(mlcanopy_type) :: mlcanopy_inst
    end function xfunc
  end interface

contains

  !-----------------------------------------------------------------------
  function hybrid (msg, p, iv, il, mlcanopy_inst, func, xa, xb, tol) result(root)
    !
    ! !DESCRIPTION:
    ! Solve for the root of a function, given initial estimates xa and xb.
    ! The root is updated until its accuracy is tol.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*) :: msg           ! String to be printed
    integer, intent(in) :: p          ! pft index for CLM g/l/c/p hierarchy
    integer, intent(in) :: iv         ! Leaf layer index
    integer, intent(in) :: il         ! Sunlit (1) or shaded (2) leaf index
    procedure (xfunc) :: func         ! Function to solve
    real(r8), intent(in) :: xa, xb    ! Initial estimates of root
    real(r8), intent(in) :: tol       ! Error tolerance
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: x0, x1                ! Estimates of root
    real(r8) :: f0, f1                ! Function value for x0 and x1
    real(r8) :: minx                  ! x0 or x1 that gives smallest function value
    real(r8) :: minf                  ! Smallest function value using x0 or x1
    real(r8) :: root                  ! Root

    integer  :: iter                  ! Iteration loop index
    real(r8) :: dx                    ! Change in root
    real(r8) :: x                     ! Updated root
    integer, parameter :: itmax = 40  ! Maximum number of iterations
    !---------------------------------------------------------------------

    x0 = xa
    f0 = func (p, iv, il, mlcanopy_inst, x0)
    if (f0 == 0._r8) then
       root = x0
       return
    end if

    x1 = xb
    f1 = func (p, iv, il, mlcanopy_inst, x1)
    if (f1 == 0._r8) then
       root = x1
       return
    end if

    if (f1 < f0) then
       minx = x1
       minf = f1
    else
       minx = x0
       minf = f0
    end if

    ! First use the secant method, and then use the brent method as a backup

    iter = 0
    do
       iter = iter + 1
       dx = -f1 * (x1 - x0) / (f1 - f0)
       x = x1 + dx
       if (abs(dx) < tol) then
          x0 = x
          exit
       end if
       x0 = x1
       f0 = f1
       x1 = x
       f1 = func (p, iv, il, mlcanopy_inst, x1)
       if (f1 < minf) then
          minx = x1
          minf = f1
       end if

       ! If a root zone is found, use the brent method for a robust backup strategy

       if (f1 * f0 < 0._r8) then
          x = zbrent (msg, p, iv, il, mlcanopy_inst, func, x0, x1, tol)
          x0 = x
          exit
       end if

       ! In case of failing to converge within itmax iterations stop at the minimum function

       if (iter > itmax) then
          f1 = func (p, iv, il, mlcanopy_inst, minx)
          x0 = minx
          exit
       end if

    end do

    root = x0

  end function hybrid

  !-----------------------------------------------------------------------
  function zbrent (msg, ip, iv, il, mlcanopy_inst, func, xa, xb, tol) result(root)
    !
    ! !DESCRIPTION:
    ! Use Brent's method to find the root of a function, which is known to exist
    ! between xa and xb. The root is updated until its accuracy is tol.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*) :: msg           ! String to be printed
    integer, intent(in) :: ip         ! pft index for CLM g/l/c/p hierarchy
    integer, intent(in) :: iv         ! Canopy layer index
    integer, intent(in) :: il         ! Sunlit (1) or shaded (2) leaf index
    procedure (xfunc) :: func         ! Function to solve
    real(r8), intent(in) :: xa, xb    ! Minimum and maximum of the variable domain to search
    real(r8), intent(in) :: tol       ! Error tolerance
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer, parameter :: itmax = 50          ! Maximum number of iterations
    real(r8), parameter :: eps = 1.e-08_r8    ! Relative error tolerance
    integer  :: iter                          ! Iteration loop index
    real(r8) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    real(r8) :: root
    !---------------------------------------------------------------------

    a = xa
    b = xb
    fa = func (ip, iv, il, mlcanopy_inst, a)
    fb = func (ip, iv, il, mlcanopy_inst, b)

    if ((fa > 0._r8 .and. fb > 0._r8) .or. (fa < 0._r8 .and. fb < 0._r8)) then
       write (iulog,*) 'zbrent: Root must be bracketed'
       write (iulog,*) 'called from: ',msg
       write (iulog,*) xa, fa
       write (iulog,*) xb, fb
       call endrun()
    end if
    c = b
    fc = fb
    iter = 0
    do
       if (iter == itmax) exit
       iter = iter + 1
       if ((fb > 0._r8 .and. fc > 0._r8) .or. (fb < 0._r8 .and. fc < 0._r8)) then
          c = a
          fc = fa
          d = b - a
          e = d
       end if
       if (abs(fc) < abs(fb)) then
          a = b
          b = c
          c = a
          fa = fb
          fb = fc
          fc = fa
       end if
       tol1 = 2._r8 * eps * abs(b) + 0.5_r8 * tol
       xm = 0.5_r8 * (c - b)
       if (abs(xm) <= tol1 .or. fb == 0._r8) exit
       if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
          s = fb / fa
          if (a == c) then
             p = 2._r8 * xm * s
             q = 1._r8 - s
          else
             q = fa / fc
             r = fb / fc
             p = s * (2._r8 * xm * q * (q - r) - (b-a) * (r - 1._r8))
             q = (q - 1._r8) * (r - 1._r8) * (s - 1._r8)
          end if
          if (p > 0._r8) q = -q
          p = abs(p)
          if (2._r8*p < min(3._r8*xm*q-abs(tol1*q),abs(e*q))) then
             e = d
             d = p / q
          else
             d = xm
             e = d
          end if
       else
          d = xm
          e = d
       end if
       a = b
       fa = fb
       if (abs(d) > tol1) then
          b = b + d
       else
          b = b + sign(tol1,xm)
       end if
       fb = func (ip, iv, il, mlcanopy_inst, b)
       if (fb == 0._r8) exit
    end do
    root = b

    if (iter == itmax) then
       write (iulog,*) 'zbrent: Maximum number of interations exceeded'
       write (iulog,*) 'called from: ',msg
       call endrun()
    end if

  end function zbrent

  !-----------------------------------------------------------------------
  subroutine quadratic (a, b, c, r1, r2)
    !
    ! !DESCRIPTION:
    ! Solve a quadratic equation for its two roots
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: a,b,c       ! Terms for quadratic equation
    real(r8), intent(out) :: r1,r2       ! Roots of quadratic equation
    !
    ! !LOCAL VARIABLES:
    real(r8) :: q                        ! Temporary term for quadratic solution
    !---------------------------------------------------------------------

    if (a == 0._r8) then
       write (iulog,*) 'Quadratic solution error: a = ',a
       call endrun()
    end if

    if (b >= 0._r8) then
       q = -0.5_r8 * (b + sqrt(b*b - 4._r8*a*c))
    else
       q = -0.5_r8 * (b - sqrt(b*b - 4._r8*a*c))
    end if

    r1 = q / a
    if (q /= 0._r8) then
       r2 = c / q
    else
       r2 = 1.e36_r8
    end if

  end subroutine quadratic

  !-----------------------------------------------------------------------
  subroutine tridiag (a, b, c, r, u, n)
    !
    ! !DESCRIPTION:
    ! Solve a tridiagonal system of equations
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    integer,  intent(in)  :: n        ! Number of soil layers
    real(r8), intent(in)  :: a(n)     ! A vector for tridiagonal solution
    real(r8), intent(in)  :: b(n)     ! B vector for tridiagonal solution
    real(r8), intent(in)  :: c(n)     ! C vector for tridiagonal solution
    real(r8), intent(in)  :: r(n)     ! R vector for tridiagonal solution
    real(r8), intent(out) :: u(n)     ! U vector for tridiagonal solution
    !
    ! !LOCAL VARIABLES:
    real(r8) :: gam(n)                ! Temporary calculation
    real(r8) :: bet                   ! Temporary calculation
    integer  :: j                     ! Soil layer index
    !---------------------------------------------------------------------

    ! Tridiagonal solution:
    !
    ! Solve for U given the set of equations F x U = R, where U is a vector
    ! of length N, R is a vector of length N, and F is an N x N tridiagonal
    ! matrix defined by the vectors A, B, C (each of length N). A(1) and
    ! C(N) are undefined and are not referenced by the subroutine.
    !
    !    | b(1) c(1)   0  ...                      |   | u(1)   |   | r(1)   |
    !    | a(2) b(2) c(2) ...                      |   | u(2)   |   | r(2)   |
    !    |                ...                      | x | ...    | = | ...    |
    !    |                ... a(n-1) b(n-1) c(n-1) |   | u(n-1) |   | r(n-1) |
    !    |                ...   0    a(n)   b(n)   |   | u(n)   |   | r(n)   |
    !

    bet = b(1)
    u(1) = r(1) / bet
    do j = 2, n
       gam(j) = c(j-1) / bet
       bet = b(j) - a(j)*gam(j)
       u(j) = (r(j) - a(j)*u(j-1)) / bet
    end do
    do j = n-1, 1, -1
       u(j) = u(j) - gam(j+1)*u(j+1)
    end do

  end subroutine tridiag

  !-----------------------------------------------------------------------
  function beta_function (p, q) result(beta)
    !
    ! !DESCRIPTION:
    ! Return the value of the beta function evaluated at p and q: B(p,q)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: p      ! Input argument
    real(r8), intent(in) :: q      ! Input argument
    real(r8) :: beta               ! Beta function: B(p,q)

    beta = exp(log_gamma_function(p) + log_gamma_function(q) - log_gamma_function(p+q))
  
  end function beta_function

  !-----------------------------------------------------------------------
  function log_gamma_function (x) result(gammaln)
    !
    ! !DESCRIPTION:
    ! Return the value of the log natural of the gamma function evaluated at x: ln(G(x)) 
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: x     ! Input argument
    real(r8) :: gammaln           ! ln(G(x))
    real(r8) :: y
    real(r8) :: tmp
    real(r8) :: ser
    integer  :: j

    real(r8), parameter :: coef(6) = (/ 76.18009172947146_r8, -86.50532032941677_r8, &
    24.01409824083091_r8, -1.231739572450155_r8, 0.1208650973866179e-02_r8, -0.5395239384953e-05_r8 /)
    real(r8), parameter :: stp = 2.5066282746310005_r8

    y = x
    tmp = x + 5.5_r8
    tmp = (x + 0.5_r8) * log(tmp) - tmp
    ser = 1.000000000190015_r8
    do j = 1, 6
       y = y + 1._r8
       ser = ser + coef(j) / y
    end do
    gammaln = tmp + log(stp*ser/x)

  end function log_gamma_function

end module MathToolsMod
