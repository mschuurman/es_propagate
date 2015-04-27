module math
!
!  Sundry mathematical definitions.
!
  use accuracy
  !$ use OMP_LIB
  implicit none
  private
  public MathFactorial, MathLogFactorial, MathLegendrePn, MathLegendrePnm
  public MathSimpleGamma, Math3J
!
  integer(ik), parameter      :: factorial_slack = 5    ! Extra factorials to produce while filling the cache
  integer(ik), save           :: factorial_max = -1     ! Largest value of the factorials cached in the table
  real(rk), allocatable, save :: factorial_table(:)
  integer(ik), save           :: log_factorial_max = -1 ! Largest value of the factorials cached in the table
  real(rk), allocatable, save :: log_factorial_table(:)
!
  contains
  !
  !  External interfaces
  !
  function MathFactorial(n) result(v)
    integer(ik), intent(in) :: n
    real(rk)                :: v
    !
    if (n<0) stop 'math%MathFactorial - domain error'
    if (n>factorial_max) call fill_factorial_table(n+factorial_slack)
    v = factorial_table(n)
  end function MathFactorial
  !
  function MathLogFactorial(n) result(v)
    integer(ik), intent(in) :: n
    real(rk)                :: v
    !
    if (n<0) stop 'math%MathLogFactorial - domain error'
    if (n>log_factorial_max) call fill_log_factorial_table(n+factorial_slack)
    v = log_factorial_table(n)
  end function MathLogFactorial
  !
  !  If all values of Pn or Pnm up to a given order are required, it is much
  !  more efficient to compute them all at once!
  !
  !  Accuracy of the recursions used to calculate Ledendre functions
  !  deteriorates in the vicinity of the polynomial roots, especially
  !  for very high orders (n,m)
  !
  function MathLegendrePn(n,x) result(pn)
    integer(ik), intent(in) :: n    ! Order of the Legendre polynomial
    real(rk), intent(in)    :: x    ! Coordinate
    real(rk)                :: pn   ! Value of the Legenre polynomial
    !
    real(rk) :: tab(0:n)
    !
    call legendrePn_table(n,x,tab)
    pn = tab(n)
  end function MathLegendrePn
  !
  function MathLegendrePnm(n,m,x) result(pnm)
    integer(ik), intent(in) :: n, m ! Order of the associated Legendre polynomial
    real(rk), intent(in)    :: x    ! Coordinate, abs(x)<=1
    real(rk)                :: pnm  ! Value of the associated Legenre polynomial
    !
    real(rk) :: tab(0:n,0:m)
    !
    call legendrePnm_table(n,m,x,tab)
    pnm = tab(n,m)
  end function MathLegendrePnm
  !
  !  Computes Wigner 3J symbols. The code below is a direct implementation
  !  of L&L 3J formulae. The accuracy of this routine is reduced relative to
  !  that is theoretically possible, due to the use of logarithms. The routine
  !  I had in MNDO99 is more accurate and can handle broader range of J values.
  !
  function Math3J(j1,j2,j3,m1,m2,m3) result(v)
    integer(ik), intent(in) :: j1, j2, j3  ! / J1 J2 J3 \
    integer(ik), intent(in) :: m1, m2, m3  ! \ M1 M2 M3 / 
    real(rk)                :: v
    !
    integer(ik) :: ij0, ij1, ij2, ij3, im1a, im1b, im2a, im2b, im3a, im3b
    integer(ik) :: t1, t2, t3
    integer(ik) :: z, minz, maxz
    real(rk)    :: logscale, logterm
    !
    !  Before we do anything, check whether this 3J symbol satisfies the
    !  vector addition constraints
    !
    ij0  =   j1 + j2 + j3 + 1
    ij1  =   j1 + j2 - j3
    ij2  =   j1 - j2 + j3
    ij3  = - j1 + j2 + j3
    im1a =   j1 - m1 ; im1b = j1 + m1
    im2a =   j2 - m2 ; im2b = j2 + m2
    im3a =   j3 - m3 ; im3b = j3 + m3
    if (ij1<0 .or. ij2<0 .or. ij3<0 .or. im1a<0 .or. im1b<0 .or. im2a<0 .or. im2b<0 .or. im3a<0 .or. im3b<0 .or. m1+m2+m3/=0) then
      v = 0
      return
    end if
    !
    logscale = MathLogFactorial(ij1)  + MathLogFactorial(ij2)  + MathLogFactorial(ij3)  &
             + MathLogFactorial(im1a) + MathLogFactorial(im1b) + MathLogFactorial(im2a) &
             + MathLogFactorial(im2b) + MathLogFactorial(im3a) + MathLogFactorial(im3b) &
             - MathLogFactorial(ij0)
    logscale = 0.5_rk * logscale
    !
    t1   = j2 - j3 - m1
    t2   = j1 + m2 - j3
    t3   = j1 - j2 - m3
    minz = max(0_ik,t1,t2)
    maxz = min(ij1,im1a,im2b)
    v = 0
    sum_terms: do z=minz,maxz,1
      logterm = logscale - MathLogFactorial(z)      - MathLogFactorial(ij1-z)  - MathLogFactorial(im1a-z) &
                         - MathLogFactorial(im2b-z) - MathLogFactorial(z-t1)   - MathLogFactorial(z-t2)
      if (abs(logterm)>=max_exp) then
        write (out,"('Math3J: Intermediate logarithm ',g12.5,' exceeds the real(rk) dynamic range.')") logterm
        write (out,"('Math3J: The 3J arguments were: ',6i10)") j1, j2, j3, m1, m2, m3
        stop 'math%Math3J - exceeded dynamic range'
      end if
      if (mod(z+t3,2)==0) then
        v = v + exp(logterm)
      else
        v = v - exp(logterm)
      end if
    end do sum_terms
    !
  end function Math3J
  !
  !  Auxiliary functions
  !
  subroutine legendrePn_table(nmax,x,pn)
    integer(ik), intent(in) :: nmax  ! Maximum order of the Legendre polynomials desired
    real(rk), intent(in)    :: x     ! Coordinate at which P values are needed
    real(rk), intent(out)   :: pn(:) ! Values of LegendreP from n=0 to n=nmax
    !
    integer(ik) :: n
    real(rk)    :: invn
    !
    if (nmax<0) stop 'math%legendreP_table - negative-order polynomial requested'
    pn(1) = 1._rk
    if (nmax<1) return
    pn(2) = x
    !
    n_recursion: do n=2,nmax
      invn = 1._rk / n
      pn(1+n) = (2._rk-invn)*x*pn(1+(n-1)) - (1._rk-invn)*pn(1+(n-2))
    end do n_recursion
  end subroutine legendrePn_table
  !
  subroutine legendrePnm_table(nmax,mmax,x,pnm)
    integer(ik), intent(in) :: nmax     ! Maximum order of the Legendre polynomials desired
    integer(ik), intent(in) :: mmax     ! Maximum order of the Legendre polynomials desired
    real(rk), intent(in)    :: x        ! Coordinate at which P values are needed, abs(x) must be <=1
    real(rk), intent(out)   :: pnm(:,:) ! Values of LegendreP from n,m=0 to n,m=nmax,mmax
                                        ! n is the first subscript; m is the second subscript
    !
    integer(ik) :: n, m
    real(rk)    :: sqfac ! sqrt(1-x**2)
    !
    sqfac = 1._rk - x**2
    if (sqfac<0) stop 'math%legendrePnm_table - domain error'
    sqfac = sqrt(sqfac)
    !
    call legendrePn_table(nmax,x,pnm(:,1))
    if (mmax<1) return
    !
    !  Special case for m=1: recursion is truncated for n=1, m=1
    !
    pnm(1+0,1+1) = 0._rk
    if (nmax>=1) pnm(1+1,1+1) = -sqfac
    m = 1
    n1_recursion: do n=2,nmax
        pnm(1+n,1+m) = pnm(1+(n-2),1+m) - (2*n-1)*sqfac*pnm(1+(n-1),1+(m-1))
    end do n1_recursion
    !
    m_recursion: do m=2,mmax
      pnm(1+0:1+min(nmax,(m-1)),1+m) = 0._rk
      nm_recursion: do n=m,nmax
        pnm(1+n,1+m) = pnm(1+(n-2),1+m) - (2*n-1)*sqfac*pnm(1+(n-1),1+(m-1))
      end do nm_recursion
    end do m_recursion
  end subroutine legendrePnm_table
  !
  subroutine fill_factorial_table(nmax)
    integer(ik), intent(in) :: nmax
    integer(ik)             :: n, alloc
    real(rk)                :: fac
    !
    !$ if (omp_in_parallel()) then
    !$   stop 'math%fill_factorial_table - unsafe call to MathFactorial'
    !$ end if
    !
    if (factorial_max>=0) then
      deallocate (factorial_table,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i10,' deallocating factorial table')") alloc
        stop 'math%fill_factorial_table - deallocate'
      end if
    end if
    !
    n   = 0
    fac = 1._rk
    !
    allocate (factorial_table(0:nmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i10,' allocating ',i10,'-element factorial table')") & 
             alloc, nmax
      stop 'math%fill_factorial_table - allocate'
    end if
    !
    fill_factorials: do while(n<=nmax-1)
      factorial_table(n) = fac
      n = n + 1
      !
      if (huge(fac)/n<=fac) then
        write (out,"(1x,i10,'! would exceed dynamic range of the chosen real kind')") n
        stop 'math%fill_factorial_table - range exceeded'
      end if
      fac = fac * n
    end do fill_factorials
    factorial_table(n) = fac
    !
    factorial_max = nmax
  end subroutine fill_factorial_table
  !
  subroutine fill_log_factorial_table(nmax)
    integer(ik), intent(in) :: nmax
    integer(ik)             :: n, alloc
    real(rk)                :: fac
    !
    !$ if (omp_in_parallel()) then
    !$   stop 'math%fill_factorial_table - unsafe call to MathLogFactorial'
    !$ end if
    !
    if (log_factorial_max>=0) then
      deallocate (log_factorial_table,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i10,' deallocating log-factorial table')") alloc
        stop 'math%fill_factorial_table - deallocate'
      end if
    end if
    !
    n   = 0
    fac = 0._rk
    !
    allocate (log_factorial_table(0:nmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i10,' allocating ',i10,'-element log-factorial table')") & 
             alloc, nmax
      stop 'math%fill_log_factorial_table - allocate'
    end if
    !
    fill_factorials: do while(n<=nmax-1)
      log_factorial_table(n) = fac
      n = n + 1
      !
      fac = fac + log(real(n,kind=rk))
    end do fill_factorials
    log_factorial_table(n) = fac
    !
    log_factorial_max = nmax
  end subroutine fill_log_factorial_table
  !
  !  Evaluation of the gamma function of complex argument seems to be a bit
  !  complicated. Let's do something very straightforward, and hopefully accurate
  !  enough:
  !    1. Use recurrence relations to _increase_ the real part of the argument
  !       until modulus is large enough to let us use asymprotic expressions
  !    2. Then apply the Stirling formula to get the result.
  !
  !  The resulting routine is not fast, but is reasonably accurate: for
  !  arguments between 0 and 10, the relative error does not exceed 5e-14
  !
  function MathSimpleGamma(z) result(v)
    complex(rk), intent(in) :: z ! Argument of the gamma function
    complex(rk)             :: v ! Gamma function
    !
    !  Stirling formula coefficients
    !
    real(rk), parameter :: c01 =               1._rk/                12._rk
    real(rk), parameter :: c02 =               1._rk/               288._rk
    real(rk), parameter :: c03 =            -139._rk/             51840._rk
    real(rk), parameter :: c04 =            -571._rk/           2488320._rk
    real(rk), parameter :: c05 =          163879._rk/         209018880._rk
    real(rk), parameter :: c06 =         5246819._rk/       75246796800._rk
    real(rk), parameter :: c07 =      -534703531._rk/      902961561600._rk
    real(rk), parameter :: c08 =     -4483131259._rk/    86684309913600._rk
    real(rk), parameter :: c09 = 432261921612371._rk/514904800886784000._rk
    !
    complex(rk) :: zr    ! Reduced argument
    complex(rk) :: logs  ! Scaling coefficient needed to reduce the result of Stirling formula
    complex(rk) :: logv  ! Logarithm of the result
    real(rk)    :: zcut  ! Minimal safe argument for the Stirling formula
    complex(rk) :: vs    ! Series part of the Stirling formula
    !
    !  To get accurate results from Stirling formula, we must make sure that
    !  the modulus of the argument is large enough to make the last contribution
    !  to the expansion small enough.
    !
    zcut = (c09/spacing(1._rk))**(1._rk/9._rk)
    ! write (out,"(' zcut = ',f25.15)") zcut
    logs = 0._rk
    zr   = z
    inflate_z: do while(abs(zr)<zcut) 
      logs = logs + log(zr*(zr+1)*(zr+2))
      zr   = zr + 3
    end do inflate_z
    ! write (out,"(' zr = ',2(1x,g25.15),' logs = ',2(1x,g25.15))") zr, logs
    !
    !  It is safe to use Stirling formula now
    !
    vs = 1._rk + (c01 + (c02 + (c03 + (c04 + (c05 + (c06 + (c07 + (c08 + c09/zr)/zr)/zr)/zr)/zr)/zr)/zr)/zr)/zr
    ! write (out,"(' vs = ',2(1x,g25.15))") vs
    logv = log(sqrt2pi) - zr + (zr-0.5_rk)*log(zr) + log(vs) - logs
    ! write (out,"(' logv = ',2(1x,g25.15))") logv
    !
    v = exp(logv)
  end function MathSimpleGamma
  !
end module math
