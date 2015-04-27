!
!  Calculation of harmonic Franck-Condon factors, using Doctorov's
!  recurrence relations. This program loosely follows the expressions
!  in CPL 138, 310 (1987), with necessary corrections for the (many) 
!  typos. The expressions have been validated against a separate
!  implementation in Mathematica (see Franck-Condon_tunneling_co2.nb),
!  which in turn have been tested against symbolic multidimensional
!  integrals.
!
!  Although integral caching suggested by Gruner and Brumer is implemented,
!  cut-offs about as efficient at terminating recursions for all (small)
!  excitation levels. As a result, caching is disabled by default. The
!  only situatiation where caching actually helps (a little) is in evaluation
!  of non-adiabatic coupling integrals.
!
!  Evaluation of non-adiabatic coupling matrix elements is OpenMP-parallel,
!  with fairly decent speedups. The (much cheaper) evaluation of the FC
!  factors is purely serial - there is no advantage in parallelizing those
!  at the integral level.
!
!  Some useful references:
!
!  1. E.V. Doktorov, I.A. Malkin, and V.I. Man'ko, "Dynamical symmetry of vibronic 
!     transitions in polyatomic molecules and the Franck-Condon principle", 
!     J. Mol. Spectroscopy 64, 302-326 (1977)
!  2. D. Gruner and P. Brumer, "Efficient evaluation of harmonic polyatomic 
!     Franck-Condon factors", Chem. Phys. Lett. 138, 310-314 (1987)
!  3. J. Dods, D. Gruner, and P. Brumer, "A genetic algorithm approach to fitting 
!     polyatomic spectra via geometry shifts", Chem. Phys. Lett. 261, 612-619 (1996).
!  
!  Version history:
!
!   December 10, 2008 - Initial version, Serguei.Patchkovskii@nrc.ca
!   December 19, 2008 - Added evaluation of non-BO transition probabilities
!   December 24, 2008 - Added somewhat sensible handling of disjoint vibrational
!                       spaces and OpenMP parallization for non-BO integrals
!   January  07, 2009 - Implemented binary-tree integral caching scheme
!
  module franck_condon
    use accuracy
    use elements
    use fc_tools
    use cache_fc
    use lapack
    use math
    use timer
    !$ use omp_lib
    implicit none
    private
    public run_fc
    !
    integer(ik)                  :: unit_scr         = 40  ! Any unused unit
    integer(ik)                  :: max_state_string = 256 ! Max length of the string describing a state
    !
    !  User-defined parameters, controlled through FC_PAR namelist
    !
    integer(ik)        :: verbose       = 0        ! Verbosity level of the output
    logical            :: use_cache     = .true.   ! Use integral cache in evaluation
    real(rk)           :: max_cache_size= 2e7_rk   ! Max. allowed cache size - about 500Mbytes
    real(rk)           :: wgt_cutoff    = 0._rk    ! Cut-off for integral summation
    real(rk)           :: cache_cutoff  = 1e-1_rk  ! Relative significance cut-off for cache lookup
    character(len=80)  :: s1_xyz_name   = ' '      ! File containing state 1 geometry
    character(len=80)  :: s2_xyz_name   = ' '      ! File containing state 2 geometry
    character(len=80)  :: s1_hess_name  = ' '      ! File containing state 1 Hessian
    character(len=80)  :: s2_hess_name  = ' '      ! File containing state 2 Hessian
    character(len=80)  :: coupling_name = ' '      ! File containing derivative coupling vector
                                                   ! for electronic wavefunctions:
                                                   !   <psi2|\Nabla|psi1>
                                                   ! Empty input means do not calculate derivative
                                                   ! coupling amplitudes
    real(rk)           :: freq_cut      = 50._rk   ! Normal modes with lower frequencies (cm^-1) are 
                                                   ! considered to be translations and rotations
    real(rk)           :: freq_dummy    = 10._rk   ! Replace very soft modes with this freqeuncy
    integer(sik)       :: nvib_max      = 4_sik    ! Highest vibrational level to consider
    real(rk)           :: evib_max      = 5000._rk ! Max. vibrational energy to consider (cm^-1)
    character(len=20)  :: task          = 'all'    ! Computational tasks. Possible choices are:
                                                   ! 'both'  - evaluate all FC factors for all modes
                                                   !           satisfying nvib_max/evib_max constraints
                                                   ! 'one'   - evaluate specific FC intergrals, as 
                                                   !           requested in the input file
                                                   ! 'left'  - evaluate all FC integrals for a specific
                                                   !           right-hand side
                                                   ! 'right' - evaluate all FC integrals for a specific
                                                   !           left-hand side
                                                   ! 'test'  - special case of "one", with the reference
                                                   !           value of the integral appearing in the
                                                   !           first column of the input stream
    real(rk)           :: report_threshold         ! Smaller integrals will not be reported
    !
    !  Everything below is computed internally
    !
    type(surfaceData)     :: s1, s2                ! Two sufraces is all we have!
    real(rk), allocatable :: coupling(:)           ! Derivative coupling vector (normal coordinates of s1 - RHS)
    real(rk)              :: zzIntegral            ! Primitive overlap integral
    real(rk), allocatable :: smat(:,:)             ! Dushinski rotation matrix
    real(rk), allocatable :: dvec(:)               ! Structural displacement vector
    real(rk), allocatable :: jmat(:,:)             
    real(rk), allocatable :: deltavec(:)           
    real(rk), allocatable :: qmat(:,:)             
    real(rk), allocatable :: pmat(:,:)             
    real(rk), allocatable :: rmat(:,:)             
    real(rk), allocatable :: rdelta(:)             
    real(rk), allocatable :: ImPdelta(:)           
    !
    real(ark)                 :: integral_count  = 0     ! Total number of integrals evaluated
    real(ark)                 :: integral_cached = 0     ! Total number of integrals looked up
    real(ark)                 :: max_cache_used  = 0._rk ! Max. size of the cache
    real(ark)                 :: cache_flushes   = 0._rk ! Number of times cache was flushed
    logical                   :: cacheSet = .false.
    integer(sik), allocatable :: n2set(:)                ! Max cache indices
    integer(sik), allocatable :: n1set(:)        
    type(fcCacheWrapper), allocatable :: fcCaches(:)     ! Per-thread integral caches
    !
    namelist /FC_PAR/                                        &
       verbose,                                              &
       task,                                                 &
       use_cache, max_cache_size,                            &
       wgt_cutoff, cache_cutoff,                             &
       freq_cut, freq_dummy,                                 &
       nvib_max, evib_max, report_threshold,                 &
       s1_xyz_name, s2_xyz_name, s1_hess_name, s2_hess_name, &
       coupling_name
    !
    contains
    !
    subroutine read_coupling_vector
      integer(ik) :: ios, id
      real(rk)    :: grad_cart(s1%nvars)
      real(rk)    :: grad_drop(s1%nvars), ovd
      !
      call TimerStart('Coupling vector')
      open (unit_scr,file=trim(coupling_name),form='formatted',action='read',status='old',iostat=ios)
      if (ios/=0) then
        write (out,"(' Error ',i6,' opening file ',a)") ios, trim(coupling_name)
        stop 'fc_tools%read_coupling_vector - open'
      end if
      read (unit_scr,*,iostat=ios) grad_cart
      if (ios/=0) then
        write (out,"(' Error ',i6,' reading file ',a)") ios, trim(coupling_name)
        stop 'fc_tools%read_coupling_vector - read'
      end if
      close (unit_scr,iostat=ios)
      if (ios/=0) then
        write (out,"(' Error ',i6,' closing file ',a)") ios, trim(coupling_name)
        stop 'fc_tools%read_coupling_vector - open'
      end if
      if (verbose>=2) then
        call print_vector(grad_cart,'Cartesian derivative coupling vector')
      end if
      !
      !  Transform Cartesian coupling vector to normal coordinates
      !
      allocate (coupling(s1%nvars))
      grad_cart = grad_cart / sqrt(s1%mass)
      coupling  = -matmul(grad_cart,s1%modes(:,1:s1%nvars))
      if (verbose>=2) then
        call print_vector(coupling,'Normal derivative coupling vector')
      end if
      !
      !  Report parts of the the derivative coupling vector which are projected
      !  onto the dummy modes.
      !
      if (verbose>=1) then
        grad_drop = 0
        scan_dropped: do id=1,s1%nvarsDummy
          ovd = -dot_product(grad_cart,s1%modes(:,id))
          write (out,"('Non-adiabatic coupling along dummy S1 mode ',i5,' (',f12.4,'/',f12.4,' cm^-1) is ',g12.6)") &
                 id, h2cm*s1%freq(id), h2cm*s1%freqTrue(id), ovd
          grad_drop = grad_drop - ovd*s1%modes(:,id)
        end do scan_dropped
        grad_drop = grad_drop * sqrt(s1%mass)
        call print_vector(grad_drop,'Neglected part of Cartesian non-adiabatic coupling vector')
      end if
      !
      call TimerStop('Coupling vector')
    end subroutine read_coupling_vector
    !
    !
    subroutine prepare_for_integral_recursion
      integer(ik) :: nv, nvt
      integer(ik) :: i, j
      real(rk)    :: sdet
      !
      call TimerStart('Prepare for recursion')
      !
      !  We'll be using these two dimensions a lot, so ...
      !
      nv  = s1%nvars
      nvt = s1%nvars
      !
      allocate (smat(nvt,nvt),dvec(nvt),jmat(nvt,nvt),deltavec(nvt), &
                qmat(nvt,nvt),pmat(nvt,nvt),rmat(nvt,nvt),rdelta(nvt), &
                ImPdelta(nvt))
      !
      smat = matmul(transpose(s2%modes(:,1:nvt)),s1%modes(:,1:nvt))
      if (verbose>=2) call print_matrix(smat,'Dushinski rotation matrix')
      !
      !  Final sanity check - Dushinski matrix must have unit determinant,
      !  or the recursions are invalid. Ideally, the "real" vibrational
      !  subspace must also have unit determinant, so that these vibrations
      !  do not project into the "dummy" space too much.
      !
      sdet = linpack_determinant(smat)
      if ( abs(abs(sdet)-1.0_rk)>spacing(100._rk) ) then
        write (out,"(/'Determinant of Dushinski rotation matrix is not unity: ',g25.15)") sdet
        write (out,"(/'Normal modes of the S1 and S2 surfaces do not span the same space,')")
        write (out,"( 'and the Franck-Condon factors are not well-defined.')")
        stop 'franck_condon%prepare_for_integral_recursion'
      end if
      if (s1%nvarsDummy==s2%nvarsDummy) then
        sdet = linpack_determinant(smat(s2%nvarsDummy+1:,s1%nvarsDummy+1:))
        if ( abs(abs(sdet)-1.0_rk)>spacing(100._rk) ) then
          write (out,"(/'Determinant of the ""true"" vibration subspace is not unity: ',g25.15)") sdet
          write (out,"( 'Absolute values of the Franck-Condon factors are not very meaningful.'/)")
        end if
      end if
      ! 
      dvec = matmul(transpose(s2%modes(:,1:nvt)),sqrt(s1%mass)*(s2%coord-s1%coord))
      if (verbose>=2) call print_vector(dvec,'Structural displacement vector')
      !
      build_jmat: do j=1,nvt
        do i=1,nvt
          jmat(i,j) = sqrt(s2%freq(i)/s1%freq(j))*smat(i,j)
        end do
      end do build_jmat
      if (verbose>=2) call print_matrix(jmat,'Doktorov J matrix')
      !
      deltavec = sqrt(s2%freq(1:nvt)) * dvec
      if (verbose>=2) call print_vector(deltavec,'Doktorov delta dector')
      !
      qmat = matmul(transpose(jmat),jmat)
      qmat_diagonal: do i=1,nvt
        qmat(i,i) = qmat(i,i) + 1._rk
      end do qmat_diagonal
      call lapack_ginverse(qmat)
      if (verbose>=2) call print_matrix(qmat,'Doktorov Q matrix')
      !
      pmat = matmul(jmat,matmul(qmat,transpose(jmat)))
      if (verbose>=2) call print_matrix(pmat,'Doktorov P matrix')
      !
      rmat = matmul(qmat,transpose(jmat))
      if (verbose>=2) call print_matrix(rmat,'Doktorov R matrix')
      !
      rdelta = matmul(rmat,deltavec)
      if (verbose>=2) call print_vector(rdelta,'R . delta')
      !
      ImPdelta = deltavec - matmul(pmat,deltavec)
      if (verbose>=2) call print_vector(ImPdelta,'(1 - P) . delta')
      !
      call TimerStop('Prepare for recursion')
    end subroutine prepare_for_integral_recursion
    !
    subroutine evaluate_primitive_integral
      real(rk)    :: exparg, det, freq
      !
      call TimerStart('Primitive integral')
      !
      freq   = product(s2%freq(1:s2%nvars)/s1%freq(1:s1%nvars))**0.25_rk
      det    = sqrt(linpack_determinant(qmat))
      exparg = sum(deltavec**2) - dot_product(deltavec,matmul(pmat,deltavec))
      zzIntegral = (sqrt(2._rk)**s1%nvars) * freq * det * exp(-0.5_rk * exparg)
      !
      if (verbose>=0) write (out,"('Primitive integral is ',g25.15)") zzIntegral
      if (verbose>=1) then
        write (out,"(t5,'=  ',g25.15,' [2**(nvars/2)]')") sqrt(2._rk)**s1%nvars
        write (out,"(t5,' x ',g25.15,' [(omega2/omega1)**0.25]')") freq
        write (out,"(t5,' x ',g25.15,' [sqrt(Det(qmat))]')") det
        write (out,"(t5,' x ',g25.15,' [exp(-0.5*...)]')") exp(-0.5_rk * exparg)
      end if
      !
      call TimerStop('Primitive integral')
    end subroutine evaluate_primitive_integral
    !
    recursive function fc_recursion_right(n2,n1,wgt,i) result(fc)
      integer(sik), intent(inout) :: n1(:), n2(:) ! Vibrational state index vectors
      real(rk), intent(in)        :: wgt          ! Weight of the requested integral
      integer(ik), intent(in)     :: i            ! Position of the largest RHS element
      real(rk)                    :: fc           ! Franck-Condon factor <n2|n1>
      !
      integer(ik) :: n1act, j, k
      real(rk)    :: term, coeff, wgtp
      !
      !$omp atomic
      integral_count = integral_count + 1
      !
      !  Recursion on the right - n1(n1max) is always decremented
      !
      n1act = n1(i)
      if (n1act==0) stop 'franck_condon%fc_recursion_right - recursion along zero index'
      n1(i) = n1(i) - 1
      fc    = 0
      !
      !  Unfortunately (?), for sensibly small numbver of quanta integral evaluation
      !  is extremely cheap. It therefore makes no sense to execute these routines
      !  in parallel internally. 
      !
      right_n2terms: do k=1,s2%nvars
        if (n2(k)<=0) cycle right_n2terms
        coeff = 2*rmat(i,k)*sqrt(real(n2(k),kind=rk)/real(n1act,kind=rk))
        wgtp  = abs(wgt*coeff)
        if (wgtp<=wgt_cutoff) cycle right_n2terms ! This integral does not contibute
        !
        n2(k) = n2(k) - 1
        term  = evaluate_fc_integral(n2,n1,wgtp)
        n2(k) = n2(k) + 1
        fc    = fc + coeff*term
      end do right_n2terms
      !
      right_n1terms: do j=1,s1%nvars
        if (n1(j)<=0) cycle right_n1terms
        coeff = 2*qmat(i,j)
        if (i==j) coeff = coeff - 1._rk
        coeff = coeff * sqrt(real(n1(j),kind=rk)/real(n1act,kind=rk))
        wgtp  = abs(wgt*coeff)
        if (wgtp<=wgt_cutoff) cycle right_n1terms ! This integral does not contibute
        n1(j) = n1(j) - 1
        term  = evaluate_fc_integral(n2,n1,wgtp)
        n1(j) = n1(j) + 1
        fc    = fc + coeff*term
      end do right_n1terms
      !
      coeff = -rdelta(i)*sqrt(2._rk/real(n1act,kind=rk))
      wgtp  = abs(wgt*coeff)
      if (wgtp>wgt_cutoff) then
        term = evaluate_fc_integral(n2,n1,wgtp)
        fc   = fc + coeff*term
      end if
      !
      !  Restore state vector
      !
      n1(i) = n1(i) + 1
    end function fc_recursion_right
    !
    recursive function fc_recursion_left(n2,n1,wgt,k) result(fc)
      integer(sik), intent(inout) :: n2(s2%nvars) ! Vibrational state index vectors
      integer(sik), intent(inout) :: n1(s1%nvars) ! Vibrational state index vectors
      real(rk), intent(in)        :: wgt          ! Weight of the requested integral in the
      integer(ik), intent(in)     :: k            ! Position of the largest LHS element
      real(rk)                    :: fc           ! Franck-Condon factor <n2|n1>
      !
      integer(ik) :: n2act, i, l
      real(rk)    :: term, coeff, wgtp
      !
      !$omp atomic
      integral_count = integral_count + 1
      !
      !  Recursion on the left - n2(k) is always decremented
      !
      n2act = n2(k)
      if (n2act==0) stop 'franck_condon%fc_recursion_left - recursion along zero index'
      n2(k) = n2(k) - 1
      fc    = 0
      !
      left_n1terms: do i=1,s1%nvars
        if (n1(i)<=0) cycle left_n1terms
        coeff = 2*rmat(i,k)*sqrt(real(n1(i),kind=rk)/real(n2act,kind=rk))
        wgtp  = abs(wgt*coeff)
        if (wgtp<=wgt_cutoff) cycle left_n1terms ! This integral does not contibute
        !
        n1(i) = n1(i) - 1
        term  = evaluate_fc_integral(n2,n1,wgtp)
        n1(i) = n1(i) + 1
        fc    = fc + coeff*term
      end do left_n1terms
      !
      left_n2terms: do l=1,s2%nvars
        if (n2(l)<=0) cycle left_n2terms
        coeff = 2*pmat(k,l)
        if (k==l) coeff = coeff - 1._rk
        coeff = coeff * sqrt(real(n2(l),kind=rk)/real(n2act,kind=rk))
        wgtp  = abs(wgt*coeff)
        if (wgtp<=wgt_cutoff) cycle left_n2terms ! This integral does not contibute
        n2(l) = n2(l) - 1
        term  = evaluate_fc_integral(n2,n1,wgtp)
        n2(l) = n2(l) + 1
        fc    = fc + coeff*term
      end do left_n2terms
      !
      coeff = ImPdelta(k)*sqrt(2._rk/real(n2act,kind=rk))
      wgtp  = abs(wgt*coeff)
      if (wgtp>wgt_cutoff) then
        term = evaluate_fc_integral(n2,n1,wgtp)
        fc   = fc + coeff*term
      end if
      !
      !  Restore state vector
      !
      n2(k) = n2(k) + 1
    end function fc_recursion_left
    !
    !  Find the location of the smallest non-zero element of an array
    !
    function lowest_excitation(n) result(imin)
      integer(sik), intent(in) :: n(:)  ! Array of integers
      integer(ik)              :: imin  ! Location of the smallest non-zero integer
                                        ! if one is present. 1 otherwise.
      !
      integer(ik)              :: i
      !
      imin = 1
      scan_array: do i=2,size(n)
        if (n(i)>0 .and. (n(i)<n(imin) .or. n(imin)==0)) imin = i
      end do scan_array
    end function lowest_excitation
    !
    !  Integral is not in the cache - evaluate it the hard way.
    !
    recursive function do_evaluate_fc_integral(n2,n1,wgt) result(fc)
      integer(sik), intent(inout) :: n1(:), n2(:) ! Vibrational state index vectors
      real(rk), intent(in)        :: wgt          ! Weight of the requested integral
      real(rk)                    :: fc           ! Franck-Condon factor <n2|n1>
      !
      integer(ik)            :: n1max, n2max
      !
      !
      !  In principle, we can perform recursion for any degree of freedom
      !  with non-zero excitation level. We will try to make the recurrences
      !  as narrow as possible, by eliminating degrees of freedom with the
      !  smallest non-zero number of quanta first. Originally, we were trying
      !  to reduce the highest excitation levels first, which leads to very
      !  broad recursion trees at high excitation levels.
      !
      n1max = lowest_excitation(n1)
      n2max = lowest_excitation(n2)
      !
      !  Terminate recursion.
      !
      if ( (n1(n1max)==0) .and. (n2(n2max)==0) ) then
        fc = zzIntegral
        return
      end if
      !
      if ( n1(n1max)==0 .or. (n1(n1max)>n2(n2max) .and. n2(n2max)/=0) ) then
        fc = fc_recursion_left(n2,n1,wgt,n2max)
      else
        fc = fc_recursion_right(n2,n1,wgt,n1max)
      end if
      !
    end function do_evaluate_fc_integral
    !
    !  Evaluate an FC integral, possibly using cache lookup
    !
    recursive function evaluate_fc_integral(n2,n1,wgt) result(fc)
      integer(sik), intent(inout) :: n1(:), n2(:) ! Vibrational state index vectors
                                                  ! Although the vectors are marked "inout",
                                                  ! the contents will remain unchanged upon
                                                  ! return.
                                                  ! Note that s2 is on the left, s1 is on the right
      real(rk), intent(in)        :: wgt          ! Weight of the requested integral in the
                                                  ! final result - useful for cut-offs!
      real(rk)                    :: fc           ! Franck-Condon factor <n2|n1>
      !
      integer(ik)                 :: ic           ! Cache index
      type(fcCache), pointer      :: cache        ! Cache location for the current integral
      !
      !  Check out integral cache first
      !
      cache => NULL()
      if (use_cache) then
        ic = 1
        !$ ic = omp_get_thread_num() + 1
        cache => cache_fc_locate(fcCaches(ic),n2,n1)
        if (associated(cache)) then
          if (cache%weight>=cache_cutoff*wgt) then
            !$omp atomic
            integral_cached = integral_cached + 1
            fc = cache%value
            !???
            !write (out,"('fetch: ',6(1x,i1),'/',6(1x,i1),' = ',2f20.12)") n2, n1, cache%value, cache%weight
            return
          end if
        end if
      end if
      !
      fc = do_evaluate_fc_integral(n2,n1,wgt)
      !???
      !write (out,"('eval: ',6(1x,i1),'/',6(1x,i1),' = ',2f20.12)") n2, n1, fc, wgt
      !
      if (associated(cache)) then
        if (wgt>cache%weight) then
          cache%value  = fc
          cache%weight = wgt
          !???
          !write (out,"('store: ',6(1x,i1),'/',6(1x,i1),' = ',2f20.12)") n2, n1, cache%value, cache%weight
        end if
      end if
    end function evaluate_fc_integral
    !
    function evaluate_derivative_coupling(n2,n1) result(amp)
      integer(sik), intent(inout) :: n2(s2%nvars) ! Vibrational state index vectors
      integer(sik), intent(inout) :: n1(s1%nvars) ! Vibrational state index vectors
                                                  ! Although the vectors are marked "inout",
                                                  ! the contents will remain unchanged upon
                                                  ! return.
                                                  ! Note that s2 is on the left, s1 is on the right
      real(rk)                    :: amp          ! Derivative coupling amplitude
      !
      integer(ik) :: s1mode
      real(rk)    :: fc_down, fc_up
      real(rk)    :: wgt_down, wgt_up
      !
      amp = 0
      !$omp parallel do default(none) schedule(guided) &
      !$omp& firstprivate(n1,n2) reduction(+:amp) &
      !$omp& private(s1mode,wgt_down,wgt_up,fc_down,fc_up) &
      !$omp& shared(s1,coupling)
      scan_s1: do s1mode=s1%nvarsDummy+1,s1%nvars
        !
        !  Use recurrence relations to calculate gradient of the RHS with respect to this mode
        !
        wgt_down =  coupling(s1mode)*sqrt(0.5_rk*s1%freq(s1mode)*n1(s1mode))
        wgt_up   = -coupling(s1mode)*sqrt(0.5_rk*s1%freq(s1mode)*(n1(s1mode)+1))
        !
        if (n1(s1mode)>=1) then
          n1(s1mode) = n1(s1mode) - 1
          fc_down = evaluate_fc_integral(n2,n1,abs(wgt_down))
          n1(s1mode) = n1(s1mode) + 1
        else
          fc_down = 0
        end if
        n1(s1mode) = n1(s1mode) + 1
        !
        !  By forcing recursion along s1mode, we guarantee that this integral
        !  will always hit a "hot" integral cache.
        !
        ! fc_up = evaluate_fc_integral(n2,n1,abs(wgt_up))
        fc_up = fc_recursion_right(n2,n1,abs(wgt_up),s1mode)
        n1(s1mode) = n1(s1mode) - 1
        amp = amp + wgt_down * fc_down + wgt_up * fc_up
      end do scan_s1
      !$omp end parallel do
    end function evaluate_derivative_coupling
    !
    function cache_size_single(n,order) result(sz)
      integer(sik), intent(in) :: n(:)     ! Excitation levels
      integer(sik), intent(in) :: order    ! Excitation order
      real(rk)                 :: sz       ! Estimated cache size
      !
      integer(ik) :: slots
      !
      if (order<=0) then
        sz = product(n+1._rk)
      else
        slots = count(n>0)
        if (order>=slots) then
          sz = (1+real(order,kind=rk)/slots)**slots
        else
          sz = 2._rk**order * MathFactorial(slots)
        end if
      end if
    end function cache_size_single
    !
    real(rk) function cache_size(n2,n1,n2order,n1order)
      integer(sik), intent(in) :: n2(:), n1(:)     ! Excitation levels
      integer(sik), intent(in) :: n2order, n1order ! Automatically generated excitation levels
      !
      cache_size = cache_size_single(n2,n2order) * cache_size_single(n1,n1order)
    end function cache_size 
    !
    subroutine do_initialize_cache(n2,n1)
      integer(sik), intent(in) :: n2(:), n1(:) ! Cached excitation levels
      !
      integer(ik)  :: max_threads, it
      !
      max_threads = 1
      !$ max_threads = omp_get_max_threads()
      allocate (fcCaches(max_threads))
      !$omp parallel do
      init_thread_cache: do it=1,max_threads
        call cache_fc_initialize(fcCaches(it),n2,n1)
      end do init_thread_cache
      !$omp end parallel do
    end subroutine do_initialize_cache
    !
    subroutine initialize_cache(n2,n1,n2order,n1order)
      integer(sik), intent(in) :: n2(:), n1(:) ! Minimal excitation levels
      integer(sik)             :: n2order      ! If n2order > 0, we will be iterating over lower
                                               ! decades of n2 - pad these positions as much as
                                               ! possible
      integer(sik)             :: n1order      ! Ditto for n1. Only one of n2order and n1order 
                                               ! can be non-zero
      !
      integer(sik) :: n2max(s2%nvars)          ! Chosen max excitation levels
      integer(sik) :: n1max(s1%nvars)
      real(rk)     :: total_size
      integer(ik)  :: i2next, i1next           ! Next variable to increment
      !
      if (.not.use_cache) return
      !
      !  Baseline indices
      !
      n2max = n2 ; n2max(:s2%nvarsDummy) = 0
      n1max = n1 ; n1max(:s1%nvarsDummy) = 0
      !
      !  If the baseline indices are compatible with the cache already in place,
      !  there is no need to do anything.
      !
      if (cacheSet) then
        if (all(n2max<=n2set) .and. all(n1max<=n1set)) then
          if (verbose>2) then
            write (out,"('Currently active cache can accommodate the new integral:')")
            write (out,"('n2: ',(t5,25(1x,i2)/))") n2max
            write (out,"('n1: ',(t5,25(1x,i2)/))") n1max
          end if
          return
        end if
        call release_cache
      end if
      !
      total_size = cache_size(n2max,n1max,n2order,n1order)
      if (total_size>max_cache_size) then
        if (verbose>=-1) then
          write (out,"(/'WARNING: Integral cache size of ',f12.0,' elements is not sufficient'" // &
                      "/'to reach the desired excitation level. Caching of FC recursions will be'" // &
                      "/'disabled from this point on.')") max_cache_size
          write (out,"(  'WARNING: Increase cache size to at least ',f12.0,' to activate caching.'/)") total_size
          use_cache = .false.
          return
        end if
      end if
      !
      if (.not.allocated(n2set)) allocate (n2set(s2%nvars))
      if (.not.allocated(n1set)) allocate (n1set(s1%nvars))
      !
      i2next = s2%nvars
      i1next = s1%nvars
      grow_cache: do while(total_size<=max_cache_size .and. i1next>s1%nvarsDummy .and. i2next>s2%nvarsDummy )
        n2set = n2max ; n1set = n1max
        if      (n2order>0) then
          n2max(i2next) = max(n2max(i2next),min(n2order,s2%maxLevel(i2next)))
          i2next        = i2next - 1
        else if (n1order>0) then
          n1max(i1next) = max(n1max(i1next),min(n1order,s1%maxLevel(i1next)))
          i1next        = i1next - 1
        else
          exit grow_cache
        end if
        total_size = cache_size(n2max,n1max,n2order,n1order)
      end do grow_cache
      !
      if (verbose>2) then
        write (out,"('Expected cache size = ',f12.0)") total_size
        write (out,"('n2: ',(t5,25(1x,i2)))") n2set
        write (out,"('n1: ',(t5,25(1x,i2)))") n1set
      end if
      !
      call do_initialize_cache(n2set,n1set)
      cacheSet = .true.
    end subroutine initialize_cache
    !
    subroutine release_cache
      integer(ik) :: it
      real(rk)    :: total_size, sz
      !
      if (.not.use_cache) return
      if (.not.cacheSet) return
      !
      total_size = 0
      !$omp parallel do reduction(+:total_size) private(sz)
      release_thread_cache: do it=1,size(fcCaches)
        call cache_fc_destroy(fcCaches(it),sz)
        total_size = total_size + sz
      end do release_thread_cache
      !$omp end parallel do
      deallocate (fcCaches)
      !
      !$omp atomic
      max_cache_used = max(total_size,max_cache_used)
      !$omp atomic
      cache_flushes  = cache_flushes + 1
      !
      if (verbose>2) then
        write (out,"('Actual size of the FC integral cache was ',f12.0,' elements')") total_size
      end if
      cacheSet = .false.
    end subroutine release_cache
    !
    recursive function increment_state(ntot,n) result(done)
      integer(sik), intent(in)    :: ntot ! Constraint on the total number of active modes
      integer(sik), intent(inout) :: n(:) ! Vibrational quanta
      logical                     :: done ! No more states satisfying the constraint
      !
      integer(ik) :: nstate
      !
      nstate = size(n)
      !
      !  Terminate recursion
      !
      if (nstate<=1) then
        done = .true.
        return
      end if
      !
      if (increment_state(ntot-n(1),n(2:))) then
        !
        !  All states to the right of the leftmost decade are explored,
        !  try to increment the left-most decade.
        !
        if (n(1)>=ntot) then
          !
          !  All combinations exhausted.
          !
          done = .true.
          return
        end if
        !
        !  Increment the leftmost decade, and reset the remaining state
        !
        n(1)      = n(1) + 1
        n(2:)     = 0
        n(nstate) = ntot - n(1)
      end if
      done = .false.
    end function increment_state
    !
    subroutine evaluate_all_fc_recursions
      integer(sik) :: n1(s1%nvars)  ! State vector - RHS
      integer(sik) :: n2(s2%nvars)  ! State vector - LHS
      !
      integer(sik) :: nvib_tot      ! Total number of active modes
      integer(sik) :: nvib1, nvib2  ! ... and for s1/s2, respectively
      real(rk)     :: evib1, evib2  ! Harmonic vibrational energies
      real(rk)     :: fc            ! Integral
      logical      :: do_dc         ! True if derivative coupling is needed
      !
      call TimerStart('Evaluate all integrals')
      !
      do_dc = (coupling_name/=' ')
      !
      if (do_dc) then
        write (out,"(/'WARNING: Evaluating derivative coupling integrals, NOT the Franck-Condon integrals')")
      end if
      write (out,"()")
      write (out,"(a18,1x,a7,1x,a7,t36,a8,t68,a8)") &
                 '# <s2|s1>', 'e2,cm-1', 'e1,cm-1', 'state2', 'state1', &
                 '# -------', '-------', '-------', '------', '------'
      excitation_level: do nvib_tot=0,nvib_max
        !
        !  Calculate all integrals at a given excitation level,
        !  for either ground or excited state
        !
        excitation_s1: do nvib1=0,nvib_tot
          nvib2=nvib_tot-nvib1
          if(verbose>=2) then
            write (out,"('# LHS quanta = ',i3,' RHS quanta = ',i3)") nvib2, nvib1
          end if
          !
          !  State enumeration starts with the maximum number of quanta in the
          !  last vibrational mode - RHS
          !
          n1 = 0 ; n1(s1%nvars) = nvib1
          states_n1: do
            evib1 = h2cm*state_energy(s1,n1,zpe=.false.)
            if (evib1<=evib_max) then
              !
              !  The left-hand-side state
              !
              n2 = 0 ; n2(s2%nvars) = nvib2
              states_n2: do
                evib2 = h2cm*state_energy(s2,n2,zpe=.false.)
                if (evib2<=evib_max) then
                  call initialize_cache(n2,n1,nvib2,0_sik)
                  if (do_dc) then
                    fc = evaluate_derivative_coupling(n2,n1)
                  else
                    fc = evaluate_fc_integral(n2,n1,1._rk)
                  end if
                  if ( (abs(fc)>report_threshold) .or. (verbose>=3) ) then
                    write (out,"(f18.13,1x,f7.1,1x,f7.1,t36,a30,t68,a30)") &
                           fc, evib2, evib1, state2ascii(n2), state2ascii(n1)
                  end if
                end if
                if (increment_state(nvib2,n2(s2%nvarsDummy+1:))) exit states_n2
              end do states_n2
            end if ! evib1<=evib_max)
            if (increment_state(nvib1,n1(s1%nvarsDummy+1:))) exit states_n1
          end do states_n1
        end do excitation_s1
      end do excitation_level
      !
      !  It is not necessary to relase caches inside the loop - initialize_cache()
      !  is smart enough to handle reallocation when needed
      !
      call release_cache
      call TimerStop('Evaluate all integrals')
    end subroutine evaluate_all_fc_recursions
    !
    !  Evaluate all FC integrals for the right-hand side, fixing the
    !  left-hand side.
    !
    subroutine evaluate_right_fc_recursions
      integer(sik)                    :: n1(s1%nvars) ! State vector - RHS
      integer(sik)                    :: n2(s2%nvars) ! State vector - LHS
      character(len=max_state_string) :: n2str        ! State names (LHS)
      !
      integer(sik)                    :: nvib1        ! Number of active modes - RHS
      real(rk)                        :: evib1, evib2 ! Harmonic vibrational energies
      real(rk)                        :: fc           ! Integral
      integer(ik)                     :: info
      logical                         :: do_dc        ! True if derivative coupling is needed
      !
      call TimerStart('Evaluate right integrals')
      !
      do_dc = (coupling_name/=' ')
      !
      if (do_dc) then
        write (out,"(/'WARNING: Evaluating derivative coupling integrals, NOT the Franck-Condon integrals')")
      end if
      write (out,"()")
      write (out,"(a18,1x,a7,1x,a7,t36,a8,t68,a8)") &
                 '# <s2|s1>', 'e2,cm-1', 'e1,cm-1', 'state2', 'state1', &
                 '# -------', '-------', '-------', '------', '------'
      !
      read_loop: do
        read(input,"(a)",iostat=info) n2str
        if (info/=0) exit read_loop
        n2    = ascii2state(trim(n2str),s2%nvars)
        evib2 = h2cm*state_energy(s2,n2,zpe=.false.)
        !
        excitation_s1: do nvib1=0,nvib_max
          if(verbose>=2) then
            write (out,"('# RHS quanta = ',i3)") nvib1
          end if
          !
          !  State enumeration starts with the maximum number of quanta in the
          !  last vibrational mode
          !
          n1 = 0 ; n1(s1%nvars) = nvib1
          states_n1: do
            evib1 = h2cm*state_energy(s1,n1,zpe=.false.)
            if (evib1<=evib_max) then
              call initialize_cache(n2,n1,0_sik,nvib1)
              if (do_dc) then
                fc = evaluate_derivative_coupling(n2,n1)
              else
                fc = evaluate_fc_integral(n2,n1,1._rk)
               end if
              if ( (abs(fc)>report_threshold) .or. (verbose>=3) ) then
                write (out,"(f18.13,1x,f7.1,1x,f7.1,t36,a30,t68,a30)") &
                       fc, evib2, evib1, state2ascii(n2), state2ascii(n1)
              end if
            end if ! evib1<=evib_max)
            if (increment_state(nvib1,n1(s1%nvarsDummy+1:))) exit states_n1
          end do states_n1
        end do excitation_s1
      end do read_loop
      !
      !  It is not necessary to relase caches inside the loop - initialize_cache()
      !  is smart enough to handle reallocation when needed
      !
      call release_cache
      call TimerStop('Evaluate right integrals')
    end subroutine evaluate_right_fc_recursions
    !
    !  Evaluate all FC integrals for the left-hand side, fixing the
    !  right-hand side.
    !
    subroutine evaluate_left_fc_recursions
      integer(sik)                    :: n1(s1%nvars) ! State vector - RHS
      integer(sik)                    :: n2(s2%nvars) ! State vector - LHS
      character(len=max_state_string) :: n1str        ! State names (LHS)
      !
      integer(sik)                    :: nvib2        ! Number of active modes - RHS
      real(rk)                        :: evib1, evib2 ! Harmonic vibrational energies
      real(rk)                        :: fc           ! Integral
      integer(ik)                     :: info
      logical                         :: do_dc        ! True if derivative coupling is needed
      !
      call TimerStart('Evaluate left integrals')
      !
      do_dc = (coupling_name/=' ')
      !
      if (do_dc) then
        write (out,"(/'WARNING: Evaluating derivative coupling integrals, NOT the Franck-Condon integrals')")
      end if
      write (out,"()")
      write (out,"(a18,1x,a7,1x,a7,t36,a8,t68,a8)") &
                 '# <s2|s1>', 'e2,cm-1', 'e1,cm-1', 'state2', 'state1', &
                 '# -------', '-------', '-------', '------', '------'
      !
      read_loop: do
        read(input,"(a)",iostat=info) n1str
        if (info/=0) exit read_loop
        n1    = ascii2state(trim(n1str),s1%nvars)
        evib1 = h2cm*state_energy(s1,n1,zpe=.false.)
        !
        excitation_s2: do nvib2=0,nvib_max
          if(verbose>=2) then
            write (out,"('# LHS quanta = ',i3)") nvib2
          end if
          !
          !  State enumeration starts with the maximum number of quanta in the
          !  last vibrational mode
          !
          n2 = 0 ; n2(s2%nvars) = nvib2
          states_n2: do
            evib2 = h2cm*state_energy(s2,n2,zpe=.false.)
            if (evib2<=evib_max) then
              call initialize_cache(n2,n1,nvib2,0_sik)
              if (do_dc) then
                fc = evaluate_derivative_coupling(n2,n1)
              else
                fc = evaluate_fc_integral(n2,n1,1._rk)
              end if
              if ( (abs(fc)>report_threshold) .or. (verbose>=3) ) then
                write (out,"(f18.13,1x,f7.1,1x,f7.1,t36,a30,t68,a30)") &
                       fc, evib2, evib1, state2ascii(n2), state2ascii(n1)
              end if
            end if ! evib2<=evib_max)
            if (increment_state(nvib2,n2(s2%nvarsDummy+1:))) exit states_n2
          end do states_n2
        end do excitation_s2
      end do read_loop
      !
      !  It is not necessary to relase caches inside the loop - initialize_cache()
      !  is smart enough to handle reallocation when needed
      !
      call release_cache
      call TimerStop('Evaluate left integrals')
    end subroutine evaluate_left_fc_recursions
    !
    !  Compute specific FC overlaps. This could be extremely inefficient!
    !
    subroutine evaluate_one_fc_integral(task)
      character(len=*),intent(in)       :: task             ! 'test' will activate special mode.
      !
      character(len=2*max_state_string) :: tmp
      character(len=max_state_string)   :: n1str, n2str     ! State names
      integer(sik)                      :: n1(s1%nvars)     ! State vectors
      integer(sik)                      :: n2(s2%nvars)     ! State vectors
      integer(ik)                       :: info, pos
      real(rk)                          :: fc, evib1, evib2
      real(rk)                          :: refval
      logical                           :: do_dc            ! True if derivative coupling is needed
      !
      call TimerStart('Evaluate one integral')
      !
      do_dc = (coupling_name/=' ')
      !
      if (do_dc) then
        write (out,"(/'WARNING: Evaluating derivative coupling integrals, NOT the Franck-Condon integrals')")
      end if
      write (out,"()")
      write (out,"(a18,1x,a7,1x,a7,t36,a8,t68,a8)") &
                 '# <s2|s1>', 'e2,cm-1', 'e1,cm-1', 'state2', 'state1', &
                 '# -------', '-------', '-------', '------', '------'
      read_loop: do
        read(input,"(a)",iostat=info) tmp
        if (info/=0) exit read_loop
        pos = index(trim(tmp),' ',back=.true.)
        n1str = tmp(pos:) ; tmp(pos:) = ' '
        pos = index(trim(tmp),' ',back=.true.)
        n2str = tmp(pos:) ; tmp(pos:) = ' '
        if (task=='test') then
          read(tmp,*,iostat=info) refval
          if (info/=0) then
            write (out,"('Error ',i8,' parsing ',a,' as the reference value')") info, trim(tmp)
            stop 'franck_condon%evaluate_one_fc_integral - bad reference'
          end if
        else
          if (tmp/=' ') then
            write (out,"('Unclaimed token: ',a,' before state specifications (',a,') and (',a,')')") &
                   trim(tmp), trim(n2str), trim(n2str)
            stop 'franck_condon%evaluate_one_fc_integral - unclaimed token'
          end if
        end if
        n2 = ascii2state(trim(n2str),s2%nvars)
        n1 = ascii2state(trim(n1str),s1%nvars)
        evib1 = h2cm*state_energy(s1,n1,zpe=.false.)
        evib2 = h2cm*state_energy(s2,n2,zpe=.false.)
        call initialize_cache(n2,n1,0_sik,0_sik)
        if (do_dc) then
          fc = evaluate_derivative_coupling(n2,n1)
        else
          fc = evaluate_fc_integral(n2,n1,1._rk)
        end if
        if (task=='test') fc = fc - refval
        write (out,"(f18.13,1x,f7.1,1x,f7.1,t36,a30,t68,a30)") &
               fc, evib2, evib1, state2ascii(n2), state2ascii(n1)
      end do read_loop
      !
      !  It is not necessary to relase caches inside the loop - initialize_cache()
      !  is smart enough to handle reallocation when needed
      !
      call release_cache
      call TimerStop('Evaluate one integral')
    end subroutine evaluate_one_fc_integral
    !
    !  Main problem driver - this is the only externally-visible routine
    !
    subroutine run_fc
      integer(ik) :: info
      !
      call TimerStart('Franck-Condon')
      call accuracyInitialize
      report_threshold = spacing(100._rk)
      wgt_cutoff       = spacing(1._rk)
      !
      read(input,nml=fc_par,iostat=info)
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=fc_par)
      write (out,"()")
      !
      !  Load geometry and Hessian for each surface, and perform normal mode
      !  analysis for each of the surfaces.
      !
      write (out,"(/'Intializing surface S1 (right-hand side)'/)")
      call initialize_surface(surf=s1,xyz_name=s1_xyz_name,hess_name=s1_hess_name, &
              verbose=verbose,freq_cut=freq_cut,freq_dummy=freq_dummy,evib_max=evib_max,nvib_max=nvib_max)
      write (out,"(/'Intializing surface S2 (left-hand side)'/)")
      call initialize_surface(surf=s2,xyz_name=s2_xyz_name,hess_name=s2_hess_name, &
              verbose=verbose,freq_cut=freq_cut,freq_dummy=freq_dummy,evib_max=evib_max,nvib_max=nvib_max)
      !
      !  A bit of sanity checking - the surfaces must be compatible
      !
      call check_surface_sanity(s1,s2)
      !
      !  Derivative coupling vector
      !
      if (coupling_name/=' ') call read_coupling_vector
      !
      !  Prepare various temporary quantities needed by the Doktorov's 
      !  recursions
      !
      call prepare_for_integral_recursion
      !
      !  Evaluate primitive overlap integral
      !
      call evaluate_primitive_integral
      !
      !  Ready to evaluate FC integrals now
      !
      select case (task)
        case default
          write (out,"('Task ""',a,'"" is not understood')") trim(task)
          stop 'franck_condon%run_fc'
        case ('all')
          call evaluate_all_fc_recursions
        case ('one','test')
          call evaluate_one_fc_integral(task)
        case ('left')
          call evaluate_left_fc_recursions
        case ('right')
          call evaluate_right_fc_recursions
      end select
      !
      write (out,"(/'Total number of integral evaluations = ',f18.0)") integral_count
      write (out,"( ' Total number of integral cache hits = ',f18.0)") integral_cached
      if (use_cache) then
        write (out,"( '      Maximum size of integral cache = ',f18.0)") max_cache_used
        write (out,"( '              Number of cache resets = ',f18.0)") cache_flushes
      end if
      write (out,"()")
      !
      call TimerStop('Franck-Condon')
      call TimerReport
    end subroutine run_fc

  end module franck_condon
  !
  program dofc
    use franck_condon
    !
    call run_fc
  end program dofc
