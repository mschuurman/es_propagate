!
!  Evaluation of short-time autocorrelation functions for wavepacket
!  evolution on a harmonic surface. The result is intended for estimation
!  of nuclear contributions to high harmonics spectra.
!
!  This is a very naive implementation - numerical efficiency was not the goal
!
!  Although the code is OpenMP-parallel, the scalability is very modest,
!  Do not expect improvements beyond 2 CPUs. Depending on the hardware
!  and the compiler, there is a possibility for a race condition leading
!  to data corruption. See comments in a_mu() below for more information.
!
!  Version history:
!
!   January 22, 2009 - Initial version
!   January 30, 2009 - OpenMP parallel version
!   February 18, 2013 - Repurposed autocorrelation.f90 code to propagate wp in presence
!                       of an external field
!
  module absorb_propagate 
    use accuracy
    use elements
    use fc_tools
    use pulse
    use cache_grx
    use lapack
    use math
    use timer
    implicit none
    private
    public run_absorb_propagate
    !
    integer(ik)                  :: unit_scr         = 40  ! Any unused unit
    integer(ik)                  :: max_state_string = 256 ! Max length of the string describing a state
    integer(ik)                  :: ntstep           = 4   ! Number of time steps to store (temporary)
    complex(rk)                  :: zeroc=(0._rk,0._rk), im=(0._rk,1._rk)
    !
    !  User-defined parameters, controlled through AC_PAR namelist
    !
    integer(ik)        :: verbose          = 0        ! Verbosity level of the output
    character(len=80)  :: xyz              = ' '      ! File containing molecular geometry
    character(len=80)  :: gs_hessian       = ' '      ! File containing Hessian of the initial state
    character(len=80)  :: es_hessian       = ' '      ! File containing Hessian of the perturbed state
    character(len=80)  :: es_gradient      = ' '      ! File containing gradient of the perturbed state
    character(len=80)  :: tdip_gradient    = ' '      ! File containing the gradient of trans. dipole (x,y,z)
    real(rk)           :: freq_cut         = 50._rk   ! Normal modes with lower frequencies (cm^-1) are 
                                                      ! considered to be translations and rotations
    real(rk)           :: freq_dummy       = 10._rk   ! Replace very soft modes with this freqeuncy
    integer(ik)        :: packet_length    = 1        ! Number of basis functions in the initial wavepacket
    real(rk)           :: gh_cutoff        = 1e-3_rk  ! G and H coefficients smaller than this fraction of the
                                                      ! largest coefficient will be set to zero.
    real(rk)           :: constant_cutoff  = 1e-12_rk ! do not proprage coefficients less than this value
    real(rk)           :: coupling_cutoff  = 1e-12_rk ! Neglect smaller couplings
    integer(ik)        :: max_nvib(100)    = 0       ! Maximum number of vibrational functions in each quanta                
    integer(ik)        :: nvib_cutoff      = -1       ! Amplitudes farther away than nvib_cutoff quanta from
                                                      ! the target set will be neglected. -1 means no cutoff
    real(rk)           :: evib_cutoff      = -1       ! States of E > evib_cutoff will be neglected, where
     
    real(rk)           :: delta_e          = 0._rk    ! Electronic energy difference between two states
    real(rk)           :: tstep            = 5._rk    ! time step for propagation (a.u.)
    real(rk)           :: tfinal           = 400._rk  ! time to propagate (a.u.)
    real(rk)           :: pulse_fwhm       = 1._rk    ! fwhm of envelope
    real(rk)           :: pulse_amp        = 1._rk    ! amplitude of pulse envelope
    real(rk)           :: pulse_freq       = 1._rk    ! pulse frequency in au
    real(rk)           :: pulse_dir(3)     = 0._rk    ! pulse polarization
    real(rk)           :: tdip0(3)         = 0._rk    ! transition dipole moment at origin


    !
    !  Everything below is computed internally
    !
    type(surfaceData)         :: gs, es               ! ground and excited state surfaces
    real(rk), allocatable     :: grad_xyz(:)          ! Cartesian gradient vector on the perturbed surface
    real(rk), allocatable     :: tdip_xyz(:,:)        ! cartesian gradient of transition dipole matrix element
    real(rk), allocatable     :: g_bar(:)             ! Gradient perturbation
    real(rk), allocatable     :: h_bar(:,:)           ! Curvature perturbation
    real(rk), allocatable     :: td_bar(:,:)          ! 1st order geometry dependent transition dipole correction
    integer(ik)               :: n_coef(2)             ! number of non-zero coefficients on each state
    integer(sik), allocatable :: n_a0(:,:)            ! List of functions in the initial wavepacket
    integer(sik), allocatable :: s_a0(:)              ! state on which to put basis function
    integer(sik), allocatable :: n_max(:)             ! maximum number of functions in each mode
    complex(rk), allocatable  :: a0(:)                ! Function weights in the initial wavepacket
    type(grCacheWrapper)      :: cache              ! Cache of the expansion coefficients for initial state
    type(grCacheWalker)       :: walker             ! Structure to walk through cache of "a" coefficients
    real(ark)                 :: cache_lookups = 0    ! Number of cache lookups for Taylor coefficients
    real(ark)                 :: cache_misses  = 0    ! Number of cache misses for Taylor coefficients
    real(ark)                 :: evaluations = 0      ! Number of evaluations for Taylor coefficients
    real(ark)                 :: neglected_amps = 0   ! Number of amplitudes neglected because of nvib_cutoff
    type(EFpulse)              :: FPulse               ! pulse information 

    !
    namelist /AP_PAR/                                                          &
       verbose,                                                                &
       freq_cut, freq_dummy, gh_cutoff, coupling_cutoff, constant_cutoff,      &
       xyz, gs_hessian, es_hessian, es_gradient, tdip_gradient,                &
       packet_length, pulse_dir, pulse_fwhm, pulse_freq, tdip0, tstep, tfinal, &
       delta_e, max_nvib,evib_cutoff, nvib_cutoff
    !
    contains
    !
    subroutine read_gradient
      integer(ik) :: ios
      !
      allocate (grad_xyz(es%nvars))
      !
      open (unit_scr,file=trim(es_gradient),form='formatted',action='read',status='old',iostat=ios)
      if (ios/=0) then
        write (out,"(' Error ',i6,' opening file ',a)") ios, trim(es_gradient)
        stop 'autocorrelation%read_gradient - open'
      end if
      read (unit_scr,*,iostat=ios) grad_xyz
      if (ios/=0) then
        write (out,"(' Error ',i6,' reading file ',a)") ios, trim(es_gradient)
        stop 'autocorrelation%read_gradient - read'
      end if
      close (unit_scr,iostat=ios)
      if (ios/=0) then
        write (out,"(' Error ',i6,' closing file ',a)") ios, trim(es_gradient)
        stop 'autocorrelation%read_gradient - open'
      end if
      if (verbose>=2) then
        call print_vector(grad_xyz,'Cartesian gradient')
      end if
    end subroutine read_gradient

    !
    subroutine read_tdip
      integer(ik)  :: i,ios
      character(1) :: label
      !
      allocate (tdip_xyz(es%nvars,3))
      !
      open (unit_scr,file=trim(tdip_gradient),form='formatted',action='read',status='old',iostat=ios)
      if (ios/=0) then
        write (out,"(' Error ',i6,' opening file ',a)") ios, trim(tdip_gradient)
        stop 'autocorrelation%read_tdip - open'
      end if
      read (unit_scr,*,iostat=ios) tdip_xyz(:,:)
      if (ios/=0) then
        write (out,"(' Error ',i6,' reading file ',a)") ios, trim(tdip_gradient)
        stop 'autocorrelation%read_tdip - read'
      end if
      close (unit_scr,iostat=ios)
      if (ios/=0) then
        write (out,"(' Error ',i6,' closing file ',a)") ios, trim(tdip_gradient)
        stop 'autocorrelation%read_tdip - open'
      end if
      if (verbose>=2) then
        print_loop: do i = 1,3
         write(label,'(i1)')i
         call print_vector(tdip_xyz(:,i),'Cartesian gradient - transition dipole, axis '//label)
        end do print_loop
      end if
    end subroutine read_tdip

    !
    !  Gradient perturbation vector
    !
    subroutine build_gbar
      integer(ik) :: im, ic
      real(rk)    :: gv
      !
      call TimerStart('Build G')
      allocate (g_bar(es%nvars))
      !
      g_bar = 0._rk
      !
      modes_loop: do im=gs%nvarsDummy+1,gs%nvars
        gv = 0._rk
        cartesian_loop: do ic=1,gs%nvars
          gv = gv + grad_xyz(ic) * gs%modes(ic,im) / sqrt(gs%mass(ic))
        end do cartesian_loop
        g_bar(im) = 0.5_rk * gv / sqrt(gs%freq(im))
      end do modes_loop
      !
      call TimerStop('Build G')
    end subroutine build_gbar

    !
    !  Transition dipole perturbation vector
    !
    subroutine build_tdbar
      integer(ik) :: im, ic, ia
      real(rk)    :: tdv
      !
      call TimerStart('Build gTD')
      allocate (td_bar(es%nvars,3))
      !
      td_bar = 0._rk
      !
      modes_loop: do im=gs%nvarsDummy+1,gs%nvars
        axis_loop: do ia=1,3
         tdv = 0._rk
         cartesian_loop: do ic=1,gs%nvars
          tdv = tdv + tdip_xyz(ic,ia) * gs%modes(ic,im) / sqrt(gs%mass(ic))
         end do cartesian_loop
         td_bar(im,ia) = 0.5_rk * tdv / sqrt(gs%freq(im))
       end do axis_loop
      end do modes_loop
      !
      call TimerStop('Build gTD')
    end subroutine build_tdbar

    !
    !  Curvature perturbation matrix
    !
    subroutine build_hbar
      integer(ik) :: im1, im2, ic1, ic2
      real(rk)    :: hv, dc
      !
      call TimerStart('Build H')
      allocate (h_bar(es%nvars,es%nvars))
      !
      h_bar = 0._rk
      !
      m1_loop: do im1=gs%nvarsDummy+1,gs%nvars
        m2_loop: do im2=gs%nvarsDummy+1,gs%nvars
          hv = 0._rk
          c1_loop: do ic1=1,gs%nvars
            c2_loop: do ic2=1,gs%nvars
              dc = (es%hess(ic2,ic1) - gs%hess(ic2,ic1)) / sqrt(gs%mass(ic2)*gs%mass(ic1))
              hv = hv + dc * gs%modes(ic2,im2) * gs%modes(ic1,im1)
            end do c2_loop
          end do c1_loop
          h_bar(im2,im1) = 0.25_rk * hv / sqrt(gs%freq(im2)*gs%freq(im1))
        end do m2_loop
      end do m1_loop
      !
      call TimerStop('Build H')
    end subroutine build_hbar

!
!
!
    subroutine screen_gbar_hbar
      real(rk)    :: eps_cut
      real(rk)    :: g_max, h_max, td_max
      integer(ik) :: ng_orig, ng_cut, nh_orig, nh_cut, nt_orig, nt_cut
      !
      ng_orig = count(abs(g_bar)>0)
      nh_orig = count(abs(h_bar)>0)
      nt_orig = count(abs(td_bar)>0)
      g_max   = maxval(abs(g_bar))
      h_max   = maxval(abs(h_bar))
      td_max  = maxval(abs(td_bar))
      eps_cut = gh_cutoff * max(g_max,h_max,td_max)
      where (abs(g_bar)<eps_cut)
        g_bar = 0
      end where
      where (abs(h_bar)<eps_cut)
        h_bar = 0
      end where
      where (abs(td_bar)<eps_cut)
        td_bar = 0
      end where
      ng_cut  = count(abs(g_bar)>0)
      nh_cut  = count(abs(h_bar)>0)
      nt_cut  = count(abs(td_bar)>0)
      !
      if (verbose>=0) then
        write (out,"()")
        write (out,"('          Largest absolute gradient coupling: ',g12.6)") g_max
        write (out,"('           Largest absolute hessian coupling: ',g12.6)") h_max
        write (out,"('  Largest absolute transition dipole element: ',g12.6)") td_max
        write (out,"('                  Coupling neglect threshold: ',g12.6)") eps_cut
        write (out,"('Non-zero gradient couplings before screening: ',i0)") ng_orig
        write (out,"('                        ...  after screening: ',i0)") ng_cut
        write (out,"(' Non-zero hessian couplings before screening: ',i0)") nh_orig
        write (out,"('                        ...  after screening: ',i0)") nh_cut
        write (out,"('  Non-zero t. dip. elements before screening: ',i0)") nt_orig
        write (out,"('                        ...  after screening: ',i0)") nt_cut
        write (out,"()")
      end if
    end subroutine screen_gbar_hbar
    !
    !
    !
    subroutine report_coupling_vectors
      integer(ik) :: im, imc, ml, hc
      real(rk)    :: freq, grad, hess, eps_g, eps_h
      logical     :: mask(gs%nvars)
      !
      write (out,"(/t5,'Significant gradient and curvature coupling vectors'/)")
      write (out,"(1x,a4,2x,a8,2x,a10,2x,a4,1x,a10)") &
                 'Mode', 'cm^-1', 'Grad', 'Mode', 'Hess', &
                 '----', '-----', '----', '----', '----'
      !
      eps_g = 1e3_rk * spacing(maxval(abs(g_bar)))
      eps_h = 1e3_rk * spacing(maxval(abs(h_bar)))
      !
      mode_loop: do im=gs%nvarsDummy+1,gs%nvars
        freq = h2cm * gs%freq(im)
        grad = g_bar(im)
        if (abs(grad)<=eps_g .and. maxval(abs(h_bar(:,im)))<=eps_h) cycle mode_loop
        !
        !  Something to report
        !
        write (out,"(1x,i4,2x,f8.2,2x,g10.4,1x)",advance='no') im, freq, grad
        mask = .true.
        hc   = 0
        hess_loop: do imc=1,gs%nvars
         ml = maxloc(abs(h_bar(:,im)),1,mask)
         hess = h_bar(ml,im)
         if (abs(hess)<=eps_h) exit hess_loop
         if (hc>=5) then
           hc = 0
           write (out,"()") 
           write (out,"(t29)",advance='no')
         end if
         hc = hc + 1
         write (out,"(1x,i4,1x,g10.4,1x)",advance='no') ml, hess
         mask(ml) = .false.
        end do hess_loop
        write (out,"()")
      end do mode_loop
      !
      write (out,"()")
      write (out,"(' Gradient and curvature coupling coefficients are in units of energy (Hartree)')")
      write (out,"()")
    end subroutine report_coupling_vectors
    !
    !  Perturbed-Hamiltonian matrix elements between two unperturbed vibrational states
    !  Matrix elements are non-zero if the two states differ by at most two quanta,
    !  and vanish identically otherwise.
    !
    function c_mu_nu(mu,nu) result(c)
      integer(sik), intent(in) :: mu(:), nu(:) ! State labels for which matrix element is needed
      real(rk)                 :: c            ! The matrix element
      !
      integer(ik) :: d_mu_nu   
      integer(ik) :: i, j
      integer(ik) :: max_i, max_j
      !
      !  Is coupling possible?
      !
      d_mu_nu = sum(abs(mu-nu))
      if (d_mu_nu>2) then
        c = 0._rk
        return ! Do not report anything if the result is zero by symmetry
      end if
      !
      !  The result -could- be non-zero; will report below
      !
      coupling: do
        !
        !  mu and nu can couple. Decide on the specific term(s)
        !
        if (d_mu_nu==0) then
          !
          !  State energy change. All modes will contribute.
          !
          c = 0._rk
          energy_loop: do i=1,gs%nvars
            c = c + 0.5_rk * h_bar(i,i) * 4._rk * (mu(i) + 0.5_rk)
          end do energy_loop
          exit coupling
        end if
        !
        scan_differences_1: do i=1,gs%nvars
          if (mu(i)/=nu(i)) exit scan_differences_1
        end do scan_differences_1
        max_i = max(mu(i),nu(i))
        !
        if (d_mu_nu==1) then
          !
          !  Gradient term
          !
          c = g_bar(i) * sqrt(2._rk*max_i)
          exit coupling
        end if
        !
        !  Off-diagonal Hessian terms
        !
        if (abs(mu(i)-nu(i))==2) then
          !
          !  Two quanta in the same mode:
          !
          c = 0.5_rk * h_bar(i,i) * 2._rk * sqrt(max_i * (max_i-1._rk))
          exit coupling
        end if
        !
        !  One quantum each in two different modes: find the remaining quantum first
        !
        scan_differences_2: do j=i+1,gs%nvars
          if (mu(j)/=nu(j)) exit scan_differences_2
        end do scan_differences_2
        max_j = max(mu(j),nu(j))
        !
        c = 0.5_rk * (h_bar(i,j)+h_bar(j,i)) * sqrt(4._rk * max_i * max_j)
        exit coupling
      end do coupling
      !
      if (verbose>=3) then
        write (out,"('C= ',g18.12,' mu= ',a30,' nu= ',a30)") c, trim(state2ascii(mu)), trim(state2ascii(nu))
      end if
      return

    end function c_mu_nu
    !
    !  Transition dipole matrix elements between two unperturbed vibrational states
    !  Matrix elements are non-zero if the two states differ by at most one quanta,
    !  and vanish identically otherwise
    !
    function d_mu_nu(mu,nu,F) result(d)
      integer(sik), intent(in) :: mu(:), nu(:) ! State labels for which matrix element is needed
      real(rk)                 :: F(3)         ! magnitude of field at current time
      real(rk)                 :: d            ! The matrix element
      !
      integer(ik) :: df_mu_nu
      integer(ik) :: i, j
      integer(ik) :: max_i, max_j
      !
      !  Is coupling possible?
      !
      df_mu_nu = sum(abs(mu-nu))
      if (df_mu_nu>1) then
        d = 0._rk
        return ! Do not report anything if the result is zero by symmetry
      end if
      !
      !  The result -could- be non-zero; will report below
      !
      coupling: do
        !
        !  mu and nu can couple. Decide on the specific term(s)
        !
        if (df_mu_nu==0) then
          !
          !  State energy change. All modes will contribute.
          !
          d = sum( tdip0 * F )
          exit coupling
        end if
        !
        scan_differences_1: do i=1,gs%nvars
          if (mu(i)/=nu(i)) exit scan_differences_1
        end do scan_differences_1
        max_i = max(mu(i),nu(i))
        !
        if (df_mu_nu==1) then
          !
          !  Gradient term
          !
          d = sqrt(2._rk*max_i) * sum( td_bar(i,:) * F )
          exit coupling
        end if

      end do coupling
      !
      if (verbose>=3) then
        write (out,"('D= ',g18.12,' mu= ',a30,' nu= ',a30)") d, trim(state2ascii(mu)), trim(state2ascii(nu))
      end if
      return

    end function d_mu_nu    

    !
    !  process_mu_nu() would be a natural candidate for a contained function
    !  of evaluate_a_mu(). Unfortunately, ifort seems to miscompile OpenMP
    !  loops with contained functions, so we have to do it the hard way.
    !
    function process_a_nu(mu,nu,cf,F,t) result(val)
      integer(sik), intent(in) :: mu(:), nu(:) ! States of interest
      complex(rk), intent(in)  :: cf
      real(rk),intent(in)      :: F(3),t
      !
      real(rk)    :: c, de
      complex(rk) :: val
      !
      val = zeroc

      c = d_mu_nu(mu,nu,F)
      if (abs(c)<coupling_cutoff) return ! No coupling - ignore mode
      de = delta_e + state_energy(gs,mu,zpe=.false.) - state_energy(gs,nu,zpe=.false.)
      val = cf * c * exp( -im * de * t)

    end function process_a_nu


    !
    !  process_mu_nu() would be a natural candidate for a contained function
    !  of evaluate_a_mu(). Unfortunately, ifort seems to miscompile OpenMP
    !  loops with contained functions, so we have to do it the hard way.
    !
    function process_b_nu(mu,nu,cf,typ,F,t) result(val)
      integer(sik), intent(in) :: mu(:), nu(:) ! States of interest
      complex(rk), intent(in)  :: cf
      integer(ik), intent(in)  :: typ          ! Order of the expansion coefficient
      real(rk),intent(in)      :: F(3),t
      !
      real(rk)    :: c, de
      complex(rk) :: val
      !
      val = zeroc
      select case(typ)

       case(0)       ! this is a second order potential perturbation term 
        c = c_mu_nu(mu,nu)
        if (abs(c)<coupling_cutoff) return ! No coupling - ignore mode
        de = state_energy(gs,nu,zpe=.false.) - state_energy(gs,mu,zpe=.false.)

       case(1)       ! this is a dipole matrix element term
        c = d_mu_nu(mu,nu,F)       
        if (abs(c)<coupling_cutoff) return ! No coupling - ignore mode
        de = delta_e + state_energy(gs,nu,zpe=.false.) - state_energy(gs,mu,zpe=.false.)

      end select

      val = cf * c * exp( im * de * t)
!      write(out,"('nu=',4(i2),' mu=',4(i2),' c=',f12.8,' cf=',f12.8,' val=',f12.8,' t=',f6.4)")nu(6:9),mu(6:9),c,abs(cf),abs(val),t

      return

    end function process_b_nu
 
    !
    !  Taylor coefficients for time evolution of the nuclear wavepacket
    !
    subroutine accumulate_at(mu,cf,F,t)
      integer(sik), intent(in) :: mu(:)  ! State of interest
      complex(rk),intent(in)   :: cf     ! Expansion coefficient
      real(rk),intent(in)      :: F(3),t
      !
      integer(sik) :: nu(gs%nvars)       ! States at lower expansion orders
      integer(ik)  :: c1, c2             ! Positions of the vibrational quanta
      integer(sik) :: n1, n1_min, n1_max ! Number of quanta at c1 position
      integer(sik) :: n2, n2_min, n2_max ! ... and c2 position
      complex(rk)  :: at
      type(grCache), pointer   :: cgr
      !
      !
      !
      !  Begin by looking at the self-coupling
      !
      nu = mu
      at = process_a_nu(mu,nu,cf,F,t)
      if(abs(at) > 0_rk) call update_cf(nu,1_ik,1_ik,at,ovrwrite=.false.)

      !
      quantum_c1: do c1=gs%nvarsDummy+1,gs%nvars
        !
        !  Only states differing by 0/1 modes will contribute to propagation of a(t) 
        !
        n1_min = max(mu(c1)-1_sik,0_sik)
        n1_max = min(mu(c1)+1_sik,n_max(c1))
        coupling_n1_only: do n1=n1_min,n1_max
          if (n1==mu(c1)) cycle coupling_n1_only ! Self was already listed
          nu(c1) = n1
          at = process_a_nu(mu,nu,cf,F,t)
         if(abs(at) > 0_rk)call update_cf(nu,1_ik,1_ik,at,ovrwrite=.false.) 
        end do coupling_n1_only
        nu(c1) = mu(c1)  ! Restore state vector
      end do quantum_c1
      return
    end subroutine accumulate_at

    !
    !
    !
    subroutine accumulate_bt(mu,cf,typ,F,t)
      integer(sik), intent(in) :: mu(:)  ! State of interest
      complex(rk),intent(in)   :: cf     ! Expansion coefficient
      integer(ik), intent(in)  :: typ     ! term type
      real(rk),intent(in)      :: F(3),t ! external field, time
      !
      integer(sik) :: nu(gs%nvars)       ! States at lower expansion orders
      integer(ik)  :: c1, c2             ! Positions of the vibrational quanta
      integer(sik) :: n1, n1_min, n1_max ! Number of quanta at c1 position
      integer(sik) :: n2, n2_min, n2_max ! ... and c2 position
      complex(rk)  :: bt,bt2
      type(grCache), pointer   :: cgr

      !
      !  We can only pick up population from states at most 2 quanta away.
      !  Let's enumerate them. We have to be careful not to list the same
      !  state more than one!
      !
      !
      !
      !  Begin by looking at the self-coupling
      !
      nu = mu
      bt = process_b_nu(mu,nu,cf,typ,F,t)
      if(abs(bt) > 0._rk)call update_cf(nu,2_ik,1_ik,bt,ovrwrite=.false.)
      !
      quantum_c1: do c1=gs%nvarsDummy+1,gs%nvars
        !
        !  Then, enumerate states where only one mode differs
        !
        n1_min = max(mu(c1)-2_sik,0_sik)
        n1_max = min(mu(c1)+2_sik,n_max(c1))
        coupling_n1_only: do n1=n1_min,n1_max
          if (n1==mu(c1)) cycle coupling_n1_only ! Self was already listed
          nu(c1) = n1
          bt = process_b_nu(mu,nu,cf,typ,F,t)
          if(abs(bt) > 0._rk) call update_cf(nu,2_ik,1_ik,bt,ovrwrite=.false.)
        end do coupling_n1_only
        !
        !
        !
        n1_min = max(mu(c1)-1_sik,0_sik)
        n1_max = min(mu(c1)+1_sik,n_max(c1))
        coupling_n1: do n1=n1_min,n1_max
          if (n1==mu(c1)) cycle coupling_n1
          nu(c1) = n1
          quantum_c2: do c2=c1+1,gs%nvars
            n2_min = max(mu(c2)-1_sik,0_sik)
            n2_max = min(mu(c2)+1_sik,n_max(c2))
            coupling_n2: do n2=n2_min,n2_max
              if (n2==mu(c2)) cycle coupling_n2
              nu(c2) = n2
              bt = process_b_nu(mu,nu,cf,typ,F,t)
              if(abs(bt) > 0._rk)call update_cf(nu,2_ik,1_ik,bt,ovrwrite=.false.)
            end do coupling_n2
            nu(c2) = mu(c2)
          end do quantum_c2
        end do coupling_n1
        nu(c1) = mu(c1)  ! Restore state vector
      end do quantum_c1

      return
    end subroutine accumulate_bt

    !
    !  Try to look up expansion coefficient in the cache, evaluating it if
    !  it was not cached before.
    !
    function retrieve_cf(mu,st,slot) result(ab)
      integer(sik), intent(in) :: mu(:)  ! State of interest
      integer(ik), intent(in)  :: st     ! state 1=gs, 2=es
      integer(ik), intent(in)  :: slot   ! 0=a(t), 1=vel, 2=adot(t), 3=adot(t+1)
      complex(rk)              :: ab      ! Expansion coefficient 
      type(grCache), pointer   :: cgr

      !
      !  Should this amplitude be neglected?
      !
      if(neglect_state(mu)) then
       ab = zeroc
       return

      else

       cgr => cache_gr_locate(cache,mu)

       ! write (out,"('Request ',a30,' ord ',i3,' found= ',l1)") trim(state2ascii(mu)), s, associated(cgr)
       if (associated(cgr)) then
        !
        !  The test below is a little convoluted, since we have no way of
        !  ensuring atomic assignments to complex cache entries. The test
        !  as written now should be safe as long as real assignment is
        !  atomic. If it is not, the race condition will occur, leading
        !  to data corruption. There is absolutely nothing we could do
        !  about it short of locking all cache updates and completely
        !  killing both the scalability and the serial performance.
        !
!        if (real (cgr%values(slot,st),kind=rk)/=cache_gr_invalid_r .and. &
!            aimag(cgr%values(slot,st))        /=cache_gr_invalid_r ) then
          cache_lookups = cache_lookups + 1
          ab = cgr%values(slot,st)
          return
!        else
!          ! take this opportunity to initialize to zero
!          cgr%values = zeroc
!          ab = cgr%values(slot,st)
!          return
!        endif
       else
         write(out,"('request out of bounds')")
         cache_misses = cache_misses + 1
         ab = zeroc
         return
       end if

      endif

    end function retrieve_cf

    !
    !  Try to look up expansion coefficient in the cache, evaluating it if
    !  it was not cached before.
    !
    subroutine update_cf(mu,st,slot,ab,ovrwrite)
      integer(sik), intent(in) :: mu(:)  ! State of interest
      integer(ik), intent(in)  :: st     ! state 1=gs, 2=es
      integer(ik), intent(in)  :: slot   ! 0=a(t), 1=vel, 2=adot(t), 3=adot(t+1)
      complex(rk), intent(in)  :: ab     ! contribution
      logical, intent(in)      :: ovrwrite    ! 0=set to value ab, 1=add value ab
      type(grCache), pointer   :: cgr

      if(.not.neglect_state(mu)) then
        
       cgr => cache_gr_locate(cache,mu)

       ! write (out,"('Request ',a30,' ord ',i3,' found= ',l1)") trim(state2ascii(mu)), s, associated(cgr)
       if (associated(cgr)) then
         !
           cache_lookups = cache_lookups + 1
           if(ovrwrite) then
            cgr%values(slot,st) = ab
           else
            cgr%values(slot,st) = cgr%values(slot,st) + ab
           endif
           return
       else
         write(out,"('should not be here...')")
         cache_misses = cache_misses + 1
         return
       endif

      endif

    end subroutine update_cf

    !
    !
    !
    function neglect_state(mu) result(neglect)
      implicit none
      integer(sik), intent(in) :: mu(:)  ! State of interest
      logical                  :: neglect 
      integer(sik)             :: nvib(packet_length)
      integer(ik)              :: is
      real(rk)                 :: evib

      neglect = .false.

      !  ...if any modes have more quanta than n_max
      if (any(n_max-mu < 0)) then
        neglected_amps = neglected_amps + 1
        neglect = .true.
        return
      endif
      ! ...if state differs more than nvib_cutoff quanta from closest
      ! reference state
      if (nvib_cutoff>=0) then
        do is = 1,packet_length
         nvib(is) = sum(abs(mu-n_a0(:,is)))
        enddo
        if (minval(nvib)>nvib_cutoff) then
          neglected_amps = neglected_amps + 1
          neglect = .true.
          return
        end if
      end if
      ! ...if state differs in energy more than evib_cutoff from
      !  nearest reference state
      if (evib_cutoff>=0) then
        evib = state_energy(gs,mu,zpe=.false.)
        if (evib>evib_cutoff) then
          neglected_amps = neglected_amps + 1
          neglect = .true.
          return
        end if
      end if

      return
    end function neglect_state

    !
    !  Read the initial vibrational wavepacket from input
    !
    subroutine read_wavepacket
      integer(ik)                       :: is    ! Vibrational state index
      character(len=2*max_state_string) :: tmp   ! I/O buffer
      character(len=max_state_string)   :: str   ! Vibrational function name
      integer(ik)                       :: pos
      integer(ik)                       :: info
      !
      allocate (n_a0(gs%nvars,packet_length),s_a0(packet_length),a0(packet_length),n_max(gs%nvars))
      !
      read_loop: do is=1,packet_length
        read(input,"(a)",iostat=info) tmp
        if (info/=0) then
          write (out,"('Run out of input reading wavepacket contribution ',i5)") is
          stop 'autocorrelation%read_wavepacket - out of states'
        end if
        ! read in the state vector
        pos = index(trim(tmp),' ',back=.true.)
        str = tmp(pos:) ; tmp(pos:) = ' '
        n_a0(:,is) = ascii2state(trim(str),gs%nvars)
        if (any(n_a0(:gs%nvarsDummy,is)/=0)) then
          write (out,"('Wavepacket component ',i0,' (',a,') requests excitation " // &
                 "of a translational/rotational mode')") is, trim(state2ascii(n_a0(:,is)))
          stop 'autocorrelation%read_wavepacket - bad excitation'
        end if
        ! read in the state coefficient
        pos = index(trim(tmp),' ',back=.true.)
        str = tmp(pos:) ; tmp(pos:) = ' ' 
        read(str,*,iostat=info) a0(is)
        if (info/=0) then
          write (out,"('Error ',i5,' converting ',a,' to a complex weight of contribution ',i5)") &
                 info, trim(tmp), is
        stop 'autocorrelation%read_wavepacket - bad weight'
        end if
        ! read in the identity of the state   
        pos = index(trim(tmp),' ',back=.true.)
        str = tmp(pos:) ; tmp(pos:) = ' ' 
        read(str,*,iostat=info) s_a0(is)
        if (info/=0) then
          write (out,"('Error ',i5,' converting ',a,' to an integer state index, packet ',i5)") &
                 info, trim(tmp), is
        stop 'autocorrelation%read_wavepacket - bad state'
        end if

      end do read_loop
      !
      !  Initialize max. excitation levels
      !
      n_max = 0_ik
      n_max(:) = max_nvib(:gs%nvars)

      !
      write (out,"(/'Norm of the initial wavepacket is ',g20.14/)") sum(abs(a0)**2)
    end subroutine read_wavepacket
    !
    !
    !
    subroutine initialize_pulse
     implicit none

     call initializePulse(FPulse,pulse_fwhm,pulse_freq)
     if(tfinal < 0) tfinal = 2._rk * FPulse%toff2

    end subroutine initialize_pulse
    !
    !
    !
    subroutine initialize_cache
      implicit none
      integer(ik)       :: is,s
      type(grCache), pointer   :: cgr
      !
      call cache_gr_initialize(cache,ntstep,n_max)

      do is = 1,packet_length
        call update_cf(n_a0(:,is),s_a0(is),0_ik,a0(is) / sqrt(dot_product(a0,a0)),ovrwrite=.true.)
      enddo 
       
    end subroutine initialize_cache
    !
    !
    !
    subroutine cache_report
      real(ark) :: entries
      !
      call cache_gr_destroy(cache,entries)
      !
      write (out,"(/' Number of amplitude evaluations  = ',f20.0 )") evaluations
      write (out,"( ' Number of amplitude cache hits   = ',f20.0 )") cache_lookups
      write (out,"( ' Number of amplitude cache misses = ',f20.0 )") cache_misses
      write (out,"( '   Number of amplitudes neglected = ',f20.0 )") neglected_amps
      write (out,"( ' Cache size (nodes)               = ',f20.0 )") entries
      write (out,"( ' Cache size (bytes)              >= ',f20.0/)") entries * (ntstep+1) * 2 * rk_bytes
    end subroutine cache_report

    !
    !
    !
    subroutine initialize_timestep
      implicit none
      integer(ik)              :: s
      integer(sik)             :: n(size(cache%n_max))
      type(grCache), pointer  :: cgr

      call walker_gr_initialize(cache,walker)
      scan_tree: do while(walker_gr_next(walker,n))
        scan_states: do s = 1,ntstep
         call update_cf(n,1_ik,s,zeroc,ovrwrite=.true.)
         call update_cf(n,2_ik,s,zeroc,ovrwrite=.true.)
        end do scan_states
      end do scan_tree
      call walker_gr_destroy(walker)     

    end subroutine initialize_timestep

    !
    ! add in contribution to all connected del_a/del_t terms
    !
    subroutine update_tderiv(F,t)
      implicit none
      real(rk),intent(in)      :: F(3),t

      integer(sik)             :: n(size(cache%n_max))
      complex(rk)              :: cf

      ! 
      ! a(t) will only contribute to propragation of b(t)
      !
      call walker_gr_initialize(cache,walker)
      scan_tree: do while(walker_gr_next(walker,n))

          ! determine contribution from a coefficients
          cf = retrieve_cf(n,1_ik,0_ik)
          if(abs(cf) > constant_cutoff ) then
           call accumulate_bt(n,cf,1_ik,F,t)
          endif

          ! determine contribution from b coefficients
          cf = retrieve_cf(n,2_ik,0_ik)
          if(abs(cf) > constant_cutoff ) then
           call accumulate_at(n,cf,F,t)
           call accumulate_bt(n,cf,0_ik,F,t)
          endif

      end do scan_tree
      call walker_gr_destroy(walker)

    end subroutine update_tderiv

    !
    ! update position using leapfrog integration
    !
    subroutine update_position
      implicit none
      integer(ik)              :: s
      integer(sik)             :: n(size(cache%n_max))
      complex(rk)              :: newcf
      complex(rk)                 :: nrma,nrmb
      type(grCache), pointer  :: cgr
      

      n_coef = 0
      call walker_gr_initialize(cache,walker)
      scan_tree: do while(walker_gr_next(walker,n))
        cgr => cache_gr_locate(cache,n)
        scan_states: do s = 1,2
         cgr%values(0,s) = cgr%values(0,s) - im * cgr%values(1,s)*tstep
         n_coef(s) = n_coef(s) + 1
        end do scan_states
      end do scan_tree
      call walker_gr_destroy(walker)

    end subroutine update_position

    !
    !
    !
    function wp_norm result(nrm)
      real(rk)                :: nrm(2)
      integer(ik)             :: s
      integer(sik)            :: n(size(cache%n_max))
      complex(rk)             :: cf

      nrm = zeroc 
      call walker_gr_initialize(cache,walker)
      scan_tree: do while(walker_gr_next(walker,n))
       scan_states: do s = 1,2
        cf = retrieve_cf(n,s,0_ik)
!        if(abs(cf) > 1e-10)write(out,"('state=',4(i2)' adds=',f12.8)")n(6:9),abs(conjg(cf)*cf)
!        if(abs(cf)**2 > 1)write(out,"('problem: ',12(i2),' |cf|=',f8.6)")n(7:18),abs(cf)
        nrm(s) = nrm(s) + abs(cf)**2 
       end do scan_states
      end do scan_tree
      call walker_gr_destroy(walker)
 
      do s = 1,2
       nrm(s) = sqrt(nrm(s))
      end do

      return
    end function wp_norm 


    !
    !  Calculate Taylor expansion for the evolution of the "interesting" part of the 
    !  wavepacket, which can overlap with the wavepacket on the unperturbed surface.
    !
    subroutine evolve_wavepacket
      integer(ik)             :: s,nstep
      real(rk)                :: t,F(3)
      complex(rk)             :: cf
      real(rk)                :: nrm(2),nrmt
      type(grCache), pointer  :: cgr
      integer(sik)            :: n(size(cache%n_max))
      !
      call TimerStart('Wavepacket evolution')

      t = 0._rk

      time_step: do
        t = t + tstep
        if(t > tfinal) exit

        ! 
        ! Update the external field
        !
        F = evaluateExternalField( FPulse,t-0.5*tstep ) * pulse_dir 

        nrm = wp_norm()
        nrmt = sqrt(nrm(1)**2 + nrm(2)**2)
        write(out,"('t= ',f9.4,' Ft= ',3(f8.5,1x),' |0|= ',f8.6,' |1|= ',f8.6,' |total|= ',f8.6)")t,F,nrm(1),nrm(2),nrmt

        ! 
        ! Set the slot in which we will accumulate right-hand side to zero
        !
        call initialize_timestep

        !
        ! Evalulate contribution of a(t) coefficients to updated left-hand side
        !  (note: will only contribute to propagation of b(t) coefficients)
        call update_tderiv(F,t-0.5_rk*tstep)
       
        ! 
        ! update a(t) and b(t)
        !
        call update_position

        !
        ! Update the velocity term
!        call update_velocity

      end do time_step


      if (verbose>=1) write (out,"()")
      call TimerStop('Wavepacket evolution')
      !
      if (verbose>=1) then
        scan_states: do s = 1,2
         call printCoefSort(s)
        end do scan_states 
      end if

      write (out,"()")
      write (out,"(' Norm Psi[0]: ',f12.8)")nrm(1)
      write (out,"(' Norm Psi[1]: ',f12.8)")nrm(2)

    end subroutine evolve_wavepacket

    ! 
    !  print sorted coefficients
    ! 
    subroutine printCoefSort(s)
     implicit none
     integer(ik),intent(in)   :: s
     integer(ik)              :: i,j,nprint
     complex(rk)              :: cf
     integer(sik)             :: n(size(cache%n_max))
     real(rk)                 :: print_cutoff = 1e-4
     type cfst
      complex(rk)              :: ab = (0._rk,0._rk)
      integer(sik),allocatable :: m(:)
     end type cfst 
     type(cfst),allocatable   :: sarr(:)

      write (out,"()")
      if(s == 1_ik) then
       write (out,"(a30)")'A coefficients at tfinal -----'
       write (out,"(1x,a21,2x,a12,2x,a30)")'a(tf)                  ','  |a(tf)|   ','harmonic basis function       '
      else
       write (out,"(a30)")'B coefficients at tfinal -----'
       write (out,"(1x,a21,2x,a12,2x,a30)")'b(tf)                  ','  |b(tf)|   ','harmonic basis function       '
      endif    

      allocate(sarr(n_coef(s)))
      do i = 1,n_coef(s)
       allocate(sarr(i)%m(size(cache%n_max)))
      end do

      nprint = 0_ik
      call walker_gr_initialize(cache,walker)
      scan_coef: do while(walker_gr_next(walker,n))
        cf = retrieve_cf(n,s,0_ik)
        if(abs(cf) > print_cutoff) then
          i = 1
          scan_list: do while( abs(cf) < abs(sarr(i)%ab) )
           i = i + 1
          end do scan_list
          ! shift existing elements down
          do j = nprint,i,-1
           sarr(j+1) = sarr(j)
          enddo
          sarr(i)%ab = cf
          sarr(i)%m = n
          nprint = nprint + 1
        endif
      end do scan_coef
      call walker_gr_destroy(walker)

      do i = 1,nprint
       write (out,"(1x,'(',f9.6,',',f9.6,')',2x,f12.8,2x,a30)")  &
                  real(sarr(i)%ab),aimag(sarr(i)%ab),abs(sarr(i)%ab),state2ascii(sarr(i)%m)
      enddo

      do i = 1,nprint
       deallocate(sarr(i)%m)
      enddo
      deallocate(sarr)
 
    end subroutine printCoefSort


   !
    !  Main problem driver - this is the only externally-visible routine
    !
    subroutine run_absorb_propagate
      integer(ik) :: info
      !
      call TimerStart('Autocorrelation')
      call TimerStart('Initialization')
      call accuracyInitialize
      !
      read(input,nml=ap_par,iostat=info)
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=ap_par)
      write (out,"()")
      !
      !  Load geometry and Hessian for each surface, and perform normal mode
      !  analysis for each of the surfaces.
      !
      write (out,"(/'Intializing neutral (unperturbed) surface'/)")
      call initialize_surface(surf=gs,xyz_name=xyz,hess_name=gs_hessian, &
              verbose=verbose,freq_cut=freq_cut,freq_dummy=freq_dummy,evib_max=1._rk,nvib_max=1_sik)
      write (out,"(/'Intializing cationic (perturbed) surface'/)")
      call initialize_surface(surf=es,xyz_name=xyz,hess_name=es_hessian, &
              verbose=verbose,freq_cut=freq_cut,freq_dummy=freq_dummy,evib_max=1._rk,nvib_max=1_sik)
      !
      !  A bit of sanity checking - the surfaces must be compatible
      !
      call check_surface_sanity(gs,es)
      !
      call read_gradient
      !
      call read_tdip
      !
      !  Determite perturbation vectors
      !
      call build_gbar
      call build_hbar
      call build_tdbar
      !
      !  Report coupling coefficients
      !
      call report_coupling_vectors
      !
      call screen_gbar_hbar
      !
      !  Read in composition of the initial wavepacket
      !
      call read_wavepacket
      !
      ! Initialize the pulse
      !
      call initialize_pulse
      !
      !  Prepare for amplitude caching
      !
      call initialize_cache
      !
      call TimerStop('Initialization')
      !
      call evolve_wavepacket
      !
      call cache_report
      !
      call TimerStop('Autocorrelation')
      call TimerReport
    end subroutine run_absorb_propagate

  end module absorb_propagate
  !
  program dofc
    use absorb_propagate
    !
    call run_absorb_propagate
  end program dofc
