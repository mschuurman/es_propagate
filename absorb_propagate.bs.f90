!
!  Simulation of short-time wavepacket evolution on  harmonic surfaces coupled
!  via an external field.
!
!  This is a very naive implementation - numerical efficiency was not the goal
!  (MSS: and yet it has mutated to take on new life)
!
!  Although the code is OpenMP-parallel, the scalability is very modest,
!  Do not expect improvements beyond 2 CPUs.
!
!  Version history:
!
!   January 22, 2009 - Initial version
!   January 30, 2009 - OpenMP parallel version
!   February 18, 2013 - Repurposed autocorrelation.f90 code to propagate wp in presence
!                       of an external field
!
  module absorb_propagate 
    !$ use OMP_LIB
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
    integer(ik)                  :: unit_scr         = 40   ! Any unused unit
    integer(ik)                  :: max_state_string = 256  ! Max length of the string describing a state
    integer(ik)                  :: ntstep           = 14   ! Number of time steps to store (temporary)
    complex(rk)                  :: zeroc            = (0._rk,0._rk)
    complex(rk)                  :: im               = (0._rk,1._rk)
    integer(ik)                  :: max_lvl          = 8    ! Maximum level of recurrence for Bulirsh-Stoer algorithm
    integer(ik)                  :: n_lvl(9)         = (/2,4,6,8,10,12,14,16,18/) ! Number of modified mid-point steps for each level of B-S
    real(rk)                     :: eps_conv         = 1e-5 ! Level of convergence for amplitudes. Setting this to an "appropriate" 
                                                            ! value requires a little more testing

    !
    !  User-defined parameters, controlled through AP_PAR namelist
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
    real(rk)           :: prop_cutoff      = 1e-12_rk ! do not proprage coefficients less than this value
    real(rk)           :: print_cutoff     = 1e-3_rk  ! cutoff for printing coefficients 
    integer(ik)        :: max_nvib(100)    = 0        ! Maximum number of vibrational functions in each quanta                
    integer(ik)        :: nvib_cutoff      = -1       ! Amplitudes farther away than nvib_cutoff quanta from
                                                      ! the target set will be neglected. -1 means no cutoff
    real(rk)           :: evib_cutoff      = -1       ! States of E > evib_cutoff will be neglected, where
    real(rk)           :: print_tstep      =-1._rk    ! Print the wavepacket coefficients every print_tstep
    real(rk)           :: delta_e          = 0._rk    ! Electronic energy difference between two states
    real(rk)           :: tstep            = 5._rk    ! time step for propagation (a.u.)
    real(rk)           :: tfinal           = 400._rk  ! time to propagate (a.u.)
    logical            :: enforce_nyquist  = .true.   ! set the maximum time step to be no larger than the nyquist period
    logical            :: use_rwa          = .false.  ! use the rotating wave approximation for the pump pulse
    real(rk)           :: pulse_fwhm       = 1._rk    ! fwhm of envelope
    real(rk)           :: pulse_freq       = 1._rk    ! pulse frequency in au
    real(rk)           :: pulse_amp        = 0.005_rk ! amplitude of pulse envelope
    real(rk)           :: pulse_phase      = 0._rk    ! pulse frequency in au
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
    real(rk)                  :: nyquist_period       ! maximum time step we can take and still sample field without aliasing
    integer(ik)               :: prop_coef(2)         ! number of coefficient to be propragated on each state (> constant_cutoff)
    integer(ik)               :: print_coef(2)        ! Number of coefficients to dump to output (> print_cutoff)
    integer(sik), allocatable :: n_a0(:,:)            ! List of functions in the initial wavepacket
    integer(sik), allocatable :: n_max(:)             ! maximum number of functions in each mode
    integer(ik), allocatable  :: s_a0(:)              ! state on which to put basis function
    complex(rk), allocatable  :: a0(:)                ! Function weights in the initial wavepacket
    type(grCacheWrapper)      :: cache                ! Cache of the expansion coefficients for initial state
    type(grCacheWalker)       :: walker               ! Structure to walk through cache of "a" coefficients
    real(ark)                 :: evaluations = 0      ! Number of evaluations for Taylor coefficients
    real(ark)                 :: neglected_amps = 0   ! Number of amplitudes neglected because of nvib_cutoff
    real(ark)                 :: extrap_t(8)          ! value of t increments for Richardson extrapolation
    type(EFpulse)             :: FPulse               ! pulse information 

    !
    namelist /AP_PAR/                                                               &
       verbose,                                                                     &
       freq_cut, freq_dummy, gh_cutoff, prop_cutoff, print_cutoff,                  &
       xyz, gs_hessian, es_hessian, es_gradient, tdip_gradient,                     &
       packet_length, pulse_dir, pulse_fwhm, pulse_freq, pulse_amp, pulse_phase,    &
       tdip0, tstep, tfinal, enforce_nyquist, use_rwa,                              &
       delta_e, max_nvib,evib_cutoff, nvib_cutoff, print_tstep
    !
    contains
    !
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
    ! zero all gbar/hbar/tdbar contributions less than eps_cut
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
        write (out,"('          Largest absolute gradient coupling: ',g13.6)") g_max
        write (out,"('           Largest absolute hessian coupling: ',g13.6)") h_max
        write (out,"('  Largest absolute transition dipole element: ',g13.6)") td_max
        write (out,"('                  Coupling neglect threshold: ',g13.6)") eps_cut
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
    ! Print predominant coupling contributions from g_bar and h_bar
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
        write (out,"(1x,i4,2x,f8.2,2x,g11.4,1x)",advance='no') im, freq, grad
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
         write (out,"(1x,i4,1x,g11.4,1x)",advance='no') ml, hess
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
    ! find stationary point and corresponding harmonic frequencies on perturbed
    ! potential energy surface
    !
    subroutine analyze_surface
     implicit none
     real,allocatable     :: H(:,:),A(:,:)
     real,allocatable     :: b(:)
     real,allocatable     :: x(:)
     integer(ik)          :: i,j     

     ! allocate
     allocate(A(gs%nvars,gs%nvars),H(gs%nvars,gs%nvars))
     allocate(b(gs%nvars),x(gs%nvars))

     ! construct coefficient matrix for potential derivative
     H = 0.
     b = 0.
     x = 0.
     loop_rows: do i=gs%nvarsDummy+1,gs%nvars 
      H(i,i) = gs%freq(i)*(gs%freq(i) + 4._rk*h_bar(i,i))
      loop_cols: do j=1,i-1
       H(i,j) = h_bar(i,j) * sqrt( gs%freq(i)*gs%freq(j) ) 
       H(j,i) = H(i,j)
      end do loop_cols
     end do loop_rows

     ! construct b vector
     loop_crds: do i=gs%nvarsDummy+1,gs%nvars
      b(i) = -2._rk * g_bar(i) * sqrt(gs%freq(i))
     end do loop_crds

     ! solve for minimum
     A = H
     call lapack_ginverse_real2(A,1.0) 
     x = matmul(A,b) 

     ! diagonalize to get harmonic frequencies
     A = H
     call lapack_ssyev(A,b)
     ! set imaginary frequencies negative
     where (b > 0._rk)
      b = sqrt(b)
     elsewhere
      b = -sqrt(abs(b))
     end where

     ! determine excited state energy relaxation
     write(out,'(" Stationary Point -------------------")')
     write(out,'(8(f10.4))')x
     write(out,'("")')
     write(out,'(" Energy Relaxation ------------------")')
     write(out,'(f15.2)')V(x)*h2cm
     write(out,'("")')
     write(out,'(" Frequencies ------------------------")')
     write(out,'(8(f10.2))')b*h2cm

     deallocate(H,A,b,x)
     return

     contains 

      ! Evaluate potential at geometry q
      function V(q) result(U)
       implicit none
       real,intent(in)  :: q(gs%nvars)
       real(rk)         :: U    
       integer(ik)      :: i,j

       U = 0._rk
       loop_diag: do i = 1,gs%nvars
        U = U + 0.5_rk * gs%freq(i) * (gs%freq(i) + 4._rk*h_bar(i,i)) * q(i)**2
        U = U + 2._rk * g_bar(i) * sqrt(gs%freq(i)) * q(i)
        loop_offdiag: do j = 1,i-1
         U = U + 4._rk * h_bar(i,j) * sqrt(gs%freq(i)*gs%freq(i)) * q(i) * q(j)
        end do loop_offdiag
       end do loop_diag

       return
      end function V 

      ! Evaluate the gradient of the potential at geometry q
      function gradV(q) result(grad)
       implicit none
       real,intent(in)        :: q(gs%nvars)
       real(rk)               :: grad(gs%nvars)
       integer(ik)               :: i,j

       do i = 1,gs%nvars
        grad(i) = grad(i) + (gs%freq(i)**2 + 4._rk*gs%freq(i)*h_bar(i,i)) * q(i)
        grad(i) = grad(i) + 2._rk * g_bar(i) * sqrt(gs%freq(i))
        loop_crds: do j=gs%nvarsDummy+1,gs%nvars
          if(j/=i) then
            grad(i) = grad(i) + 4._rk * h_bar(i,j) * sqrt(gs%freq(i)*gs%freq(j)) * q(j)
          endif
        end do loop_crds
       enddo

       return
      end function gradV

    end subroutine analyze_surface

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
        write (out,"('C= ',g19.12,' mu= ',a30,' nu= ',a30)") c, trim(state2ascii(mu)), trim(state2ascii(nu))
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
      integer(ik) :: i
      integer(ik) :: max_i
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
        write (out,"('D= ',g19.12,' mu= ',a30,' nu= ',a30)") d, trim(state2ascii(mu)), trim(state2ascii(nu))
      end if
      return

    end function d_mu_nu    

    !
    !  compute contribution of state nu to coefficient for vibrational state mu
    !  on the ground electronic state.
    !  process_a_nu() would be a natural candidate for a contained function
    !  of accumulate_at(). Unfortunately, ifort seems to miscompile OpenMP
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
!      if (abs(c)<coupling_cutoff) return ! No coupling - ignore mode
      de = delta_e + state_energy(gs,mu,zpe=.false.) - state_energy(gs,nu,zpe=.false.)
      val = cf * c * exp( -im * de * t)
      return

    end function process_a_nu

    !
    !  compute contribution of state nu to coefficient for vibrational state mu
    !  on the upper electronic state.
    !  process_b_nu() would be a natural candidate for a contained function
    !  of accumulate_bt(). Unfortunately, ifort seems to miscompile OpenMP
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
!        if (abs(c)<coupling_cutoff) return ! No coupling - ignore mode
        de = state_energy(gs,nu,zpe=.false.) - state_energy(gs,mu,zpe=.false.)

       case(1)       ! this is a dipole matrix element term
        c = d_mu_nu(mu,nu,F)       
!        if (abs(c)<coupling_cutoff) return ! No coupling - ignore mode
        de = delta_e + state_energy(gs,nu,zpe=.false.) - state_energy(gs,mu,zpe=.false.)

      end select

      val = cf * c * exp( im * de * t)
!      write(out,"('nu=',4(i2),' mu=',4(i2),' c=',f12.8,' cf=',f12.8,' val=',f12.8,' t=',f6.4)")nu(6:9),mu(6:9),c,abs(cf),abs(val),t

      return

    end function process_b_nu
 
    !
    !  Taylor coefficients for time evolution of the nuclear wavepacket on lower state
    !
    subroutine accumulate_at(mu,cf,F,t,dxdt)
      integer(sik), intent(in) :: mu(:)  ! State of interest
      complex(rk),intent(in)   :: cf     ! Expansion coefficient
      real(rk),intent(in)      :: F(3),t
      integer(ik),intent(in)   :: dxdt   ! slot to accumulate
      !
      integer(sik) :: nu(gs%nvars)       ! States at lower expansion orders
      integer(ik)  :: c1                ! Positions of the vibrational quanta
      integer(sik) :: n1, n1_min, n1_max ! Number of quanta at c1 position
      complex(rk)  :: at
      !
      !  Begin by looking at the self-coupling
      !
      !$omp parallel default(none) &
      !$omp& private(nu,at,c1,n1_min,n1_max,n1) &
      !$omp& shared(mu,gs,n_max,cf,F,t,dxdt)
      !
      nu = mu
      !$omp single      
      at = process_a_nu(mu,nu,cf,F,t)
      call update_cf(nu,1_ik,dxdt,at,ovrwrite=.false.)
      !$omp end single nowait
      !
      !$omp do schedule(dynamic)
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
          if(abs(at) > 0_rk)call update_cf(nu,1_ik,dxdt,at,ovrwrite=.false.) 
        end do coupling_n1_only
        nu(c1) = mu(c1)  ! Restore state vector
      end do quantum_c1
      !$omp end do
      !$omp end parallel
      return
    end subroutine accumulate_at

    !
    !  Taylor coefficients for time evolution of the nuclear wavepacket on upper state
    !
    subroutine accumulate_bt(mu,cf,typ,F,t,dxdt)
      integer(sik), intent(in) :: mu(:)  ! State of interest
      complex(rk),intent(in)   :: cf     ! Expansion coefficient
      integer(ik), intent(in)  :: typ     ! term type
      real(rk),intent(in)      :: F(3),t ! external field, time
      integer(ik),intent(in)   :: dxdt   ! slot to accumulate
      !
      integer(sik) :: nu(gs%nvars)       ! States at lower expansion orders
      integer(ik)  :: c1, c2             ! Positions of the vibrational quanta
      integer(sik) :: n1, n1_min, n1_max ! Number of quanta at c1 position
      integer(sik) :: n2, n2_min, n2_max ! ... and c2 position
      complex(rk)  :: bt

      !
      !  Begin by looking at the self-coupling
      !$omp parallel default(none) &
      !$omp& private(nu,bt,c1,c2,n1_min,n1_max,n1,n2_min,n2_max,n2) &
      !$omp& shared(mu,gs,n_max,cf,typ,F,t,dxdt)

      nu = mu
      !$omp single
      bt = process_b_nu(mu,nu,cf,typ,F,t)
      call update_cf(nu,2_ik,dxdt,bt,ovrwrite=.false.)
      !$omp end single nowait

      !
      !$omp do schedule(dynamic)
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
          call update_cf(nu,2_ik,dxdt,bt,ovrwrite=.false.)
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
              call update_cf(nu,2_ik,dxdt,bt,ovrwrite=.false.)
            end do coupling_n2
            nu(c2) = mu(c2)
          end do quantum_c2
        end do coupling_n1
        nu(c1) = mu(c1)  ! Restore state vector
      end do quantum_c1
      !$omp end do
      !$omp end parallel

      return
    end subroutine accumulate_bt

    !
    !  Look up expansion coefficient in the cache
    !
    function retrieve_cf(mu,st,slot) result(ab)
      integer(sik), intent(in) :: mu(:)  ! State of interest
      integer(ik), intent(in)  :: st     ! state 1=gs, 2=es
      integer(ik), intent(in)  :: slot   ! index in "values" array 
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
       if(.not.associated(cgr))stop 'tried to access non-existent coefficient...'
       ab = cgr%values(slot,st)
       return

      endif

    end function retrieve_cf

    !
    !  Update expansion coefficient in the cache
    !
    subroutine update_cf(mu,st,slot,ab,ovrwrite)
      integer(sik), intent(in) :: mu(:)  ! State of interest
      integer(ik), intent(in)  :: st     ! state 1=gs, 2=es
      integer(ik), intent(in)  :: slot   ! index in "values" array 
      complex(rk), intent(in)  :: ab     ! contribution
      logical, intent(in)      :: ovrwrite    ! 0=set to value ab, 1=add value ab
      type(grCache), pointer   :: cgr

      if(neglect_state(mu))return
        
       cgr => cache_gr_locate(cache,mu)

       ! write (out,"('Request ',a30,' ord ',i3,' found= ',l1)") trim(state2ascii(mu)), s, associated(cgr)
       if(.not.associated(cgr))stop 'tried to update non-existent cf...'
       if(ovrwrite) then
         cgr%values(slot,st) = ab            
        else
         cgr%values(slot,st) = cgr%values(slot,st) + ab
       endif

       return
    end subroutine update_cf

    !
    !  Take linear combinations of elements in values array (given by pos), using
    !   coefficients in coef and place result at index "targ"
    !   Currently this is done for both state slots.
    !
    subroutine combine_cf(mu,targ,n_sum,pos,coef)
      integer(sik), intent(in) :: mu(:)  ! State of interest
      integer(ik), intent(in)  :: targ         ! target slot to put sum
      integer(ik), intent(in)  :: n_sum  ! number of terms in sum
      integer(ik), intent(in)  :: pos(n_sum)   ! position of terms in sum
      complex(rk), intent(in)  :: coef(n_sum)  ! coefficient of terms in sum
      ! 
      integer(ik)              :: i,s
      complex(rk)              :: cf
      type(grCache), pointer   :: cgr

      if(neglect_state(mu))return

       cgr => cache_gr_locate(cache,mu)

       ! write (out,"('Request ',a30,' ord ',i3,' found= ',l1)") trim(state2ascii(mu)), s, associated(cgr)
       if(.not.associated(cgr)) stop 'tried to access non-existent cf..'
       !
       loop_states: do s = 1,2
        cf = zeroc
        loop_sum: do i = 1,n_sum
         cf = cf + coef(i)*cgr%values(pos(i),s)
        end do loop_sum
        cgr%values(targ,s) = cf
       end do loop_states

       return
    end subroutine combine_cf

    !
    ! Return true if passed state is outside one of a number of tolerances
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
      write (out,"(/'Norm of the initial wavepacket is ',g21.14/)") sum(abs(a0)**2)
    end subroutine read_wavepacket
    !
    ! Not much happens here, need to set gaussian envelope width mostly
    !
    subroutine initialize_pulse
     implicit none

     call initializePulse(FPulse,pulse_fwhm,pulse_freq,pulse_amp,pulse_phase)
     if(tfinal < 0) tfinal = 2._rk * FPulse%toff2
     nyquist_period = 1. / (Fpulse%omega)
     if(use_rwa) enforce_nyquist = .false.

     return
     
    end subroutine initialize_pulse
    !
    ! Initialize coefficients at t = 0
    !
    subroutine initialize_cache
      implicit none
      integer(ik)       :: is
      !
      call cache_gr_initialize(cache,ntstep,n_max)

      prop_coef = 0_ik
      print_coef = 0_ik
      do is = 1,packet_length
        call update_cf(n_a0(:,is),s_a0(is),0_ik,a0(is) / sqrt(dot_product(a0,a0)),ovrwrite=.true.)
        prop_coef(s_a0(is)) = prop_coef(s_a0(is)) + 1_ik
        print_coef(s_a0(is)) = print_coef(s_a0(is)) + 1_ik
      enddo

      do is = 1,2
       call printCoefSort(is,0._rk,.true.,.false.)
       write (out,"(' ')")
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
    ! add in contribution to all connected del_a/del_t and del_b/del_t terms
    !
    subroutine update_tderiv(t,x,dxdt)
      implicit none
      real(rk),intent(in)      :: t
      integer(ik),intent(in)   :: x,dxdt !x=slot to pull coefs, dxdt=slot to put derivs

      integer(sik)             :: n(size(cache%n_max))
      complex(rk)              :: cf
      real(rk)                 :: F(3)

      F = evaluateExternalField( FPulse,t,use_rwa ) * pulse_dir 
      !F = evaluateSinPulse(FPulse,t)* pulse_dir
      ! 
      ! a(t) will only contribute to propragation of b(t)
      !
      call walker_gr_initialize(cache,walker)
      scan_tree: do while(walker_gr_next(walker,n))

          ! determine contribution from a coefficients
          cf = retrieve_cf(n,1_ik,x)
          if(abs(cf) > prop_cutoff ) then
           call accumulate_bt(n,cf,1_ik,F,t,dxdt)
          endif

          ! determine contribution from b coefficients
          cf = retrieve_cf(n,2_ik,x)
          if(abs(cf) > prop_cutoff ) then
           call accumulate_at(n,cf,F,t,dxdt)
           call accumulate_bt(n,cf,0_ik,F,t,dxdt)
          endif

      end do scan_tree
      call walker_gr_destroy(walker)

    end subroutine update_tderiv

    !
    ! take the coefficients determined from integration procedure, put them
    !  into the "current time" slot (i.e. slot=0)
    !
    subroutine update_coefficients(nrm)
     implicit none
     real(rk),intent(inout)   :: nrm(2)
     integer(ik)            :: s
     integer(sik)           :: n(size(cache%n_max))
     complex(rk)            :: cf

     nrm = 0._rk 
     print_coef = 0
     prop_coef = 0
     call walker_gr_initialize(cache,walker)
     scan_tree: do while(walker_gr_next(walker,n))
       scan_states: do s = 1,2
        cf = retrieve_cf(n,s,2_ik)
        if(abs(cf) > prop_cutoff)then
         call update_cf(n,s,0_ik,cf,ovrwrite=.true.)
         prop_coef(s) = prop_coef(s) + 1
         if(abs(cf) > print_cutoff)print_coef(s) = print_coef(s) + 1
        else
         call update_cf(n,s,0_ik,zeroc,ovrwrite=.true.)
        endif
        nrm(s) = nrm(s) + abs(cf)**2
       end do scan_states
     end do scan_tree
     call walker_gr_destroy(walker)

     do s = 1,2
      nrm(s) = sqrt(nrm(s))
     end do

    end subroutine update_coefficients

    !
    !  Calculate Taylor expansion for the evolution of the "interesting" part of the 
    !  wavepacket, which can overlap with the wavepacket on the unperturbed surface.
    !
    subroutine evolve_wavepacket
      integer(ik)             :: s,t_cnt
      integer(ik)             :: row_max,row_active
      real(rk)                :: F,nrm(2)
      real(rk)                :: t,dt_attempt,dt_conv,dt_next
      real(rk)                :: a(max_lvl+1),alpha(max_lvl,max_lvl)
      logical                 :: first
      !
      call TimerStart('Wavepacket evolution')

      t = 0._rk
      t_cnt = 0_ik
      dt_attempt = tstep
      call bs_init(eps_conv,row_max,row_active,a,alpha)

      first = .true.
      time_step: do

        ! 
        ! Print coefficients if at print_tstep
        !        
        if(t >= t_cnt*print_tstep .and. print_tstep /= -1._rk) then
         loop_states: do s = 1,2
          call printCoefSort(s,t,.false.,.true.)
         end do loop_states
         t_cnt = t_cnt + 1
        endif

        ! 
        ! Set the slot in which we will accumulate right-hand side to zero
        !
        call initialize_timestep

        !
        ! Evalulate contribution of a(t) coefficients to updated left-hand side
        !  (note: will only contribute to propagation of b(t) coefficients)
        ! first two slots correspond to x,dx/dt of current time step
        call update_tderiv(t,0_ik,1_ik)
       
        ! 
        ! update a(t) and b(t)
        !
        call bs_step(first,t,dt_attempt,dt_conv,dt_next,a,alpha,row_max,row_active)

        !
        ! move updated a(t) and b(t) to appropriate location
        !
        call update_coefficients(nrm)

        ! 
        ! print out current status
        !
        F = evaluateExternalField( FPulse,t,use_rwa )
        write(out,"('t= ',f9.3,' dt=',f9.5,' Ft= ',f8.5,' |0|= ',f9.6,' |1|= ',f9.6,' |total|= ',f9.6,' N=',i10)") &
                   t,dt_attempt,F,nrm(1),nrm(2),sqrt(sum(nrm**2)),sum(prop_coef)

        ! 
        ! update time, and timestep
        dt_attempt = min(dt_next,tfinal-t)
        if(dt_attempt <= 0._rk)exit
 
        t = t + dt_conv
        first = .false.
       
      end do time_step

      call TimerStop('Wavepacket evolution')
      !
      if (verbose>=1 ) then
        print_states: do s = 1,2
         call printCoefSort(s,t,.true.,.false.)
         write (out,"(' ')")   
        end do print_states 
      end if

      write (out,"(' ')")
      write (out,"(' Norm Psi[0]: ',f12.8)")nrm(1)
      write (out,"(' Norm Psi[1]: ',f12.8)")nrm(2)

    end subroutine evolve_wavepacket

    ! 
    !  print sorted coefficients
    ! 
    subroutine printCoefSort(s,t,print_summary,print_file)
     implicit none
     integer(ik),intent(in)   :: s
     real(rk),intent(in)      :: t
     logical,intent(in)       :: print_summary,print_file
     integer(ik)              :: i,j,nbat,ncol,nprint,ios
     complex(rk)              :: cf
     integer(sik)             :: n(size(cache%n_max))
     real(rk)                 :: nwgts(maxval(n_max),gs%nvars)
     real(rk)                 :: wtsum
     character(len=7)         :: fname(2) = (/'a_coef.','b_coef.'/)
     character(len=10)        :: t_str
     type cfst
      complex(rk)              :: ab = (0._rk,0._rk)
      integer(sik),allocatable :: m(:)
     end type cfst 
     type(cfst),allocatable   :: sarr(:)

      nwgts = 0._rk

      ! print summary output 
      if(print_summary) then
       write (out,"()")
       if(s == 1_ik) then
        write (out,"(a21,f9.2,a6)")'A coefficients at t= ',t,' -----'
        write (out,"(1x,a21,2x,a12,2x,a30,a11)")'a(tf)                  ','  |a(tf)|   ','harmonic basis function       ','energy (eV)'
       else
        write (out,"(a21,f9.2,a6)")'B coefficients at t= ',t,' -----'
        write (out,"(1x,a21,2x,a12,2x,a30,a11)")'b(tf)                  ','  |b(tf)|   ','harmonic basis function       ','energy (eV)'
       endif    
      endif

      allocate(sarr(print_coef(s)))
      do i = 1,print_coef(s)
       allocate(sarr(i)%m(size(cache%n_max)))
      end do

      ! print all coefficients greater than print_cutoff in descending order
      nprint = 0_ik
      call walker_gr_initialize(cache,walker)
      scan_coef: do while(walker_gr_next(walker,n))
        cf = retrieve_cf(n,s,0_ik)
        do i =gs%nvarsDummy+1,gs%nvars
         if(n(i)/=0)nwgts(n(i),i) = nwgts(n(i),i) + abs(cf)
        enddo
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

      if(nprint > 0) then
       if(print_file) then
        write(t_str,'(i10)')floor(t+0.5_rk)
        open (unit_scr,file=fname(s)//trim(adjustl(t_str)),form='formatted',action='write',status='replace',iostat=ios) 
        if (ios/=0) then
          write (out,"(' Error ',i6,' opening file ',a)") ios, trim(fname(s))
          stop 'absorb_propagate%printCoefSort - open'
        end if
       endif
       do i = 1,nprint
         if(print_summary) write (out,1000) real(sarr(i)%ab),aimag(sarr(i)%ab),abs(sarr(i)%ab),state2ascii(sarr(i)%m),state_energy(gs,sarr(i)%m,zpe=.false.)*h2ev
         if(print_file) write(unit_scr,1001)real(sarr(i)%ab),aimag(sarr(i)%ab),abs(sarr(i)%ab),state2ascii(sarr(i)%m),state_energy(gs,sarr(i)%m,zpe=.false.)*h2ev
       enddo
       close(unit_scr)
      endif

      ! print basis analsyis
      if(print_summary.and.t.gt.0._rk) then
       ncol = 10
       write(out,"()")
       write(out,"('  Basis set analysis ---------------------')")
       write(out,"(2x,a4,1x,a20,1x,a23)")'Mode','| Total contribution','| Weighted contribution'
       do i = gs%nvarsDummy+1,gs%nvars
        nbat = ceiling(1.*n_max(i)/ncol)
        wtsum = sum(nwgts(:,i))
        if(wtsum == 0._rk)wtsum = 1._rk
        write(out,"()")
        write(out,1002)i,sum(nwgts(:,i)),nwgts(1:min(ncol,n_max(i)),i)/wtsum
        do j = 2,nbat 
         write(out,1003)nwgts(ncol*(j-1)+1:min(j*ncol,n_max(i)),i)/wtsum
        enddo
       enddo
      endif

      do i = 1,print_coef(s)
       deallocate(sarr(i)%m)
      enddo
      deallocate(sarr)
 
1000  format(1x,'(',f9.6,',',f9.6,')',2x,f12.8,2x,a30,f12.4)
1001  format(1x,'(',f9.6,',',f9.6,')',2x,f12.8,2x,a30,f12.4)
1002  format(1x,i5,1x,'|',f19.8,1x,'|',1x,10(f7.4))
1003  format(29x,10(f7.4))
    end subroutine printCoefSort

    !------------------- Start Numerical Integration routines ------------------!

    !
    ! Initialize Bulirsch-Stoer propagation
    !
    subroutine bs_init(eps,kmax,kopt,a,alfa) 
     implicit none
     real(rk),intent(in)       :: eps
     integer(ik),intent(out)   :: kmax
     integer(ik),intent(inout) :: kopt
     real(rk),intent(out)      :: a(max_lvl+1)
     real(rk),intent(out)      :: alfa(max_lvl,max_lvl)
     !
     integer(ik)               :: k,q
     real(rk)                  :: eps_safe
     real(rk)                  :: safe = 0.25_rk

     eps_safe = eps * safe

     ! compute the Ak coefficients
     a(1) = n_lvl(1) + 1
     init_ak: do k = 1,max_lvl
      a(k+1) = a(k) + n_lvl(k+1)
     end do init_ak

     ! compute the alpha matrix
     loop_col: do q = 2,max_lvl
      loop_row: do k = 1,q-1
       alfa(k,q) = eps_safe ** ((a(k+1)-a(q+1)) / ((a(q+1)-a(1)+1._rk)*(2._rk*k+1._rk)))
      end do loop_row
     end do loop_col

     ! determine what level of convergence is required
     do kopt = 2,max_lvl-1     
      if(a(kopt+1) > a(kopt)*alfa(kopt-1,kopt))exit
     enddo
     kmax = kopt

     return
    end subroutine bs_init

    !
    ! Perform a Bulirsch-Stoer Step, modified from "Numerical Recipes".
    !
    subroutine bs_step(bs_reinit,t,dt_attempt,dt_actual,dt_next,a,alfa,kmax,kopt)
      implicit none
      logical,intent(in)          :: bs_reinit    ! .true. if we are initializing the order window
      real(rk),intent(in)         :: t            ! current time
      real(rk),intent(in)         :: dt_attempt   ! time step to attempt
      real(rk),intent(out)        :: dt_actual    ! time step converged
      real(rk),intent(out)        :: dt_next      ! size of next time step
      real(rk),intent(in)         :: a(max_lvl+1) ! ak array
      real(rk),intent(in)         :: alfa(max_lvl,max_lvl) ! alpha matrix
      integer(ik),intent(in)      :: kmax
      integer(ik),intent(inout)   :: kopt
      !
      integer(ik)                 :: k,l
      real(rk)                    :: t_current, t_estimate
      real(rk)                    :: max_error,err(max_lvl)
      real(rk)                    :: red,work,fact,work_min,scal
      logical                     :: converged
      logical                     :: reduce_dt    
      !
      real(rk)                    :: safe1 = 0.25_rk
      real(rk)                    :: safe2 = 0.70_rk
      real(rk)                    :: scale_max = 0.8_rk
      real(rk)                    :: red_min = 0.70
      real(rk)                    :: red_max = 0.00001_rk 

      dt_actual = dt_attempt
      converged = .false.
      reduce_dt = .false.

      try_stepsize: do while( .not.converged ) 

       loop_level: do k = 1,kmax
         t_current = t + dt_actual
         if(t_current == t)stop 'step size underflow in bs_step'
         ! Do modified midpoint integration from t to t+dt
         call modified_midpoint(t,dt_actual,n_lvl(k),k)
         ! esimate error 
         t_estimate = (dt_actual / n_lvl(k))**2
         ! extrapolate results
         call rational_extrap(k,t_estimate,max_error)
 
         if(k /= 1) then
           max_error = max_error / eps_conv
           err(k-1) = (max_error/safe1)**(1._rk / (2._rk*k-1._rk))

           if(k >= (kopt-1) .or. bs_reinit) then
            
             ! converged, exit loop
             if(max_error < 1._rk) then
               converged = .true.
               exit 
             endif 

             ! Not
             if (k == kmax .or. k == (kopt+1)) then   
               red = safe2/err(k-1)
               exit
             else if (k == kopt .and. alfa(kopt-1,kopt) < err(k-1)) then
               red = 1._rk / err(k-1) 
               exit
             else if (k == max_lvl .and. alfa(kopt,kopt-1) < err(k-1)) then
               red = alfa(k-1,kopt-1) * safe2 / err(k-1)
               exit
             else if (alfa(k-1,kopt) < err(k-1)) then
               red = alfa(k-1,kopt-1) / err(k-1) 
               exit
             endif

           endif ! end if (k>=(krow-1) || bs_reinit)
         endif ! if (k/=1)

       end do loop_level

       if( .not.converged ) then
        red = min(red,red_min)
        red = max(red,red_max)
        dt_actual = dt_actual * red
        reduce_dt = .true.
       endif 

      end do try_stepsize

      ! Determine convergence properties for new step size
      work_min = safe_max
      update_kopt: do l = 1,k-1
        fact = max(err(l),scale_max)
        work = fact * a(l+1)
        if(work < work_min) then
          scal     = fact
          work_min = work
          kopt     = l+1
        endif
      end do update_kopt

     ! determine size of next time step
     dt_next = dt_actual / scal

     ! Check for increase in step size 
     if(kopt >= k .and. kopt /= kmax .and. .not.reduce_dt) then
       fact = max(scal/alfa(kopt-1,kopt),scale_max)
       if(a(kopt+1)*fact <= work_min) then
         dt_next = dt_actual / fact
         kopt = kopt + 1
       endif
     endif

     ! make nyquist sampling of external field sets maximum step size 
     if(dt_next > nyquist_period .and. enforce_nyquist ) dt_next = nyquist_period
     return
    end subroutine bs_step

    !
    ! Perform a rational function extrapolation of successive modified
    ! midpoint integrations
    !
    subroutine rational_extrap(ilvl,t_estimate,error_max)
      implicit none
      integer(ik),intent(in)  :: ilvl
      real(rk),intent(in)     :: t_estimate
      real(rk),intent(out)    :: error_max
      !
      integer(ik)             :: k,s
      real(rk)                :: fx(max_lvl)
      complex(rk)             :: yy,c,v,b1,b,ddy
      !
      integer(sik)            :: n(size(cache%n_max))
      type(grCache), pointer  :: cgr

      extrap_t(ilvl) = t_estimate

      if(ilvl > 1) then
       error_max      = safe_max**(-1.) 

       do k = 1,ilvl-1
        fx(k+1) = extrap_t(ilvl-k) / t_estimate 
       end do

       call walker_gr_initialize(cache,walker)
       scan_tree: do while(walker_gr_next(walker,n))
        scan_states: do s = 1,2
          cgr => cache_gr_locate(cache,n)
          yy = cgr%values(2_ik,s)
          v  = cgr%values(5_ik,s)
          c  = yy
          cgr%values(5_ik,s) = yy
          scan_orders: do k = 2,ilvl
            b1 = fx(k) * v
            b  = b1 - c
            if(b /= zeroc) then
              b   = (c-v)/b
              ddy = c*b
              c   = b1*b
            else
              ddy = v
            endif
            if(k /= ilvl)v = cgr%values(4_ik+k,s)
            cgr%values(4_ik+k,s) = ddy
                              yy = yy + ddy
          end do scan_orders 
          cgr%values(2_ik,s) = yy
          error_max = max(error_max,abs(ddy/1._rk)) ! assume scale factor is 1. for now

        end do scan_states 
       end do scan_tree
       call walker_gr_destroy(walker)

      else
       error_max = 1._rk
      endif

      return
    end subroutine rational_extrap

    !
    ! Modified Mid-point integration
    !
    subroutine modified_midpoint(t0,tstep,nstep,k_lvl)
      implicit none
      real(rk),intent(in)     :: t0,tstep
      integer(ik),intent(in)  :: nstep
      integer(ik),intent(in)  :: k_lvl
      !
      integer(ik)             :: s,i
      integer(ik)             :: active_indices(3) = (/2_ik, 3_ik, 4_ik/)
      real(rk)                :: ti,h
      integer(sik)            :: n(size(cache%n_max))
      complex(rk)             :: re
       
      re = (1._rk,0._rk)
      h = tstep / nstep

      ! Take first step
      call walker_gr_initialize(cache,walker)
      scan_initial: do while(walker_gr_next(walker,n))
         call combine_cf(n, 2_ik, 1_ik, (/ 0 /), (/ re /) )
         call combine_cf(n, 3_ik, 2_ik, (/0, 1/), (/ re , -im*h /) )
         zero_tderiv: do s = 1,2
          call update_cf(n,s,4_ik,zeroc,ovrwrite=.true.)
         end do zero_tderiv
      end do scan_initial
      call walker_gr_destroy(walker)

      ti = t0 + h
      ! pull coefficients from 3, accumulate in 4
      call update_tderiv(ti,3_ik,4_ik)
      
      loop_steps: do i = 2,nstep

        call walker_gr_initialize(cache,walker)
        scan_steps: do while(walker_gr_next(walker,n))
          call iterate_midpoint(n, h, active_indices)
        end do scan_steps
        call walker_gr_destroy(walker)

        ti = ti + h 
        call update_tderiv(ti,3_ik,4_ik)
      
       end do loop_steps

       ! last step
        call walker_gr_initialize(cache,walker)
        scan_final: do while(walker_gr_next(walker,n))
          call combine_cf(n, 2_ik, 3_ik, active_indices , (/0.5_rk*re, 0.5_rk*re, -im * 0.5_rk*h /) )          
          if(k_lvl == 1) call combine_cf(n, 5_ik, 1_ik, (/ 2_ik /), (/1._rk*re/) )
          end do scan_final
        call walker_gr_destroy(walker)

    end subroutine modified_midpoint

    !
    !  Iterate the modified midpoint routine, moving the y(n+1) to the y(n). We assuming
    !  y'(n) has already been evaluated and is sitting at indices(3). 
    !
    subroutine iterate_midpoint(mu,step,indices)
      integer(sik), intent(in) :: mu(:)       ! State of interest
      real(rk),intent(in)      :: step        ! step size
      integer(ik), intent(in)  :: indices(3)  ! coefficient of terms in sum
      !
      integer(ik)              :: s
      complex(rk)              :: cf
      type(grCache), pointer   :: cgr

      if(neglect_state(mu))return

       cgr => cache_gr_locate(cache,mu)

       if (associated(cgr)) then
          loop_states: do s = 1,2
           cf = cgr%values(indices(1),s) - im * cgr%values(indices(3),s)* 2._rk * step
           cgr%values(indices(1),s) = cgr%values(indices(2),s)
           cgr%values(indices(2),s) = cf
           cgr%values(indices(3),s) = zeroc
          end do loop_states
       else
         write(out,"('should not be here...')")
         return
       endif

    end subroutine iterate_midpoint

    !--------------- End Numerical Integration routines ------------------!

    !
    !  Main problem driver - this is the only externally-visible routine
    !
    subroutine run_absorb_propagate
      integer(ik) :: info
      !
      call TimerStart('Propagation')
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
      !  Report coupling coefficients, and screen by magnitude
      !
      call report_coupling_vectors
      call screen_gbar_hbar
      call analyze_surface
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
      call TimerStop('Propagation')
      call TimerReport
    end subroutine run_absorb_propagate

  end module absorb_propagate
  !
  program dofc
    use absorb_propagate
    !
    call run_absorb_propagate
  end program dofc
