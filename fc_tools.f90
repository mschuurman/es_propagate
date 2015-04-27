!
  module fc_tools
    use accuracy
    use elements
    use lapack
    use math
    use timer
    implicit none
    private
    public surfaceData
    public initialize_surface, check_surface_sanity, print_vector, print_matrix
    public state2ascii, ascii2state, state_energy
    !
    integer(ik)                  :: unit_scr         = 40  ! Any unused unit
    integer(ik)                  :: max_state_string = 256 ! Max length of the string describing a state
    !
    !  Collection of all information related to each harmonic surface
    !
    type surfaceData
      integer(ik)                    :: natoms        ! Number of atoms
      integer(ik)                    :: nvars         ! Number of degrees of freedom
      integer(ik)                    :: nvarsDummy    ! Number of dummy degrees of freedom. Dummy degrees
                                                      ! of freedom always have lower numbers.
      character(len=20), allocatable :: label(:)      ! Atom names
      real(rk), allocatable          :: xyz_ang(:,:)  ! Cartesian coordiates of each atom - Angstrom
      real(rk), allocatable          :: mass_amu(:)   ! Mass of each atom - AMU
      real(rk), allocatable          :: coord(:)      ! Cartesian coordinates - atomic units
      real(rk), allocatable          :: mass(:)       ! Mass for each coordinate - atomic units
      real(rk), allocatable          :: hess(:,:)     ! Hessian - atomic units
      real(rk), allocatable          :: modes(:,:)    ! Normal modes
      real(rk), allocatable          :: freq (:)      ! Frequencies of the normal modes, with soft modes
                                                      ! replaced by dummy vibrations
      real(rk), allocatable          :: freqTrue(:)   ! Frequencing of the normal nodes
      integer(sik), allocatable      :: maxLevel(:)   ! Max. automatic excitation level possible
    end type ! surfaceData
    !
    contains
    !
    !  Structure is expected in a mostly-standard XYZ format. Optionally, we'll
    !  accept isotopic mass (atomic units) after the coordinates.
    !
    subroutine read_xyz(surf,xyz_name)
      type(surfaceData), intent(out) :: surf
      character(len=*), intent(in)   :: xyz_name
      !
      character(len=256) :: line_buf
      character(len=80)  :: action
      integer(ik)        :: ios, line, iat
      !
      line = 0
      input_block: do 
        action = 'opening'
        open(unit_scr,file=trim(xyz_name),form='formatted',action='read',status='old',iostat=ios)
        if (ios/=0) exit input_block
        action = 'reading'
        !
        line = line + 1
        read(unit_scr,*,iostat=ios) surf%natoms
        if (ios/=0) exit input_block
        allocate (surf%label(surf%natoms),surf%xyz_ang(3,surf%natoms), &
                  surf%mass_amu(surf%natoms),surf%hess(3*surf%natoms,3*surf%natoms))
        !
        line = line + 1
        read(unit_scr,"()",iostat=ios)
        if (ios/=0) exit input_block
        !
        read_atoms: do iat=1,surf%natoms
          line = line + 1
          read(unit_scr,"(a)",iostat=ios) line_buf
          if (ios/=0) exit input_block
          !
          !  Try to parse it as element, X, Y, Z, and mass
          !
          read(line_buf,*,iostat=ios) surf%label(iat), surf%xyz_ang(:,iat), surf%mass_amu(iat)
          if (ios/=0) then
            !
            !  Did not work, try to pass as element, X, Y, and Z
            !
            read(line_buf,*,iostat=ios) surf%label(iat), surf%xyz_ang(:,iat)
            if (ios/=0 ) exit input_block
            surf%mass_amu(iat) = get_element_mass(trim(surf%label(iat)))
            if (surf%mass_amu(iat)<0) then
              write (out,"('Element ',a,' is not recognized at line ',i8,' of file ',a)") &
                     trim(surf%label(iat)), line, trim(xyz_name)
              stop 'fc_tools%read_xyz'
            end if
          end if
        end do read_atoms
        !
        action = 'closing'
        close(unit_scr,iostat=ios)
        return
      end do input_block
      write (out,"('Error ',i5,' while ',a,' file ',a,' at line ',i8)") ios, trim(action), trim(xyz_name), line
      stop 'fc_tools%read_xyz'
    end subroutine read_xyz
    !
    !  Read Hessian, in GAMESS $HESS format (sans $HESS and header lines)
    !
    subroutine read_hessian(surf,hess_name)
      type(surfaceData), intent(inout) :: surf
      character(len=*), intent(in)     :: hess_name
      !
      character(len=80)  :: action
      integer(ik)        :: ios, line
      integer(ik)        :: irow, icol
      integer(ik)        :: ichk, icont, ic
      real(rk)           :: values(5), err
      !
      line = 0
      input_block: do 
        action = 'opening'
        open(unit_scr,file=trim(hess_name),form='formatted',action='read',status='old',iostat=ios)
        if (ios/=0) exit input_block
        action = 'reading'
        !
        irow = 1 ; icol = 1
        read_lines: do while(irow<=3*surf%natoms)
          line = line + 1
          action = 'reading'
          read(unit_scr,"(i2,i3,5e15.8)",iostat=ios) ichk, icont, values
          if (ios/=0) exit input_block
          action = 'parsing'
          if (ichk/=mod(irow,100)) exit input_block
          stuff_values: do ic=1,5
            surf%hess(irow,icol) = values(ic)
            icol = icol + 1
            if (icol>3*surf%natoms) then
              irow = irow + 1
              icol = 1
              exit stuff_values
            end if
          end do stuff_values
        end do read_lines
        !
        action = 'closing'
        close(unit_scr,iostat=ios)
        !
        err = maxval(abs(surf%hess-transpose(surf%hess)))
        write (out,"('Hessian matrix in ',a,' deviates from symmetry by at most ',e12.5)") &
               trim(hess_name), err
        return
      end do input_block
      write (out,"('Error ',i5,' while ',a,' file ',a,' at line ',i8)") ios, trim(action), trim(hess_name), line
      stop 'fc_tools%read_hessian'
    end subroutine read_hessian
    !
    !  Perform harmonic frequency analysis
    !
    subroutine harmonic_frequency_analysis(surf,verbose,freq_cut,freq_dummy,evib_max,nvib_max)
      type(surfaceData), intent(inout) :: surf       ! Harmonic vibrational surface
      integer(ik), intent(in)          :: verbose    ! Verbosity level
      real(rk), intent(in)             :: freq_cut   ! Low cut-off frequency
      real(rk), intent(in)             :: freq_dummy ! Replacement for low (and imaginary) frequencies
      real(rk), intent(in)             :: evib_max   ! Max. intended vibrational energy, in cm^-1
      integer(sik), intent(in)         :: nvib_max   ! Max. intended number of quanta
      !
      integer(ik)       :: iat, i, j, imod, loc(1)
      character(len=1)  :: s_imag
      character(len=4)  :: s_skip
      character(len=50) :: buf
      real(rk)          :: norm
      !
      !  Convert units to something sensible and rearrange coordinates to avoid
      !  making distinction between Cartesian directions
      !
      surf%nvars = 3*surf%natoms
      allocate (surf%coord(surf%nvars),surf%mass(surf%nvars),surf%modes(surf%nvars,surf%nvars), &
                surf%freq(surf%nvars),surf%freqTrue(surf%nvars),surf%maxLevel(surf%nvars))
      stuff_atoms: do iat=1,surf%natoms
        surf%coord(3*iat-2:3*iat) = surf%xyz_ang(:,iat)/abohr
        surf%mass (3*iat-2:3*iat) = surf%mass_amu(iat)/electron_mass_in_amu
      end do stuff_atoms
      !
      !  Prepare mass-weighted Hessian and diagonalize it
      !
      do j=1,surf%nvars
        do i=1,surf%nvars
          surf%modes(i,j) = surf%hess(i,j)/sqrt(surf%mass(i)*surf%mass(j))
        end do
      end do
      call lapack_syev(surf%modes,surf%freq)
      !
      !  Fix phases of the normal modes such that the largest element is
      !  positive. Phase of the normal modes appears to be important!
      !
      fix_mode_phase: do imod=1,surf%nvars
        loc = maxloc(abs(surf%modes(:,imod)))
        if (surf%modes(loc(1),imod)<0._rk) then
          surf%modes(:,imod) = -surf%modes(:,imod)
        end if
      end do fix_mode_phase
      !
      !  Lapack routine returns squared normal-mode frequencies in ascending
      !  order. Count the number of "bad" vibrations, which will have to be
      !  replaced by dummies.
      !
      if (verbose>-1) then
        write (out,"(/'Frequencies of the harmonic normal vibrational modes in cm^-1:')")
      end if
      surf%nvarsDummy = 0
      screen_modes: do imod=1,surf%nvars
        if (surf%freq(imod)<=0) then
          s_imag          = 'I'
          surf%freq(imod) = -sqrt(-surf%freq(imod))
        else
          s_imag          = ' '
          surf%freq(imod) =  sqrt( surf%freq(imod))
        end if
        surf%freqTrue(imod) = surf%freq(imod)
        s_skip = ' '
        if (surf%freq(imod)<freq_cut/h2cm) then
          surf%nvarsDummy     = surf%nvarsDummy + 1
          surf%freq(imod)     = freq_dummy/h2cm
          surf%maxLevel(imod) = 0
          s_skip              = 'skip'
        else
          if (surf%freq(imod)>(evib_max/nvib_max)/h2cm) then
            surf%maxLevel(imod) = (evib_max/h2cm) / surf%freq(imod)
          else
            surf%maxLevel(imod) = nvib_max
          end if
        end if
        !
        if (verbose>-1) then
          write (out,"(i4,1x,f12.4,1x,a1,1x,a4,2x,' <= (x ',i4,')')") imod,h2cm*abs(surf%freqTrue(imod)), &
                 s_imag, s_skip, surf%maxLevel(imod)
        end if
      end do screen_modes
      !
      if (surf%nvarsDummy>0 .and. verbose>-1) then
        write (out,"(/'Modes marked as ""skip"" have been replaced with dummy vibrations at ',f6.1,' cm^-1')") freq_dummy
        write (out,"( 'These modes will not be vibrationally excited in the Franck-Condon factor calculations')")
      end if
      if (verbose>-1) write (out,"()")
      !
      !  Report normal modes
      !
      if (verbose>=1) then
        print_modes: do imod=1,surf%nvars
          write (buf,"('Mode ',i0,' (',f12.4,'/',f12.4,')')") imod, h2cm*surf%freq(imod), h2cm*surf%freqTrue(imod)
          call print_vector(surf%modes(:,imod),trim(buf))
        end do print_modes
      end if
      !
      !  Do a bit of sanity testing - the eigenvectors must be orthonormal
      !
      sanity: do i=1,surf%nvars
        do j=1,i
          norm = dot_product(surf%modes(:,i),surf%modes(:,j))
          if ( (i/=j .and. abs(norm)>10*spacing(1._rk)) .or. &
               (i==j .and. abs(norm-1.0_rk)>100*spacing(1._rk)) ) then
            write (out,"('Normal modes ',i5,' and ',i5,' are not orthonormal enough: <i|j> = ',g25.18)") i, j, norm
            stop 'fc_tools%harmonic_frequency_analysis - bad modes'
          end if
        end do
      end do sanity
    end subroutine harmonic_frequency_analysis
    !
    subroutine initialize_surface(surf,xyz_name,hess_name,verbose,freq_cut,freq_dummy,evib_max,nvib_max)
      type(surfaceData), intent(out) :: surf        ! Harmonic vibrational surface
      character(len=*), intent(in)   :: xyz_name    ! Name of the file containing XYZ coordinates
      character(len=*), intent(in)   :: hess_name   ! Name of the file containing the Hessian
      integer(ik), intent(in)        :: verbose     ! Verbosity level
      real(rk), intent(in)           :: freq_cut    ! Low cut-off frequency
      real(rk), intent(in)           :: freq_dummy  ! Replacement for low (and imaginary) frequencies
      real(rk), intent(in)           :: evib_max    ! Max. intended vibrational energy, in cm^-1
      integer(sik), intent(in)       :: nvib_max    ! Max. intended number of quanta
      !
      call TimerStart('Initialize surface')
      call read_xyz(surf,xyz_name)
      call read_hessian(surf,hess_name)
      call harmonic_frequency_analysis(surf,verbose,freq_cut,freq_dummy,evib_max,nvib_max)
      call TimerStop('Initialize surface')
    end subroutine initialize_surface
    !
    !  Make sure our harmonic surfaces are compatible.
    !
    subroutine check_surface_sanity(s1,s2)
      type(surfaceData), intent(in) :: s1, s2  ! Two harmonic surfaces involved in the calculation
      integer(ik)                   :: ideg
      !
      call TimerStart('Check sanity')
      if (s1%natoms /= s2%natoms) then
        write (out,"('Number of atoms does not match: ',i4,' vs ',i4)") s1%natoms, s2%natoms
        stop 'fc_tools%check_surface_sanity - number of atoms'
      end if
      if (s1%nvars /= s2%nvars) then
        write (out,"('Number of vibrations does not match: ',i4,' vs ',i4)") s1%nvars, s2%nvars
        stop 'fc_tools%check_surface_sanity - number of vibrations'
      end if
      if (s1%nvarsDummy /= s2%nvarsDummy) then
        write (out,"(/'Number of dummy vibrations is different on S1 and S2: ',i5,' vs. ',i5)") &
               s1%nvarsDummy, s2%nvarsDummy
        write (out,"( 'The absolute values of the Franck-Condon overlaps will not be meaningful')")
        write (out,"( 'The calculation will proceed nonetheless - handle with care.'/)")
      end if
      if (any(abs(s1%mass-s2%mass)>10*spacing(max(abs(s1%mass),abs(s2%mass))))) then
        write (out,"('Atomic masses do not match:')")
        write (out,"((5(1x,i5,1x,e8.3,2x)))") (ideg,s1%mass(ideg)-s2%mass(ideg),ideg=1,s1%nvars)
        stop 'fc_tools%check_surface_sanity - masses'
      end if
      !
      !  It makes sense to do the calculation if we get this far
      !
      call TimerStop('Check sanity')
    end subroutine check_surface_sanity
    !
    function state_energy(s,n,zpe) result(e)
      type(surfaceData), intent(in) :: s    ! Harmonic surface
      integer(sik), intent(in)      :: n(:) ! State vector
      logical, intent(in)           :: zpe  ! Include the zero-point energy?
      real(rk)                      :: e    ! Energy of the state
      !
      e = sum(n*s%freqTrue(1:s%nvars))
      if (zpe) then
        e = e + 0.5_rk*sum(s%freqTrue(1:s%nvars))
      end if
    end function state_energy
    !
    subroutine print_matrix(mat,name)
      real(rk), intent(in)         :: mat(:,:)
      character(len=*), intent(in) :: name
      !
      integer(ik), parameter :: cb = 5  ! Batch size
      integer(ik)            :: nr, nc  ! Number of rows and columns
      integer(ik)            :: r,  c   ! Row and column
      integer(ik)            :: c1, c2  ! Current batch of columns
      !
      write (out,"(/t5,a/)") trim(name)
      !
      nr = size(mat,dim=1)
      nc = size(mat,dim=2)
      col_batch: do c1=1,nc,cb
        c2 = min(nc,c1+cb-1)
        write (out,"()")
        write (out,"(1x,5x,16(1x,i16))") (c,c=c1,c2)
        !
        rows: do r=1,nr
          write (out,"(1x,i5,16(1x,e16.10))") r, (mat(r,c),c=c1,c2)
        end do rows
        !
        write (out,"()")
      end do col_batch
    end subroutine print_matrix
    !
    subroutine print_vector(vec,name)
      real(rk), intent(in)         :: vec(:)
      character(len=*), intent(in) :: name
      !
      integer(ik), parameter :: cb = 5  ! Batch size
      integer(ik)            :: nc      ! Number of columns
      integer(ik)            :: c       ! Row and column
      integer(ik)            :: c1, c2  ! Current batch of columns
      !
      write (out,"(/t5,a/)") trim(name)
      !
      nc = size(vec)
      write (out,"(1x,5x,sp,16(14x,i3))") (c,c=1,min(cb,nc))
      batch: do c1=1,nc,cb
        c2 = min(nc,c1+cb-1)
        write (out,"(1x,i5.5,16(1x,e16.10))") c1-1, (vec(c),c=c1,c2)
      end do batch
      write (out,"()")
    end subroutine print_vector
    !
    function state2ascii(n) result(s)
      integer(sik), intent(in)        :: n(:) ! Vibrational state vector
      character(len=max_state_string) :: s    ! "Compact" text desribing the state
      !
      integer(ik)       :: is, s_next, len_buf
      character(len=20) :: buf
      !
      s = ' '
      s_next = 1
      scan_states: do is=1,size(n)
        if (n(is)==0) cycle scan_states
        ! write statements below use Fortran-95 formatting feature (i0)
        if (n(is)==1) then
          write (buf,"(',',i0)") is
        else
          write (buf,"(',',i0,'^',i0)") is, n(is)
        end if
        len_buf = len_trim(buf)
        if (s_next+len_buf>max_state_string) then
          write (out,"('Compile-time limit on the lenhth of a vibrational state name exceeded.')")
          write (out,"('Increase parameter ""max_state_string"" in franck_condon.f90 and recompile.')")
          stop 'franck_condon%state2ascii - static limit exceeded'
        end if
        s(s_next:s_next+len_buf-1) = buf(:len_buf)
        s_next = s_next + len_buf
      end do scan_states
      !
      if (s_next==1) then
        !
        !  Special case - ground state
        !
        s = "0"
      else
        !
        !  Delete the leading comma - it's not needed
        !
        s = s(2:) 
      end if
    end function state2ascii
    !
    function ascii2state(s,ns) result(n)
      character(len=*),intent(in) :: s     ! "Compact" text desribing the state
      integer(ik), intent(in)     :: ns    ! Number of vibrational modes
      integer(sik)                :: n(ns) ! Vibrational state vector
      !
      integer(ik)       :: s_next, s_last, caret, info, is, io
      character(len=50) :: fmt_buf
      !
      n = 0
      !
      !  Skip leading blanks
      !
      s_next = index(trim(s),' ',back=.true.)+1
      !
      if (s(s_next:)=='0') return ! Special case: ground state
      !
      scan_string: do
        if (s(s_next:)==' ') return
        s_last = index(s(s_next:),',')
        if (s_last<=0) then
          s_last = len_trim(s(s_next:))+1
        endif
        !
        caret = index(s(s_next:s_next+s_last-2),'^')
        if (caret/=0) then
          read(s(s_next+caret:s_next+s_last-2),"(i10)",iostat=info) io
          if (info/=0) then
            write (out,"('Error parsing excitation level specification:')") 
            exit scan_string
          end if
        else
          io = 1
          caret = s_last
        end if
        read(s(s_next:s_next+caret-2),"(i10)",iostat=info) is
        if (info/=0) then
          write (out,"('Error parsing mode specification:')") 
          exit scan_string
        end if
        if (is<1 .or. is>ns) then
          write (out,"('State number ',i5,' is outside of the valid range')") is
          exit scan_string
        end if
        if (io<0 .or. io>huge(1_sik)) then
          write (out,"('State occupation ',i5,' is outside of the valid range')") is
          exit scan_string
        end if
        n(is) = io
        !
        s_next = s_next + s_last 
      end do scan_string
      !
      write (out,"(t1,'""',a,'""')") trim(s)
      write (fmt_buf,"('(T',i0,',''>'',A,T',i0,',''<'')')") s_next, s_next+s_last
      write (out,fmt_buf) s(s_next:s_next+s_last-2)
      stop 'franck_condon%ascii2state'
    end function ascii2state
    !
  end module fc_tools
