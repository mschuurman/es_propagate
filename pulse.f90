!
!  Electric field pulse shapes
!
 module pulse
   use accuracy
   implicit none
   private initializeAlpha
   public initializePulse,evaluateExternalField,EFPulse,evaluateSinPulse

!
!  Data pulse
!
   type EFPulse 
     real(rk)          :: fwhm            ! fwhm of envelope
     real(rk)          :: amplitude       ! amplitude of field
     real(rk)          :: omega           ! Frequency 
     real(rk)          :: alpha           ! gaussian exponent (needs to be initialized)
     real(rk)          :: phi             ! Laser carrier phase (radians)
     real(rk)          :: t0              ! center of the envelope
     real(rk)          :: toff1           ! start of decay to 0
     real(rk)          :: toff2           ! decay at 0
   end type EFPulse

  contains

    ! 
    ! initialize a pulse to sensible defaults, given a 
    ! fwhm and frequency
    !
    subroutine initializePulse(Pulse, fwhm, freq, amp, phase)
      type(EFPulse), intent(inout) :: Pulse
      real(rk),intent(in)          :: fwhm            ! fwhm of envelope
      real(rk),intent(in)          :: freq           ! Frequency 
      real(rk),intent(in)          :: amp
      real(rk),intent(in)          :: phase

      Pulse%fwhm = fwhm
      Pulse%omega = freq * twopi
      Pulse%amplitude = amp
      Pulse%phi = phase
      Pulse%toff1 = 0.5*fwhm*sqrt(log(100._rk)/log(2._rk)) ! gaussian up to 0.01 of maximum
      Pulse%toff2 = Pulse%toff1 + 0.1_rk*Pulse%toff1    ! stretch time to zero rapidly (by default)
      Pulse%t0    = Pulse%toff2

      write(out,"('t0=',f12.4,' 0.5*fwhm=',f12.4,' toff1=',f12.4,' toff2=',f12.4)")Pulse%t0,0.5*fwhm,Pulse%toff1,Pulse%toff2

      Pulse%alpha = initializeAlpha(Pulse)

    end subroutine initializePulse

    !
    !
    !
    function initializeAlpha( Pulse ) result(ialpha)
      type(EFPulse), intent(in)  :: Pulse
      real(ark)                  :: ialpha
      real(ark)                  :: alpha0,alpha
      !
      alpha0     = alpha_iteration(0._ark)
      alpha      = alpha_iteration(alpha0)
      ialpha = alpha0
      fixed_point: do while(abs(ialpha-alpha)>=100._ark*spacing(alpha))
        ialpha = alpha
        alpha  = alpha_iteration(alpha)
      end do fixed_point
      !
      write (out,"(/'Gaussian vector-potential envelope exponent = ',g19.11)") ialpha
      write (out,"( '    The long-pulse approximate GVP exponent = ',g19.11/)") alpha0
    
      return

      contains 

      real(ark) function alpha_iteration(alpha)
        real(ark), intent(in) :: alpha
        !
        alpha_iteration = (2._ark/Pulse%fwhm**2) * (log(2._ark) + log(1._ark + (Pulse%fwhm*alpha/Pulse%omega)**2))
      end function alpha_iteration

    end function initializeAlpha

    !
    ! Evaluate pulse of the form: cos(w1*t)*sin(Pi*t/carrier)^2
    !
    function evaluateSinPulse(Pulse,t,rwa) result(ElectricField)
     type(EFPulse),intent(in)   :: Pulse
     real(ark),intent(in)       :: t
     logical,intent(in)         :: rwa
     real(ark)                  :: ElectricField
     real(ark)                  :: costerm

     costerm = 1._ark
     if(.not.rwa)costerm = cos(Pulse%omega*t)
     ElectricField = Pulse%amplitude* costerm * (sin(t*Pi/(2.*Pulse%fwhm))**2)

     return
    end function evaluateSinPulse

    !
    ! Evaluate the external field at time ft
    !
    function evaluateExternalField( Pulse,ft,rwa ) result(ElectricField)
      type(EFPulse) , intent(in) :: Pulse      ! Object containing pulse parameterization
      real(ark), intent(in)      :: ft         ! Current time
      logical,intent(in)         :: rwa        ! employ rotating wave approximation
      real(ark)                  :: ElectricField
      real(ark)                  :: dt = 0.0020_rk
      integer                    :: out=6
      !
      ElectricField = 0._ark
      !
      ! Numerical differentiation; don't forget the overall minus sign
      !
     
      if(rwa)then
       ElectricField = Pulse%omega * GaussianVP(ft,rwa)
      else
       ElectricField = -(GaussianVP(ft+0.1_ark*dt,rwa) - GaussianVP(ft-0.1_ark*dt,rwa))/(0.2_ark*dt)
      endif
    
      return
      !
      contains
        !
        !  Gaussian envelope of the vector-potential, with guaranteed switch-off
        !
        real(ark) function GaussianVP(time,rwa)
          real(ark), intent(in) :: time       ! Time
          logical,intent(in)    :: rwa        ! employ rotating wave
          real(ark)             :: t_carrier  ! Time used for the carrier part
          real(ark)             :: t_envelope ! Time used for the envelope part
          real(ark)             :: t_red      ! Time reduced into the [0:1] range for stretching
          real(ark)             :: exp_arg    ! Argument for the exponential envelope
          real(ark)             :: ph         ! Carrier phase at this time
          !
          GaussianVP = 0._ark
          t_carrier  = time - Pulse%t0
          if (Pulse%toff2<=0._ark .or. abs(t_carrier)<=Pulse%toff2) then
            !
            !  The field is not zero; decide the effective envelope time
            !
            t_envelope = t_carrier
            if (Pulse%toff1<Pulse%toff2 .and. abs(t_envelope)>Pulse%toff1) then
              !
              !  We are in the envelope stretch zone; scale the envelope time.
              !  We do not need the sign for the envelope part, so do not bother with that
              !
              t_red      = (abs(t_envelope) - Pulse%toff1)/(Pulse%toff2-Pulse%toff1)
              t_red      = min(1._ark-spacing(10._ark),t_red)
              t_envelope = Pulse%toff1 + (2._ark/pi) * (Pulse%toff2-Pulse%toff1) * tan((pi/2._ark) * t_red)
            end if
            exp_arg = Pulse%alpha*t_envelope**2
            if (exp_arg<=max_exp) then
              ph = Pulse%omega * t_carrier + Pulse%phi
              if(rwa) then
               GaussianVP = (Pulse%amplitude/Pulse%omega) * exp(-exp_arg)
              else
               GaussianVP = (Pulse%amplitude/Pulse%omega) * exp(-exp_arg) * cos(ph)
              endif
            end if
          end if
        end function GaussianVP

    end function evaluateExternalField

 end module pulse
