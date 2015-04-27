module accuracy
!
!  Sundry mathematical and physical definitions.
!
  implicit none
  private
  public sik, ik, hik, rk, ark, input, out, safe_max, max_exp, pi, twopi, fourpi, sqrtpi
  public sqrt2pi, femtosecond, abohr, h2ev, h2cm, k_Boltzmann, Hartree, bar, vlight
  public R_gas, N_Avogadro, unified_atomic_mass_unit, electron_mass_in_amu
  public sik_bytes, ik_bytes, hik_bytes, rk_bytes, ark_bytes
  public accuracyInitialize
!
  integer, parameter :: sik         = selected_int_kind(4)       ! Small integers
  integer, parameter :: ik          = selected_int_kind(8)       ! "Normal" integers. This must map on
                                                                 ! C "int" type, or the DX interface won't
                                                                 ! work correctly.
  integer, parameter :: hik         = selected_int_kind(15)      ! "Pointer" integers - sufficient to store
                                                                 ! memory address
  integer, parameter :: rk          = selected_real_kind(14,17)  ! "Normal" reals and complex (complexi? :-)
  integer, parameter :: ark         = selected_real_kind(14,17)  ! "Accurate" reals and complex (complexi? :-)
  integer, parameter :: input       = 5                          ! Standard input I/O channel
  integer, parameter :: out         = 6                          ! Output I/O channel

  real(rk)           :: safe_max                                 ! Largest number we want to work with
  real(rk)           :: max_exp                                  ! Largest number OK for exponentiating
  real(rk)           :: pi, twopi, fourpi, sqrtpi, sqrt2pi       ! Pi, 2*Pi, 4*Pi, SQRT(2*Pi)
  integer(ik)        :: sik_bytes, ik_bytes, hik_bytes, rk_bytes, ark_bytes
  !
  !  Physical constants - must be here to guarantee consistency
  !
  real(ark), parameter :: femtosecond = 41.341373336561364_ark   ! Conversion factor: from femtoseconds to au[t]
  real(ark), parameter :: abohr       = 0.5291772083_ark         ! Conversion factor: from au[l] (Bohrs) to Angstrom
  real(ark), parameter :: h2ev        = 27.2113845_ark           ! Conversion factor: from Hartree to electron-volts
  real(ark), parameter :: h2cm        = 219474.6313705_ark       ! Conversion factor: from Hartree to cm^-1 (wavenumbers)
  real(ark), parameter :: k_Boltzmann = 1.0_rk/315774.65_ark     ! Conversion factor: from Kelvin to Hartree (aka Boltzmann constant)
  real(ark), parameter :: Hartree     = 4.35974417e-18_ark       ! Conversion factor: from au[e] (Hartree) to Joules
  real(ark), parameter :: bar         = 101325._ark              ! Conversion factor: from bars to Pascal (aka standard pressure)
  real(ark), parameter :: vlight      = 137.0359991_ark          ! Speed of light in atomic units; also the
                                                                 ! inverse of the fine structure constant
  real(ark), parameter :: R_gas       = 8.314472_ark             ! Molar gas constant, in J/mol-K
  real(ark), parameter :: N_Avogadro  = 6.0221415e23_ark         ! Avogadro number, in particles/mole
  real(ark), parameter :: unified_atomic_mass_unit = 1.660538782e-27_ark ! Unified atomic mass unit, in kg
  real(ark), parameter :: electron_mass_in_amu     = 5.4857990943e-4_ark ! Electron mass in u-amu
  !
  contains

  subroutine accuracyInitialize
    safe_max = huge(1.0_rk)**(0.25_rk)
    max_exp  = log(safe_max)
    pi       = 4.0_rk * atan2(1.0_rk,1.0_rk)
    twopi    = 2.0_rk * pi
    fourpi   = 4.0_rk * pi
    sqrt2pi  = sqrt(twopi)
    sqrtpi   = sqrt(pi)
    !
    !  Type sizes
    !
    sik_bytes = int_bytes(radix(1_sik),digits(1_sik))
    ik_bytes  = int_bytes(radix(1_ik ),digits(1_ik ))
    hik_bytes = int_bytes(radix(1_hik),digits(1_hik))
    rk_bytes  = real_bytes(radix(1._rk) ,digits(1._rk ),maxexponent(1._rk )-minexponent(1._rk ))
    ark_bytes = real_bytes(radix(1._ark),digits(1._ark),maxexponent(1._ark)-minexponent(1._ark))
    write (out,"('Compiled integer sizes: ',3(i4,1x),' bytes')") sik_bytes, ik_bytes, hik_bytes
    write (out,"('Compiled real sizes: ',2(i4,1x),' bytes')") rk_bytes, ark_bytes
  end subroutine accuracyInitialize

  integer(ik) function int_bytes(radix,digits)
    integer(ik), intent(in) :: radix, digits
    !
    real(ark) :: bits
    !
    bits = digits*(log(real(radix,kind=ark))/log(2._ark))
    int_bytes = ceiling((1+bits)/8._ark)
  end function int_bytes

  integer(ik) function real_bytes(radix,digits,exp_range)
    integer(ik), intent(in) :: radix, digits, exp_range
    !
    real(ark) :: exp_bits, mant_bits
    !
    exp_bits   = log(real(exp_range,kind=ark))/log(2._ark)
    mant_bits  = digits*(log(real(radix,kind=ark))/log(2._ark))
    real_bytes = ceiling((exp_bits+mant_bits)/8._ark)
  end function real_bytes

end module accuracy
