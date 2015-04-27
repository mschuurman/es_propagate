!
!  Element properties
!
  module elements
    use accuracy
    implicit none
    private
    public get_element_mass
    !
    contains
    !
    function get_element_mass(name) result(mass)
      character(len=*), intent(in) :: name
      real(rk)                     :: mass
      !
      select case (name)
        case default;
          write (out,"('Element ',a,' is not in the database')") trim(name)
          stop 'elements%get_element_mass'
        case ('H','h');        mass =  1.0079_rk
        case ('He','he','HE'); mass =  4.00260_rk
        case ('Li','li','LI'); mass =  6.941_rk
        case ('Be','be','BE'); mass =  9.01218_rk
        case ('B','b');        mass = 10.811_rk
        case ('C','c');        mass = 12.011_rk
        case ('N','n');        mass = 14.0067_rk
        case ('O','o');        mass = 15.9994_rk
        case ('F','f');        mass = 18.9984_rk
        case ('Ne','ne','NE'); mass = 20.1797_rk
        case ('Na','na','NA'); mass = 22.98977_rk
        case ('P','p');        mass = 30.9738_rk
        case ('S','s');        mass = 32.065_rk
        case ('Cl','cl','CL'); mass = 35.453_rk
        case ('I','i');        mass = 126.90447_rk
      end select
    end function get_element_mass
  end module elements
