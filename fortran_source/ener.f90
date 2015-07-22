double precision function ener( num, q, r, kmax, boxlen, eta, rksqmax, fullhi, fourpicell, hlab2, rsqmax )

implicit none

integer, intent( in ) :: num
double precision, dimension( num, 3 ), intent(in) :: r
double precision, dimension( num ), intent(in) :: q
double precision, dimension( 3 ), intent(in) :: boxlen
integer, dimension( 3 ), intent(in) :: kmax
double precision, intent(in) :: eta, fourpicell, rksqmax, rsqmax
double precision, dimension( 3, 3 ), intent(in) :: fullhi, hlab2

interface

    function rgdrealE( num, q, r, boxlen, hlab2, eta, rsqmax )
        implicit none
        double precision :: rgdrealE
        integer, intent(in) :: num
        double precision, dimension( num ), intent(in) :: q
        double precision, dimension( num, 3 ), intent(in) :: r
        double precision, dimension( 3 ), intent(in) :: boxlen
        double precision, dimension( 3, 3 ), intent(in) :: hlab2
        double precision, intent(in) :: eta, rsqmax
    end function rgdrealE

    function rgdrecipE( num, q, r, kmax, boxlen, eta, rksqmax, fullhi, fourpicell )
        implicit none
        double precision rgdrecipE
        integer, intent(in) :: num
        double precision, dimension(num), intent(in) :: q
        double precision, dimension(num,3), intent(in) :: r
        double precision, dimension(3), intent(in) :: boxlen
        double precision, dimension(3,3), intent(in) :: fullhi
        double precision, intent(in) :: eta
        double precision, intent(in) :: rksqmax
        integer, dimension(3), intent(in) :: kmax
        double precision, intent(in) :: fourpicell
    end function rgdrecipE

end interface

ener = rgdrecipE( num, q, r, kmax, boxlen, eta, rksqmax, fullhi, fourpicell ) + &
       rgdrealE( num, q, r, boxlen, hlab2, eta, rsqmax )

return
end function ener
