program main

!use ewald_mod

implicit none 

double precision, dimension(3) :: boxlen
integer, dimension(3) :: kmax
double precision, dimension(:,:), allocatable :: r
double precision, dimension(:), allocatable :: q
integer :: num
double precision :: eta, fourpicell, rsqmax, rksqmax
double precision, dimension( 3, 3 ) :: fullhi, hlab2

interface
    subroutine setup( boxlen, kmax, num, q, r, eta, rksqmax, rsqmax, fullhi, hlab2, fourpicell )
        implicit none
        integer, dimension(3), intent(out) :: kmax
        integer, intent(out) :: num
        double precision, dimension( :, : ), allocatable, intent(out) :: r
        double precision, dimension( : ), allocatable, intent(out) :: q
        double precision, dimension(3,3), intent(out) :: hlab2
        double precision, dimension(3) :: boxlen
        double precision :: eta   !Ewald smearing parameter
        double precision :: rksqmax, rsqmax, fourpicell
        double precision, dimension(3,3) :: fullhi
    end subroutine setup

    function ener( num, q, r, kmax, boxlen, eta, rksqmax, fullhi, fourpicell, hlab2, rsqmax )
        implicit none
        double precision :: ener
        integer, intent(in) :: num
        double precision, dimension( num ), intent(in) :: q
        double precision, dimension( num, 3 ), intent(in) :: r
        integer, dimension(3), intent(in) :: kmax
        double precision, dimension(3), intent(in) :: boxlen
        double precision, intent(in) :: eta, rksqmax, fourpicell, rsqmax
        double precision, dimension(3,3), intent(in) :: fullhi, hlab2
    end function ener
end interface

call setup( boxlen, kmax, num, q, r, eta, rksqmax, rsqmax, fullhi, hlab2, fourpicell )

write(6,*) num
write(6,*) q(1)
write(6,*) kmax
write(6,*) boxlen
write(6,*) eta
write(6,*) rksqmax
write(6,*) fullhi
write(6,*) fourpicell
write(6,*) hlab2
write(6,*) rsqmax

write(6,*) ener( num, q, r, kmax, boxlen, eta, rksqmax, fullhi, fourpicell, hlab2, rsqmax )

write(6,*) "-1.1470399097117384E+02 <<"
stop
end program
