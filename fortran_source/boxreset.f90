!***********************************************************
!
! Sets/resets all box length dependent variables.
!
!***********************************************************
subroutine boxreset( h, fullh, fullhi, boxlen, fourpicell )

implicit none

double precision, intent(out) :: fourpicell
double precision, dimension(3,3), intent(in) :: h
double precision, dimension(3,3), intent(out) :: fullh
double precision, dimension(3), intent(in) :: boxlen
integer :: i
double precision, dimension(3) :: lengths
double precision,  dimension(10) :: dcellinfo
double precision :: cellvol
double precision, dimension(3,3), intent(out) :: fullhi 
double precision, parameter :: fourpi = 12.56637061d0

interface
    subroutine invert( h, hi )
        implicit none
        double precision, dimension(3,3), intent(in) :: h
        double precision, dimension(3,3), intent(out) :: hi
    end subroutine invert
end interface

lengths(1)=boxlen(1)
lengths(2)=boxlen(2)
lengths(3)=boxlen(3)

call dcell(h,dcellinfo)

cellvol=dcellinfo(10)*boxlen(1)*boxlen(2)*boxlen(3)
fourpicell=fourpi/cellvol
! Invert the cell matrix for use with the force calculations
! and the Ewald summation.
do i=1,3
   fullh(i,:)=h(i,:)*lengths(:)
enddo
call invert(fullh,fullhi)

return
end subroutine
