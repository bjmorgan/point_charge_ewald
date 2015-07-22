subroutine shake

use boxdata, only: h3,hi3,h2

implicit none

integer :: j,k,l,it
double precision :: numer,detmh,lambda,denom,r

double precision, dimension(3,3) :: hih

do it=1,10
   hi3(1,1)=h3(2,2)*h3(3,3)-h3(3,2)*h3(2,3)
   hi3(2,1)=h3(3,1)*h3(2,3)-h3(2,1)*h3(3,3)
   hi3(3,1)=h3(2,1)*h3(3,2)-h3(3,1)*h3(2,2)
   hi3(1,2)=h3(3,2)*h3(1,3)-h3(1,2)*h3(3,3)
   hi3(2,2)=h3(1,1)*h3(3,3)-h3(3,1)*h3(1,3)
   hi3(3,2)=h3(3,1)*h3(1,2)-h3(1,1)*h3(3,2)
   hi3(1,3)=h3(1,2)*h3(2,3)-h3(2,2)*h3(1,3)
   hi3(2,3)=h3(2,1)*h3(1,3)-h3(1,1)*h3(2,3)
   hi3(3,3)=h3(1,1)*h3(2,2)-h3(2,1)*h3(1,2)

   detmh=h3(1,1)*hi3(1,1)+h3(1,2)*hi3(2,1)+h3(1,3)*hi3(3,1)
   r=0.d0
   if (abs(detmh).gt.0.d0) then
       r = 1.d0/detmh
   endif

   hi3=r*hi3

   hih=matmul(hi3,h2)

   numer=detmh-1.0d0
   denom=detmh*(hih(1,1)+hih(2,2)+hih(3,3))
   lambda=-numer/denom

   h3=h3+(lambda*h2)
enddo   

return
end subroutine
