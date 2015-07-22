!************************************************************
!   Sets up the ion velocities on a Gaussian distribution.
!************************************************************

SUBROUTINE velset

USE commondata, ONLY: nspec,nsp,vx,vy,vz,vartrans,dummy

IMPLICIT NONE

INTEGER :: n,i,j
DOUBLE PRECISION :: zero,gauss

n=1
do i=1,nspec
   do j=1,nsp(i)
      vx(n)=gauss(dummy)*vartrans(i)
      vy(n)=gauss(dummy)*vartrans(i)
      vz(n)=gauss(dummy)*vartrans(i)
      n=n+1
   enddo   
enddo   
!
! Rescale the velocities.
!
call rescale

write(*,*)'**** Velocities set up ****'

return
END SUBROUTINE
