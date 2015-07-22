subroutine kset( kmax, conv1, convfac, eta, rksqmax, fullhi )

implicit none

integer :: kmaxt
double precision :: rkmin
double precision, parameter :: twopi = 6.283185307
integer, dimension(3), intent(out) :: kmax
double precision, dimension(3,3) :: fullhit !Several cell matrices
double precision, dimension(3,3), intent(in) :: fullhi !Several cell matrices
double precision, intent(out) :: rksqmax     !Reciprocal space cutoff
double precision, dimension(10) :: bh
double precision, intent(in) :: eta    !Ewald smearing parameters
double precision, intent(in) :: conv1, convfac  !Convergence factors for reciprocal space summations...

fullhit=transpose(fullhi)

call dcell(fullhit,bh)

rkmin=(twopi*dmin1(dabs(bh(7)),dabs(bh(8)),dabs(bh(9))))
rksqmax=-log(conv1)*4.0d0*eta*eta/(rkmin*rkmin)
kmaxt=int(dsqrt(rksqmax))

kmax(1)=int((rkmin/(twopi*bh(7)))*float(kmaxt))+1
kmax(2)=int((rkmin/(twopi*bh(8)))*float(kmaxt))+1
kmax(3)=int((rkmin/(twopi*bh(9)))*float(kmaxt))+1

rksqmax=rksqmax*convfac

return
end subroutine
