subroutine setup( boxlen, kmax, num, q, r, eta, rksqmax, rsqmax, fullhi, hlab2, fourpicell )

implicit none

double precision, dimension(3), intent(out) :: boxlen
integer, dimension(3), intent(out) :: kmax
integer, intent(out) :: num
double precision, dimension(:,:), allocatable, intent(out) :: r
double precision, dimension(:), allocatable, intent(out) :: q
double precision, intent(out) :: eta
integer :: nspec,nunitcellx,nunitcelly,nunitcellz
double precision :: rsqmax,rcut,rksqmax       !Several cutoffs
double precision, dimension(3,3) :: h,fullhi,fullhit,fullh !Several cell matrices
double precision, dimension(10) :: bh, b
integer :: i,j,l,m,nx,ny,nz
double precision :: chgsum, fourpicell
double precision, dimension(:,:), allocatable :: bx, by, bz   
double precision, dimension(:), allocatable :: chg
integer, allocatable, dimension(:) :: nsp
double precision, dimension(3,3) :: h3, hi3
double precision, dimension(3,3), intent(out) :: hlab2
integer, allocatable, dimension(:) :: nunitcellpos !Number of ions of each species in the unit cell
double precision :: etainpt   !Ewald smearing parameters
double precision :: conv1,convfac  !Convergence factors for reciprocal space summations...
character(len=20), allocatable, dimension(:) :: cellcoordfile   !X.mat, etc... 
double precision :: a0, b0, c0
double precision :: vol3

interface

    subroutine dcell( h, b )
        implicit none
        double precision, dimension(3,3), intent(in) :: h
        double precision, dimension(10), intent(out) :: b
    end subroutine dcell

    subroutine invert( h, hi )
        implicit none
        double precision, dimension(3,3), intent(in) :: h
        double precision, dimension(3,3), intent(out) :: hi
    end subroutine invert

    subroutine boxreset( h, fullh, fullhi, boxlen, fourpicell )
        implicit none
        double precision, dimension(3,3), intent(in) :: h
        double precision, dimension(3,3), intent(out) :: fullh, fullhi
        double precision, dimension(3), intent(in) :: boxlen
        double precision, intent(out) :: fourpicell
    end subroutine boxreset

    subroutine kset( kmax, conv1, convfac, eta, rksqmax, fullhi )
        implicit none
        integer, dimension(3), intent(out) :: kmax
        double precision, dimension(3,3), intent(in) :: fullhi !Several cell matrices
        double precision, intent(out) :: rksqmax     !Reciprocal space cutoff
        double precision, intent(in) :: eta    !Ewald smearing parameters
        double precision, intent(in) :: conv1, convfac  !Convergence factors for reciprocal space
    end subroutine kset

end interface

!  Read the file containing the simulation details

open(10,file='new_runtime.inpt',status='old')
read(10,*) nspec
!  Allocate all nspec-dimensional arrays

allocate ( nsp(nspec),chg(nspec) )
allocate ( nunitcellpos(nspec), cellcoordfile(nspec) )
          
read(10,*) ( nsp(i), i=1, nspec )
read(10,*) ( chg(i), i=1, nspec )
read(10,*) etainpt
read(10,*) rcut
read(10,*) conv1
read(10,*) convfac

open (10,file='crystal_cell.inpt',status='old')
read(10,*) h(1,1),h(2,1),h(3,1)
read(10,*) h(1,2),h(2,2),h(3,2)
read(10,*) h(1,3),h(2,3),h(3,3)
read(10,*) boxlen(1)
read(10,*) boxlen(2)
read(10,*) boxlen(3)
read(10,*) nunitcellx
read(10,*) nunitcelly
read(10,*) nunitcellz
do i=1,nspec
  read(10,*) nunitcellpos(i)
  read(10,'(a)') cellcoordfile(i)
enddo
read(10,*) a0
read(10,*) b0
read(10,*) c0
close(10)

call dcell( h, b )

vol3 = boxlen(1)*boxlen(2)*boxlen(3)*b(10)
hlab2 = h
h3(:,1) = h(:,1)*boxlen(1)/(vol3**(1.0d0/3.0d0))
h3(:,2) = h(:,2)*boxlen(2)/(vol3**(1.0d0/3.0d0))
h3(:,3) = h(:,3)*boxlen(3)/(vol3**(1.0d0/3.0d0))

call invert( h3, hi3 )

! setup
l=0
chgsum=0.d0
num=sum(nsp)

allocate( r( num, 3 ) )
allocate( q( num ) )

do i=1,nspec
   do j=1,nsp(i)
      l=l+1
      q(l)=chg(i)
      chgsum=chgsum+(chg(i)*chg(i))
   enddo   
enddo    

allocate( bx( nspec, 7000), by( nspec, 7000), bz( nspec, 7000 ) )

do i=1,nspec
   open (7,file=cellcoordfile(i),status='old')
   do j=1,nunitcellpos(i)
      read (7,*) bx(i,j),by(i,j),bz(i,j)
   enddo    
   close(7)
enddo   

! the .mat files that contain the bx,y,z arrays are in the
! cell-frame and to be visualised, should be transformed into
! the lab-frame.
!
m=1
do i=1,nspec
   do j=1,nunitcellpos(i)
      do nx=1,nunitcellx
         do ny=1,nunitcelly
            do nz=1,nunitcellz
               r(m,1)=a0*(nx-1)+bx(i,j)*a0
               r(m,2)=b0*(ny-1)+by(i,j)*b0
               r(m,3)=c0*(nz-1)+bz(i,j)*c0
               m=m+1
            enddo   
         enddo   
      enddo   
   enddo   
enddo   
!
! Enforce periodic boundary conditions.
!
do i=1,num
    if(r(i,1).lt.0.0d0) r(i,1)=r(i,1)+boxlen(1)
    if(r(i,2).lt.0.0d0) r(i,2)=r(i,2)+boxlen(2)
    if(r(i,3).lt.0.0d0) r(i,3)=r(i,3)+boxlen(3)
    if(r(i,1).gt.boxlen(1)) r(i,1)=r(i,1)-boxlen(1)
    if(r(i,2).gt.boxlen(2)) r(i,2)=r(i,2)-boxlen(2)
    if(r(i,3).gt.boxlen(3)) r(i,3)=r(i,3)-boxlen(3)
enddo

call boxreset( h, fullh, fullhi, boxlen, fourpicell )
call dcell( fullh, bh )
rsqmax=rcut**2.0d0
fullhit=transpose(fullhi)
call dcell( fullhit, bh )

eta=etainpt/(2.0d0*rcut)

call kset( kmax, conv1, convfac, eta, rksqmax, fullhi )

return
end subroutine setup
