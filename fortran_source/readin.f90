!************************************************************
!   Reads in the data from the .inpt file.
!************************************************************

subroutine readin( boxlen, kmax, a0, b0, c0 )

!use commondata
!use boxdata
!use recipdata, only: convfac,conv1
!use fshift

implicit none

integer :: i, j, k, ii
integer :: num,nanion,ncation,nspec,nummon,nunitcellx,nunitcelly,nunitcellz,nionunit
integer, allocatable, dimension(:) :: nsp
double precision, allocatable, dimension(:) :: chg,q                 !charges
double precision :: relax,free,eps1,eps2,eps3,veps1,veps2,veps3,vol3,relaxb,relaxb2
double precision, dimension(3,3) :: h2, hi2, h3, hi3, vg2, vg3, hlab2, hlab3
integer, allocatable, dimension(:) :: nunitcellpos !Number of ions of each species in the unit cell
double precision :: etainpt,eta,etapi,etaconst,etasq    !Ewald smearing parameters
double precision :: rsqmax,rsqmaxsr,rcut      !Several cutoffs
double precision :: conv1,convfac  !Convergence factors for reciprocal space summations...
double precision, dimension(3) :: boxlen
integer, dimension(3) :: kmax
double precision, dimension(3,3) :: h,hi,fullhi,fullhit,fullh !Several cell matrices
double precision, dimension(10) :: bh, b
CHARACTER(len=20), allocatable, dimension(:) :: cellcoordfile   !X.mat, etc... Will this work??
double precision :: a0,b0,c0,cellvol,fourpicell,cellvol3rec,fac

!  Read the file containing the simulation details

open(10,file='new_runtime.inpt',status='old')
read(10,*) nionunit  
read(10,*) nspec
!  Allocate all nspec-dimensional arrays

allocate ( nsp(0:nspec),chg(nspec) )
allocate ( nunitcellpos(nspec), cellcoordfile(nspec) )
          
read(10,*) ( nsp(i), i=1, nspec)
read(10,*) ( chg(i), i=1, nspec)
read(10,*)etainpt
read(10,*)rcut
read(10,*)conv1
read(10,*)convfac

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

call dcell(h,b)
vol3=boxlen(1)*boxlen(2)*boxlen(3)*b(10)
hlab3=h
hlab2=h 
h3(:,1)=h(:,1)*boxlen(1)/(vol3**(1.0d0/3.0d0))
h3(:,2)=h(:,2)*boxlen(2)/(vol3**(1.0d0/3.0d0))
h3(:,3)=h(:,3)*boxlen(3)/(vol3**(1.0d0/3.0d0))

call invert(h3,hi3)

return
end subroutine
