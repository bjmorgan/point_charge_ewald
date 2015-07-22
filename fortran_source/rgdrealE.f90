double precision function rgdrealE( num, q, r, boxlen, hlab2, eta, rsqmax )

!use commondata, only: engpetot

!implicit none

integer, intent(in) :: num
double precision, dimension( num ), intent(in) :: q
double precision, dimension( num, 3 ), intent(in) :: r
double precision, dimension( 3 ), intent(in) :: boxlen
double precision, dimension( 3) :: halfboxrec
double precision, dimension( 3, 3 ), intent(in) :: hlab2
double precision, intent(in) :: eta, rsqmax
integer :: i, j
double precision :: erfunc
double precision :: drsq, dr, errfuncr
double precision, dimension(3) :: drsav
double precision, dimension(3) :: drcf
                    
rgdrealE = 0.d0
halfboxrec = 2.0 / boxlen

do j=2, num
   do i=1, j-1

      drcf = r(i,:) - r(j,:)

      drcf(1)=drcf(1)-boxlen(1)*int( drcf(1)*halfboxrec(1) )
      drcf(2)=drcf(2)-boxlen(2)*int( drcf(2)*halfboxrec(2) )
      drcf(3)=drcf(3)-boxlen(3)*int( drcf(3)*halfboxrec(3) )

      drsav(1)=hlab2(1,1)*drcf(1)+hlab2(1,2)*drcf(2)+hlab2(1,3)*drcf(3)
      drsav(2)=hlab2(2,1)*drcf(1)+hlab2(2,2)*drcf(2)+hlab2(2,3)*drcf(3)
      drsav(3)=hlab2(3,1)*drcf(1)+hlab2(3,2)*drcf(2)+hlab2(3,3)*drcf(3)

      drsq=dot_product( drsav, drsav )

      if (drsq.ge.rsqmax) cycle

      dr=dsqrt(drsq)

      errfuncr = q(i) * q(j) * erfunc(eta*dr) / dr
      rgdrealE = rgdrealE + errfuncr
   enddo   
enddo   

end function

