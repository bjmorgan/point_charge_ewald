double precision function rgdrecipE( num, q, r, kmax, boxlen, eta, rksqmax, fullhi, fourpicell )

implicit none

integer, intent(in) :: num
double precision, dimension(num), intent(in) :: q
double precision, dimension(num,3), intent(in) :: r
double precision, dimension(3), intent(in) :: boxlen
double precision, dimension(3,3), intent(in) :: fullhi
double precision, intent(in) :: eta
double precision, intent(in) :: rksqmax
integer, dimension(3), intent(in) :: kmax
double precision, intent(in) :: fourpicell

integer :: i,l,m,n,ll,mm,nn,nmin,mmin
double precision :: rl,rkx1,rky1,rkz1,rm,rkx2,rky2, &
                    rkz2,rn,rkx3,rky3,rkz3,xkk,ckcs,ckss, &
                    akk, efac
double precision, dimension(num) :: clmall,slmall,ckcnoqall,cksnoqall
double precision, dimension( num, 0:kmax(1)+1 ) :: elcall, elsall
double precision, dimension( num, 0:kmax(2)+1 ) :: emcall, emsall
double precision, dimension( num, 0:kmax(3)+1 ) :: encall, ensall
double precision, dimension(3) :: twopibox
double precision, parameter :: twopi = 6.283185307
double precision :: chgsum, chgcorrec, etapi, etaconst
double precision :: engpetot 

twopibox = twopi / boxlen
etapi = 2.0d0 * eta / 1.7724538509055160249d0
etaconst = -1.0d0 / ( 4.0d0 * eta * eta )
engpetot = 0.0d0

chgsum = dot_product( q, q )
chgcorrec=chgsum*etapi*0.5d0

do i=1,num
    elcall(i,0) = 1.0d0
    emcall(i,0) = 1.0d0
    encall(i,0) = 1.0d0
    elsall(i,0) = 0.0d0
    emsall(i,0) = 0.0d0
    ensall(i,0) = 0.0d0

    elcall(i,1) = cos(twopibox(1)*r(i,1))
    emcall(i,1) = cos(twopibox(2)*r(i,2))
    encall(i,1) = cos(twopibox(3)*r(i,3))
    elsall(i,1) = sin(twopibox(1)*r(i,1))
    emsall(i,1) = sin(twopibox(2)*r(i,2))
    ensall(i,1) = sin(twopibox(3)*r(i,3))
enddo    

do l=2,kmax(1)
   do i=1,num
      elcall(i,l)=elcall(i,l-1)*elcall(i,1)-elsall(i,l-1)*elsall(i,1)
      elsall(i,l)=elsall(i,l-1)*elcall(i,1)+elcall(i,l-1)*elsall(i,1)
   enddo    
enddo    

do l=2,kmax(2)
   do i=1,num
      emcall(i,l)=emcall(i,l-1)*emcall(i,1)-emsall(i,l-1)*emsall(i,1)
      emsall(i,l)=emsall(i,l-1)*emcall(i,1)+emcall(i,l-1)*emsall(i,1)
   enddo    
enddo    

do l=2,kmax(3)
   do i=1,num
      encall(i,l)=encall(i,l-1)*encall(i,1)-ensall(i,l-1)*ensall(i,1)
      ensall(i,l)=ensall(i,l-1)*encall(i,1)+encall(i,l-1)*ensall(i,1)
   enddo    
enddo    

mmin=0
nmin=1

do ll=0,kmax(1)
   l=iabs(ll)
   rl=dble(ll)*twopi

   rkx1=rl*fullhi(1,1)
   rky1=rl*fullhi(1,2)
   rkz1=rl*fullhi(1,3)

   do mm=mmin,kmax(2)
      m=iabs(mm)
      rm=dble(mm)*twopi

      rkx2=rkx1+rm*fullhi(2,1)
      rky2=rky1+rm*fullhi(2,2)
      rkz2=rkz1+rm*fullhi(2,3)

      if(mm.ge.0) then
         clmall(:)=elcall(:,l)*emcall(:,m)-elsall(:,l)*emsall(:,m)
         slmall(:)=elsall(:,l)*emcall(:,m)+emsall(:,m)*elcall(:,l)
      else
         clmall(:)=elcall(:,l)*emcall(:,m)+elsall(:,l)*emsall(:,m)
         slmall(:)=elsall(:,l)*emcall(:,m)-emsall(:,m)*elcall(:,l)
      endif

      do nn=nmin,kmax(3)
         n=iabs(nn)
         rn=dble(nn)*twopi

         rkx3=rkx2+rn*fullhi(3,1)
         rky3=rky2+rn*fullhi(3,2)
         rkz3=rkz2+rn*fullhi(3,3)

         xkk=rkx3*rkx3+rky3*rky3+rkz3*rkz3

         if(xkk.le.rksqmax) then
            if(nn.ge.0) then
               ckcnoqall(:)=clmall(:)*encall(:,n)-slmall(:)*ensall(:,n)
               cksnoqall(:)=slmall(:)*encall(:,n)+clmall(:)*ensall(:,n)
            else
               ckcnoqall(:)=clmall(:)*encall(:,n)+slmall(:)*ensall(:,n)
               cksnoqall(:)=slmall(:)*encall(:,n)-clmall(:)*ensall(:,n)
            endif

            ckcs=sum(ckcnoqall*q)
            ckss=sum(cksnoqall*q)

            akk=exp(etaconst*xkk)/xkk

            efac=akk*(ckss*ckss+ckcs*ckcs)
            engpetot = engpetot + efac

         endif
      enddo    
      nmin=-kmax(3)
   enddo    
   mmin=-kmax(2)
enddo    
rgdrecipE = fourpicell * engpetot - chgcorrec

return
end 
