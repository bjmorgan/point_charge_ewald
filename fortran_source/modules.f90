module commondata

implicit none
save

! Constants

double precision, parameter :: PI=3.14159265358979323d0
double precision, parameter :: PISQ=PI*PI
double precision, parameter :: SQRPI=1.7724538509055160249d0
double precision, parameter :: PITHREEHALF=SQRPI*PI
double precision, parameter :: TWOPI=2.D0*PI 
double precision, parameter :: FOURPI=4.D0*PI
double precision, parameter :: EIGHTPI=8.D0*PI
double precision, parameter :: FOURPISQ=4.D0*PISQ
double precision, parameter :: ONETHIRD=1.D0/3.D0

! Conversion of units 

double precision, parameter ::  &
        emass=9.109534d-28,     &      !Electron mass in grams
        boltz=3.166829689d-6,   &      !(Inverse Kelvins)
        avo=6.022045d023               !Avogadro's number

! Functions

double precision, allocatable, dimension(:) :: x, y, z     ! particle coordinates
double precision, allocatable, dimension(:) :: vx, vy, vz  ! particle velocities
double precision, allocatable, dimension(:) :: frrx, frry, frrz  ! forces on particles 
double precision, allocatable, dimension(:) :: frrx3, frry3, frrz3  ! forces on particles 
! in fact they are the internal coordinates, ranging from 0 to 1
! ie   r=x a1 + y a2 + z a3  NOT r=x i + y j + z k
double precision, allocatable, dimension(:,:) :: dxsav, dysav, dzsav !Interatomic distances
double precision, allocatable, dimension(:,:) :: errfunc !Error function 
double precision :: stsrxx, stsrxy, stsrxz, stsryy, stsryz, stsrzz
double precision :: stcxx, stcxy, stcxz, stcyy, stcyz, stczz
double precision :: stpxx, stpxy, stpxz, stpyy, stpyz, stpzz
double precision :: stp2xx, stp2xy, stp2xz, stp2yy, stp2yz, stp2zz
double precision :: stpsrxx, stpsrxy, stpsrxz, stpsryy, stpsryz, stpsrzz
double precision :: stewxx, stewxy, stewxz, stewyy, stewyz, stewzz, stewyx, stewzx, stewzy
double precision :: stpqquadxx, stpqquadyy, stpqquadzz, stpqquadxy, stpqquadxz, stpqquadyz
double precision :: stpdipquadxx, stpdipquadyy, stpdipquadzz, stpdipquadxy, stpdipquadxz, stpdipquadyz
double precision :: stpquadquadxx, stpquadquadyy, stpquadquadzz, stpquadquadxy, stpquadquadxz, stpquadquadyz
double precision :: stqsrxx, stqsryy, stqsrzz, stqsrxy, stqsrxz, stqsryz
double precision :: pext, pextstruct
                         !Components of the stress tensor (Change in future?)
double precision, dimension(3,3) :: pint2   !Internal pressure tensor

double precision :: dtime   !Time step
logical, allocatable, dimension(:) :: fixed_posx, fixed_posy, fixed_posz

integer, allocatable, dimension(:) :: ntype !1 for 1st comp, 2 for 2nd comp ...
integer, allocatable, dimension(:) :: nsp   !Number of ions of each species...
integer, allocatable, dimension(:) :: nunitcellpos !Number of ions of each species in the unit cell
integer, allocatable, dimension(:) :: nmon   !targets ions for which to show detailed information...
integer :: num,nanion,ncation,nspec,nummon,nunitcellx,nunitcelly,nunitcellz,nionunit
integer :: nconj

double precision, allocatable, dimension(:) :: elecx,elecy,elecz !Electric fields
double precision, allocatable, dimension(:) :: elecxq,elecyq,eleczq !Electric fields
double precision, allocatable, dimension(:) :: elecxsr,elecysr,eleczsr !Electric fields (short-range contribution)
double precision, allocatable, dimension(:) :: exx,eyy,ezz,exy,exz,eyz !Electric field gradients
double precision, allocatable, dimension(:) :: exxq,eyyq,ezzq,exyq,exzq,eyzq !Electric field gradients
double precision, allocatable, dimension(:) :: exxsr,eyysr,ezzsr,exysr,exzsr,eyzsr !Electric field gradients (short-range contribution)

double precision, allocatable, dimension(:) :: amass,hmass,recamass  !masses
double precision, allocatable, dimension(:) :: chg,q                 !charges
double precision, allocatable, dimension(:) :: xmu,ymu,zmu     !dipoles
double precision, allocatable, dimension(:) :: srdipx,srdipy,srdipz    !SR dipoles
double precision, allocatable, dimension(:) :: asdipx,asdipy,asdipz    !LR dipoles
double precision, allocatable, dimension(:) :: quadxx,quadyy,quadzz,quadxy,quadxz,quadyz !quadrupoles
double precision, allocatable, dimension(:) :: srquadxx,srquadyy,srquadzz,srquadxy,srquadxz,srquadyz    !SR quadrupoles
double precision, allocatable, dimension(:) :: asquadxx,asquadyy,asquadzz,asquadxy,asquadxz,asquadyz    !LR quadrupoles

double precision, allocatable, dimension(:) :: xk1,xk2,xk3,xk4,engeff !polarizabilites
double precision, allocatable, dimension(:) :: alppolar,Bpolar,Cpolar,gammapolar
double precision, allocatable, dimension(:) :: alppolar1,alppolardel,alppolareff, alphadelta
double precision, allocatable, dimension(:) :: Bpolar1,Cpolar1,gammapolar1
double precision                            :: deltalpha

double precision :: engpetot, tranke, tke, trantemp, tranke_old !Potential and kinetic energies, and temperatures
double precision, allocatable, dimension(:) :: reseng,quadeng,dipquadeng,dipsqeng !PIM
double precision :: selfengtot,epsselfengtot,quaimselfengtot !CIM/AIM selfenergies

double precision, allocatable, dimension(:) :: delta           !CIM distortion
double precision, allocatable, dimension(:) :: epsilonx,epsilony,epsilonz !DAIM distortion
double precision, allocatable, dimension(:) :: quaimxx,quaimyy,quaimzz,quaimxy,quaimxz,quaimyz     !QUAIM distortion

double precision, allocatable, dimension(:,:) :: ftalp,ftbeta,ftgamma,ftb,ftb2,ftb3
double precision, allocatable, dimension(:,:) :: rvph,rvpn,rvpr4,ftalpx,ftbx
integer         , allocatable, dimension(:,:) :: nrpower
                                              ! Short-range repulsion potential
double precision, allocatable, dimension(:,:) :: engft1,engft2,engft3
                                !CIM CG "Forces"
double precision, allocatable, dimension(:,:) :: engft1dotx,engft1doty,engft1dotz
double precision, allocatable, dimension(:,:) :: engft2dotx,engft2doty,engft2dotz
double precision, allocatable, dimension(:,:) :: engft3dotx,engft3doty,engft3dotz
                                !DAIM CG "forces"
double precision, allocatable, dimension(:,:) :: engft1dotxx,engft1dotyy,engft1dotzz
double precision, allocatable, dimension(:,:) :: engft1dotxy,engft1dotxz,engft1dotyz
double precision, allocatable, dimension(:,:) :: engft2dotxx,engft2dotyy,engft2dotzz
double precision, allocatable, dimension(:,:) :: engft2dotxy,engft2dotxz,engft2dotyz
double precision, allocatable, dimension(:,:) :: engft3dotxx,engft3dotyy,engft3dotzz
double precision, allocatable, dimension(:,:) :: engft3dotxy,engft3dotxz,engft3dotyz
                                !QUAIM CG "Forces"

double precision, allocatable, dimension(:,:) :: ftc,ftd         !(dispersion potential)
double precision, allocatable, dimension(:,:) :: dddamp,dddamp2,dddamp3,dddamp4,dddamp5,dddamp6    !Dipole-dipole dispersion damping
double precision, allocatable, dimension(:,:) :: dqdamp,dqdamp2,dqdamp3,dqdamp4,dqdamp5,dqdamp6,dqdamp7,dqdamp8    !Dipole-quadrupole dispersion damping

double precision, allocatable, dimension(:,:) :: dampa,dampfac,fgb,fgc !polarization damping
integer, allocatable, dimension(:,:) :: nkdamp,nkfg            !polarization damping
double precision, allocatable, dimension(:,:) :: gausA, gausB, gausR ! gaussian shift parameters

double precision, allocatable, dimension(:) :: selfgam,selfB !CIM/AIM selfenergy parameters
double precision :: extraalpha,extrab
double precision, allocatable, dimension(:) :: selfeps,selfC,selfeps1,selfeps2
double precision, allocatable, dimension(:) :: selfquaim,selfH,selfquaim1

integer, dimension(100) :: nstpbrk
integer :: nmat1,nmat2
integer :: nmsdcalltime, nperrdf, nrdfcall, nrscale, npereng, npervel, nperfri, npercell, nperquench
integer :: nsofar,ntotstp,nrun,nstep
integer, dimension(6) :: cellconstraints !for structural relaxation of cell
double precision :: trantkb, gtrankin, dummy, amovefac
double precision, dimension(3,3) :: amove
double precision, allocatable, dimension(:) :: vartrans
double precision :: chgcorrec,dispcorrec1,dispcorrec2   !Ewald selfenergies
double precision :: etainpt,eta,etapi,etaconst,etasq    !Ewald smearing parameters
double precision :: rsqmax,rsqmaxsr,rcut,rcutsr       !Several cutoffs
double precision :: dispalpsq,dispalp,dispalp3,dispalp4,dispalp6  !For Dispersion Ewald summation
logical :: static
logical :: constrainedlog

CHARACTER(len=20), allocatable, dimension(:) :: cellcoordfile   !X.mat, etc... Will this work??
CHARACTER(len=20) :: fileout

double precision, dimension(0:1000) :: rdftot      !Total and partial g(r)'s
double precision, allocatable, dimension(:,:,:) :: rdfpart

integer :: nsteptramp, nsteppramp ! temperature / pressure ramp
double precision :: deltatramp, deltapramp, ewcorr

double precision :: relax,free,eps1,eps2,eps3,veps1,veps2,veps3,vol3,relaxb,relaxb2 !Thermostat-Barostat
double precision, dimension(5) :: vpzeta1,vpzeta2,vpzeta3,pzeta2,pzeta3    !Thermostat
double precision, dimension(5) :: vbzeta1,vbzeta2,vbzeta3,bzeta2,bzeta3    !barostat
double precision :: tcell,tvol,tbzeta,tpzeta,PeeVee,pzeta,bzeta,Pcomx,Pcomy,Pcomz
double precision :: CUEp,CUEp2,CUEprec,CUEp2rec,dom,CUEb,CUEb2,CUEbrec,CUEb2rec,Wgo,W,Wrec,Wgorec
double precision :: totmz0 ! target total dipole moment along z
double precision :: totmz  ! actual total dipole moment along z
double precision :: tol,ftol,tolaim,ftolaim,tolstruc,ftolstruc,engtol !Tolerances

logical :: rim, dippimlog, quadpimlog, cimlog, daimlog, quaimlog, rimlog, epplog, ftlog, xftlog, rvplog, gxftlog !Potential function
logical :: ooaimlog,oodamplog,cluslog,environmentalaimlog,environmentalpimlog !Levels of theory
logical, allocatable, dimension(:) :: polarizablelog,deformablelog !Level of theory for each species
logical :: displace, restart, rescalelog, velinit, relaxconfig, dynam, moveions, relaxcell, harmonic_quench     !Simulation options
logical :: conjgradlog,msdcalllog,ewlog,conjgradaimlog 
logical :: chargequadlog,dipquadlog,fullewaldlog    !Level of Ewald summation
logical :: veldumplog,crddumplog,chgdumplog,fulldumplog,rdfcall   !Write to disk
logical :: pereng,perfric,pervel,percell,nrscalelog             !Write to disk
logical :: forfl,nth,nib,nab,ortho      !Variable cell
logical :: verbose  !to control the amount of information written to disk
logical :: firstiter !To optimize the number of calls to cgener...
logical :: tramplog, pramplog ! T/p ramp
logical :: lscalc    ! Do the light scattering calculation
logical :: stoponerr !if true calculation will halt if multipoles do not converge

end module commondata

!----------------------------------------!

module boxdata

IMPLICIT NONE
SAVE

double precision, dimension(3,3) :: h,hi,fullhi,fullhit,fullh !Several cell matrices
double precision, dimension(3,3) :: h2, hi2, h3, hi3, vg2, vg3, hlab2, hlab3       
double precision, dimension(10) :: bh, b, bee

double precision :: boxlenx,boxleny,boxlenz,halfboxx,halfboxy,halfboxz
double precision :: halfboxxrec,halfboxyrec,halfboxzrec,halfboxminsq
double precision :: twopiboxx,twopiboxy,twopiboxz

double precision :: a0,b0,c0,cellvol,fourpicell,cellvol3rec,fac

END MODULE BOXDATA

!----------------------------------------!

MODULE RECIPDATA

IMPLICIT NONE
SAVE

double precision, allocatable, dimension(:,:) :: elcall,emcall,encall,elsall,emsall,ensall
                                     !sines and cosines for recipE (Change in future??)
double precision, allocatable, dimension(:,:) :: bdisp !for recipr. dispersion energy

double precision, allocatable, dimension(:,:) :: sk_ds !Structure factor
integer :: nbin !Structure factor related
integer, dimension(0:1000) :: norm_ds !Structure factor related

double precision :: qmueng,xmumueng,qquadeng,dipquadengrec,quadquadengrec
                                !Reciprocal space energy contributions
integer :: kmaxx,kmaxy,kmaxz,ksqmax    !Reciprocal space maximum vectors...
double precision :: conv1,convfac  !Convergence factors for reciprocal space summations...
double precision :: rksqmax     !Reciprocal space cutoff
double precision :: slens(10) !Scattering lens !Added by D. Marrocchelli 05/03/2008

END MODULE RECIPDATA
!----------------------------------------!

MODULE CGDATA    

IMPLICIT NONE
SAVE

double precision, allocatable, dimension(:) :: PCOMAIM,XICOMAIM
double precision, allocatable, dimension(:) :: PCOMSTRUCT,XICOMSTRUCT

END MODULE CGDATA

!--------------------------------------!

module fshift !used to apply a fermi-function derived energy shift to species 1

implicit none
save

logical :: fshiftlog
double precision :: fshiftmag
double precision :: fshiftlambda
double precision :: fshift_x0, fshift_x1

end module fshift

