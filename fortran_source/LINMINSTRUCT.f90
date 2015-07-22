SUBROUTINE LINMINSTRUCT(P,XI,FRET)

USE commondata, ONLY: num,tolstruc
USE cgdata, ONLY: PCOMSTRUCT,XICOMSTRUCT

IMPLICIT NONE

INTEGER :: N,J
DOUBLE PRECISION, INTENT(OUT) :: FRET
DOUBLE PRECISION, DIMENSION(3*num+6), INTENT(INOUT) :: P,XI
DOUBLE PRECISION :: AXX,BXX,XX,FA,FX,FB,XMIN,DBRENT,F1DIMSTRUCT

EXTERNAL F1DIMSTRUCT,DF1DIMSTRUCT

pcomstruct = p
xicomstruct = xi
AXX=0.00D0
XX=1.00D0
BXX=2.00D0
CALL MNBRAKSTRUCT(AXX,XX,BXX,FA,FX,FB)
FRET=DBRENT(AXX,XX,BXX,F1DIMSTRUCT,DF1DIMSTRUCT,TOLSTRUC,XMIN)
! FRET=BRENT(AXX,XX,BXX,F1DIMSTRUCT,TOLSTRUC,XMIN)

XI=XMIN*XI
P=P+XI

RETURN
END SUBROUTINE
