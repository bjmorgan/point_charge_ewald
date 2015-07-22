subroutine fshiftener

use commondata, only:  nsp, x, frrx, engpetot
use fshift

implicit none

integer :: i
double precision :: xf, exp1, exp2
double precision :: fdenom1, fdenom2, fshift1, fshift2, eng1, eng2
double precision :: fshifteng, fshiftforce
double precision :: fshifteng_tot

fshifteng_tot = 0.0d0

do i=1, nsp(1) ! apply shifted energy potential to species 1
    ! at the moment only applied to x direction
    xf = x(i)
    exp1 = exp((xf - fshift_x0) / fshiftlambda)
    exp2 = exp((fshift_x1 - xf) / fshiftlambda)
    fdenom1 = 1.0 + exp1
    fdenom2 = 1.0 + exp2

    if (abs(exp1) > huge(exp1)) then
        eng1 = 0.0d0
        fshift1 = 0.0d0
    else
        eng1 = 1.0 / fdenom1
        fshift1 = exp1 / (fdenom1 * fdenom1)
    endif

    if (abs(exp2) > huge(exp2)) then
        eng2 = 0.0d0
        fshift2 = 0.0d0
    else
        eng2 = 1.0 / fdenom2
        fshift2 = exp2 / (fdenom2 * fdenom2)
    endif

    fshifteng = fshiftmag * (eng1 + eng2)
    fshiftforce = fshiftmag / fshiftlambda * (fshift1 - fshift2)

    engpetot = engpetot + fshifteng
    fshifteng_tot = fshifteng_tot + fshifteng
    frrx(i) = frrx(i) + fshiftforce
enddo

return
end subroutine
