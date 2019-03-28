program main      

use pentadiagonal
use precision, only: dp 
use variables

implicit none 
 
open(1, file = "M.dat")
!Iniciar variables


do i = 0,n
 x(i) = xo + dfloat(i)*dx; y(i) =yo + dfloat(i)*dy
 do j =0,n
  u(i,j) = uo(x(i),y(j),0.0_dp)
  v(i,j) = vo(x(i),y(j),0.0_dp)
 end do
end do 

s=0.0_dp

call Rellenar(u,v,n,up)
!call thomas(M,uord,n,zres)
print*, up
 

contains
real function uo(x,y,t)
 implicit none 
 real(dp), intent(in) :: x,y,t
 uo = 0.25_dp*(3.0_dp - 1.0_dp/(1.0_dp+exp(-0.03125_dp*Re*(t+4.0_dp*x+4.0_dp*y))))
end function uo

real function vo(x,y,t)
 implicit none 
 real(dp), intent(in) :: x,y,t
 vo = 0.25_dp*(3.0_dp + 1.0_dp/(1.0_dp+exp(-0.03125_dp*Re*(t+4.0_dp*x+4.0_dp*y))))
end function vo


 
end program main