module pentadiagonal  
use precision, only: dp
use variables 
public :: nind


contains
 
 integer function nind(ii,jj,nn)
  implicit none 
  integer,     intent(in) :: ii,jj
  integer(dp), intent(in) :: nn
  nind = jj + (ii-1)*nn
 end function nind



 subroutine Rellenar(uu,vv,nn,u)
  implicit none 
  save
  real(dp),     intent(in), dimension(0:nn,0:nn) :: uu,vv
  integer(dp),  intent(in) :: nn
  integer(dp)              :: kk
  real(dp),     intent(out),dimension(0:nn,0:nn) :: u
  

  do i = 1,n-1
   do j = 1,n-1
    b1(i,j) = -(1/dy)*(0.25_dp*(v(i,j-1)+v(i,j))+1.0_dp/(Re*dy))
    b2(i,j) = -(1/dx)*(0.25_dp*(u(i-1,j)+u(i,j))+1.0_dp/(Re*dx))
    b3(i,j) =  (1/dx)*(0.25_dp*(u(i+1,j)-u(i-1,j))+2.0_dp/(Re*dx))
    b4(i,j) =  (1/dy)*(0.25_dp*(v(i,j+1)-v(i,j-1))+2.0_dp/(Re*dy))
    b5(i,j) =  (1/dx)*(0.25_dp*(u(i,j)+v(i+1,j))-1.0_dp/(Re*dx))
    b6(i,j) =  (1/dy)*(0.25_dp*(v(i,j)+v(i,j+1))-1.0_dp/(Re*dy))
    !termino fuente conocido 
    s(i,j) = -b2(i,j)*v(i-1,j) +(0.5_dp*dt -b3(i,j))*v(i,j) -b4(i,j)*v(i+1,j)
    print*, "1"
   end do
  end do
  
  !direccion y
  do i= 1,n-1
   do j = 1, n-1
    a(j) = b1(i,j)
    b(j) = (b4(i,j)-0.5_dp*dt)
    c(j) = b5(i,j)
    r(j) = s(i,j)
    print*, "2"
   end do
  

   ! Aca se considera la conicion de frontera, esto falta   
   print*, "estoy en thomas"
   call thomas(a,b,c,r,n,q)

   do j=1,n-1
    z(i,j) = q(j)
   end do

  end do 

  !direccion x
  do j = 1,n-1
   do i = 1,n-1
    a(i) = b2(i,j) 
    b(i) = b4(i,j)+0.5_dp*dt
    c(i) = b4(i,j)
    r(i) = z(i,j)
   end do 

   ! Aca para la otra direccion las condiciones de frontera

   call thomas(a,b,c,r,n,q)

   do i = 1, n-1 
    u(i,j) = q(i)
   end do 
  end do
  return
 end subroutine

 subroutine thomas(a,b,c,r,n,x)
  implicit none
  real(dp),     dimension(1:n-1) :: a,b,c,r
  integer(dp),  intent(in) :: n
  real(dp),     intent(out),dimension(1:n-1) :: x
  
  do i=1,n 
   b(i) = b(i) - a(i)/b(i-1)*c(i-1)
   r(i) = r(i) - a(i)/b(i-1)*r(i-1)
  end do
  ! backward substitution phase 
  x(n) = r(n)/b(n)
  do i=n-1,1,-1
   x(i) = (r(i)-c(i)*x(i+1))/b(i)
  end do
 return 
 end subroutine

end module pentadiagonal