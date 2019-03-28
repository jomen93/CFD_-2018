program compacto
implicit none

integer, parameter :: dp=selected_real_kind(15, 307)

integer :: n,i
real(kind=dp),allocatable:: y(:), dydx(:), dy2dx(:)
real(kind=dp) :: x, h, error
real(kind=dp), parameter :: pi=4.0D0*atan(1.0D0)
integer, parameter :: xmin=90.0d0
double precision :: k0,k1,k2, x1

open(1, file = "c1.dat")

 k0=1.1250d0
 k1=4.D0*DATAN(1.D0)/(8.0*xmin)
 k2=(3.5d0*k1)/4.0

do n=xmin, 10000, 20 !BUCLE PARA EL CAMBIO DE TAMAÑO DE MALLA

allocate (y(n), dydx(n)) !asignar matrices de malla

h= (2.0d0*pi)/(n-1.0d0) !espacio de malla, asumiendo x de 0->2PI

do i=1,n
   x= -pi + h*(DBLE(i)- 1.0d0)
   y(i)=exp(-k0*(x**2))*(1.1d0*COS(k1*x)-1.3d0*SIN(k2*x)) !DEFINICIÓN DE FUNCIÓN f(x)
end do

call dpadex (y, n, h, dydx) !calcula dydx


do i=2, n
if (i>20 .AND. i <(n-20)) then
x= -pi + h*(DBLE(i)- 1.0d0)
error= error + (((dydx(i) - exactder(x))**2)*h)       !   DEFINICIÓN DE ERROR
end if
!write(*, "(2ES14.7)") dydx(i), exactder(x)
end do
write(1, *) h, SQRT(error)/n  !DATOS QUE IMPRIME
 
deallocate(y, dydx) !terminar

end do

close(1)
contains

subroutine dpadex(f, nx, h, s)
real(kind=dp), intent(in):: f(nx), h
real(kind=dp),intent(out) :: s(nx)
real(kind=dp) :: hx(nx)
real(kind=dp) :: b0, c0, d(nx)
integer :: i, nx

b0= 14.0/9.0
c0= 1.0/9.0

do i=1,nx
   hx(i)= pi
end do

d(1)= 0.0
d(2)= 0.5* b0*(f(3)-f(1))/h + 0.25*c0*(f(4)-f(2))/h

do i=3, nx-2
d(i)= 0.5* b0*(f(i+1)-f(i-1))/h + 0.25*c0*(f(i+2)-f(i-2))/h
end do
d(nx)=d(1)
d(nx-1)=d(2)

call thomas(nx,d,s)

do i=1, nx
s(i)= s(i)/hx(i)
end do
return
end

double precision function exactder(xp) !DERIVADA EXACTA
  double precision, INTENT(IN) :: xp
  double precision, parameter :: pi=4.D0*DATAN(1.D0)
  double precision, parameter :: xmin=90.0d0
  double precision :: k0,k1,k2
    k0=1.1250d0
    k1=pi/(8.0*xmin)
    k2=(3.5d0*k1)/4.0 

    exactder=exp(-k0*(xp**2))*(-2.2d0*k0*xp*COS(k1*xp)-1.1d0*k1*SIN(k1*xp)+2.6d0*k0*xp*SIN(k2*xp)-1.3d0*k2*COS(k2*xp))
end function exactder


subroutine thomas(n,d,s)

real(kind=dp), INTENT(IN):: d(n)
real(kind=dp) :: s(n) 
real(kind=dp) :: alpha(n), beta(n), y(n), a(n), b(n), c(n)
integer :: i, n

do i=1, n
b(i)= 1.0
a(i)= 1.0/3.0
c(i)= 1.0/3.0
end do

do i=1, n
if (n .EQ. 1) then
beta(1)= c(1)/b(1)
y(1)= 0.0
return
end if
alpha(i) = b(i) - a(i) * beta(i-1)
beta(i) = c(i)/ alpha(i)
y(i) = (d(i) - a(i)*y(i-1))/alpha(i)
end do
do i=1, n-1
s(i) = y(i) - beta(i)*s(i+1)
end do
s(n)= y(n)
return
end


end program compacto
