program compacto2
implicit none

integer, parameter :: dp=selected_real_kind(15, 307)

integer :: n,i
real(kind=dp),allocatable:: y(:), dy2dx(:)
real(kind=dp) :: x, h, error
real(kind=dp), parameter :: pi=4.0D0*atan(1.0D0)
integer, parameter :: xmin=80.0d0
double precision :: k0,k1,k2

open(1,file = "c2.dat")

 k0=1.1250d0
 k1=4.D0*DATAN(1.D0)/(8.0*xmin)
 k2=(3.5d0*k1)/4.0

do n=xmin, 100000, 20 !BUCLE PARA EL CAMBIO DE TAMAÑO DE MALLA

allocate (y(n), dy2dx(n)) !asignar matrices de malla

h= (2.0d0*pi)/(n-1.0d0) !espacio de malla, asumiendo x de 0->2PI

do i=1,n
   x= -pi + h*(DBLE(i)- 1.0d0)
   y(i)=exp(-k0*(x**2))*(1.1d0*COS(k1*x)-1.3d0*SIN(k2*x)) !DEFINICIÓN DE FUNCIÓN f(x)
end do

call d2padx (y, n, h, dy2dx)

do i=1, n
if (i>20 .AND. i <(n-20)) then
x= -pi + h*(DBLE(i)- 1.0d0)
error= error + (((dy2dx(i) - exact2der(x))**2)*h)
end if
!write(*, "(2ES14.7)") dy2dx(i), exact2der(x)
end do
write(1,*) h, SQRT(error)/n !DATOS QUE IMPRIME
 
deallocate(y, dy2dx) !terminar

end do
close(1)
contains

subroutine d2padx(f, nx, h, s)
real(kind=dp), intent(in):: f(nx), h
real(kind=dp),intent(out) :: s(nx)
!real(kind=dp) :: x(nx), hx(nx), h2x(nx)
real(kind=dp) :: b0, c0, d(nx), m
integer :: i, nx

b0= 12.0/11.0
c0= 3.0/11.0
m= h**2

d(1)= b0*(2*f(2)-2*f(1))/m + 0.25*c0*(2*f(3)-2*f(1))/m
d(2)= b0*(f(3)-2*f(2)+f(1))/m + 0.25*c0*(f(4)-2*f(2)+f(2))/m
 
do i=3, nx-2
d(i)= b0*(f(i+1)-2*f(i)+f(i-1))/m + 0.25*c0*(f(i+2)-2*f(i)+f(i-2))/m
end do
d(nx)=d(1)
d(nx-1)=d(2)

call thomas(nx,d,s)

do i=1, nx
s(i)= s(i)/(2*pi)
end do
return
end

real(kind=dp) function exact2der(xp) !SEGUNDA DERIVADA EXACTA
  double precision, INTENT(IN) :: xp
  double precision, parameter :: pi=4.D0*DATAN(1.D0)
  double precision, parameter :: xmin=80.0d0
  double precision :: k0,k1,k2
    k0=1.1250d0
    k1=pi/(8.0*xmin)
    k2=(3.5d0*k1)/4.0 

    exact2der= 4.4*exp(-k0*(xp**2))*(((k0**2)*(xp**2)*COS(k1*xp)) -1.18182*(k0**2)*(xp**2)*SIN(k2*xp) &
+(k0*k1*xp)*SIN(k1*xp)-0.5*k0*COS(k1*xp)+0.590909*k0*SIN(k2*xp)+1.18182*k0*k2*xp*COS(k2*xp)-0.25*(k1**2)*COS(k1*xp) &
+0.295455*(k2**2)*SIN(k2*xp))
end function exact2der

subroutine thomas(n,d,s)

real(kind=dp), INTENT(IN):: d(n)
real(kind=dp) :: s(n) 
real(kind=dp) :: alpha(n), beta(n), y(n), a(n), b(n), c(n)
integer :: i, n

do i=1, n
b(i)= 1.0
a(i)= 2.0/11.0
c(i)= 2.0/11.0
end do

do i=1, n
if (n .EQ. 1) then
beta(1)= c(1)/b(1)
y(1)= d(1)/b(1)
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


end program compacto2
