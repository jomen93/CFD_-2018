program derivada1
implicit none

integer, parameter :: dp=selected_real_kind(15, 307)

integer :: n,i
real(kind=dp),allocatable:: y(:), dydx(:)
real(kind=dp) :: x, dx, error
real(kind=dp), parameter :: pi=4.0D0*atan(1.0D0)
integer, parameter :: xmin=90.0d0
double precision :: k0,k1,k2

open(1,file = "d1_4to.dat")

 k0=1.1250d0
 k1=4.D0*DATAN(1.D0)/(8.0*xmin)
 k2=(3.5d0*k1)/4.0

do n=xmin, 10000, 20 !BUCLE PARA EL CAMBIO DE TAMAÑO DE MALLA

allocate (y(n), dydx(n)) !asignar matrices de malla

dx= (2.0d0*pi)/(n-1.0d0) !espacio de malla, asumiendo x de 0->2PI

do i=1,n
   x= -pi + dx*(DBLE(i)- 1.0d0)
   y(i)=exp(-k0*(x**2))*(1.1d0*COS(k1*x)-1.3d0*SIN(k2*x)) !DEFINICIÓN DE FUNCIÓN f(x)
end do

call derivative (y,n,dx,dydx) !calcula dydx

do i=1, n
if (i>20 .AND. i <(n-20)) then
x= -pi + dx*(DBLE(i)- 1.0d0)
error= error + ((dydx(i) - exactder(x))**2)*dx !DEFINICIÓN DE ERROR
end if
end do

write(1, *) dx, SQRT(error)/n  !DATOS QUE IMPRIME
 
deallocate(y, dydx) !terminar

end do
close(1)

contains

subroutine derivative (a,np, h, aprime) !DERIVADA NUMÉRICA
 real(kind=dp), intent(in):: a(np), h
 real(kind=dp), intent(out):: aprime(np)
 integer :: i, np !variable local

do i= 1, np-1
 !aprime(i)= (a(i+1)-a(i-1))/(2*h) !1ra derivada 2do Orden
 aprime(i)= (-a(i+2)+ 8*a(i+1)- 8*a(i-1)+a(i-2))/(12*h) !1ra derivada 4to Orden
end do
aprime(np)= aprime(1)  !trato de hacerlo periódico

end subroutine derivative

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

end program derivada1
