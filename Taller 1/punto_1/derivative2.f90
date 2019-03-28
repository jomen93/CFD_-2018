program derivada2
implicit none

integer, parameter :: dp=selected_real_kind(15, 307)

integer :: n,i
real(kind=dp),allocatable:: y(:), dy2dx(:)
real(kind=dp) :: x, dx, error
real(kind=dp), parameter :: pi=4.0D0*atan(1.0D0)
integer, parameter :: xmin=80.0d0
double precision :: k0,k1,k2

open(1,file = "d2_6to.dat")

 k0=1.1250d0
 k1=4.D0*DATAN(1.D0)/(8.0*xmin)
 k2=(3.5d0*k1)/4.0

do n=xmin, 100000, 20 !BUCLE TAMAÑO DE MALLA

allocate (y(n), dy2dx(n)) !asignar matrices de malla

dx= (2.0d0*pi)/(n-1.0d0) !espacio de malla, asumiendo x de 0->2PI

do i=1,n
   x= -pi + dx*(DBLE(i)- 1.0d0)
   y(i)=exp(-k0*(x**2))*(1.1d0*COS(k1*x)-1.3d0*SIN(k2*x))
end do

call derivative (y,n,dx,dy2dx) !calcula dydx2

do i=1, n
if (i>20 .AND. i <(n-20)) then
x= -pi + dx*(DBLE(i)- 1.0d0)
error= error + (((dy2dx(i) - exact2der(x))**2)*dx)
end if
!write(*, "(2ES14.7)") dy2dx(i), exact2der(x)
end do

write(1,*) dx, SQRT(error)/n


deallocate(y, dy2dx) !terminar

end do

close(1)
contains

subroutine derivative (a,np, h, aprime) !SEGUNDA DERIVADA NUMÉRICA
 real(kind=dp), intent(in):: a(np), h
 real(kind=dp), intent(out):: aprime(np)
 integer :: i, np !variable local

 do i= 3, np-2
 !aprime(i)= (a(i+1)- 2*a(i)+a(i-1))/(h*h) !2da derivada 2do Orden
 aprime(i)= (-a(i+2)+16*a(i+1)-30*a(i)+16*a(i-1)-a(i-2))/(12*(h*h)) !2da derivada 4to Orden
 end do
 aprime(np)= aprime(1)  !trato de hacerlo periódico

end subroutine derivative

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

end program derivada2
