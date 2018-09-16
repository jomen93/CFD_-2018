program DiffFin
 implicit none

 integer, parameter :: dp = selected_real_kind(15,307)
 integer :: i 
 real :: dt, rho, Cd, R,M,V,u,x,t,usig,w,Ao,Aref
 real(kind = dp), parameter :: pi = 4.0D0*atan(1.0D0)
 real :: dxdt, dudt, dx

 write(*,*) "Ingrese el radio, coeficiente de arrastre , la masa "
 read(*,*) R,Cd,M

 u = 0.0
 x = 0.0
 rho = 2.0
 dt = 0.00001
 V = 10.0
 t = 0
 Aref = 1.0
 w = 0.1
 Ao= 1.0


 open(10,file = "Resultado.dat")
 do i = 1, 100
  dxdt = u
  !dudt = ((rho*Cd*pi*R*R)/(2*M))*(V*V - 2*u*V * u*u) 
  dudt = ((rho*Cd*pi*R*R)/(2*M))*(Aref + Ao*sin(w*t)*Aref + Ao*sin(w*t) - 2*u*Aref + Ao*sin(w*t) * u*u) 
  u = u + t*dudt
  x = x + u*t
  t = t + dt
  write(10,*) t," ", u ," ",x," ", dudt
  write(*,*) t," ", x ," ",u," ", dudt
 end do 
 write(*,*) "Done!"

end program DiffFin