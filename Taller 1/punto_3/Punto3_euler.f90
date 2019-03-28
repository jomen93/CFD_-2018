
program main
 implicit none
 
 integer, parameter :: dp = selected_real_kind(15)
 integer :: i,n,j,k,tf
 real(dp), parameter :: l=1.0,rho=1.0,u=1.0,gamma=0.025
 real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
 real(dp) :: dx,dt,d,c,t,r,h
 REAL(dp) :: a1,a2,a3,a4
 real(dp), allocatable :: phi(:),phi1(:),x(:),phi2(:)
 
 t = 0.0;tf = 1000000
 dx = 1/41.0!0.001
 n = int(l/dx); r = 0.99
 allocate(phi(n),phi1(n),x(n),phi2(n))
 
 h = (l*(1-r))/(1-r**(n-1))
 x(1) = 0.0
 x(2) = h
 
 do i=2, n-1
  x(i+1) = h + r*x(i)
 enddo

 do k = 1,n
  phi(k) = 1.0*sin(k*0.5*pi/n)
  phi2(k) = 1.0*sin(x(k)*0.5*pi)
 end do 
 phi(1) = 0.0
 phi1(1) = 0.0
 
 dt = 0.000001
 d = gamma*dt/(rho*dx*dx)
 c = u*dt/dx
 
 open(1,file = "DCT_CDS.dat")
 open(2,file = "DCT_CDS_t.dat")
 open(4,file = "DCT_CDS_NU.dat")
 open(5,file = "malla_CDS.dat")

 if (d < 0.5 .and. c < 2*d) then
 do j = 1, tf
  do i = 2,n-1
   a1 = (2.0*gamma*dt)/( rho*(x(i+1)-x(i))*(x(i)-x(i-1))  )
   a2 = (2.0*gamma*dt)/( rho*(x(i+1)-x(i-1))*(x(i+1)-x(i)) )
   a3 = u*dt/(x(i+1)-x(i-1))
   a4 = (2.0*gamma*dt)/(rho*(x(i+1)-x(i-1))*(x(i)-x(i-1)) )
   phi2(i) = (1-a1)*phi(i) + (a2-a3)*phi(i+1) + (a4+a3)*phi(i-1)

   phi(i)= (1.0-2.0*d)*phi(i)+(d-0.5*c)*phi(i+1)+(d+0.5*c)*phi(i-1)
  end do
   
  phi(n)= 1.0
  phi2(n) = 1.0


  if (j == 1) then
   write(1,*) phi
   write(2,*) t
   write(4,*) phi2
  endif
  if (j == 1*tf/5) then
   write(1,*) phi
   write(2,*) t
   write(4,*) phi2
  endif
  if (j == 2*tf/5 ) then
   write(1,*) phi
   write(2,*) t
   write(4,*) phi2
  endif
  if (j == 3*tf/5 ) then
   write(1,*) phi
   write(2,*) t
   write(4,*) phi2
  endif
  if (j == 4*tf/5) then
   write(1,*) phi
   write(2,*) t
   write(4,*) phi2
  endif
  if (j == tf ) then
   write(1,*) phi
   write(2,*) t
   write(4,*) phi2
  endif
 end do 

 write(*,*) "Done Euler!"
 else
  print*, "Para que sea estable el programa se debe asegurar"
  print*, "valores adecuados para dt y dx "
  print*, "d = ",d,"c = ", c
 endif

 write(5,*) x
 close(1)
 close(2)
 close(4)
 close(5)
 deallocate(phi,phi1,phi2)

end program main
