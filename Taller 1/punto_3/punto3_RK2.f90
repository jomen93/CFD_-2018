program main
 implicit none
 integer, parameter :: dp = selected_real_kind(15)

 integer :: i,n,j,k,tf
 real(dp), parameter :: l=1.0,rho=1.0,u=1.0,gamma=0.025
 real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
 real(dp) :: dx,dt,t,k1,k2,a1,a2,a3,a4,l1,l2,r,Pe,h
 real(dp), allocatable :: phi(:),x(:),phi1(:)
 
 !read*, n,name # 
 !Nombra archivos para cuando se utilzia varios nodos
 t = 0.0;tf = 50000
 dx = 1/41.0!0.001
 n = int(l/dx)
 Pe = rho*u/gamma
 r=0.9
 allocate(phi(n),x(n),phi1(n))
 
 
 !call malla(n-1,l,r,x)
 h = (l*(1-r))/(1-r**(n-1))
 x(1) = 0.0
 x(2) = h
 do i=2, n-1
  x(i+1) = h + r*x(i)
 enddo

 
 do k = 1,n
  phi(k) = 1.0*sin(k*0.5*pi/n)
  phi1(k) = 1.0*sin(x(k)*0.5*pi)
 end do 
 phi(1) =0.0
 phi1(1) = 0.0
 dt = 0.00001

 open(2,file = "DCT_RK2_t.dat")
 open(1,file = "DCT_RK2.dat")
 open(4,file = "DCT_RK2_NU.dat")
 open(3,file = "malla_RK2.dat")
 !open(2,file = "DCT_RK2_t_"//name//".dat")
 !open(1,file = "DCT_RK2_"//name//".dat")
 !open(4,file = "DCT_RK2_NU_"//name//".dat")
 !open(3,file = "malla_RK2_"//name//".dat")

 do j = 1, tf
  do i = 2,n-1
   

   k1 = (0.5*(-phi(i+1)+phi(i-1))/dx) + (1.0/Pe)*( (phi(i+1)-2.0*phi(i)+phi(i-1) )/(dx*dx) )
   k2 = (0.5*(-phi(i+1)+phi(i-1))/dx) + (1.0/Pe)*( (phi(i+1)-2.0*phi(i)+phi(i-1) )/(dx*dx) ) + 0.5*dt*k1
   phi(i) = phi(i) + dt*(k2+k1)

   a1 = (2.0*gamma)/( rho*(x(i+1)-x(i))*(x(i)-x(i-1))  )
   a2 = (2.0*gamma)/( rho*(x(i+1)-x(i-1))*(x(i+1)-x(i)) )
   a3 = u/(x(i+1)-x(i-1))
   a4 = (2.0*gamma)/(rho*(x(i+1)-x(i-1))*(x(i)-x(i-1)) )

   l1 = (a3+a4)*phi1(i-1) -a1*phi1(i) + (a2-a3)*phi1(i+1)
   l2 = (a3+a4)*phi1(i-1) -a1*phi1(i) + (a2-a3)*phi1(i+1) + 0.5*dt*l1
   phi1(i) = phi1(i) + dt*(l1+l2)

  end do
  
  phi(n) = 1.0
  phi1(n) = 1.0

  if (j == 1) then
   write(1,*) phi
   write(2,*) t
   write(4,*) phi1
  endif
  if (j == 1*tf/5) then
   write(1,*) phi
   write(2,*) t
   write(4,*) phi1
  endif
  if (j == 2*tf/5 ) then
   write(1,*) phi
   write(2,*) t
   write(4,*) phi1
  endif
  if (j == 3*tf/5 ) then
   write(1,*) phi
   write(2,*) t
   write(4,*) phi1
  endif
  if (j == 4*tf/5) then
   write(1,*) phi
   write(2,*) t
   write(4,*) phi1
  endif
  if (j == tf ) then
   write(1,*) phi
   write(2,*) t
   write(4,*) phi1
  endif
 end do 
 write(3,*) x
 write(*,*) "Done RK2!" 
 close(1)
 close(2)
 close(3)
 close(4)
 deallocate(phi)

end program main





