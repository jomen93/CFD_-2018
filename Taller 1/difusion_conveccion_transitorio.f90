program main
 implicit none
 integer, parameter :: dp = selected_real_kind(15)

 integer :: i,n,j,k,tf
 real(dp), parameter :: l=1.0,rho=1.0,u=1.0,gamma=0.025
 real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
 real(dp) :: dx,dt,beta,Pe,d,c,t,t1,k1,k2,r
 REAL(dp) :: a1,a2,a3,a4
 real(dp), allocatable :: phi(:),phi1(:),x(:),phi2(:)
 
 t = 0.0;tf = 1000000
 t1 =0.0 
 dx = 1/31.0!0.001
 n = l/dx; r = 0.98
 allocate(phi(n),phi1(n),x(n),phi2(n))
 x(1) = 0.0
 
 call malla(n-1,l,r,x)

 do k = 1,n
  phi(k) = 1.0*sin(k*0.5*pi/n)
  phi1(k) = 1.0*sin(k*0.5*pi/n)
  phi2(k) = 1.0*sin(x(k)*0.5*pi)
 end do 
 
 print*, (x(n)-x(n-1))*(x(n-1)-x(n-2))
 print*, (x(n)-x(n-1))*(x(n)-x(n-2))
 print*, (x(n)-x(n-2))*(x(n-1)-x(n-2))

 Pe = rho*u*l/gamma
 beta = l/(u)
 dt = 0.000001
 d = gamma*dt/(rho*dx*dx)
 c = u*dt/dx
 print*, "Condiciones de frontera: ",phi(1),phi(n)



 ! ---------------- Forward Euler -------------------------
 open(1,file = "DCT_CDS.dat")
 open(2,file = "DCT_CDS_t.dat")
 open(3,file = "DCT_RK2.dat")
 open(4,file = "DCT_RK2_NU_E.dat")
 open(5,file = "malla.dat")

 if (d < 0.5 .and. c < 2*d) then
 do j = 1, tf
  do i = 2,n-1
   a1 = (2.0*gamma*dt)/( rho*(x(i+1)-x(i))*(x(i)-x(i-1))  )
   a2 = (2.0*gamma*dt)/( rho*(x(i+1)-x(i-1))*(x(i+1)-x(i)) )
   a3 = u*dt/(x(i+1)-x(i-1))
   a4 = (2.0*gamma*dt)/(rho*(x(i+1)-x(i-1))*(x(i)-x(i-1)) )
   phi2(i) = (1-a1)*phi(i) + (a2-a3)*phi(i+1) + (a4+a3)*phi(i-1)


   phi(i)= (1.0-2.0*d)*phi(i)+(d-0.5*c)*phi(i+1)+(d+0.5*c)*phi(i-1)
   k1 = (0.5*(-phi(i+1)+phi(i-1))/dx) + (1.0/Pe)*( (phi(i+1)-2.0*phi(i)+phi(i-1) )/(dx*dx) )
   k2 = (0.5*(-phi(i+1)+phi(i-1))/dx) + (1.0/Pe)*( (phi(i+1)-2.0*phi(i)+phi(i-1) )/(dx*dx) ) + 0.5*dt*k1
   phi1(i) = phi1(i) + dt*(k2+k1) 
  end do
   
  
  phi(n)= 1.0
  phi1(n) = 1.0
  phi2(n) = 1.0

  !print*, phi2(j)
  !t = t + dt
  !t1 = t1 + dt1
  if (j == 1) then
   write(1,*) phi
   write(2,*) t
   write(3,*) phi1
   write(4,*) phi2
  endif
  if (j == 1*tf/5) then
   write(1,*) phi
   write(2,*) t
   write(3,*) phi1
   write(4,*) phi2
  endif
  if (j == 2*tf/5 ) then
   write(1,*) phi
   write(2,*) t
   write(3,*) phi1
   write(4,*) phi2
  endif
  if (j == 3*tf/5 ) then
   write(1,*) phi
   write(2,*) t
   write(3,*) phi1
   write(4,*) phi2
  endif
  if (j == 4*tf/5) then
   write(1,*) phi
   write(2,*) t
   write(3,*) phi1
   write(4,*) phi2
  endif
  if (j == tf ) then
   write(1,*) phi
   write(2,*) t
   write(3,*) phi1
   write(4,*) phi2
  endif

  
 end do 
 !write(*,*) t,phi


 write(*,*) "Done Euler!"
 else
  print*, "Para que sea estable el programa se debe asegurar"
  print*, "valores adecuados para dt y dx "
  print*, "d = ",d,"c = ", c
 endif

 write(5,*) x
 close(1)
 close(2)
 close(3)
 close(4)
 close(5)
 deallocate(phi,phi1,phi2)



 contains 
 subroutine malla(n,L,r,x)
 implicit none
 integer , intent(in) :: n
 real(dp), intent(in) :: L,r
 real(dp), dimension(n+1), intent(out) :: x
 real(dp) :: s, dx
 if(abs(r-1)<1e-9) then
  s = n
 else
  s = (1-(r**n))/(1-r)
 endif

 dx = L/s
 x(1) = 0.0

 do i=2, n+1
  x(i) = x(i-1)+dx
  dx = r*dx
 end do 
 end subroutine malla

end program main
