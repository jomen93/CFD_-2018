program B1D
 implicit none
 integer, parameter :: dp = selected_real_kind(15)
 real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
 real(dp) alpha, beta,dx,dt,t,nu,l,sigma
 real(dp), allocatable :: u(:)
 integer(dp) i,n,j,tf,mu
 
 dt = 0.000001
 dx = 0.001
 nu = 1.0
 t = 0.0; tf =1000
 nu = 1.0
 l = 1.0

 alpha = nu*dt/(dx*dx)
 beta = 0.5*dt/dx
 n = int(l/dx)
 mu = n/2
 sigma = 100

allocate(u(n))
u = 0.0_dp
do i= 1,n
 u(i) = 1000.0*(1/(sigma*sqrt(2*pi)))*exp(-0.5*((i-mu)/(sigma))**2) 
end do


open(1, file = "B1D.dat")

do j=2,tf
 u(1) = 0
 do i=2,n-1
  u(i) = u(i)*(1-2*alpha-beta*(u(i+1)-u(i-1))) + alpha*(u(i+1)+u(i-1))
 end do
 t = t + dt
 u(n) = (4.0_dp*u(n-1) - u(n-2))/3.0_dp
 print*, u(2),u(n/2),u(n-1)

 if (j == 1) then
 write(1,*) u 
 endif  
 if (j == 1*tf/5) then
 write(1,*) u
 endif  
 if (j == 2*tf/5) then
 write(1,*) u
 endif  
 if (j == 3*tf/5) then
 write(1,*) u
 endif  
 if (j == 4*tf/5) then
 write(1,*) u
 endif  
 if (j == tf) then
 write(1,*) u
 endif  

end do 

close(1)
    
end program B1D