program main
 implicit none 
 integer, parameter :: dp = selected_real_kind(15)
 real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
 integer :: i,j,n,t,tf,mu
 real(dp) :: nu,dt,dx,l,r,sigma,Re
 real(dp), allocatable :: u(:)
 
 dx = 0.001!l/n
 dt = 0.000001
 t = 0.0; tf = 80000
 nu = 1.0 
 l =1.0
 r = nu*dt/(dx*dx)
 n=int(l/dx)
 mu = n/5
 Re = 0.01
 sigma = 100

 allocate(u(n))
 !u = 0.0_dp
 do i= 1,n
 !u(i) = 1.0*(1/(sigma*sqrt(2*pi)))*exp(-0.5*((i-mu)/(sigma))**2) 
 u(i) = 0.25*(3.0 - 1.0/(1+exp(Re*real(-4*i)/32)))
end do
 !(4.0_dp*u(n+1) - u(n+2))/3.0_dp

 open(1,file = "FV.dat")
 
 do j=1,tf
 u(1) = (4.0_dp*u(n+1) + u(n+2))/3.0_dp
  do i = 2,n-1
   u(i) = u(i) + 0.5_dp*nu*(u(i+1)-2.0_dp*u(i)+u(i-1))*dt/(dx**2) - 0.125_dp*(dt/dx)*((u(i+1)+u(i))**2 - (u(i-1)+u(i))**2)
  end do
 u(n) = (4.0_dp*u(n-1) - u(n-2))/3.0_dp

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

 !write(1,*) u
 !close(1)

end program main

