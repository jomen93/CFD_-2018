program aleta_trasitoria
 implicit none
 integer, parameter :: dp = selected_real_kind(15)
 integer :: i,j,n
 real(dp) :: thetao, r1,r2,dx,dt,A,l
 real(dp), allocatable :: theta(:), t(:)

 A = 100.0_dp
 l = 1.0_dp
 dx = 0.01_dp
 dt = ((2*dx*dx)/(A*dx*dx+4))*0.1_dp
 n = l/dx
 r1 = dt/(dx*dx)
 r2 = 1.0_dp-A*dt -2.0_dp*r1
 allocate(theta(n),t(n))

 thetao = 1.0_dp
 theta(1) = thetao
    
open(1,file = "theta_transitorio.dat") 
 do j = 1,1000000
  !aux = theta
  do i = 2,n-1
   theta(i) = r1*theta(i-1) + r2*theta(i) +r1*theta(i+1)
  end do 
  theta(n) = theta(n-1)
  print*, theta(n/4),theta(n/2),theta(3*n/4)
 end do
 write(1,*) theta
end program aleta_trasitoria