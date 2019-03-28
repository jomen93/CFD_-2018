program CDS_NU
 implicit none
 integer, parameter :: dp = selected_real_kind(15)
 integer :: i,j,nodos
 real(dp), parameter :: l = 1.0, rho=1.0,u=1.0,gamma=0.01667
 real(dp) :: dx,Pe,d_x,r1,r2
 real(dp), allocatable :: phi(:),x(:)

 nodos = 11.0
 dx = l/nodos
 allocate(phi(nodos),x(nodos))
 phi(1) = 0.0_dp 
 x(1) = 0.0_dp
 x(2) = 0.31_dp
 Pe = rho*u*l/gamma

 !construccion de la grilla no uniforme 
 open(1, file = "CDS_NU.dat")
 do j = 1,6
  do i = 2,nodos-1
   x(i+1) = 0.7_dp*(x(i)-x(i-1)) + x(i)
   r1 = x(i+1)-x(i) ; r2 = x(i)-x(i-1)
   phi(i) = (r1*r2/(r1+r2))*((1.0_dp/r1 - 0.5_dp*Pe)*phi(i+1) + (1.0_dp/r2 + 0.5_dp*Pe)*phi(i-1))
  end do
 !x(nodos) = 0.0125
 phi(1) = 0.0_dp
 phi(nodos) = 1.0_dp
 end do 
 write(*,*) x!x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11)
 write(1,*) x
 write(1,*) phi

end program CDS_NU