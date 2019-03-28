
!*******************************************************************************
program main

 implicit none

 integer, parameter :: dp =  selected_real_kind(15)
 integer :: i,j=0
 integer :: n
 real(dp) :: thetao
 real(dp) :: Beta,dx,l, long
 real(dp), allocatable :: theta(:), t(:),aux(:)

 print*, "Ingrese el dx del sistema"
 read*, dx
 print*, "Ingrese el numero caracteristico del sistema l: "
 read*, l

 n = 2**7_dp
 allocate(theta(n), t(n),aux(n))
 thetao = 1.0_dp 
 Beta = 2.0_dp + dx*dx*l
 theta(1) = thetao
 
 
 open(1,file = "theta.dat") 
 do j = 1,1000000
  aux = theta
  do i = 2,n-1
   theta(i) = (1.0_dp/Beta)*(theta(i+1)+theta(i-1))
  end do 
  theta(n) = (2.0_dp/Beta)*theta(n-1)
 end do
 write(1,*) theta

end program main

!***************************************************








