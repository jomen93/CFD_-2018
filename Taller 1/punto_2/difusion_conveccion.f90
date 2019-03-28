program main
 implicit none

 integer, parameter :: dp = selected_real_kind(15)
 integer :: i,n,j,k,e,nodos,numnodos
 real(dp), parameter :: l=1.0,rho=1.0,u=1.0,gamma=0.01667
 real(dp) :: dx,beta,b,Pe,r1,r2
 real(dp), allocatable :: phi(:),x(:),phi1(:)


 print*, "Upwind ------------------------------------>1"
 print*, "Diferencias centradas no uniformes -------->2"
 print*, "Diferencias centradas uniformes ----------->3"
 read*, nodos, e
 
 dx = 1.0/nodos
 !dx = 1.0/101.0
 n = l/dx
 !nodos = 11.0
 allocate(phi(n),phi1(nodos),x(nodos))
 phi(1) = 0.0_dp
 phi1(1) = 0.0_dp
 x(1) = 0.0_dp
 x(2) = 0.31_dp
 Pe = rho*u*l/gamma
 b = 2*gamma/(rho*u*l)
 beta =b/dx
 

 select case(e)

 case(1)

  open(1, file  = "DC_UDS.dat")
  do j = 1,1000000
   do i = 2,n-1
    phi(i) = (2.0_dp/(dx*Pe+4.0_dp))*phi(i+1) + ((dx*Pe+2.0_dp)/(dx*Pe+4.0_dp))*phi(i-1)
    end do 
   phi(n) = 1.0_dp
  end do 
  write(*,*) phi
  write(1,*) phi
 



 case(2)
  open(2, file = "CDS_NU.dat")
  do j = 1,10
   do i = 2,nodos-1
    x(i+1) = 0.7_dp*(x(i)-x(i-1)) + x(i)
    r1 = x(i+1)-x(i) ; r2 = x(i)-x(i-1)
    phi1(i) = (r1*r2/(r1+r2))*((1.0_dp/r1 - 0.5_dp*Pe)*phi1(i+1) + (1.0_dp/r2 + 0.5_dp*Pe)*phi1(i-1))
    write(2,*) x(i),phi1(i)
   end do
  phi1(1) = 0.0_dp
  phi1(nodos) = 1.0_dp
  end do 
  write(*,*) phi1
  !write(2,*) x,phi1
  !write(2,*) phi1
 
 case(3)
  open(1, file  = "DC_CDS.dat")
  do j = 1,1000000
   do i = 2,n-1
    phi(i) = ( (beta+1.0_dp)*phi(i-1)+(beta-1.0_dp)*phi(i+1) )/(2.0_dp*beta)
    end do
 
   phi(n) = 1.0_dp
  end do 
  write(*,*) phi
  write(1,*) phi


 end select



end program main


! Verificar el algoritmo 