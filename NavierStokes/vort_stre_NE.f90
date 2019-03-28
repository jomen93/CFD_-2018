program main
 implicit none
 integer, parameter :: dp = selected_real_kind(15)
 integer, parameter :: NumFile_u = 10, NumFile_v =11
 integer :: i,j,n,k,kend,m
 character(200) FileNameu,FileNamev
 real(dp), parameter :: dt= 0.0001
 real(dp) :: t1,t2,t3,a1,a2,a3,a4,a5,nu, vo
 real(dp) :: c1,c2,c3,c4,c5,d,rho,dx,dy,l
 real(dp), allocatable :: u(:,:), v(:,:),w(:,:),psi(:,:),p(:,:),x(:),y(:)

 n = 128
 l =10.0
 dx = l/n
 dy = l/n
 nu = 1.0_dp
 vo =1.0_dp
 rho = 1.0_dp
 d = dx**2 + dy**2
 
 allocate(u(0:n,0:n),v(0:n,0:n),w(0:n,0:n),psi(0:n,0:n),p(0:n,0:n),x(0:n),y(0:n))

 open(1, file = "u.dat")
 open(2, file = "v.dat")
 open(3, file = "n.dat")
 !open(4, file = "w.dat")
 open(50, file = "p.dat")
 !open(6, file = "psi.dat")
 open(7, file = "x.dat")
 open(8, file = "y.dat")


 
 do i = 1,n
   x(i) = i
   y(i) = i
 end do

 u = 0.0_dp
 v = 0.0_dp 
 p = 0.0_dp

read(*,"(1i10)") kend 

do k = 1,kend
 write(*,fmt="(a3,a,t21,f6.2,a)",advance="no") achar(13), "progreso completo:    ", (real(k)/real(kend))*100.0," %"

! condiciones de frontera

 do i = 2,n-1
 !lado izquierdo
  
  psi(1,i) = vo
  w(1,i) = 0.0_dp
  p(1,i) = (-4.0_dp*p(2,i) + p(3,i))/3.0_dp
  
 !lado derecho
  
  psi(n,i) = (4.0_dp*psi(n-1,i)-psi(n-2,i))/3.0_dp
  w(n,i) = (4.0_dp*w(n-1,i)-w(n-2,i))/3.0_dp
  p(n,i) = (4.0_dp*p(n-1,i)-p(n-2,i))/3.0_dp

 
 !Arriba
  
  psi(i,n) = (4.0_dp*psi(i,n-1)-psi(i,n-2))/3.0_dp
  w(i,n) = 0.0_dp
  p(i,n) = (4.0_dp*p(i,n-1)-p(i,n-2))/3.0_dp

  
 !Abajo
  
  psi(i,1) = (4.0_dp*psi(i,2)-psi(i,3))/3.0_dp
  w(i,1) = 2.0_dp*(psi(i,1)-psi(i,2))/dy**2
  p(i,1) = (4.0_dp*p(i,2) - p(i,3))/3.0_dp

  

 !esquinas

  w(1,1) = (w(1,2)+w(2,1)+w(2,2))*0.333_dp
  w(1,n) = (w(2,n)+w(1,n-1)+w(2,n-1))*0.333_dp
  w(n,1) = (w(n,2)+w(n-1,1)+w(n-1,2))*0.333_dp
  w(n,n) = (w(n-1,n)+w(n-1,n-1)+w(n,n-1))
  p(1,1) = 0.0!(p(1,2)+p(2,1)+p(2,2))*0.333_dp
  p(1,n) = (p(2,n)+p(1,n-1)+p(2,n-1))*0.333_dp
  p(n,1) = (p(n,2)+p(n-1,1)+p(n-1,2))*0.333_dp
  p(n,n) = (p(n-1,n)+p(n-1,n-1)+p(n,n-1))
  psi(1,1) = 0.0!(psi(1,2)+psi(2,1)+psi(2,2))*0.333_dp
  psi(1,n) = (psi(2,n)+psi(1,n-1)+psi(2,n-1))*0.333_dp
  psi(n,1) = (psi(n,2)+psi(n-1,1)+psi(n-1,2))*0.333_dp
  psi(n,n) = (psi(n-1,n)+psi(n-1,n-1)+psi(n,n-1))*0.333_dp
 
 end do 

 !obstaculo
 !do i = int(n*0.15) - 5, int(n*0.15) + 5
 ! do j = 0*n/11, 1*n/11
 !  psi(i,j) = 0.0_dp
 !  w(i,j) = 0.0_dp
 ! end do
 !end do

 do i = int(n*0.15) - 5, int(n*0.15) + 5
  do j = 5*n/11, 6*n/11
   psi(i,j) = 0.0_dp
   p(i,j) = 0.0_dp
  end do
 end do 



 !vorticidad
 do i = int(n*0.15) - 5, int(n*0.15) + 5
  w(i,5*n/11) = 2.0_dp*(psi(i,5*n/11) - psi(i-1,5*n/11))/d
  w(i,6*n/11) = 2.0_dp*(psi(i,6*n/11) - psi(i+1,6*n/11))/d 
 end do 

 !do i = 5*n/12, 7*n/12
 ! w(int(n*0.8) - 5,i) = 2.0_dp*(psi(int(n*0.8) - 5,i-1)-psi(int(n*0.8) - 5,i-1))/d 
 ! w(int(n*0.8) + 5,i) = 2.0_dp*(psi(int(n*0.8) + 5,i+1)-psi(int(n*0.8) + 5,i+1))/d
 !end do


 ! calculo de la vorticidad a cada paso

 do i = 2,n-1
  do j = 2,n-1
   w(i,j) = 0.5_dp*(v(i+1,j)-v(i-1,j))/dx - 0.5_dp*(u(i,j+1)-u(i,j-1))/dy 
  end do
 end do 



 do i = 2,n-1
  do j = 2,n-1
   t1 = dx**2 * dy**2 *w(i,j) 
   t2 = dy**2 *(psi(i+1,j)+psi(i-1,j))
   t3 = dx**2 *(psi(i,j+1)+psi(i,j-1))
   psi(i,j) = 0.5_dp*(t1+t2+t3)/(d)

   a1 = 1.0_dp-2.0_dp*nu*dt*((1.0_dp/dx**2) + (1.0_dp/dy**2))
   a2 = (dt/dx)*((nu/dx) - 0.25_dp*(psi(i,j+1) - psi(i,j-1))/dy)
   a3 = (dt/dx)*((nu/dx) + 0.25_dp*(psi(i,j+1) - psi(i,j-1))/dy)
   a4 = (dt/dy)*((nu/dy) - 0.25_dp*(psi(i+1,j) - psi(i-1,j))/dx)
   a5 = (dt/dy)*((nu/dy) + 0.25_dp*(psi(i+1,j) - psi(i-1,j))/dx)
   w(i,j) = a1*w(i,j) + a2*w(i+1,j)+a3*w(i-1,j) + a4*w(i,j+1)+a5*w(i,j-1)

   t1 = dx**2 * dy**2 *w(i,j) 
   t2 = dy**2 *(psi(i+1,j)+psi(i-1,j))
   t3 = dx**2 *(psi(i,j+1)+psi(i,j-1))
   psi(i,j) = 0.5_dp*(t1+t2+t3)/(d)

   c1 = p(i+1,j)+p(i-1,j)
   c2 = p(i,j+1)+p(i,j-1)
   c3 = psi(i+1,j+1)-psi(i-1,j+1)-psi(i+1,j-1)+psi(i-1,j-1)
   c4 = psi(i+1,j)-2.0_dp*psi(i,j)+psi(i-1,j)
   c5 = psi(i,j+1)-2.0_dp*psi(i,j)+psi(i,j+1)
   p(i,j) = 0.5_dp*(dy**2 * t1 + dx**2 * t2)/d + (rho/d)*(0.0625_dp*c3 - c4*c5)
  end do
 end do 

 !do i = 2,n-1
  !do j = 2,n-1
  ! a1 = 1.0_dp-2.0_dp*nu*dt*((1.0_dp/dx**2) + (1.0_dp/dy**2))
  ! a2 = (dt/dx)*((nu/dx) - 0.25_dp*(psi(i,j+1) - psi(i,j-1))/dy)
  ! a3 = (dt/dx)*((nu/dx) + 0.25_dp*(psi(i,j+1) - psi(i,j-1))/dy)
  ! a4 = (dt/dy)*((nu/dy) - 0.25_dp*(psi(i+1,j) - psi(i-1,j))/dx)
  ! a5 = (dt/dy)*((nu/dy) + 0.25_dp*(psi(i+1,j) - psi(i-1,j))/dx)
  ! w(i,j) = a1*w(i,j) + a2*w(i+1,j)+a3*w(i-1,j) + a4*w(i,j+1)+a5*w(i,j-1)
  !end do
 !end do 

 !do i = 2,n-1
 !  do j = 2,n-1
 !   t1 = dx**2 * dy**2 *w(i,j) 
 !   t2 = dy**2 *(psi(i+1,j)+psi(i-1,j))
 !   t3 = dx**2 *(psi(i,j+1)+psi(i,j-1))
 !   psi(i,j) = 0.5_dp*(t1+t2+t3)/(d)
 !  end do
 !end do 

! Calculo de la presion

! do i = 2, n-1
!  do j = 2, n-1
!   c1 = p(i+1,j)+p(i-1,j)
!   c2 = p(i,j+1)+p(i,j-1)
!   c3 = psi(i+1,j+1)-psi(i-1,j+1)-psi(i+1,j-1)+psi(i-1,j-1)
!   c4 = psi(i+1,j)-2.0_dp*psi(i,j)+psi(i-1,j)
!   c5 = psi(i,j+1)-2.0_dp*psi(i,j)+psi(i,j+1)
!   p(i,j) = 0.5_dp*(dy**2 * t1 + dx**2 * t2)/d + (rho/d)*(0.0625_dp*c3 - c4*c5)           
!  end do 
! end do


! Calculo del nuevo campo de velocidades
 

 do m = 1,1000

  if(k == m*kend/1000 )then

   write(FileNameu, '(a,i0,a)') "u_",m,".txt"
   write(FileNamev, '(a,i0,a)') "v_",m,".txt"
   open(NumFile_u,file=FileNameu)
   open(NumFile_v,file=FileNamev)
 
   do i = 2,n-1
    do j = 2,n-1
     u(i,j) = 0.5_dp*(psi(i,j+1)-psi(i,j-1))/dx
     v(i,j) = 0.5_dp*(psi(i-1,j)-psi(i+1,j))/dy
    end do 
   end do

   write(NumFile_u,*) u
   write(NumFile_v,*) v
   close(NumFile_u)
   close(NumFile_v)
  end if
 end do 




end do  


 write(1,*) u
 write(2,*) v
 write(3,*) n
 !write(4,*) w 
 write(50,*) p 
 !write(6,*) psi 
 write(7,*) x 
 write(8,*) y 

 close(1)
 close(2)
 close(3)
 close(50)
 close(7)
 close(8)

end program main






