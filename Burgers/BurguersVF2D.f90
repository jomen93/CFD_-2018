program Burguers2D
 implicit none

 integer,         parameter :: dp = selected_real_kind(15)
 real(dp),        parameter :: pi = 4.0_dp*atan(1.0_dp)
 real(dp),        parameter :: dt = 0.0000001, dx = 0.05, dy = 0.05, t=0.0
 real(dp),        parameter :: l = 1.0   
 real(dp),        parameter :: Re = 80  
 integer(dp),     parameter :: n = l/dx ,to = 0, tf = 10     
 integer                    :: i,j,k
 real(dp)                   :: xi,yi,ko,no,u1,u2
 
 real(dp),        dimension(:,:), allocatable :: u, v
 real(dp),        dimension(:,:), allocatable :: up, vp
 real(dp),        dimension(:)  , allocatable :: x,y
 


 !print *, "Tamaño del sistema", n

 allocate(u(0,n),v(0,n))
 allocate(up(0,n),vp(0,n))
 allocate(x(n),y(n))

 open(1, file = "u.dat")
 open(2, file = "v.dat")
 open(3, file = "x.dat")
 open(4, file = "y.dat")
 open(5, file = "n.dat")

 print *, "Tamano del sistema ",n
 ! Inicializar la variables 
 do i = 0,n
   x(i) = i*dx
   y(i) = i*dy
 end do


do i=1,n
 do j=1,n
  u1 = 0.25_dp*((u(i+1,j))**2 + (u(i-1,j)+u(i,j))**2)-(1/Re)*(0.5_dp*(u(i+1,j)-2.0_dp*u(i,j)+u(i-1,j)/dx))
  u2 = 0.25_dp*((u(i,j+1)+u(i,j))*(v(i,j+1)+v(i,j)) - (u(i-1,j)+u(i,j))*(v(j-i,j)+v(i,j))) & 
       - (1/Re)*(2.0_dp*(u(i,j+1)-2.0_dp*u(i,j)+u(i,j-1)/dy))
  up(i,j) = u1*dt/dx + u2*dt/dy 
 enddo
enddo


do i=0,n
  do j=0,1
   write(1,*) u(i,j)
   write(2,*) v(i,j)
  end do
end do



!open(unit=2, file='graph1.txt', ACTION="write", STATUS="replace")
!do i=1,N
!     write(2, '(1000F14.7)')( real(Vec(i,j)) ,j=1,M)
!end do


 
 !do i = 0,n
   !u(i,j) = (uo(x(i),y(j),0.0_dp), j = 1,n)
   !v(i,j) = vo(x(i),y(j),0.0_dp)
 !end do

 !do i=0,n
  !do j=0,1
   !write(1,*) u(i,j)
   !write(2,*) v(i,j)
  !end do
 !end do
 
 

 !do k = to,tf
  !write(*,fmt="(a3,a,t21,f6.2,a)",advance="no") achar(13), "progreso completado:    ", (real(k)/real(tf))*100.0," %"
  
  !u(1,1) = uo(dx,dy,dt)
  !u(1,n) = uo(0.0_dp,no,t+dt)
  !u(n,1) = uo(no,0.0_dp,t+dt)
  !u(n,n) = uo(no,no,t+dt)
  
  !v(1,1) = (v(1,2)+v(2,1)+v(2,2))/3!vo(0.0_dp,0.0_dp,0.0_dp)
  !v(1,n) = vo(0.0_dp,no,t+dt)
  !v(n,1) = vo(no,0.0_dp,t+dt)
  !v(n,n) = vo(no,no,t+dt)

  !do i = 0,n
  !xi = real(i)
  !lado izquierdo
   !u(0,i) = uo(0.0_dp,y(i),dt)
   !v(0,i) = vo(0.0_dp,y(i),dt)
  !lado derecho
   !u(n,i) = uo(no,xi,t+dt)
   !v(n,i) = vo(no,xi,t+dt)
  !abajo
   !u(i,1) = uo(xi,0.0_dp,t+dt)
   !v(i,1) = vo(xi,0.0_dp,t+dt)
  !arriba
   !u(i,n) = uo(xi,no,t+dt)
   !v(i,n) = vo(xi,no,t+dt)
  !end do 
  
  !do i = 2,n-1
  ! do j = 2,n-1
  !  up(i,j) = u(i,j)+(0.5_dp*dt/(dx**2 *Re))*(u(i+1,j)-u(i-1,j)+v(i+1,j)-v(i-1,j)) &
  !            -0.5_dp*dt*u(i,j)*(u(i+1,j)-u(i-1,j))/dx &
  !            -0.5_dp*dt*v(i,j)*(u(i,j+1)-u(i,j-1))/dy
  !  vp(i,j) = v(i,j)+(0.5_dp*dt/(dy**2 *Re))*(u(i,j+1)-u(i,j-1)+v(i,j+1)-v(i,j-1)) &
  !            -0.5_dp*dt*u(i,j)*(v(i+1,j)-v(i-1,j))/dx &
  !            -0.5_dp*dt*v(i,j)*(v(i,j+1)-v(i,j-1))/dy
  ! end do 
  !end do
  
  !do i = 1,n
  ! do j = 1,n
  !  u(i,j) = up(i,j)
  !  v(i,j) = vp(i,j)
  ! end do
  !end do

 !end do 

 

 write(5,*) n
 
 close(1)
 close(2)
 close(3)
 close(4)
 close(5)

 !print *, (up(j,0), j = 0,3)
 



contains
real function uo(x,y,t)
 implicit none 
 real(dp), intent(in) :: x,y,t
 uo = 0.25_dp*(3.0_dp - 1.0_dp/(1.0_dp+exp(-0.03125_dp*Re*(t+4.0_dp*x+4.0_dp*y))))
end function uo

real function vo(x,y,t)
 implicit none 
 real(dp), intent(in) :: x,y,t
 vo = 0.25_dp*(3.0_dp + 1.0_dp/(1.0_dp+exp(-0.03125_dp*Re*(t+4.0_dp*x+4.0_dp*y))))
end function vo

end program Burguers2D







