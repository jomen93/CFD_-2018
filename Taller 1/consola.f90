program rungekutta4
 implicit none

 integer, parameter :: dp = selected_real_kind(15)
 real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
 integer :: i,e
 real(dp) :: Ao,kx,ky,w,alpha,dt,g,beta,R,tau,M,C,rho,rhoe,fe
 real(dp) :: t,y,x,u,v,dxdt,dudt,dydt,dvdt
 real(dp) :: unew,vnew,xnew,ynew
 real(dp) :: alphau1,alphau2,alphav1,alphav2
 real(dp) :: kx1,kx2,kx3,kx4,lx1,lx2,lx3,lx4
 real(dp) :: ky1,ky2,ky3,ky4,ly1,ly2,ly3,ly4

 u = 0.01; v = 0.01
 x = 0.0; y = 0.0
 dt = 0.001; t = 0.0
 Ao = 10.0; w = 15.0
 kx = 1.0; ky = 1.0

 rho = 1.0
 rhoe = 10.0
 C = 1.0
 R = 0.5
 M = 1.0
 tau = 0.01
 g = 10.0

 alpha = (rho*C*pi*R**3)/(2*M)
 beta = tau/R
 fe = tau*beta*g*((rho/rhoe)-1)
 


 open(14,file = "RKEsferax.dat")
  open(15, file = "RKEsferay.dat") 
  open(14,file = "RKEsferax.dat")
  open(15, file = "RKEsferay.dat") 
  do i = 1,20000

   kx1 = u
   lx1 = fx(x,y,u,v,t)
   
   ky1 = v
   ly1 = fy(x,y,u,v,t)
  
   kx2 = (u + 0.5*lx1)
   lx2 = fx(x+0.5*kx1,y+0.5*ky1,u+0.5*lx1,v+0.5*ly1,t+0.5*dt)
   
   ky2 = (v + 0.5*ly1)
   ly2 = fy(x+0.5*kx1,y+0.5*ky1,u+0.5*lx1,v+0.5*ly1,t+0.5*dt)

   kx3 = (u + 0.5*lx2)
   lx3 = fx(x+0.5*kx2,y+0.5*ky2,u+0.5*lx2,v+0.5*ly2,t+0.5*dt)
   
   ky3 = (v + 0.5*ly2)
   ly3 = fy(x+0.5*kx2,y+0.5*ky2,u+0.5*lx2,v+0.5*ly2,t+0.5*dt)

   kx4 = (u + lx3)
   lx4 = fx(x+kx3,y+ky3,u+lx3,v+ly3,t+dt)

   x = x + (1.0/6.0)*dt*(kx1 + 2.0*kx2 + 2.0*kx3 + kx4)
   y = y + (1.0/6.0)*dt*(ky1 + 2.0*ky2 + 2.0*ky3 + ky4)
   u = u + (1.0/6.0)*dt*(lx1 + 2.0*lx2 + 2.0*lx3 + lx4)
   v = v + (1.0/6.0)*dt*(ly1 + 2.0*ly2 + 2.0*ly3 + ly4)
   t = t + dt
  
   write(14,*) t," ",x," ",u
   write(15,*) t," ",y," ",v 
   !write(*,*) t," ",x," ",y," ",u," ",v
  end do 
  close(14)
  close(15)
  !write(*,*) "Done Runge Kutta-4!"


 contains 
 real function fx(x,y,u,v,t)
  implicit none 
  real(dp) , intent(in) :: x,y,u.v,t
  fx = alpha*(beta*vx(x,y,t)-u)**2
 end function fx

 real function fy(x,y,u,v,t)
  implicit none 
  real(dp) , intent(in) :: x,y,v,u,t
  fy = alpha*(beta*vy(x,y,t)-v)**2 + fe
 end function fy

 real function gx(x,y,u,v,t)
  implicit none 
  real(dp) , intent(in) :: x,y,u.v,t
  gx = u
 end function gx

 real function gy(x,y,u,v,t)
  implicit none 
  real(dp) , intent(in) :: x,y,v,u,t
  gy = v
 end function gy

 real function vx(x,y,t)
  implicit none 
  real(dp) , intent(in) :: x,y,t
  vx = (Ao/kx)*sin(w*t*tau)*sin(kx*x*R)*sin(ky*y*R)
 end function vx

 real function vy(x,y,t)
  implicit none 
  real(dp) , intent(in) :: x,y,t
  vy = (Ao/kx)*sin(w*t*tau)*cos(kx*x*R)*cos(ky*y*R)
 end function vy
end program rungekutta4