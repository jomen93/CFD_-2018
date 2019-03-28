program FEesfera
 implicit none
 
 integer, parameter :: dp = selected_real_kind(15)
 real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
 integer :: i, iteraciones
 real(dp) :: Ao,kx,ky,w,alpha,dt,g,beta,R,tau,M,C,rho,rhoe,fe
 real(dp) :: t,dxdt,dudt,dydt,dvdt
 real(dp) :: xFE,yFE,xCN,yCN,xRK,yRK,uFE, vFE,uCN,vCN,uRK,vRK
 real(dp) :: unew,vnew,xnew,ynew
 real(dp) :: alphau1,alphau2,alphav1,alphav2
 real(dp) :: kx1,kx2,kx3,kx4,lx1,lx2,lx3,lx4
 real(dp) :: ky1,ky2,ky3,ky4,ly1,ly2,ly3,ly4
 
 iteraciones = 10000

 uFE = 0.01; vFE = 0.01
 xFE = 0.0; yFE = 0.0

 uCN = 0.01; vCN = 0.01
 xCN = 0.0; yCN = 0.0

 uRK = 0.01; vRK = 0.01
 xRK = 0.0; yRK = 0.0

 dt = 0.0001; t = 0.0
 Ao = 10.0; w = 1.0
 kx = 1.0; ky = 1.0

 rho = 997.0 !kg/m3
 rhoe = 7874.0 !kg/m3
 C = 0.47
 R = 0.1     !m
 M = 4*pi*R**3*rhoe/3.0 ! kg
 g = 9.81
 
 tau = 80*sqrt(R/g)!s
 alpha = (rho*C*pi*R**3)/(2*M)
 beta = tau/R !s/m
 fe = tau*beta*g*((rho/rhoe)-1)
 
 print*, alpha

 open(1,file = "FEesfera.dat")
 open(2,file = "posiciones")
 open(3,file = "velocidades")
 open(4,file = "Reporte general")

 

 do i = 1,iteraciones


  !-----------------forward euler ---------------------------------
  dxdt = uFE
  dydt = vFE
  dudt = fx(xFE,yFE,uFE,vFE,t)
  dvdt = fy(xFE,yFE,uFE,vFE,t)
  uFE = uFE + dt*dudt
  vFE = vFE + dt*dvdt
  xFE = xFE + uFE*dt
  yFE = yFE + vFE*dt
  !------------------ termina forward euler ----------------------




  !-----------------Cranck Nicholson ---------------------------------
  alphau1 = 1.0 + alpha*dt*beta*vx(xCN,yCN,t+dt)-alpha*dt*0.5*uCN
  alphau2 = alpha*dt*0.5*beta*beta*vx(xCN,yCN,t+dt)*vx(xCN,yCN,t+dt)
  alphav1 = 1.0 + alpha*dt*beta*vy(xCN,yCN,t+dt)-alpha*dt*0.5*vCN
  alphav2 = alpha*dt*0.5*beta*beta*vy(xCN,yCN,t+dt)*vy(xCN,yCN,t+dt)
  unew = (uCN + 0.5*dt*alpha*(beta*beta*vx(xCN,yCN,t)*vx(xCN,yCN,t)-2.0*beta*vx(xCN,yCN,t)*uCN + uCN*uCN)+alphau2)/alphau1
  vnew = (vCN + 0.5*dt*alpha*(beta*beta*vy(xCN,yCN,t)*vy(xCN,yCN,t)-2.0*beta*vy(xCN,yCN,t)*vCN + vCN*vCN)+alphav2+dt*Fe)/alphav1
  
  xnew = xCN + dt*(unew)
  ynew = yCN + dt*(vnew)

  uCN = unew
  vCN = vnew
  xCN = xnew
  yCN = ynew

  !-----------------termina Cranck Nicholson ---------------------------------



  !-----------------Runge kutta 4 ---------------------------------
  kx1 = dt*uRK
  lx1 = dt*fx(xRK,yRK,uRK,vRK,t)
   
  ky1 = dt*vRK
  ly1 = dt*fy(xRK,yRK,uRK,vRK,t)
  
  kx2 = dt*(uRK + 0.5*lx1)
  lx2 = dt*fx(xRK+0.5*kx1,yRK+0.5*ky1,uRK+0.5*lx1,vRK+0.5*ly1,t+0.5*dt)
   
  ky2 = dt*(vRK + 0.5*ly1)
  ly2 = dt*fy(xRK+0.5*kx1,yRK+0.5*ky1,uRK+0.5*lx1,vRK+0.5*ly1,t+0.5*dt)

  kx3 = dt*(uRK + 0.5*lx2)
  lx3 = dt*fx(xRK+0.5*kx2,yRK+0.5*ky2,uRK+0.5*lx2,vRK+0.5*ly2,t+0.5*dt)
   
  ky3 = dt*(vRK + 0.5*ly2)
  ly3 = dt*fy(xRK+0.5*kx2,yRK+0.5*ky2,uRK+0.5*lx2,vRK+0.5*ly2,t+0.5*dt)

  kx4 = dt*(uRK + lx3)
  lx4 = dt*fx(xRK+kx3,yRK+ky3,uRK+lx3,vRK+ly3,t+dt)

  xRK = xRK+ (1.0/6.0)*(kx1 + 2.0*kx2 + 2.0*kx3 + kx4)
  yRK = yRK + (1.0/6.0)*(ky1 + 2.0*ky2 + 2.0*ky3 + ky4)
  uRK = uRK + (1.0/6.0)*(lx1 + 2.0*lx2 + 2.0*lx3 + lx4)
  vRK = vRK + (1.0/6.0)*(ly1 + 2.0*ly2 + 2.0*ly3 + ly4)
  
  t = t + dt

  write(*,*) error(R*xRK,R*xFE),error(R*xRK,R*xCN)
   
  write(1,*) t*tau," ",R*xFE," ",R*yFE," ",R*xCN," ",R*yCN," ",R*xRK," ", &
             R*yRK," ",uFE/beta," ",vFE/beta," ",uCN/beta," ",vRK/beta," ", &
             uRK/beta," ",vRK/beta, error(xRK,xFE),error(xRK,xCN),&
             error(yRK,yFE),error(yRK,yCN),error(uRK,uFE),error(uRK,uCN),error(vRK,vFE),error(vRK,vCN)

  write(2,*) R*xFE," ",R*yFE,R*xCN," ",R*yCN," ",R*xRK," ",R*yRK
  write(3,*) uFE/beta," ",vFE/beta," ",uCN/beta," ",vRK/beta,uRK/beta," ",vRK/beta
  write(4,*) tau," ", A," ",w," ", kx," ",ky 
  
  end do 

  close(1)
  close(2)
  close(3)
  close(4)
 
 contains 
 real function fx(x,y,u,v,t)
  implicit none 
  real(dp) , intent(in) :: x,y,u,v,t
  fx = alpha*(beta*vx(x,y,t)-u)**2
 end function fx

 real function fy(x,y,u,v,t)
  implicit none 
  real(dp) , intent(in) :: x,y,u,v,t
  fy = alpha*(beta*vy(x,y,t)-v)**2 + fe
 end function fy

 real function gx(x,y,u,v,t)
  implicit none 
  real(dp) , intent(in) :: x,y,u,v,t
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

 real function error(xref,x)
  implicit none
  real(dp) , intent(in) :: xref,x
  error = abs((xref-x)/xref)
 end function



end program FEesfera






   










