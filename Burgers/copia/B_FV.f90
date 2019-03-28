program uBurgersCD
implicit none
integer::i,j,k,nx,ny,nt,ns,nf
real*8 ::dx,dy,dt,x0,xL,y0,yL,t,Tmax,R,nu
real*8,allocatable ::u(:,:),x(:),y(:),v(:,:)

!Domain
x0 = 0.0d0 !left
xL = 1.0d0 !right

y0 = 0.0d0 !bottom
yL = 1.0d0 !up

!number of points
nx = 25
ny = 25

!grid spacing (spatial)
dx = (xL-x0)/dfloat(nx)
dy = (yL-y0)/dfloat(ny)

!spatial coordinates 
allocate(x(0:nx))
do i=0,nx
x(i) = x0 + dfloat(i)*dx
end do

allocate(y(0:ny))
do j=0,ny
y(j) = y0 + dfloat(j)*dy
end do


!maximum time desired
Tmax = 1.0d0

!diffusion constant:
R = 1000.0d0
nu= 1/R

!time step
dt = 0.001d0

!number of points in time
nt = nint(Tmax/dt)

!number of snapshots to plot
ns = nt

!frequency for plotting
nf = nint(dfloat(nt)/dfloat(ns))

!u, v: convected variables
allocate(u(0:nx,0:ny))
allocate(v(0:nx,0:ny))

!initial condition
t = 0.0d0
do j=0,ny
do i=0,nx
u(i,j) = 0.75d0-(1/(4.0d0*(1.0d0+exp(0.03125*R*(-4*i+4*j)))))
end do
end do

do j=0,ny
do i=0,nx
v(i,j) = 0.75d0+(1/(4.0d0*(1.0d0+exp(0.03125*R*(-4*i+4*j)))))
end do
end do

!Plot initial condition
open(18,file='uCD.plt')
open(1,file = 't1.dat')
open(2,file = 't2.dat')
open(3,file = 't3.dat')
open(4,file = 't4.dat')
open(5,file = 't5.dat')
open(6,file = 't6.dat')
open(7,file = 'tiempos.plt')
do j=0,ny
do i=0,nx
write(18,*) x(i),y(j),u(i,j)
end do
end do


!time integration
do k=1,nt
    
    call CN(nx,ny,dx,dy,dt,t,nu,x,y,u,v)
      
    !update t
    t = t+dt 
if (mod(k,nf).eq.0) then
!write(18,100)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
do j=0,ny
do i=0,nx
if (k == 1) then
   write(1,*) x(i),y(j),u(i,j),v(i,j)
   write(7,*) t
  endif
  if (k == 1*nt/5) then
   write(2,*) x(i),y(j),u(i,j),v(i,j)
   write(7,*) t
  endif
  if (k == 1*nt/4 ) then
   write(3,*) x(i),y(j),u(i,j),v(i,j)
   write(7,*) t
  endif
  if (k == 2*nt/5 ) then
   write(4,*) x(i),y(j),u(i,j),v(i,j)
   write(7,*) t
  endif
  if (k == 2*nt/4 ) then
   write(5,*) x(i),y(j),u(i,j),v(i,j)
   write(7,*) t
  endif
  if (k == 3*nt/5) then
   write(6,*) x(i),y(j),u(i,j),v(i,j)
   write(7,*) t
  endif
 write(18,*) x(i),y(j),u(i,j),v(i,j)
end do
end do
end if

    
	print*,k,t,maxval(u)


end do
  
close(18)


100 format(a16,i8,a4,i8,a10,f10.4,a3)
end

!-----------------------------------------------------------------------------!
!Crank-Nicolson(CN) scheme with ADI
!-----------------------------------------------------------------------------!
subroutine CN(nx,ny,dx,dy,dt,t,nu,x,y,u,v)
implicit none
integer::nx,ny,i,j,m,n,d,l
real*8 ::dx,dy,dt,t,nu,bx,by,ae,an, ua, ub
real*8 ::u(0:nx,0:ny),x(0:nx),y(0:ny), v(0:nx,0:ny)
real*8,allocatable ::a(:),b(:),c(:),r(:),q(:)
real*8 ::s(0:nx,0:ny),p(0:nx,0:ny),z(0:nx,0:ny),be(0:nx,0:ny),bw(0:nx,0:ny),bn(0:nx,0:ny),bs(0:nx,0:ny)

do j=0,ny
do i=0,nx
s(i,j) = 0.0d0
p(i,j) = 0.0d0
z(i,j) = 0.0d0
end do
end do

ae = 0.5d0*nu*dt/(dx*dx)
an = 0.5d0*nu*dt/(dy*dy)


do i=0,nx
l=nint(i+0.5)
m=nint(i-0.5)
end do
do j=0,ny
n= nint(j+0.5)
d= nint(j-0.5)
end do


do j=0,ny
do l=0,nx
be(l,j)= 0.5d0*u(l,j)*dt/(dx)
end do
do m=0,nx
bw(m,j)= 0.5d0*u(m,j)*dt/(dx)
end do
end do

do i=0,nx
do n=0,ny
bn(i,n)= 0.5d0*v(i,n)*dt/(dy)
end do
do d=0,ny
bs(i,d)= 0.5d0*v(i,d)*dt/(dy)
end do
end do

!compute source term in known step
do j=1,ny-1
do i=0,nx
do n=0, ny
do d=0, ny
s(i,j) = u(i,j) + (an-bn(i,n))*(u(i,j+1))+(bs(i,d)- bn(i,n)-2.0d0*an)*u(i,j)+(bs(i,d)+an)*(u(i,j-1))
end do
end do
end do
end do


do j=1,ny-1
do i=1,nx-1
do l=1,nx-1
do m=1,nx-1
p(i,j) = s(i,j) +(ae-be(l,j))*s(i+1,j) + (bw(m,j)-be(l,j)-2.0d0*ae)*s(i,j)+(bw(m,j)+ae)*s(i-1,j) 
end do
end do
end do
end do

!x-sweep to compute intermediate values:
do j=1,ny-1
do n=1,ny-1
do d=1,ny-1  

	!Build coefficient matrix:
	allocate(a(1:nx-1),b(1:nx-1),c(1:nx-1),r(1:nx-1),q(1:nx-1))

	do i=1,nx-1
	do m=1,nx-1
	do l=1,nx-1
	a(i) = -(bw(m,j)+ae)
	b(i) = (1.0d0-bw(m,j)+be(l,j) +2.0d0*ae)
	c(i) = (be(l,j)-ae)
	r(i) = p(i,j) 
	end do
	end do
	end do
    
	!apply boundary conditions
    ua= u(0,j) + (an-bn(0,n))*(u(0,j+1))+(bs(0,d)- bn(0,n)-2.0d0*an)*u(0,j)+(bs(0,d)+an)*(u(0,j-1))
	! ua = u(0,j)  
 	 !ub = u(nx,j)

    ub= u(nx,j) + (an-bn(nx,n))*(u(nx,j+1))+(bs(nx,d)- bn(nx,n)-2.0d0*an)*u(nx,j)+(bs(nx,d)+an)*(u(nx,j-1))
    
	r(1)   = r(1) - a(1)*ua        !b.c.
	r(nx-1) = r(nx-1) - c(nx-1)*ub !b.c.
    
	call tdma(a,b,c,r,q,1,nx-1)

	!assign solutions for as z
	do i=1,nx-1
	z(i,j)=q(i)
	end do
    z(0,j) =ua
    z(nx,j)=ub

    deallocate(a,b,c,r,q)

end do
end do
end do


!y-sweep to compute final solution:
do i=1,nx-1
do l=1,nx-1
do m=1,nx-1
  
	!Build coefficient matrix:
	allocate(a(1:ny-1),b(1:ny-1),c(1:ny-1),r(1:ny-1),q(1:ny-1))

	do j=1,ny-1
	do n=1,ny-1
	do d=1,ny-1
	a(j) = -(bs(i,d)+an)
	b(j) = (1.0d0-bs(i,d)+bn(i,n)+2.0d0*an)
	c(j) = (bn(i,n)-an)
	r(j) = z(i,j) 
	end do
	end do
	end do
    
	!apply boundary conditions
    ua = u(i,0)  
    ub = u(i,ny) 
    
	r(1)   = r(1) - a(1)*ua        !b.c.
	r(ny-1) = r(ny-1) - c(ny-1)*ub !b.c.
    
	call tdma(a,b,c,r,q,1,ny-1)

	!assign solutions for as z
	do j=1,ny-1
	u(i,j)=q(j)
	end do

    deallocate(a,b,c,r,q)
end do
end do
end do


end 

!------------------------------------------------------------------!
!Tridiagonal matrix algorithm (TDMA)
!Thomas algorithm
!solution tridiagonal systems
!a: lower diagonal
!b: main diagonal
!c: upper diagonal
!r: source vector
!x: solution vector
!   for indices s(start) to e(end)
!   i: s,s+1,s+2, ....,i,....,e 
!
!Note: a(s) and c(e) are dummy coefficients, not used.
!------------------------------------------------------------------!

subroutine tdma(a,b,c,r,x,s,e)
implicit none
integer s,e,i
real*8, dimension(s:e) ::a,b,c,r,x    

! forward elimination phase
do i=s+1,e
b(i) = b(i) - a(i)/b(i-1)*c(i-1)
r(i) = r(i) - a(i)/b(i-1)*r(i-1)
end do
! backward substitution phase 
x(e) = r(e)/b(e)
do i=e-1,s,-1
x(i) = (r(i)-c(i)*x(i+1))/b(i)
end do

return
end


