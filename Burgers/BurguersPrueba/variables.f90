module variables
 use precision 
 implicit none 
 save
 real(dp),    parameter          :: pi = 4.0_dp*atan(1.0_dp)
 real(dp),    parameter          :: l = 1.0   
 real(dp),    parameter          :: Re = 80  
 integer(dp), parameter          :: n = 5 ,too = 0, tf = 10     
 integer                         :: i,j,k
 real(dp),    parameter          :: xo = 0.0_dp,xl =1.0_dp,yo=0.0_dp,yl=1.0_dp
 real(dp),    parameter          :: dt = 0.001, dx = (xl-xo)/dfloat(n), dy = (yl-yo)/dfloat(n), t=0.0
 real(dp),    dimension(1:n-1)      :: a,b,c,r,q
 real(dp),    dimension(0:n,0:n) :: u,v,up,vp
 real(dp),    dimension(0:n) :: x,y
 real(dp),    dimension(1:n**2) :: uord,vord
 real(dp),    dimension(1:n**2,1:n**2) :: M
 real(dp),    dimension(0:n,0:n) :: b1,b2,b3,b4,b5,b6,S,z
 
 public 
end module variables 