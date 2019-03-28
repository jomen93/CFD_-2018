import numpy as np

dt = 0.001 ; dx = 0.05 ; dy = 0.05
l = 1.0 
Re = 0.01
n = int(l/dx)
tf = 1; t = 0

x = np.linspace(0,l,n)
y = np.linspace(0,l,n)

def uo(x,y,t,Re):
    return (3.0 + 1.0*(1.0+np.exp(Re*(t+4.0*x+4.0*y)/32)))/4

def vo(x,y,t,Re):
    return (3.0 - 1.0*(1.0+np.exp(Re*(t+4.0*x+4.0*y)/32)))/4

du = open("u.txt","w")
dv = open("v.txt","w")

# Inicializacion de las matrices
u = v = np.zeros((n,n))
for i in range(n):
	for j in range(n):
		u[i,j] = uo(x[i],y[j],0,Re)
		v[i,j] = vo(x[i],y[j],0,Re)


# Ciclo del metodo
up = vp = np.zeros((n,n))
u[:,0][1:-1] = 0.5*u[:,1][1:-1] + 0.3*u[:,2][1:-1] +0.2*u[:,3][1:-1]#uo(x,y[0],dt,Re)[1:-1] # lado izquierdo en x
u[:,-1][1:-1] = 0.5*u[:,-2][1:-1] + 0.3*u[:,-3][1:-1] +0.2*u[:,-4][1:-1]#uo(x,y[-1],dt,Re)[1:-1] # lado derecho en x
u[0][1:-1] = 0.5*u[1][1:-1] + 0.3*u[2][1:-1] +0.2*u[3][1:-1]#uo(x[0],y,dt,Re)[1:-1] # Abajo 
u[-1][1:-1] = 0.5*u[-2][1:-1] + 0.3*u[-3][1:-1] +0.2*u[-4][1:-1]#uo(x[-1],y,dt,Re)[1:-1] # Arriba

v[:,0][1:-1] = 0.5*v[:,1][1:-1] + 0.3*v[:,2][1:-1] +0.2*v[:,3][1:-1]#uo(x,y[0],dt,Re)[1:-1] # lado izquierdo en x
v[:,-1][1:-1] = 0.5*v[:,-2][1:-1] + 0.3*v[:,-3][1:-1] +0.2*v[:,-4][1:-1]#uo(x,y[-1],dt,Re)[1:-1] # lado derecho en x
v[0][1:-1] = 0.5*v[1][1:-1] + 0.3*v[2][1:-1] +0.2*v[3][1:-1]#uo(x[0],y,dt,Re)[1:-1] # Abajo 
v[-1][1:-1] = 0.5*v[-2][1:-1] + 0.3*v[-3][1:-1] +0.2*v[-4][1:-1]#uo(x[-1],y,dt,Re)[1:-1] # Arriba

u[0,0] = (u[1,0]+u[0,1]+u[1,1])/3.
u[0,-1] = (u[0,-2]+u[1,-1]+u[1,-2])/3.
u[-1,0] = (u[-2,0]+u[-1,1]+u[-2,1])/3.
u[-1,-1] = (u[-2,-1]+u[-1,-2]+u[-2,-1])/3.

v[0,0] = (v[1,0]+v[0,1]+v[1,1])/3.
v[0,-1] = (v[0,-2]+v[1,-1]+v[1,-2])/3.
v[-1,0] = (v[-2,0]+v[-1,1]+v[-2,1])/3.
v[-1,-1] = (v[-2,-1]+v[-1,-2]+v[-2,-1])/3.



for i in range(0,n-1):
	for j in range(0,n-1):
		up[i,j] = u[i,j] + 0.25*((u[i+1,j]+u[i,j])**2-(u[i-1,j]+u[i,j])**2)*dt/dx
		-(0.5/Re)*((u[i+1,j]-2*u[i,j]+u[i-1,j])/dx)*dt/dx
		+0.25*((u[i,j+1]+u[i,j])*(v[i,j+1]+v[i,j]) - (u[i-1,j]+u[i,j])*(v[i-1,j]+v[i,j]))*dt/dy
		-(0.5/Re)*((u[i,j+1]-2*u[i,j]+u[i,j-1])/dy)*dt/dy
		vp[i,j] = v[i,j] + 0.25*((u[i+1,j]+u[i,j])*(v[i+1,j]+v[i,j])-(u[i-1,j]+u[i,j])*(v[i-1,j]+v[i,j]))*dt/dx
		-(0.5/Re)*((v[i+1,j]-2*v[i,j]+v[i-1,j])/dx)*dt/dx
		+0.25*((v[i,j+1]+v[i,j])**2 - (v[i,j-1]+v[i,j])**2)*dt/dy
		-(0.5/Re)*((v[i,j+1]-2*v[i,j]+v[i,j-1])/dy)*dt/dy
		
u = up
v = vp
t = t + dt





np.savetxt(du,u)
np.savetxt(dv,v)

