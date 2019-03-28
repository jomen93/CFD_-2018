import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
from matplotlib.patches import Rectangle

k = raw_input()
x = np.loadtxt("x.dat")
y = np.loadtxt("y.dat")
w = np.loadtxt("w.dat")
p = np.loadtxt("p.dat")
#psi = np.loadtxt("psi.dat")
dim = int(np.loadtxt("n.dat"))+1
n = int(np.loadtxt("n.dat"))



#u = np.loadtxt("FilesAnimation/u_"+str(k)+".txt")
#v = np.loadtxt("FilesAnimation/v_"+str(k)+".txt")
u = np.loadtxt("Files/u_100.txt")
v = np.loadtxt("Files/v_100.txt")
U = u.reshape((dim,dim))
V = v.reshape((dim,dim))
P = p.reshape((dim,dim))
#W = w.reshape((dim,dim))
#Psi = psi.reshape((dim,dim))

f,ax = plt.subplots(1,2,figsize=(14,5))
st = f.suptitle("Flujo a traves de un objeto", fontsize="x-large")
st.set_y(1.0)
#currentAxis = plt.gca()
ax[0].add_patch(Rectangle((int(n*0.15) - 5, 5.0*n/11.0), 10.0, 0.8*n/11.0,alpha=1,facecolor="k"))
#currentAxis.add_patch(Rectangle((int(n*0.8) - 5, 5.0*n/12.0), 10.0, n/6.0,alpha=1,facecolor="k"))
#currentAxis.add_patch(Rectangle((int(n*0.8) - 5, 8.0*n/12.0), 10.0, n/6.0,alpha=1,facecolor="k"))
M = np.hypot(U,V)
im1 = ax[0].quiver(x, y, U,V, M , cmap=plt.cm.jet, width=0.025,scale=100.0)
im1.set_clim(0, 1)
#ax[0].streamplot(x,y,U,V,color="w",linewidth=1.2,density=1.0,arrowstyle="->",arrowsize=1.5)
ax[0].set_title('Magnitud de la velocidad')
ax[0].set_xlim(0,dim-1)
ax[0].set_ylim(0,dim-1)
ax[0].set_xlabel("$x$")
ax[0].set_ylabel("$y$")
f.colorbar(im1,ax=ax[0])
#ax[0].set_aspect("auto")
#ax[0].axis('equal',shrink = 1.0)


#f.colorbar(im,ax =ax, shrink = 1.0)
#Q = plt.quiver(x,y,U,V,M,cmap = cm.jet,width=0.02,scale = 100.0)
#plt.streamplot(x,y,U,V,color="w",linewidth=1.2,density=1.0,arrowstyle="->",arrowsize=1.5)
#m = cm.ScalarMappable(cmap = cm.jet)
#m.set_array(M)
#plt.colorbar(m)
#plt.colorbar(m, extend = "both")
#plt.savefig("CampoVel"+str(i)+".png")


#im2=ax[1].contourf(x,y,P,cmap="hot")
im2=ax[1].imshow(P,cmap="gist_yarg")
im2.set_clim(0,0.1)
ax[1].add_patch(Rectangle((int(n*0.15) - 5, 5.0*n/11.0), 10.0, 0.8*n/11.0,alpha=1,facecolor="k"))

#v = np.linspace(0, 1, 10, endpoint=True)
ax[1].set_xlim(0,dim-1)
ax[1].set_ylim(0,dim-1)
ax[1].set_title("Presion")
ax[1].set_xlabel("$x$")
f.colorbar(im2,ax=ax[1])
#f.colorbar(im2,ax=ax[1],ticks =v)
#ax[1].set_aspect("auto")
#ax[1].axis('equal')



plt.savefig("Evolucion.png")
#plt.show()
"""
plt.imshow(W)
plt.title("Vorticidad")
plt.colorbar()
plt.show()

plt.imshow(Psi)
plt.title("Lineas de corriente")
plt.colorbar()
plt.show()

"""

