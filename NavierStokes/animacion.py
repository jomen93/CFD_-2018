import numpy as np 
import matplotlib 
import time 
matplotlib.use("TKAgg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import ArtistAnimation
from matplotlib.patches import Rectangle

k = raw_input()
x = np.loadtxt("x.dat")
y = np.loadtxt("y.dat")
dim = int(np.loadtxt("n.dat"))+1
n = int(np.loadtxt("n.dat"))

#Lectura archivos
fig = plt.figure(figsize=(10,10),facecolor="w")
ims = [] 
start = time.time()
for i in range(1,1000):
	u = np.loadtxt("Files/u_"+str(i)+".txt")
	v = np.loadtxt("Files/v_"+str(i)+".txt")
	U = u.reshape((dim,dim))
	V = v.reshape((dim,dim))
	M = np.hypot(U,V)
	im = plt.contourf(x,y,U**2+V**2,cmap=cm.jet,vamx=6,vmin=0)
	ims.append(im.collections)
cbar = plt.colorbar(im)
cbar.set_clim(0,1)
#cbar.set_ticks(np.linspace(0,6,7))
cbar.set_label("Magnitud de la velocidad")
cA = plt.gca()
cA.add_patch(Rectangle((int(n*0.15) - 5, 5.0*n/11.0), 10.0, 0.8*n/11.0,alpha=1,facecolor="k"))
plt.xlim(0,dim-1)
plt.ylim(0,dim-1)
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.axes().set_aspect("equal")
ani = ArtistAnimation(fig,ims,interval=100,repeat=True)
ani.save('flujo.gif',writer='imagemagick')
ani.save('flujo.mp4',writer='ffmpeg')
#plt.show()





#fig = plt.figure(figsize=(6.1,5),facecolor = "w")
#currentAxis = plt.gca()
#currentAxis.add_patch(Rectangle((int(n*0.8) - 5, 2.0*n/12.0), 10.0, n/6.0,alpha=1,facecolor="k"))
#currentAxis.add_patch(Rectangle((int(n*0.8) - 5, 5.0*n/12.0), 10.0, n/6.0,alpha=1,facecolor="k"))
#currentAxis.add_patch(Rectangle((int(n*0.8) - 5, 8.0*n/12.0), 10.0, n/6.0,alpha=1,facecolor="k"))
#plt.title("Magnitud velocidades")
#plt.xlabel("$x$")
#plt.ylabel("$y$")
#M = np.hypot(uvel[0],vvel[0])
#M = np.hypot(uvel[0],vvel[0])
#plt.xlim(0,dim-1)
#plt.ylim(0,dim-1)
#Q = plt.quiver(x,y,uvel[0],vvel[0],M,cmap = cm.jet,width=0.02,scale = 100.0)
#plt.streamplot(x,y,U,V,color="w",linewidth=1.2,density=1.0,arrowstyle="->",arrowsize=1.5)
#m = cm.ScalarMappable(cmap = cm.jet)
#m.set_array(M)
#plt.colorbar(m)
#plt.savefig("CampoVel"+str(i)+".png")
plt.show()



"""
fig = plt.figure(facecolor="w")
ims=[]

plt.imshow(I)
plt.show()

"""

