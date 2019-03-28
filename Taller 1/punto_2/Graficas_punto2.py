import numpy as np 
import matplotlib.pyplot as plt

l = 1.0; u=1.0;rho = 1.0; gamma = 0.025

markers1 = ["k.","g.","r.","b.","y.","m."]
markers2 = ["k-","g-","r-","b-","y-","m-"]

def grafica(x,y,xnon,ynon,nombre):
	for i in range(np.shape(ynon)[1]):
		plt.plot(xnon,ynon[:,i],markers1[i],label = "$t_{"+str(i)+"} \\approx "+str(i/5.)+" $")
		plt.plot(x,y[:,i],markers1[i])

	plt.title("$ \\phi \quad Vs \quad x$ ")
	plt.ylabel("$\\phi$")
	plt.xlabel("$x$")
	plt.grid(True)
	plt.savefig(nombre,format ="eps",dpi = 400)

i = raw_input()

if(i==1):
	phiRK2NU = np.loadtxt('DCT_CDS_NU.dat', unpack = True)
	phiT = np.loadtxt('DCT_CDS.dat', unpack = True)
	xRK2NU = np.loadtxt('malla_CDS.dat', unpack = True)
	xT = np.linspace(0,1,len(phiT))
	grafica(xT,phiT,xRK2NU,phiRK2NU,"punto3_euler.eps")

##if(i==2):
#	phiRKNU = np.loadtxt('DCT_RK2_NU.dat', unpack = True)
#	phiRK = np.loadtxt('DCT_RK2.dat', unpack = True)
#	xRKNU = np.loadtxt('malla_RK2.dat', unpack = True)
#	xRK = np.linspace(0,1,len(phiRK))
#	grafica(xT,phiT,xRK2NU,phiRK2NU,"punto3_RK2.eps")





#for k in range(np.shape(phiRK2NU)[1]):
#	plt.plot(xRK2NU,phiRK2NU[:,k],markers1[k],label = "$t_{"+str(k)+"} \\approx "+str(k/5.)+" $")
#	plt.plot(xT,phiT[:,k],markers1[k])
#plt.title("$ \\phi \quad Vs \quad x$ ")
#plt.ylabel("$\\phi$")
#plt.xlabel("$x$")
#plt.grid(True)
#plt.savefig("punto3_euler.eps",format ="eps",dpi = 400)
	


#if(i==2):

	#phiRK = np.loadtxt('DCT_RK2.dat', unpack = True)
    #phiRKNU = np.loadtxt('DCT_RK2_NU.dat', unpack = True)
    #T = np.loadtxt('DCT_CDS_t.dat', unpack = True)
    ##xRK = np.linspace(0,1,len(phiRK))
	#grafica(xRKNU,phiRKNU,xRK,phiRK,"punto3_RK2.eps")