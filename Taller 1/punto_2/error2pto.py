import numpy as np 
import matplotlib.pyplot as plt 

l = 1.0; u=1.0;rho = 1.0; gamma = 0.025
Pe = rho*u*l/gamma

def phiT(x):
	return (np.exp(x*Pe/l)-1.0)/(np.exp(Pe)-1.0)


phi = np.loadtxt("DC_CDS.dat", unpack = True)
phiuds = np.loadtxt("DC_UDS.dat", unpack = True)
phicds = np.loadtxt("CDS_NU.dat", unpack = True)

phix = np.linspace(0,l,len(phi))
xuds = np.linspace(0,l,len(phiuds))
xcds = np.linspace(0,l,len(phicds))
x = np.linspace(0,l,1000)

plt.plot(x,phiT(x),"b-",label = "Teorica")
plt.plot(phix,phi,"r+",label = "CDS")
plt.plot(xuds,phiuds,"k*",label = "UDS")
plt.plot(phicds[:,0],phicds[:,1],"g.",label = "CDS NU")
plt.title("$\\phi \quad Vs \quad  x$")
plt.xlabel("x")
plt.ylabel("$\\phi$")
plt.xlim(0.8,1.01)
plt.ylim(-0.1,1.05)
plt.grid(True)
plt.legend(loc = 2)
plt.show()

def Error(phiT,phiE):
	suma = 0.0
	argument = np.linspace(0,l,len(phiT(phiE)))
	for i in range(0,len(argument)):
		suma += abs(phiT(argument)[i]-phiE[i])
	return suma/len(argument)






print Error(phiT,phi), Error(phiT,phiuds)







