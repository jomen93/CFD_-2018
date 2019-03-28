#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np 
import matplotlib.pyplot as plt

u = np.loadtxt("B1D.dat")
nombre_archivo = "Burgers1D.png"

plt.plot(u[0],"b--",label = "$t_{1}$")
plt.plot(u[1],"k--",label = "$t_{2}$")
plt.plot(u[2],"r--",label = "$t_{3}$")
plt.plot(u[3],"y--",label = "$t_{4}$")
plt.plot(u[4],"m--",label = "$t_{5}$")
plt.xlabel("$x$")
plt.ylabel("$u$")
plt.ticklabel_format(useOffset=False)
plt.title("Perfil ecuacion de Burgers 1D")
plt.grid(True)
plt.legend()
plt.savefig(nombre_archivo)
plt.show()
