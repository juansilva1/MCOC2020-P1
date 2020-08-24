# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 14:16:43 2020

@author: jpsil
"""


import scipy as sp
from scipy.integrate import odeint
import matplotlib.pylab as plt

# Unidades base
cm = 0.01 # m
inch = 2.54*cm
g = 9.81 #m/s^2

#coeficiente de arrastre
ro = 1.2225 #kg/m^3
cd = 0.47
D = 8.5*inch
r = D/2
A = sp.pi * r**2
CD = 0.5*ro*cd*A

# masa
m = 15 #kg

Viento=[0,10,20]
def bala(z,t):
    zp = sp.zeros(4)
    zp[0] = z[2]
    zp[1] = z[3]
    v = z[2:4] # saca las ultimas 2 componentes
    v[0]=v[0]-Viento1
    v2=sp.dot(v,v)
    vnorm = sp.sqrt(v2)
    FD = -CD*v2*(v/vnorm)
    zp[2] = FD[0]/m
    zp[3] = FD[1]/m -g
    
    return zp

plt.figure(1)
for i in Viento:
    Viento1=i
    # vector de tiempo
    t = sp.linspace(0,6,1001)
    vi = 100*1000./3600.  
    # parte en el origen y tiene vx = vy = 2. m/s
    z0 = sp.array([0,0,vi,vi])
    sol = odeint(bala,z0,t)
    
    x = sol[:,0]
    y = sol[:,1]



    plt.plot(x,y)
plt.axis([0,150,0,50])
plt.title('Trayectoria para distintos vientos')
plt.legend(['V = 0 m/s','V = 10.0 m/s','V = 20.0 m/s'])
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
plt.grid(True)
nombre='Grafico 1'
plt.savefig(f"{nombre}.png")
plt.close()
