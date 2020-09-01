# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 09:14:36 2020

@author: jpsil
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import math as mt
from scipy.integrate import odeint
# Datos

m = 1.0 #kg
f = 1.0 #Hz
ξ = 0.2
ω = 2*mt.pi*f
k = m*(ω**2)
c = 2*ξ*ω*m

# Datos Reales obtenidos tras resolver la EDO y ocupar la formula:
# https://es.wikipedia.org/wiki/Oscilador_arm%C3%B3nico 
# "Oscilador con amortiguamiento débil"
# Estos valores se obtienen al reemplazar las condiciones iniciales

A=1.06413
phi=2.7968
def real(t):
    lista=[]
    for i in t:
        eje_y=(-A*mt.exp(-((c*i)/(2*m))))*mt.cos(ω*i+phi)
        lista.append(eje_y)
    return lista

def f_odeint(z,t):

    y = z[0]
    yp = z[1]
    # ypp = -A*np.cos(ω*t+fi)*(ω**2)
    zp=np.zeros(2)
    zp[0]=z[1]
    zp[1]=-(c*yp+k*y)/m

    return zp


def eulerint(zp,z0,t,Nsubdivisiones=1):
    Nt=len(t)
    Ndim=len(z0)
    z=np.zeros((Nt,Ndim))
    z[0,:]=z0
    for i in range(1,Nt):
        t_anterior=t[i-1]
        dt=(t[i] - t[i-1])/Nsubdivisiones
        z_temp=z[i-1,:].copy()
        for k in range(Nsubdivisiones):
            z_temp+= dt*zp(z_temp,t_anterior+k*dt)
        z[i,:]=z_temp
    return z

"""
Graficando
"""

t=sp.linspace(0,4.,100)
z0=np.array([1,1])
# Real o analitica
plt.plot(t,real(t),'k',linewidth=2)

# odeint
z_odeint=odeint(f_odeint,z0,t)
plt.plot(t,z_odeint[:,0])

# euler Nsubdivisiones=1
z_eulerint=eulerint(f_odeint , z0, t,Nsubdivisiones=1)
plt.plot(t,z_eulerint[:,0],'g--')

# euler Nsubdivisiones=10
z_eulerint=eulerint(f_odeint , z0, t,Nsubdivisiones=10)
plt.plot(t,z_eulerint[:,0],'r--')

# euler Nsubdivisiones=1
z_eulerint=eulerint(f_odeint , z0, t,Nsubdivisiones=100)
plt.plot(t,z_eulerint[:,0],'--',color='orange')

plt.xlabel('t')
plt.ylabel('x(t)')
plt.legend(['Analítica','odeint','eulerint_1','eulerint_10','eulerint_100'])
nombre='Oscilador_Armónico'
plt.savefig(f"Grafico {nombre}.png")
plt.show()
