# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 19:55:25 2020

@author: jpsil
"""

import numpy as np
import scipy as sp
from scipy.integrate import odeint
import math as mt

import warnings
warnings.filterwarnings("ignore")

from matplotlib.pylab import *
from leer_eof import leer_eof
from matplotlib import pyplot as plt
from sys import argv
fname = argv[1]
t, x, y, z, vx, vy, vz = leer_eof(fname)

# Punto inicial
vector_inicial_real=np.array([x[0],y[0],z[0],vx[0],vy[0],vz[0]])

# Punto final
vector_final_real=np.array([x[-1],y[-1],z[-1],vx[-1],vy[-1],vz[-1]])

# Datos
radio_tierra=6378136.3 #m
masa_tierra=(5.97219*(10**24)) #kg
G=(6.67392*(10**-11))# m3/s2kg
omega=7.2921150e-5# rad/s

km3=(1000.)**3
km5=(1000.)**5
km6=(1000.)**6


mu=398600.4415*km3

# print(mu)
# mu=G*masa_tierra
# print(mu)

J2=1.7555280486257942e+25

J3=-2.6191328602512207e+29


# algunas funciones utiles para el futuro:


# Esta funcion se encarga de relacionar las variables del problema
def bala(z,t):
    
    c=np.cos(omega*t)
    s=np.sin(omega*t)
    R=np.array([[c,  -s,  0],
                [s, c,  0],
                [0,  0,  1]])
    Rp=omega*np.array([[-s,  -c,  0],
                       [c, -s,  0],
                       [0,  0,  0]])
    Rpp=(omega**2)*np.array([[-c,  s,  0],
                             [-s, -c,  0],
                             [0,  0,  0]])
    
    
    x=z[0:3]
    xp=z[3:6]
    
    r=np.sqrt(np.dot(x,x))
    
    xstill=R@x
    rnorm=xstill/r
    Fg=-(mu/r**2)*rnorm
    
    z2=xstill[2]**2
    rflat=xstill[0]**2+xstill[1]**2
    
    FJ2=(J2*xstill)/(r**7)
    FJ2[0]=FJ2[0]*(6*z2-1.5*rflat)
    FJ2[1]=FJ2[1]*(6*z2-1.5*rflat)
    FJ2[2]=FJ2[2]*(3*z2-4.5*rflat)
    
    # Fx=J2 * (x/(r**7)) * (6*(z**2) - 1.5*((x**2) + (y**2)))
    # Fy=J2 * (y/(r**7)) * (6*(z**2) - 1.5*((x**2) + (y**2)))
    # Fz=J2 * (z/(r**7)) * (3*(z**2) - 4.5*((x**2) + (y**2)))
    
    
    FJ3=np.zeros(3)
    FJ3[0]=((J3*xstill[0]*xstill[2])/(r**9)) *(10*z2-7.5*rflat)
    FJ3[1]=((J3*xstill[1]*xstill[2])/(r**9)) *(10*z2-7.5*rflat)
    FJ3[2]=J3*                    (1/(r**9)) *(4*z2*(z2-3*rflat)+1.5*(rflat**2))  
    
    # Fx=J3 * ((x*z)/(r**9)) * (10*(z**2) - 7.5*((x**2) + (y**2)))
    # Fy=J3 * ((y*z)/(r**9)) * (10*(z**2) - 7.5*((x**2) + (y**2)))
    # Fz=J3 * ((1)/(r**9)) * (4*(z**2)*((z**2)-3*((x**2) + (y**2))) + 1.5*(((x**2) + (y**2))**2))
    
    zp=np.zeros(6)
    zp[0:3]=xp
    
    zp[3:6]=R.T@(Fg+FJ2+FJ3-(2*Rp@xp+Rpp@x))
    return zp

from matplotlib.pylab import *
from scipy.integrate import odeint
import datetime as dt
from leer_eof import leer_eof
import numpy as np

#Leer desde comando
"""
#leer de la linea de comando
from sys import argv
eofname=argv[1]
"""

import xml
import xml.etree.ElementTree as ET
from numpy import zeros
import datetime as dt

#Leer desde terminal
eofname=fname

sat_t,sat_x,sat_y,sat_z,sat_vx,sat_vy,sat_vz = leer_eof(eofname)
z0=[sat_x[0],sat_y[0],sat_z[0],sat_vx[0],sat_vy[0],sat_vz[0]]
vector = odeint(bala,z0,t)
# print(sol)




x = vector[:,0]
y = vector[:,1]
z = vector[:,2]
vx = vector[:,3]
vy = vector[:,4]
vz = vector[:,5]










"""
FUNCIONA
"""

from matplotlib.pylab import *
from scipy.integrate import odeint
import datetime as dt
from leer_eof import leer_eof
import numpy as np



#leer de la linea de comando
from sys import argv
archivo=argv[1]



archivo=open(f'{fname}')

string_name=str(fname)
archivo_final=string_name.replace('.EOF','.PRED')

archivo1=open(archivo_final,'w')

contador=0
lista_datos=[]
for line in archivo:
    cambio=line
    if '<X unit="m"' in line:
        cambio=line.replace(line,(' '*6+'<X unit="m">'+str(x[contador])+'</X>' + "\n" ))
        
    elif '<Y unit="m"' in line:
        cambio=line.replace(line,(' '*6+'<Y unit="m">'+str(y[contador])+'</Y>' + "\n" ))
        
    elif '<Z unit="m"' in line:
        cambio=line.replace(line,(' '*6+'<Z unit="m">'+str(z[contador])+'</Z>' + "\n" ))
    
    elif '<VX unit="m/s"' in line:
        cambio=line.replace(line,(' '*6+'<VX unit="m">'+str(vx[contador])+'</VX>' + "\n" ))
        
    elif '<VY unit="m/s"' in line: 
        cambio=line.replace(line,(' '*6+'<VY unit="m">'+str(vy[contador])+'</VY>' + "\n" ))
        
    elif '<VZ unit="m/s"' in line:
        cambio=line.replace(line,(' '*6+'<VZ unit="m">'+str(vz[contador])+'</VZ>' + "\n" ))
        contador+=1
    archivo1.write(cambio)
archivo.close()


# for i in lista_datos:
    
# archivo1.close()