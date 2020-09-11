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
    Rs=np.array([[c,  -s,  0],
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
    
    xstill=Rs@x
    rnorm=xstill/r
    Fg=-(mu/r**2)*rnorm
    
    z2=xstill[2]**2
    rflat=xstill[0]**2+xstill[1]**2
    
    
    Fx= -5.26658414587738e+25*xstill[0]*xstill[2]**2*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-3.5) - 3.0*xstill[0]*(2.63329207293869e+25*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2) - 8.77764024312897e+24)*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-2.5) - 4.0*xstill[0]*(-6.54783215062805e+29*xstill[2]**3/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(3/2) + 3.92869929037683e+29*xstill[2]/sqrt(xstill[0]**2 + xstill[1]**2 + xstill[2]**2))*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-3.0) - 9.0*xstill[0]*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-5.5)*(-1.1114550082219e+64*xstill[2]**8/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**4 + 2.07471601534756e+64*xstill[2]**6/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**3 - 1.1969515473159e+64*xstill[2]**4/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**2 + 2.17627554057436e+63*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2) - 6.04520983492878e+61) - 8.0*xstill[0]*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-5.0)*(-1.61248026818891e+57*xstill[2]**7/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(7/2) + 2.60477581784362e+57*xstill[2]**5/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(5/2) - 1.18398900811074e+57*xstill[2]**3/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(3/2) + 1.31554334234526e+56*xstill[2]/sqrt(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)) - 7.0*xstill[0]*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-4.5)*(2.09075678222103e+50*xstill[2]**6/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**3 - 2.85103197575595e+50*xstill[2]**4/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**2 + 9.50343991918651e+49*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2) - 4.525447580565e+48) - 6.0*xstill[0]*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-4.0)*(-7.54485777124404e+42*xstill[2]**5/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(5/2) + 8.38317530138227e+42*xstill[2]**3/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(3/2) - 1.79639470743906e+42*xstill[2]/sqrt(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)) - 5.0*xstill[0]*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-3.5)*(-4.67333271026623e+36*xstill[2]**4/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**2 + 4.00571375165677e+36*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2) - 4.00571375165677e+35) + (1.96434964518842e+30*xstill[0]*xstill[2]**3/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(5/2) - 3.92869929037683e+29*xstill[0]*xstill[2]/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(3/2))*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-2.0) + (1.86933308410649e+37*xstill[0]*xstill[2]**4/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**3 - 8.01142750331354e+36*xstill[0]*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**2)*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-2.5) + (xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-4.5)*(8.89164006577524e+64*xstill[0]*xstill[2]**8/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**5 - 1.24482960920853e+65*xstill[0]*xstill[2]**6/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**4 + 4.78780618926359e+64*xstill[0]*xstill[2]**4/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**3 - 4.35255108114872e+63*xstill[0]*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**2) + (xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-4.0)*(1.12873618773224e+58*xstill[0]*xstill[2]**7/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(9/2) - 1.30238790892181e+58*xstill[0]*xstill[2]**5/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(7/2) + 3.55196702433221e+57*xstill[0]*xstill[2]**3/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(5/2) - 1.31554334234526e+56*xstill[0]*xstill[2]/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(3/2)) + (xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-3.5)*(-1.25445406933262e+51*xstill[0]*xstill[2]**6/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**4 + 1.14041279030238e+51*xstill[0]*xstill[2]**4/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**3 - 1.9006879838373e+50*xstill[0]*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**2) + (xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-3.0)*(3.77242888562202e+43*xstill[0]*xstill[2]**5/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(7/2) - 2.51495259041468e+43*xstill[0]*xstill[2]**3/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(5/2) + 1.79639470743906e+42*xstill[0]*xstill[2]/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(3/2))
    Fy= -5.26658414587738e+25*xstill[1]*xstill[2]**2*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-3.5) - 3.0*xstill[1]*(2.63329207293869e+25*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2) - 8.77764024312897e+24)*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-2.5) - 4.0*xstill[1]*(-6.54783215062805e+29*xstill[2]**3/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(3/2) + 3.92869929037683e+29*xstill[2]/sqrt(xstill[0]**2 + xstill[1]**2 + xstill[2]**2))*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-3.0) - 9.0*xstill[1]*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-5.5)*(-1.1114550082219e+64*xstill[2]**8/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**4 + 2.07471601534756e+64*xstill[2]**6/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**3 - 1.1969515473159e+64*xstill[2]**4/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**2 + 2.17627554057436e+63*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2) - 6.04520983492878e+61) - 8.0*xstill[1]*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-5.0)*(-1.61248026818891e+57*xstill[2]**7/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(7/2) + 2.60477581784362e+57*xstill[2]**5/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(5/2) - 1.18398900811074e+57*xstill[2]**3/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(3/2) + 1.31554334234526e+56*xstill[2]/sqrt(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)) - 7.0*xstill[1]*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-4.5)*(2.09075678222103e+50*xstill[2]**6/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**3 - 2.85103197575595e+50*xstill[2]**4/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**2 + 9.50343991918651e+49*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2) - 4.525447580565e+48) - 6.0*xstill[1]*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-4.0)*(-7.54485777124404e+42*xstill[2]**5/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(5/2) + 8.38317530138227e+42*xstill[2]**3/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(3/2) - 1.79639470743906e+42*xstill[2]/sqrt(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)) - 5.0*xstill[1]*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-3.5)*(-4.67333271026623e+36*xstill[2]**4/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**2 + 4.00571375165677e+36*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2) - 4.00571375165677e+35) + (1.96434964518842e+30*xstill[1]*xstill[2]**3/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(5/2) - 3.92869929037683e+29*xstill[1]*xstill[2]/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(3/2))*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-2.0) + (1.86933308410649e+37*xstill[1]*xstill[2]**4/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**3 - 8.01142750331354e+36*xstill[1]*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**2)*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-2.5) + (xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-4.5)*(8.89164006577524e+64*xstill[1]*xstill[2]**8/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**5 - 1.24482960920853e+65*xstill[1]*xstill[2]**6/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**4 + 4.78780618926359e+64*xstill[1]*xstill[2]**4/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**3 - 4.35255108114872e+63*xstill[1]*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**2) + (xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-4.0)*(1.12873618773224e+58*xstill[1]*xstill[2]**7/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(9/2) - 1.30238790892181e+58*xstill[1]*xstill[2]**5/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(7/2) + 3.55196702433221e+57*xstill[1]*xstill[2]**3/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(5/2) - 1.31554334234526e+56*xstill[1]*xstill[2]/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(3/2)) + (xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-3.5)*(-1.25445406933262e+51*xstill[1]*xstill[2]**6/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**4 + 1.14041279030238e+51*xstill[1]*xstill[2]**4/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**3 - 1.9006879838373e+50*xstill[1]*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**2) + (xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-3.0)*(3.77242888562202e+43*xstill[1]*xstill[2]**5/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(7/2) - 2.51495259041468e+43*xstill[1]*xstill[2]**3/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(5/2) + 1.79639470743906e+42*xstill[1]*xstill[2]/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(3/2))
    Fz= -3.0*xstill[2]*(2.63329207293869e+25*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2) - 8.77764024312897e+24)*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-2.5) - 4.0*xstill[2]*(-6.54783215062805e+29*xstill[2]**3/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(3/2) + 3.92869929037683e+29*xstill[2]/sqrt(xstill[0]**2 + xstill[1]**2 + xstill[2]**2))*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-3.0) - 9.0*xstill[2]*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-5.5)*(-1.1114550082219e+64*xstill[2]**8/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**4 + 2.07471601534756e+64*xstill[2]**6/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**3 - 1.1969515473159e+64*xstill[2]**4/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**2 + 2.17627554057436e+63*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2) - 6.04520983492878e+61) - 8.0*xstill[2]*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-5.0)*(-1.61248026818891e+57*xstill[2]**7/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(7/2) + 2.60477581784362e+57*xstill[2]**5/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(5/2) - 1.18398900811074e+57*xstill[2]**3/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(3/2) + 1.31554334234526e+56*xstill[2]/sqrt(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)) - 7.0*xstill[2]*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-4.5)*(2.09075678222103e+50*xstill[2]**6/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**3 - 2.85103197575595e+50*xstill[2]**4/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**2 + 9.50343991918651e+49*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2) - 4.525447580565e+48) - 6.0*xstill[2]*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-4.0)*(-7.54485777124404e+42*xstill[2]**5/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(5/2) + 8.38317530138227e+42*xstill[2]**3/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(3/2) - 1.79639470743906e+42*xstill[2]/sqrt(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)) - 5.0*xstill[2]*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-3.5)*(-4.67333271026623e+36*xstill[2]**4/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**2 + 4.00571375165677e+36*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2) - 4.00571375165677e+35) + (-5.26658414587738e+25*xstill[2]**3/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**2 + 5.26658414587738e+25*xstill[2]/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2))*(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-1.5) + (xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-4.5)*(8.89164006577524e+64*xstill[2]**9/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**5 - 2.13399361578606e+65*xstill[2]**7/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**4 + 1.72361022813489e+65*xstill[2]**5/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**3 - 5.22306129737846e+64*xstill[2]**3/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**2 + 4.35255108114872e+63*xstill[2]/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)) + (xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-4.0)*(1.12873618773224e+58*xstill[2]**8/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(9/2) - 2.43112409665405e+58*xstill[2]**6/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(7/2) + 1.65758461135503e+58*xstill[2]**4/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(5/2) - 3.68352135856674e+57*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(3/2) + 1.31554334234526e+56/sqrt(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)) + (xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-3.5)*(-1.25445406933262e+51*xstill[2]**7/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**4 + 2.394866859635e+51*xstill[2]**5/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**3 - 1.33048158868611e+51*xstill[2]**3/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**2 + 1.9006879838373e+50*xstill[2]/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)) + (xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-3.0)*(3.77242888562202e+43*xstill[2]**6/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(7/2) - 6.2873814760367e+43*xstill[2]**4/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(5/2) + 2.69459206115859e+43*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(3/2) - 1.79639470743906e+42/sqrt(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)) + (xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-2.5)*(1.86933308410649e+37*xstill[2]**5/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**3 - 2.67047583443785e+37*xstill[2]**3/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**2 + 8.01142750331354e+36*xstill[2]/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)) + (xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(-2.0)*(1.96434964518842e+30*xstill[2]**4/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(5/2) - 2.3572195742261e+30*xstill[2]**2/(xstill[0]**2 + xstill[1]**2 + xstill[2]**2)**(3/2) + 3.92869929037683e+29/sqrt(xstill[0]**2 + xstill[1]**2 + xstill[2]**2))
    
    
    # xr=xstill
    
    JNN=-np.array([Fx,Fy,Fz])
    zp=np.zeros(6)
    zp[0:3]=xp
    
    zp[3:6]=Rs.T@(Fg+JNN-(2*Rp@xp+Rpp@x))
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