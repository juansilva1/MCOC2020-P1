# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 07:45:03 2020

@author: jpsil
"""

import numpy as np
import scipy as sp
from scipy.integrate import odeint
import math as mt


# Datos
radio_tierra=6371000 #m
masa_tierra=(5.97219*(10**24)) #kg
G=(6.67392*(10**-11))# m3/s2kg
omega=(7.27*(0.00001))# rad/s


# Esta funcion me crea 3 matrices, las cuales son matrices de rotacion
# para pasar de un sistema de coordenadas inersiales a un sistema
# de coordenadas solidario a la tierra en rotacion.
def matriz_R(omega,t):
    teta=omega*t
    R_0=np.array([[np.cos(teta),-np.sin(teta),0],
                  [np.sin(teta), np.cos(teta),0],
                  [     0,             0,     1]])
    
    R_1=omega*np.array([[-np.sin(teta),-np.cos(teta),0],
                        [np.cos(teta), -np.sin(teta),0],
                        [     0,             0,      0]])
    
    R_2=(omega**2)*np.array([[-np.cos(teta),  np.sin(teta),0],
                             [-np.sin(teta), -np.cos(teta),0],
                             [      0,              0,     0]])
    if omega==0:
        vector_R=[R_0,np.zeros((3,3)),np.zeros((3,3))]
        return vector_R
    vector_R=(R_0,R_1,R_2)
    return vector_R


# Esta funcion recibe un vector posicion y retorna la distancia r.
def distancia_satelite_tierra(vector_posicion):
    r=(np.dot(vector_posicion,vector_posicion))**0.5
    return r


# Esta funcion se encarga de relacionar las variables del problema
def bala(z,t):
    
    # Creo un vector zp de 6x1
    zp = sp.zeros(6)
    
    # los primeros 3 "espacios" de zp son los ultimos de z.
    zp[0] = z[3]
    zp[1] = z[4]
    zp[2] = z[5]
    
    # Parto a escribir la ecuacion de estado.
    """
    ¿Ecuacion de estado?
    """
    
    # Parte 1
    parte1_a=(-G*masa_tierra)
    r3=((distancia_satelite_tierra(z[0:3]))**3)
    parte1_final=np.dot((parte1_a/r3),z[0:3])
    
    vector_R=matriz_R(omega, t)
    
    R_0=vector_R[0]
    R_1=vector_R[1]
    R_2=vector_R[2]
    
    # Parte 2
    parte2_a=-np.transpose(R_0)
    parte2_b=np.dot(R_2,z[0:3])
    parte2_c=2*np.dot(R_1,z[3:6])
    parte2_d=parte2_b+parte2_c
    parte2_final=np.dot(parte2_a,parte2_d)
    
    # Parte final = Parte 1 + Parte 2
    parte3_final=parte1_final+parte2_final
    
    # Relaciono los 3 ultimos "espacios" de zp con el vector de 3x1 resultante de la Parte final.
    zp[3:6]=parte3_final
    
    return zp


# Creo un Vector de tiempo en segundos.
segundos=3600*4
t = sp.linspace(0,segundos,segundos+1)

# Por el momento me doy la velocidad que yo desee.
v_km=24551 # km/h
vi = (v_km*1000)/3600

# Defino mis coordenadas iniciales [xi,yi=0,z0=0,vx=0,vy=?,vz=0].
z0 = sp.array([radio_tierra+700000 ,0,0,0,vi,0])

# Utilizo la funcion odeint para resolver el sistema EDO.
sol = odeint(bala,z0,t)

"""
Grafica: Recorrido satelite
"""

import matplotlib.pylab as plt

# Con esto obtengo listas de las posiciones x,y,z por cada segundo que pasa.
x = sol[:,0]
y = sol[:,1]
z = sol[:,2]

# Grafico aparte esta figura, la cual corresponde a las posiciones x,y
# en el espacio, en resumen esto me indica el recorrido del satelite.
plt.figure(1)
plt.plot(x,y)

# Cambio medidas de los ejes, para que sea más entendible.
ytks=[6000000,4000000,2000000,0,-2000000,-4000000,-6000000]
ytks_labels=[str(i/1000) for i in ytks]
xtks=[-7500000,-5000000,-2500000,0,2500000,5000000,7500000]
xtks_labels=[str(j/1000) for j in xtks]

# Le coloco nombre a los ejes.
plt.yticks(ytks,ytks_labels)
plt.xticks(xtks,xtks_labels,rotation=45)

plt.ylabel('(km)')
plt.xlabel('(km)')

"""
Graficas: Ayudan a visualisar 

1) El punto inicial y final. 
2) La atmosfera terrestre.
"""
# 1)
# Aqui realizo algo opcional y vizual, esto me permite ver que tan
# serca o lejos esta el satelite de los 80 km de atmosfera terrestre.

# Punto inicial: Cuadrado Verde
plt.plot([radio_tierra+700000,radio_tierra+700000],[0,0],'gs')

# Punto final: Triangulo Rojo
plt.plot([sol[len(sol)-1][0],sol[len(sol)-1][0]],[sol[len(sol)-1][1],sol[len(sol)-1][1]],'r^')


# 2) 

from matplotlib import pyplot as plt

# Creo un circulo que reprecenta la tierra (aproximacion)
num_segmentos = 20000
rad = 6371000+80000

# Me permite cambiar el centro del circulo
cx = 0
cy = 0

angulo = np.linspace(0, 2*np.pi, num_segmentos+1)
x1 = rad * np.cos(angulo) + cx
y1= rad * np.sin(angulo) + cy

# Grafico la "Tierra" como un circulo
plt.plot(x1, y1)
plt.legend(['Recorrido Satelite','Inicio','Fin','Tierra'],loc='upper right')

# Ordena y permite que las graficas no se deformen
plt.tight_layout()
plt.axis('equal')

nombre='Recorrido_Satelite_xy'
plt.savefig(f"Grafico {nombre}.png")

"""
Graficas: 
1) x(t) 
2) y(t) 
3) z(t) 
4) 
   4.1) r(t)
   4.2) Superficie de la Tierra
   4.3) Atmosfera
"""


plt.figure(2)

# 1)
plt.subplot(3,1,1)
plt.plot(t,x)

ytks_1=[5000000,0,-5000000]
ytks_labels_1=[str(i/1000) for i in ytks_1]
xtks_1=[0,2000,4000,6000,8000,10000,12000,14000]

plt.yticks(ytks_1,ytks_labels_1)
plt.xticks(xtks_1,[])

plt.title('x(t)')
plt.ylabel('(km)')

# 2)
plt.subplot(3,1,2)
plt.plot(t,y)

ytks_2=[5000000,0,-5000000]
ytks_labels_2=[str(i/1000) for i in ytks_2]
xtks_2=[0,2000,4000,6000,8000,10000,12000,14000]

plt.yticks(ytks_2,ytks_labels_2)
plt.xticks(xtks_2,[])

plt.title('y(t)')
plt.ylabel('(km)')

# 3)
plt.subplot(3,1,3)
plt.plot(t,z)

ytks_3=[0.05,0,-0.05]
ytks_labels_3=[str(i) for i in ytks_3]
xtks_3=[0,2000,4000,6000,8000,10000,12000,14000]

plt.yticks(ytks_3,ytks_labels_3)
plt.xticks(xtks_3,xtks_3,rotation=45)

plt.title('z(t)')
plt.ylabel('(m)')
plt.xlabel('Tiempo transcurrido (segs)')

nombre='Historias_xyz_Orbitas_Tiempo'
plt.savefig(f"Grafico {nombre}.png")

plt.tight_layout()

# 4)
plt.figure(3)

# 4.1)
plt.plot(t,np.sqrt(x**2+y**2+z**2))

ytks_4=[7000000,6500000]
ytks_labels_4=[str(i/1000) for i in ytks_4]
xtks_4=[0,2000,4000,6000,8000,10000,12000,14000]

plt.yticks(ytks_4,ytks_labels_4)
plt.xticks(xtks_4,xtks_4,rotation=45)

plt.title('r(t)')
plt.ylabel('(km)')
plt.xlabel('Tiempo transcurrido (segs)')

# 4.2)
plt.hlines(radio_tierra,0,segundos,color='blue')

# 4.3)
plt.hlines(rad,0,segundos,color='red')
plt.title('r(t)')

# Ordena y permite que las graficas no se deformen
plt.tight_layout()
nombre='Distancia_TS_vs_Tiempo'
plt.savefig(f"Grafico {nombre}.png")


"""
Preguntas a realizar:
    
1) Velocidad minima?
2) Como calcular el tiempo justo para las 2 vueltas???
"""


"""
pruebas para ver la velocidad limite
"""

segundos=3600*3
t = sp.linspace(0,segundos,segundos+1)

for j in range(24500,25000):
    vi = (j*1000)/3600
    # Defino mis coordenadas iniciales [xi,yi=0,z0=0,vx=0,vy=?,vz=0].
    z0 = sp.array([radio_tierra+700000 ,0,0,0,vi,0])
    
    # Utilizo la funcion odeint para resolver el sistema EDO.
    sol_p = odeint(bala,z0,t)
    
    xp = sol_p[:,0]
    yp = sol_p[:,1]
    zp = sol_p[:,2]
    
    lista=[np.sqrt((xp[int(i)])**2+(yp[int(i)])**2+(zp[int(i)])**2) for i in t]

    if min(lista)>=rad:
        print('La velocidad necesaria minima para rozar la atmosfera es:',j)
        break
    # if np.sqrt(xp[i]**2+y[i]**2+z[i]**2) <rad:
    #     print()
        