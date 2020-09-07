
# Entrega6

## Código del programa utilizado, "Modo concurso".

```

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
radio_tierra=6371000 #m
masa_tierra=(5.97219*(10**24)) #kg
G=(6.67392*(10**-11))# m3/s2kg
omega=(7.27*(0.00001))# rad/s
r=10


# algunas funciones utiles para el futuro:

# Esta funcion recibe un vector posicion y retorna la distancia r.
def distancia_satelite_tierra(vector_posicion):
    r=(np.dot(vector_posicion,vector_posicion))**0.5
    return r

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

# J2
def T_J2(x,y,z,r):
   
    J2=1.75553*(10**10)*(1000**5)
    
    Fx=J2 * (x/(r**7)) * (6*(z**2) - 1.5*((x**2) + (y**2)))
    Fy=J2 * (y/(r**7)) * (6*(z**2) - 1.5*((x**2) + (y**2)))
    Fz=J2 * (z/(r**7)) * (3*(z**2) - 4.5*((x**2) + (y**2)))
    
    J2_final=[Fx,Fy,Fz]
    
    return np.array(J2_final)

# J3
def T_J3(x,y,z,r):

    J3= -2.61913*(10**11)*(1000**6)
    
    Fx=J3 * ((x*z)/(r**9)) * (10*(z**2) - 7.5*((x**2) + (y**2)))
    Fy=J3 * ((y*z)/(r**9)) * (10*(z**2) - 7.5*((x**2) + (y**2)))
    Fz=J3 * ((1)/(r**9)) * (4*(z**2)*((z**2)-3*((x**2) + (y**2))) + 1.5*(((x**2) + (y**2))**2))
    
    J3_final=[Fx,Fy,Fz]
    
    return np.array(J3_final)   


# Esta función me entrega el potencial más exacto
def potencial(x,y,z,r,Ignorar_J2=0,Ignorar_J3=0):

    J2=T_J2(x,y,z,r)*Ignorar_J2
    J3=T_J3(x,y,z,r)*Ignorar_J3

    potencial_f=-1*G*(masa_tierra/r)*(np.ones(3)+J2+J3)
    potencial_f=(J2+J3)

    return potencial_f

# Esta funcion se encarga de relacionar las variables del problema
def bala(z,t,Ignorar_mejoras=0,Ignorar_J2=0,Ignorar_J3=0):
    r=distancia_satelite_tierra(z[0:3])
    
    # Creo un vector zp de 6x1
    zp = np.zeros(6)

    # Los primeros 3 "espacios" de zp son los ultimos de z.
    zp[0] = z[3]
    zp[1] = z[4]
    zp[2] = z[5]
    
    # Parte potencial
    # Aplicando correcciones a los calculos
    potenciall=potencial(z[0],z[1],z[2],r,Ignorar_J2,Ignorar_J3)*Ignorar_mejoras
    
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
    
    # Parte final = Parte 1 + Parte 2 + Parte potencial
    parte3_final=parte1_final+parte2_final+potenciall

    # Relaciono los 3 ultimos "espacios" de zp con el vector de 3x1 resultante de la Parte final.
    zp[3:6]=parte3_final 
    
    return zp

# Me entrega la distancia entre el vector final real y el de la predicción.
def diferencia_real_predicha(sol,v_real):
    
    v_p_xyz=sol[-1][:3]
    v_r_xyz=v_real[:3]
    
    distancia=np.sqrt((v_p_xyz[0]-v_r_xyz[0])**2+(v_p_xyz[1]-v_r_xyz[1])**2+(v_p_xyz[2]-v_r_xyz[2])**2)
    
    return distancia


"""
PARTE b) Graficas de deriva sin J2 ni J3
"""

# z0 = vector_inicial_real
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
    lista_datos.append(cambio)
archivo.close()

string_name=str(fname)
archivo_final=string_name.replace('.EOF','.PRED')

archivo1=open(archivo_final,'w')
for i in lista_datos:
    archivo1.write(i)
archivo1.close()

```
