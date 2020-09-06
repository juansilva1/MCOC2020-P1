# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 07:04:34 2020

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

t, x, y, z, vx, vy, vz = leer_eof("S1A_OPER_AUX_POEORB_OPOD_20200814T120815_V20200724T225942_20200726T005942 (2).EOF")

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
PARTE a) Graficas de x(t), y(t), z(t)
"""

# Gráfico x(t)
plt.figure(1)
plt.subplot(3,1,1)
plt.plot(t/3600,x)
ytks_1=[5000000,-5000000]
ytks_labels_1=[str(i/1000) for i in ytks_1]
plt.yticks(ytks_1,ytks_labels_1)
plt.ylabel('X (KM)')
plt.title('Posición')

# Gráfico y(t)
plt.subplot(3,1,2)
plt.plot(t/3600,y)
ytks_2=[5000000,-5000000]
ytks_labels_2=[str(i/1000) for i in ytks_2]
plt.yticks(ytks_2,ytks_labels_2)
plt.ylabel('Y (KM)')

# Gráfico z(t)
plt.subplot(3,1,3)
plt.plot(t/3600,z)
ytks_3=[5000000,-5000000]
ytks_labels_3=[str(i/1000) for i in ytks_3]
plt.yticks(ytks_3,ytks_labels_3)
plt.xlabel('Tiempo, t(horas)')
plt.ylabel('Z (KM)')


plt.tight_layout()

nombre='Pregunta1_Gráficos_posición_vs_tiempo'
plt.savefig(f"{nombre}.png")


"""
PARTE b) Graficas de deriva sin J2 ni J3
"""

z0 = vector_inicial_real

# Me defino la función eulerint
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

# la utilizaré para contar el tiempo que se demora cada sistema (odeint y eulerint)
from time import perf_counter

# Tiempo de odeint sin mejoras
t1_odeint = perf_counter()
sol = odeint(bala,z0,t)
t2_odeint = perf_counter()

dt_odeint = t2_odeint - t1_odeint

# Tiempo de eulerint sin mejoras
t1_eulerint = perf_counter()
sol_1 = eulerint(bala, z0, t,Nsubdivisiones=1)
t2_eulerint = perf_counter()

dt_eulerint = t2_eulerint - t1_eulerint


def diferencia_vectores(sol,v_real):
    
    v_p_xyz=sol[:3]
    v_r_xyz=v_real[:3]
    
    distancia=np.sqrt((v_p_xyz[0]-v_r_xyz[0])**2+(v_p_xyz[1]-v_r_xyz[1])**2+(v_p_xyz[2]-v_r_xyz[2])**2)
    
    return distancia

#Caso comparacion de distancias entre odeint y real
lista_deriva=[]
c=0
for i in sol[:,0]:
    te=t[0:c+1]
    aa=diferencia_vectores([i,sol[:,1][c],sol[:,2][c]],[x[c],y[c],z[c],vx[c],vy[c],vz[c]])
    lista_deriva.append(aa)
    c+=1
plt.figure(2)
plt.plot(t/3600,lista_deriva)
ytks_3=[0,100000,200000,300000,400000,500000]
ytks_labels_3=[str(i/1000) for i in ytks_3]
plt.yticks(ytks_3,ytks_labels_3)
plt.xlabel('Tiempo, t(horas)')
plt.ylabel('Deriva, d (KM)')
plt.title(f'Distancia entre posicion real y predicha por odeint dmax = {lista_deriva[-1]/1000}')
plt.tight_layout()
nombre='Pregunta2_Gráfica_deriva_odeint_vs_real'
plt.savefig(f"{nombre}.png")

#Caso comparacion de distancias entre eulerint y real
lista_deriva=[]
c=0
for i in sol_1[:,0]:
    te=t[0:c+1]
    aa=diferencia_vectores([i,sol_1[:,1][c],sol_1[:,2][c]],[x[c],y[c],z[c],vx[c],vy[c],vz[c]])
    lista_deriva.append(aa)
    c+=1

plt.figure(3)
plt.plot(t/3600,lista_deriva)
ytks_3=[0,0.5*10000000,1*10000000,1.5*10000000,2.0*10000000]
ytks_labels_3=[str(i/1000) for i in ytks_3]
plt.yticks(ytks_3,ytks_labels_3)
plt.xlabel('Tiempo, t(horas)')
plt.ylabel('Deriva, d (KM)')
plt.title(f'Distancia entre posicion real y predicha por eulerint dmax = {lista_deriva[-1]/1000}')
plt.tight_layout()
nombre='Pregunta2_Gráfica_deriva_eulerint_vs_real'
plt.savefig(f"{nombre}.png")

#Caso comparacion de distancias entre odeint y eulerint
lista_deriva=[]
c=0
for i in sol_1[:,0]:
    te=t[0:c+1]
    aa=diferencia_vectores([i,sol_1[:,1][c],sol_1[:,2][c]],[sol[:,0][c],sol[:,1][c],sol[:,2][c],sol[:,3][c],sol[:,4][c],sol[:,5][c]])
    lista_deriva.append(aa)
    c+=1

plt.figure(4)
plt.plot(t/3600,lista_deriva)
ytks_3=[0,0.5*10000000,1*10000000,1.5*10000000,2.0*10000000]
ytks_labels_3=[str(i/1000) for i in ytks_3]
plt.yticks(ytks_3,ytks_labels_3)
plt.xlabel('Tiempo, t(horas)')
plt.ylabel('Deriva, d (KM)')
plt.title(f'Distancia entre posicion predicha por odeint y eulerint dmax = {lista_deriva[-1]/1000}')
plt.tight_layout()
nombre='Pregunta2_Gráfica_deriva_odeint_vs_eulerint'
plt.savefig(f"{nombre}.png")


# Muy extra, entrega un gráfico de barras para indicar la diferencia entre los
# tiempos obtenidos por cada función.
plt.figure(5)
plt.bar(['odeint','eulerint'],[dt_odeint,dt_eulerint],width=0.3,align='center')
plt.legend([f'odeint = {np.round(dt_odeint,3)}(s) \neulerint = {np.round(dt_eulerint,3)}(s)'])

plt.ylabel('Tiempo de ejecución (s)')
plt.title('Tiempos de ejecución odeint vs eulerint')
plt.tight_layout()

nombre='Pregunta2_Gráfica_demora_odeint_vs_eulerint'
plt.savefig(f"{nombre}.png")
plt.show()


"""
PARTE c) Numero de Nsubdivisiones para que el error de euler int sea menor al 1%
"""
# Si se quiere ver el procedimiento, descomentar las lineas 329 a la 379
"""
#Caso comparacion de distancias entre eulerint y real 

# sol_error = eulerint(bala, z0, t,Nsubdivisiones=1)

# Error:
def error(sol_error,si):
    # Error en porcentaje compara la distancia final con 0, la distancia exacta
    error_x=((sol_error[:,0][-1]-vector_final_real[0])/x[si])*100 
    error_y=((sol_error[:,1][-1]-vector_final_real[1])/y[si])*100 
    error_z=((sol_error[:,2][-1]-vector_final_real[2])/z[si])*100 
    promedio_errores=(error_x+error_y+error_z)/3
    if  promedio_errores>1:
        return False
    else:
        return True 
print("tiempo",len(t))
for i in range(51,101):
    print("paso numero",i*100)
    sol_error_prueba=eulerint(bala, z0, t[:10],i*100)
    if not error(sol_error_prueba,9):
        print('sol_error_prueba =',i*100)
        continue
    
    sol_error=eulerint(bala, z0, t,i*100)
    if error(sol_error,-1):
        print(f"Lo lograste, necesitas {i*1000} Nsubdiviciones para que el error sea menor al 1%")
        lista_deriva=[]
        c=0
        for i in sol_error[:,0]:
            te=t[0:c+1]
            aa=diferencia_vectores([i,sol_error[:,1][c],sol_error[:,2][c]],[x[c],y[c],z[c],vx[c],vy[c],vz[c]])
            lista_deriva.append(aa)
            c+=1
        plt.figure(6)
        plt.plot(t/3600,lista_deriva)
        ytks_3=[0,0.5*10000000,1*10000000,1.5*10000000,2.0*10000000]
        ytks_labels_3=[str(i/1000) for i in ytks_3]
        plt.yticks(ytks_3,ytks_labels_3)
        plt.xlabel('Tiempo, t(horas)')
        plt.ylabel('Deriva, d (KM)')
        plt.title(f'Distancia entre posicion real y predicha por eulerint dmax = {lista_deriva[-1]/1000} < 1%')
        plt.tight_layout()
        nombre='Pregunta3_Gráfica_deriva_eulerint_vs_real_con error_menor_al_1_porciento'
        plt.savefig(f"{nombre}.png")
        break
        
    else:
        continue

"""
"""
PARTE d) Implementando J2 y J3
"""

"""
Implementando solamente J2
"""

# Esta funcion se encarga de relacionar las variables del problema
def bala_J2(z,t,Ignorar_mejoras=1,Ignorar_J2=1,Ignorar_J3=0):
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

# la utilizaré para contar el tiempo que se demora cada sistema (odeint y eulerint)
from time import perf_counter

# Tiempo de odeint sin mejoras
t1_odeint = perf_counter()
sol = odeint(bala_J2,z0,t)
t2_odeint = perf_counter()

dt_odeint = t2_odeint - t1_odeint

# Tiempo de eulerint sin mejoras
t1_eulerint = perf_counter()
sol_1 = eulerint(bala_J2, z0, t,Nsubdivisiones=1)
t2_eulerint = perf_counter()

dt_eulerint = t2_eulerint - t1_eulerint


#Caso comparacion de distancias entre odeint y real
lista_deriva=[]
c=0
for i in sol[:,0]:
    te=t[0:c+1]
    aa=diferencia_vectores([i,sol[:,1][c],sol[:,2][c]],[x[c],y[c],z[c],vx[c],vy[c],vz[c]])
    lista_deriva.append(aa)
    c+=1
plt.figure(7)
plt.plot(t/3600,lista_deriva)
ytks_3=[0,100000,200000,300000,400000,500000]
ytks_labels_3=[str(i/1000) for i in ytks_3]
plt.yticks(ytks_3,ytks_labels_3)
plt.xlabel('Tiempo, t(horas)')
plt.ylabel('Deriva, d (KM)')
plt.title(f'Distancia entre posicion real y predicha por odeint dmax = {lista_deriva[-1]/1000}')
plt.tight_layout()
nombre='Pregunta4_Gráfica_J2_deriva_odeint_vs_real'
plt.savefig(f"{nombre}.png")

#Caso comparacion de distancias entre eulerint y real
lista_deriva=[]
c=0
for i in sol_1[:,0]:
    te=t[0:c+1]
    aa=diferencia_vectores([i,sol_1[:,1][c],sol_1[:,2][c]],[x[c],y[c],z[c],vx[c],vy[c],vz[c]])
    lista_deriva.append(aa)
    c+=1

plt.figure(8)
plt.plot(t/3600,lista_deriva)
ytks_3=[0,0.5*10000000,1*10000000,1.5*10000000,2.0*10000000]
ytks_labels_3=[str(i/1000) for i in ytks_3]
plt.yticks(ytks_3,ytks_labels_3)
plt.xlabel('Tiempo, t(horas)')
plt.ylabel('Deriva, d (KM)')
plt.title(f'Distancia entre posicion real y predicha por eulerint dmax = {lista_deriva[-1]/1000}')
plt.tight_layout()
nombre='Pregunta4_Gráfica_J2_deriva_eulerint_vs_real'
plt.savefig(f"{nombre}.png")

#Caso comparacion de distancias entre odeint y eulerint
lista_deriva=[]
c=0
for i in sol_1[:,0]:
    te=t[0:c+1]
    aa=diferencia_vectores([i,sol_1[:,1][c],sol_1[:,2][c]],[sol[:,0][c],sol[:,1][c],sol[:,2][c],sol[:,3][c],sol[:,4][c],sol[:,5][c]])
    lista_deriva.append(aa)
    c+=1

plt.figure(9)
plt.plot(t/3600,lista_deriva)
ytks_3=[0,0.5*10000000,1*10000000,1.5*10000000,2.0*10000000]
ytks_labels_3=[str(i/1000) for i in ytks_3]
plt.yticks(ytks_3,ytks_labels_3)
plt.xlabel('Tiempo, t(horas)')
plt.ylabel('Deriva, d (KM)')
plt.title(f'Distancia entre posicion predicha por odeint y eulerint dmax = {lista_deriva[-1]/1000}')
plt.tight_layout()
nombre='Pregunta4_Gráfica_J2_deriva_odeint_vs_eulerint'
plt.savefig(f"{nombre}.png")


# Muy extra, entrega un gráfico de barras para indicar la diferencia entre los
# tiempos obtenidos por cada función.
plt.figure(10)
plt.bar(['odeint','eulerint'],[dt_odeint,dt_eulerint],width=0.3,align='center')
plt.legend([f'odeint = {np.round(dt_odeint,3)}(s) \neulerint = {np.round(dt_eulerint,3)}(s)'])

plt.ylabel('Tiempo de ejecución (s)')
plt.title('Tiempos de ejecución odeint vs eulerint')
plt.tight_layout()

nombre='Pregunta4_Gráfica_J2_demora_odeint_vs_eulerint'
plt.savefig(f"{nombre}.png")
plt.show()



"""
Implementando solamente J3
"""

# Esta funcion se encarga de relacionar las variables del problema
def bala_J3(z,t,Ignorar_mejoras=1,Ignorar_J2=0,Ignorar_J3=1):
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

# la utilizaré para contar el tiempo que se demora cada sistema (odeint y eulerint)
from time import perf_counter

# Tiempo de odeint sin mejoras
t1_odeint = perf_counter()
sol = odeint(bala_J3,z0,t)
t2_odeint = perf_counter()

dt_odeint = t2_odeint - t1_odeint

# Tiempo de eulerint sin mejoras
t1_eulerint = perf_counter()
sol_1 = eulerint(bala_J3, z0, t,Nsubdivisiones=1)
t2_eulerint = perf_counter()

dt_eulerint = t2_eulerint - t1_eulerint


#Caso comparacion de distancias entre odeint y real
lista_deriva=[]
c=0
for i in sol[:,0]:
    te=t[0:c+1]
    aa=diferencia_vectores([i,sol[:,1][c],sol[:,2][c]],[x[c],y[c],z[c],vx[c],vy[c],vz[c]])
    lista_deriva.append(aa)
    c+=1
plt.figure(11)
plt.plot(t/3600,lista_deriva)
ytks_3=[0,100000,200000,300000,400000,500000]
ytks_labels_3=[str(i/1000) for i in ytks_3]
plt.yticks(ytks_3,ytks_labels_3)
plt.xlabel('Tiempo, t(horas)')
plt.ylabel('Deriva, d (KM)')
plt.title(f'Distancia entre posicion real y predicha por odeint dmax = {lista_deriva[-1]/1000}')
plt.tight_layout()
nombre='Pregunta4_Gráfica_J3_deriva_odeint_vs_real'
plt.savefig(f"{nombre}.png")

#Caso comparacion de distancias entre eulerint y real
lista_deriva=[]
c=0
for i in sol_1[:,0]:
    te=t[0:c+1]
    aa=diferencia_vectores([i,sol_1[:,1][c],sol_1[:,2][c]],[x[c],y[c],z[c],vx[c],vy[c],vz[c]])
    lista_deriva.append(aa)
    c+=1

plt.figure(12)
plt.plot(t/3600,lista_deriva)
ytks_3=[0,0.5*10000000,1*10000000,1.5*10000000,2.0*10000000]
ytks_labels_3=[str(i/1000) for i in ytks_3]
plt.yticks(ytks_3,ytks_labels_3)
plt.xlabel('Tiempo, t(horas)')
plt.ylabel('Deriva, d (KM)')
plt.title(f'Distancia entre posicion real y predicha por eulerint dmax = {lista_deriva[-1]/1000}')
plt.tight_layout()
nombre='Pregunta4_Gráfica_J3_deriva_eulerint_vs_real'
plt.savefig(f"{nombre}.png")

#Caso comparacion de distancias entre odeint y eulerint
lista_deriva=[]
c=0
for i in sol_1[:,0]:
    te=t[0:c+1]
    aa=diferencia_vectores([i,sol_1[:,1][c],sol_1[:,2][c]],[sol[:,0][c],sol[:,1][c],sol[:,2][c],sol[:,3][c],sol[:,4][c],sol[:,5][c]])
    lista_deriva.append(aa)
    c+=1

plt.figure(13)
plt.plot(t/3600,lista_deriva)
ytks_3=[0,0.5*10000000,1*10000000,1.5*10000000,2.0*10000000]
ytks_labels_3=[str(i/1000) for i in ytks_3]
plt.yticks(ytks_3,ytks_labels_3)
plt.xlabel('Tiempo, t(horas)')
plt.ylabel('Deriva, d (KM)')
plt.title(f'Distancia entre posicion predicha por odeint y eulerint dmax = {lista_deriva[-1]/1000}')
plt.tight_layout()
nombre='Pregunta4_Gráfica_J3_deriva_odeint_vs_eulerint'
plt.savefig(f"{nombre}.png")


# Muy extra, entrega un gráfico de barras para indicar la diferencia entre los
# tiempos obtenidos por cada función.
plt.figure(14)
plt.bar(['odeint','eulerint'],[dt_odeint,dt_eulerint],width=0.3,align='center')
plt.legend([f'odeint = {np.round(dt_odeint,3)}(s) \neulerint = {np.round(dt_eulerint,3)}(s)'])

plt.ylabel('Tiempo de ejecución (s)')
plt.title('Tiempos de ejecución odeint vs eulerint')
plt.tight_layout()

nombre='Pregunta4_Gráfica_J3_demora_odeint_vs_eulerint'
plt.savefig(f"{nombre}.png")
plt.show()


"""
Implementando J2 y J3
"""

# Esta funcion se encarga de relacionar las variables del problema
def bala_J2_y_J3(z,t,Ignorar_mejoras=1,Ignorar_J2=1,Ignorar_J3=1):
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

# la utilizaré para contar el tiempo que se demora cada sistema (odeint y eulerint)
from time import perf_counter

# Tiempo de odeint sin mejoras
t1_odeint = perf_counter()
sol = odeint(bala_J2_y_J3,z0,t)
t2_odeint = perf_counter()

dt_odeint = t2_odeint - t1_odeint

# Tiempo de eulerint sin mejoras
t1_eulerint = perf_counter()
sol_1 = eulerint(bala_J2_y_J3, z0, t,Nsubdivisiones=1)
t2_eulerint = perf_counter()

dt_eulerint = t2_eulerint - t1_eulerint


#Caso comparacion de distancias entre odeint y real
lista_deriva=[]
c=0
for i in sol[:,0]:
    te=t[0:c+1]
    aa=diferencia_vectores([i,sol[:,1][c],sol[:,2][c]],[x[c],y[c],z[c],vx[c],vy[c],vz[c]])
    lista_deriva.append(aa)
    c+=1
plt.figure(15)
plt.plot(t/3600,lista_deriva)
ytks_3=[0,100000,200000,300000,400000,500000]
ytks_labels_3=[str(i/1000) for i in ytks_3]
plt.yticks(ytks_3,ytks_labels_3)
plt.xlabel('Tiempo, t(horas)')
plt.ylabel('Deriva, d (KM)')
plt.title(f'Distancia entre posicion real y predicha por odeint dmax = {lista_deriva[-1]/1000}')
plt.tight_layout()
nombre='Pregunta4_Gráfica_J2_y_J3_deriva_odeint_vs_real'
plt.savefig(f"{nombre}.png")

#Caso comparacion de distancias entre eulerint y real
lista_deriva=[]
c=0
for i in sol_1[:,0]:
    te=t[0:c+1]
    aa=diferencia_vectores([i,sol_1[:,1][c],sol_1[:,2][c]],[x[c],y[c],z[c],vx[c],vy[c],vz[c]])
    lista_deriva.append(aa)
    c+=1

plt.figure(16)
plt.plot(t/3600,lista_deriva)
ytks_3=[0,0.5*10000000,1*10000000,1.5*10000000,2.0*10000000]
ytks_labels_3=[str(i/1000) for i in ytks_3]
plt.yticks(ytks_3,ytks_labels_3)
plt.xlabel('Tiempo, t(horas)')
plt.ylabel('Deriva, d (KM)')
plt.title(f'Distancia entre posicion real y predicha por eulerint dmax = {lista_deriva[-1]/1000}')
plt.tight_layout()
nombre='Pregunta4_Gráfica_J2_y_J3_deriva_eulerint_vs_real'
plt.savefig(f"{nombre}.png")

#Caso comparacion de distancias entre odeint y eulerint
lista_deriva=[]
c=0
for i in sol_1[:,0]:
    te=t[0:c+1]
    aa=diferencia_vectores([i,sol_1[:,1][c],sol_1[:,2][c]],[sol[:,0][c],sol[:,1][c],sol[:,2][c],sol[:,3][c],sol[:,4][c],sol[:,5][c]])
    lista_deriva.append(aa)
    c+=1

plt.figure(17)
plt.plot(t/3600,lista_deriva)
ytks_3=[0,0.5*10000000,1*10000000,1.5*10000000,2.0*10000000]
ytks_labels_3=[str(i/1000) for i in ytks_3]
plt.yticks(ytks_3,ytks_labels_3)
plt.xlabel('Tiempo, t(horas)')
plt.ylabel('Deriva, d (KM)')
plt.title(f'Distancia entre posicion predicha por odeint y eulerint dmax = {lista_deriva[-1]/1000}')
plt.tight_layout()
nombre='Pregunta4_Gráfica_J2_y_J3_deriva_odeint_vs_eulerint'
plt.savefig(f"{nombre}.png")


# Muy extra, entrega un gráfico de barras para indicar la diferencia entre los
# tiempos obtenidos por cada función.
plt.figure(18)
plt.bar(['odeint','eulerint'],[dt_odeint,dt_eulerint],width=0.3,align='center')
plt.legend([f'odeint = {np.round(dt_odeint,3)}(s) \neulerint = {np.round(dt_eulerint,3)}(s)'])

plt.ylabel('Tiempo de ejecución (s)')
plt.title('Tiempos de ejecución odeint vs eulerint')
plt.tight_layout()

nombre='Pregunta4_Gráfica_J2_y_J3_demora_odeint_vs_eulerint'
plt.savefig(f"{nombre}.png")
plt.show()


# Si se quiere ver el procedimiento de la letra c), descomentar las lineas 329 a la 379 /(tarda muchisimo))