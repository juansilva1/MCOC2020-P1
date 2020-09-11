![Grafico 1](https://user-images.githubusercontent.com/69159364/91100102-932a5600-e632-11ea-84e0-177b4dc84e90.png)


# Entrega 2

*Nota1: los archivos pedidos se encuentran en la carpeta con el mismo nombre de esta entrega*

*Nota2: hay preguntas que estan constestadas en canvas ("ecuación en espacio de estado del sistema (escribala con el editor de ecuaciones)")*

### Contestar preguntas

1) Si se coloca el satélite a una altura de 700km (en la dirección x, ojo con el radio de la tierra) y se le da una velocidad inicial tangencial y(0) = vt ¿Cuanto debe valer vt de modo que el satélite efectivamente orbite sin caer de vuelta dentro de la atmosfera (asuma que esta comienza a una altura de 80km)?

**R:** Tomando en cuenta lo anterior, y asumiendo que el satelite da como minimo 2 vueltas a la Tierra (2 orbitas), se necesita una velosidad de **24551 [km/h] (6819.72 [m/s]).** Esto permite que el satelite de las 2 vueltas, sin caer de vuelta dentro de la atmosfera terrestre.

2) ¿Como encontró vt?

**R:** Primero al tanteo vi que con 25000 km/h las 2 órbitas ya no tocaban la Tierra (también grafiqué la Tierra con su atmósfera de 80 km), del mismo modo vi que con 24500 km/h todavía las órbitas tocaban la Tierra, por lo que teniendo este intervalo, proseguí a hacer un `for i in range (24500,25000)`, con esto formaba una lista que se formaba tomando en cuenta las 2 rotaciones y me entregaba la distancia del satélite a la tierra, luego obtenía el mínimo de la lista. Si este mínimo era menor al radio de la tierra más la atmósfera, entonces seguía adelante. así cuando obtenía la velocidad "justa", el ciclo de detenía y me arrojaba la ultima velocidad, en este caso fue **24551 [km/h].**

## Gráficos

### *Grafique y pegue en su respuesta las historias de tiempo de x(t), y(t), z(t) para dos órbitas completas.*

![Grafico Historias_xyz_Orbitas_Tiempo](https://user-images.githubusercontent.com/69159364/91487552-e80cdd00-e87b-11ea-9815-a21b75c3e127.png)

### *Grafique también la distancia al centro de la tierra del satélite vs. el tiempo, r(t), indicando la superficie de la tierra y la atmósfera como una línea de distancia constante.*

*Nota: la linea horizontal roja es la atmosfera, la linea horizontal azul es la Tierra*

![Grafico Distancia_TS_vs_Tiempo](https://user-images.githubusercontent.com/69159364/91487992-a3ce0c80-e87c-11ea-8f2b-4f38704748bd.png)

### *Extra: Representación de la tierra con su atmósfera y más de 2 órbitas*

![Grafico Recorrido_Satelite_xy](https://user-images.githubusercontent.com/69159364/91488016-ad577480-e87c-11ea-969b-7c6040493a80.png)

**Nota:** ***Estos graficos toman en cuenta la velocidad limite, es por esto que a simple vista perece como si se estuvieran tocando, si se hace zoom se puede apreciar que el satelite no llega a tocar a la atmosfera. Si se desea apresiar esta diferencia a simple vista recomiendo escoger una velocidad mayor, por ejemplo 25000 [km/h]***

*El codigo imprime al final la velocidad de **24551 [km/h]***

# Entrega 3

No se pidio nada de comentarios en esta entrega, pero aquí van unos extra:

Diferencia en metros tomando en cuenta tiempo UTC es de: **525005.8656269278 [m]**

*Nota: el programa se encuentra en la carpeta "Entrega3".*


# Entrega 4

## Análisis

Se puede apreciar que la función analítica, es muy similar a la función odeint, esto ocurre pues esta función "automáticamente" toma subdivisiones, lo que hace que el resultado aproximado sea muy cercano al real.

Por otro lado, las funciones euleint son aproximaciones que dependen del número de subdivisiones que posea el intervalo, si este es muy pequeño, el resultado es "desastroso", no se asemeja a la función analítica, pero si el número de subdivisiones aumenta los valores son cada vez más cercanos al original, por lo que el error se corrige significativamente. Se puede observar que para un Nsubdivisiones=100 el resultado al menos a simple vista resulta satisfactorio.

![Grafico Oscilador_Armónico](https://user-images.githubusercontent.com/69159364/91896708-4ca6ae00-ec67-11ea-8221-c7c645fe4e93.png)

*Nota: el programa se encuentra en la carpeta "Entrega4".*




# Entrega5

_Nota: los títulos de algunos gráficos se me cortaron por el largo de los mismos._

## Parte 1
**_En este gráfico se presentan las 3 coordenadas x,y,z con respecto al tiempo. Están consideradas solamente las coordenadas reales del satélite, no las obtenidas con odeint ni con eulerint._**

![Pregunta1_Gráficos_posición_vs_tiempo](https://user-images.githubusercontent.com/69159364/92931278-9c306b00-f419-11ea-8321-e03f8827177b.png)

## Pate 2


**_Para esta parte realice varios gráficos extras que ayudan a comprender de mejor forma el funcionamiento de las funciones odeint y eulerint._**

### Gráfico odeint vs real

![Pregunta2_Gráfica_deriva_odeint_vs_real](https://user-images.githubusercontent.com/69159364/92931516-f2051300-f419-11ea-8a59-f30e93fac24a.png)

### Gráfico eulerint vs real

![Pregunta2_Gráfica_deriva_eulerint_vs_real](https://user-images.githubusercontent.com/69159364/92931511-f0d3e600-f419-11ea-81c3-5dc92ecfb941.png)

### Gráfico odeint vs eulerint

![Pregunta2_Gráfica_deriva_odeint_vs_eulerint](https://user-images.githubusercontent.com/69159364/92931514-f2051300-f419-11ea-99d8-d57edee32230.png)

**R:** La deriva entre odeint y eulerint en el punto final es de 19723.989 kilómetros

### Gráfico de barras, me indica el tiempo de odeint y eulerint

**_Este gráfico sí que es extra, pero resulta más cómodo de leer y revisar que solo ver números._**

![Pregunta2_Gráfica_demora_odeint_vs_eulerint](https://user-images.githubusercontent.com/69159364/92932580-760bca80-f41b-11ea-9215-dddd7ed99d9f.png)

**R:** en todas las corridas realizadas odeint es más rápido que eulerint y además arroja resultados más cercanos a la realidad. Tal como aparece en el gráfico de barras anterior, el tiempo de ejecución de **odeint es de 0.57 (s)** y el de **eulerint es de 1.117 (s).** // obviamente varia si se vuelve a correr, pero siempre se ve esta relación de mayor y menor (tiempo de eulerint>tiempo de odeint).

## Parte 3

**R:** Se tienen que utilizar más de 10000 Nsuvdiviciones para que el error sea menor al 1%.
El tiempo de ejecución de eulerint es extremadamente lento para un número de Nsubdivisiones muy grandes, yo creé un método que revisaba los primeros 10 intervalos de tiempo y si este ya era diferente de 1% lo omitía y seguía con el siguiente, así fue del 100 al 10000 (de 100 en 100) y el tiempo de ejecución solo para los 10 intervalos era ya de unos 5-10 segs en los números más grandes, así que si en vez de 10 intervalos utilizo los 9361, el programa se demoraría demasiado tiempo. 

**A continuación coloco el gráfico para Nsubdiviciones=1000**

_Nota: esté lo tuve que correr aparte, a, solo lo corrí para 1000, no hice desde 100 hasta 10000, si se desea obtener este gráfico se tiene que modificar la parte 3 del código para que solo arroje un gráfico y además se tendrían que arreglar sus límites (del gráfico)._

![Pregunta3_Gráfica_deriva_eulerint_vs_real_con error_menor_al_1_porciento](https://user-images.githubusercontent.com/69159364/92932582-773cf780-f41b-11ea-8b87-f6450fb302ca.png)

Resumen: eulerint es bastante ineficiente. A pesar de lo anterior igual redujo significativamente el error, con 1000 subdivisiones, probablemente con un número mucho mayor de seguro se lograría el objetivo, pero la idea no es dejar corriendo el computador todo un día.

## Parte 4

**_Aquí también realice varios gráficos de más, considere los casos solo con J2, solo con J3 y con ambos_** 

Como resumen J2 es el que hace una mayor corrección de los 2.

**Colocaré solo los 4 gráficos de J2 + J3**


### Gráfico odeint vs real utilizando J2 y J3

![Pregunta4_Gráfica_J2_y_J3_deriva_odeint_vs_real](https://user-images.githubusercontent.com/69159364/92937469-bec68200-f421-11ea-9e05-74e4cb9acd8a.png)

### Gráfico eulerint vs real utilizando J2 y J3

![Pregunta4_Gráfica_J2_y_J3_deriva_eulerint_vs_real](https://user-images.githubusercontent.com/69159364/92933166-4a3d1480-f41c-11ea-8f43-a372c8064e8d.png)

### Gráfico odeint vs eulerint utilizando J2 y J3

![Pregunta4_Gráfica_J2_y_J3_deriva_odeint_vs_eulerint](https://user-images.githubusercontent.com/69159364/92933168-4ad5ab00-f41c-11ea-9d46-22219cd2e4f7.png)

**R:** La deriva entre odeint y eulerint en el punto final es de 19408.0368 kilómetros. Experimento una reducción de más de 300 kilómetros.

### Gráfico de barras, tiempo de ejecución de odeint y eulerint utilizando J2 y J3

![Pregunta4_Gráfica_J2_y_J3_demora_odeint_vs_eulerint](https://user-images.githubusercontent.com/69159364/92933163-49a47e00-f41c-11ea-9537-46b047f0f384.png)

**R:** Nuevamente odeint es significativamente más rápido que eulerint. en este gráfico el odeint se demora **0.742(s)** y el eulerint se demora **1.751(s).**

==> A pesar de que los tiempos aumentan, la precicion tambien lo hace aumenta.

_Nota: se se desea ver el resto de gráficos, ingresar en a la carpeta "Entrega5", en esta se encuentra todo lo relacionado con esta entrega (archivos .py, .EOF y .png)_




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


**_Nota: si se desea ver el resto de archivos (.EOF y "leer_eof.py"), estos se encuentran en la carpeta "Entrega6."_**



# Entrega 7

## Mejoras realizadas

1) la primera fue corregir la entrega anterior, en la cual estaba utilizando mal los términos J2 y J3.

2) mejore la precisión de los términos J2 y J3.

3) Aumenté la velocidad al eliminar uno de los ciclos.

4) Implemente hasta el J8, lo que redujo significativamente la distancia entre los puntos reales y finales del satelite.

### Triste historia

*_Si desea continuar con la historia ver la carpeta "Entrega7,"_*

### Historía semi feliz

*_Si desea continuar con la historia ver la carpeta "Entrega7,"_*

*_Nota: Codigo se encuentra tambien en la carpeta "Entrega7."_*

