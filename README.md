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

![Pregunta1_Gráficos_posición_vs_tiempo](https://user-images.githubusercontent.com/69159364/92327729-37989900-f032-11ea-97c8-44d88956a955.png)

## Pate 2

**_Para esta parte realice varios gráficos extras que ayudan a comprender de mejor forma el funcionamiento de las funciones odeint y eulerint._**

### Gráfico odeint vs real

![Pregunta2_Gráfica_deriva_odeint_vs_real](https://user-images.githubusercontent.com/69159364/92327906-521f4200-f033-11ea-8b3c-1f3ad02c16dc.png)

### Gráfico eulerint vs real

![Pregunta2_Gráfica_deriva_eulerint_vs_real](https://user-images.githubusercontent.com/69159364/92327904-5186ab80-f033-11ea-8aa8-266e30cb1257.png)

### Gráfico odeint vs eulerint

![Pregunta2_Gráfica_deriva_odeint_vs_eulerint](https://user-images.githubusercontent.com/69159364/92327905-521f4200-f033-11ea-8079-0dd868a96e86.png)

**R:** La deriva entre odeint y eulerint en el punto final es de 19666.822 kilómetros

### Gráfico de barras, me indica el tiempo de odeint y eulerint

**_Este gráfico sí que es extra, pero resulta más cómodo de leer y revisar que solo ver números._**

![Pregunta2_Gráfica_demora_odeint_vs_eulerint](https://user-images.githubusercontent.com/69159364/92327903-50ee1500-f033-11ea-8525-d3a1e95ad601.png)

**R:** en todas las corridas realizadas odeint es más rápido que eulerint y además arroja resultados más cercanos a la realidad. Tal como aparece en el gráfico de barras anterior, el tiempo de ejecución de **odeint es de 0.362 (s)** y el de **eulerint es de 0,914 (s).** // obviamente varia si se vuelve a correr, pero siempre se ve esta relación de mayor y menor (tiempo de eulerint>tiempo de odeint).

## Parte 3

**R:** Se tienen que utilizar más de 10000 Nsuvdiviciones para que el error sea menor al 1%.
El tiempo de ejecución de eulerint es extremadamente lento para un número de Nsubdivisiones muy grandes, yo creé un método que revisaba los primeros 10 intervalos de tiempo y si este ya era diferente de 1% lo omitía y seguía con el siguiente, así fue del 100 al 10000 (de 100 en 100) y el tiempo de ejecución solo para los 10 intervalos era ya de unos 5-10 segs en los números más grandes, así que si en vez de 10 intervalos utilizo los 9361, el programa se demoraría demasiado tiempo. 

**A continuación coloco el gráfico para Nsubdiviciones=1000**

_Nota: esté lo tuve que correr aparte, a, solo lo corrí para 1000, no hice desde 100 hasta 10000, si se desea obtener este gráfico se tiene que modificar la parte 3 del código para que solo arroje un gráfico y además se tendrían que arreglar sus límites (del gráfico)._

![Pregunta3_Gráfica_deriva_eulerint_vs_real_con error_menor_al_1_porciento](https://user-images.githubusercontent.com/69159364/92335466-89aae000-f06d-11ea-9162-a2a0d4d34a02.png)

Resumen: eulerint es bastante ineficiente. A pesar de lo anterior igual redujo significativamente el error, con 1000 subdivisiones, probablemente con un número mucho mayor de seguro se lograría el objetivo, pero la idea no es dejar corriendo el computador todo un día.

## Parte 4

**_Aquí también realice varios gráficos de más, considere los casos solo con J2, solo con J3 y con ambos_** 

Como resumen J2 es el que hace una mayor corrección de los 2.

**Colocaré solo los 4 gráficos de J2 + J3**


### Gráfico odeint vs real utilizando J2 y J3

![Pregunta4_Gráfica_J2_y_J3_deriva_odeint_vs_real](https://user-images.githubusercontent.com/69159364/92328441-361d9f80-f037-11ea-8289-c7b9fd334f4b.png)

### Gráfico eulerint vs real utilizando J2 y J3

![Pregunta4_Gráfica_J2_y_J3_deriva_eulerint_vs_real](https://user-images.githubusercontent.com/69159364/92328439-35850900-f037-11ea-9525-c0708c7a2c11.png)

### Gráfico odeint vs eulerint utilizando J2 y J3

![Pregunta4_Gráfica_J2_y_J3_deriva_odeint_vs_eulerint](https://user-images.githubusercontent.com/69159364/92328440-361d9f80-f037-11ea-890a-2523e2e4f6a5.png)

**R:** La deriva entre odeint y eulerint en el punto final es de 19343.5190 kilómetros. Experimento una reducción de más de 300 kilómetros.

### Gráfico de barras, tiempo de ejecución de odeint y eulerint utilizando J2 y J3

![Pregunta4_Gráfica_J2_y_J3_demora_odeint_vs_eulerint](https://user-images.githubusercontent.com/69159364/92328438-34ec7280-f037-11ea-8f7c-7720596cf766.png)

**R:** Nuevamente odeint es significativamente más rápido que eulerint. en este gráfico el odeint se demora 0.327(s) y el eulerint se demora 0.813(s).

_Nota: se se desea ver el resto de gráficos, ingresar en a la carpeta "Entrega5", en esta se encuentra todo lo relacionado con esta entrega (archivos .py, .EOF y .png)_


