
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

