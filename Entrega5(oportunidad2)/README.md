
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

