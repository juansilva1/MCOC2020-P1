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

*El codigo imprime al final la velocidad de **24551 [km/h]** *
