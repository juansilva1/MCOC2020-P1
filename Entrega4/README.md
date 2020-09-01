# Entrega 4

## Análisis

Se puede apreciar que la función analítica, es muy similar a la función odeint, esto ocurre pues esta función "automáticamente" toma subdivisiones, lo que hace que el resultado aproximado sea muy cercano al real.

Por otro lado, las funciones euleint son aproximaciones que dependen del número de subdivisiones que posea el intervalo, si este es muy pequeño, el resultado es "desastroso", no se asemeja a la función analítica, pero si el número de subdivisiones aumenta los valores son cada vez más cercanos al original, por lo que el error se corrige significativamente. Se puede observar que para un Nsubdivisiones=100 el resultado al menos a simple vista resulta satisfactorio.

![Grafico Oscilador_Armónico](https://user-images.githubusercontent.com/69159364/91896708-4ca6ae00-ec67-11ea-8221-c7c645fe4e93.png)
