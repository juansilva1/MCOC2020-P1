# Entrega 4

## Analisis

Se puede apreciar que la función analítica, es muy similar a la función odeint, esto ocurre pues esta funcion "automaticamente" toma subdibiciones, lo que hace que el resultado aproximado sea muy sercano al real.

por otro lado las funciones euleint son aproximaciones que tependen del numero de subdiviciones que posea el intervalo, si este es muy pequeño, el resultado es "desastrozo", no se asemeja a la función analitica, pero si el numero de subdiviciones aumenta los valores son cada vez más cercanos al original, por lo que el error se corrige significativamente. Se puede observar que para un Nsubdivisiones=100 el resultado al menos a simplevista resulta satisfactorio.

