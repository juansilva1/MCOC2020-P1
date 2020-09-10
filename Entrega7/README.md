# Entrega 7

## Mejoras realizadas

1) la primera fue corregir la entrega anterior, en la cual estaba utilizando mal los términos J2 y J3.

2) mejore la precisión de los términos J2 y J3.

3) Aumenté la velocidad al eliminar uno de los ciclos.

### Triste historia (muy extra)


Desde el lunes Empecé a realizar mejoras en el código (partí arreglando los errores) y me propuse a obtener Jn hasta J8 y los términos que acompañaban C15 y S15. para esto tras 3 días de código logré estas funciones:

```
zonal_coefficients=[-0.1082635854e-2, 0.2532435346e-5, 0.1619331205e-5,
                    0.2277161016e-6, -0.5396484906e-6, 0.3513684422e-6,
                    0.2025187152e-6]
km3=(1000.)**3
mu=398600.4415*km3
R=6378136.3
def creador_J(zonal_coefficients):
    contador=2
    lista=[]
    for i in zonal_coefficients:
        resultado=-(i*mu*(R**contador))
        contador+=1
        lista.append(resultado)
    return lista

"""
Lista hasta J8 en [m y s]
"""
lista_Jn=[1.7555280486257942e+25, -2.6191328602512207e+29, -1.0681903337751385e+36, -9.580771773008308e+41, 1.4481432257808013e+49, -6.013912422149775e+55, -2.210819596773952e+62]

"""
Creacion de lista para Pn
"""
def factorial(n):
    if n<2:
        return 1
    else:
        return n*factorial(n-1) 

def combinatoria(n,k):
    return(factorial(n)/(factorial(k)*factorial(n-k)))

def Pn(n):
    resultado=0
    
    for k in range(n+1):
        resultado+=combinatoria(n,k)*combinatoria(n+k,k)*(((sin(teta)-1)/2)**k)
    return expand((resultado))

def J(n):
    a=lista_Jn[n-2]
    paso1=(Jnn*Pn(n))/(r**(n+1))
    paso2=paso1.subs(sin(teta),z/r) # retornar paso2 lo deja como en wiki
    paso3=paso2.subs(r,sqrt(x**2 + y**2 + z**2))

    return paso3
    
# paso a los P_nm y la sumatoria doble

#lista que pocee hasta C_S (n=8,m=8)
lista1=[['2', '1', '-0.24140000522221e-09', '0.15430999737844e-08'],
        ['3', '1', '0.21927988018965e-05', '0.26801189379726e-06'],
        ['4', '1', '-0.50872530365024e-06', '-0.44945993508117e-06'], 
        ['5', '1', '-0.53716510187662e-07', '-0.80663463828530e-07'], 
        ['6', '1', '-0.59877976856303e-07', '0.21164664354382e-07'], 
        ['7', '1', '0.20514872797672e-06', '0.69369893525908e-07'], 
        ['8', '1', '0.16034587141379e-07', '0.40199781599510e-07'], 
        ['2', '2', '0.15745360427672e-05', '-0.90386807301869e-06'], 
        ['3', '2', '0.30901604455583e-06', '-0.21140239785975e-06'], 
        ['4', '2', '0.78412230752366e-07', '0.14815545694714e-06'], 
        ['5', '2', '0.10559053538674e-06', '-0.52326723987632e-07'], 
        ['6', '2', '0.60120988437373e-08', '-0.46503948132217e-07'], 
        ['7', '2', '0.32844904836492e-07', '0.92823143885084e-08'], 
        ['8', '2', '0.65765423316743e-08', '0.53813164055056e-08'], 
        ['3', '3', '0.10055885741455e-06', '0.19720132389889e-06'], 
        ['4', '3', '0.59215743214072e-07', '-0.12011291831397e-07'], 
        ['5', '3', '-0.14926153867389e-07', '-0.71008771406986e-08'], 
        ['6', '3', '0.11822664115915e-08', '0.18431336880625e-09'], 
        ['7', '3', '0.35285405191512e-08', '-0.30611502382788e-08'], 
        ['8', '3', '-0.19463581555399e-09', '-0.87235195047605e-09'], 
        ['4', '4', '-0.39823957404129e-08', '0.65256058113396e-08'], 
        ['5', '4', '-0.22979123502681e-08', '0.38730050770804e-09'], 
        ['6', '4', '-0.32641389117891e-09', '-0.17844913348882e-08'], 
        ['7', '4', '-0.58511949148624e-09', '-0.26361822157867e-09'], 
        ['8', '4', '-0.31893580211856e-09', '0.91177355887255e-10'], 
        ['5', '5', '0.43047675045029e-09', '-0.16482039468636e-08'], 
        ['6', '5', '-0.21557711513900e-09', '-0.43291816989540e-09'], 
        ['7', '5', '0.58184856030873e-12', '0.63972526639235e-11'], 
        ['8', '5', '-0.46151734306628e-11', '0.16125208346784e-10'], 
        ['6', '6', '0.22136925556741e-11', '-0.55277122205966e-10'], 
        ['7', '6', '-0.24907176820596e-10', '0.10534878629266e-10'], 
        ['8', '6', '-0.18393642697634e-11', '0.86277431674150e-11'], 
        ['7', '7', '0.25590780149873e-13', '0.44759834144751e-12'], 
        ['8', '7', '0.34297618184624e-12', '0.38147656686685e-12'], 
        ['8', '8', '-0.15803322891725e-12', '0.15353381397148e-12']]





# print(10e02)
archivo=open('C_S.txt')
lista=[]
diccionario={}
for line in archivo:
    line=line.split()
    line[2]=line[2].replace('D','e')
    line[3]=line[3].replace('D','e')
    lista.append(line)
# print(lista)

# ==> desde aca se puede copiar (incluir lista superior y datos)
for i in lista1:
    diccionario[f'{i[0]}_{i[1]}']=[float32(i[2]),float32(i[3])]
# print(diccionario)
# hasta acá tengo un diccionario que como llave tiene n_m y entrega una lista [C,S]
"""
CS_nm
"""
def CS_nm(n,m):
    n_m=f'{n}_{m}'
    C=diccionario[n_m][0]
    S=diccionario[n_m][1]
    
    C_f=-C*mu*(R**n)
    S_f=-S*mu*(R**n)

    return [C_f,S_f]
    


# hasta este punto ya tengo los C y S correctos hasta el n=8,m=8

    
"""
P_nm
"""

def P_nm(n,m):
    funcion_Pn=Pn(n)
    resultado=((-1)**m) * ((1-((sin(teta))**2))**(m/2)) * diff(funcion_Pn,sin(teta),m)
    return resultado


"""
Unir todo Pnm
"""

def Unir_P_nm (n,m):
    
    elemento_suma=0
    
    Nt=linspace(2,n,n-1)
    for i in Nt:
        n1=linspace(1,i,i)
        for j in n1:
            C=CS_nm(i,j)[0]
            S=CS_nm(i,j)[1]
            elemento_suma += ( P_nm(i,j) * (C*cos(j*fi) + S*sin(m*fi)) )/(r**(i+1))
    return elemento_suma

# hasta aca tengo las 2 sumatorias completas

# posterior a esto cree una funcion derivadas

def derivadas(ex):
    ex1=ex.subs(r,sqrt(x**2 + y**2 + z**2))
    dx=-diff(ex1,x)
    dx=factor(dx)
    
    dy=-diff(ex,y)
    dy=factor(dy)
    
    dz=-diff(ex,z)
    dz=factor(dz)
    
    return np.array([dx,dy,dz])
 
```

**Con todo esto solo quedaba unirlo y probarlo: resultado ==> no funciona.**


```
El profe manda una versión mucho más simple de lo que estaba intentando, me decido a probarlo obtengo las derivadas hasta J15
```
**lo remplazo: ==> no me funciona**

*_Dato: Fin de la historia ==> tras 4 días programando las mejoras no valieron realmente todo el tiempo invertido._*
