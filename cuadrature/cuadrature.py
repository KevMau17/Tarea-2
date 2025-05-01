#!/usr/bin/env python3


"""Cálculo de sin(x²) mediante el método de cuadratura Gaussiana.

Éste módulo calcula la integral de 0 a pi de la función sin(x²) con el método de cuadratura Gaussiana.

El modulo contiene las siguientes funciones:

- `gaussxw(N)` - Retorna una tupla con los valores 'x' y 'w' sin escalar.
- `gaussxwab(a, b, x, w)` - Retorna la tupla de valores 'x' y 'w' pero escalados al intervalo [a,b].
- `integrando(x)` - Expresa la función sin(x²) como una función que recibe un número y retorna otro.
- `integral(a,b,funcion,N)` - Retorna la integral calculada de 'funcion' en el intervalo [a,b] con N segmentos.
"""

from scipy.special import legendre
import matplotlib.pyplot as plt
import numpy as np

def gaussxw(N):
    """Determina los ceros de un polinomio de Legredre de orden N (x) así como el peso de los términos del método de cuadratura Gaussiana (w). Estos valores solo funcionan para una integral en el intervalo de -1 a 1, para otro intervalo ver función 'gaussxwab'.

    Examples:
        >>> gaussxw(3)
        (array([-0.77459667, 0.0, 0.77459667]), array([0.55555556, 0.88888889, 0.55555556]))

    Args:
        N (int): Número de segmentos para la cuadratur Gaussiana

    Returns:
        (tuple): Retorna una tupla de ndarrays ('x' y 'w') con N elementos cada uno.

    """

    x, w = np.polynomial.legendre.leggauss(N)    
    return x, w

def gaussxwab(a, b, x, w):
    """Escala los valores de 'x' y 'w' para los límites de la integral que se calculará.

    Examples:
        >>> gaussxwab(0,3,array([-0.77459667, 0.0, 0.77459667]), array([0.55555556, 0.88888889, 0.55555556]))
        (array([0.338105, 1.5     , 2.661895]), array([0.83333334, 1.33333334, 0.83333334]))

    Args:
        a (float): Límite inferior de la integral a realizar.
        b (float): Límite superior de la integral a realizar.
        x (ndarray): Valores de los ceros de un polinomio de Legendre de un orden N.
        w (ndarray): Valores de los pesos de los términos de la suma de cuadratura Gaussiana.

    Returns:
        (tuple): Retorna una tupla de ndarrays con N elementos, ahora escalados al intervalo [a,b].

    """

    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w

def integrando(x):
    """Aquí se escribe la función que se desea integrar. En éste caso: sin(x²).

    Examples:
        >>> integrando(2**-0.5*np.pi**0.5)
        1.0

    Args:
        x (float): Un valor numérico.

    Returns:
        (float): Retorna el resultado de aplicarle a 'x' la función.
        
    """

    return np.sin(x**2)

def integral(a, b, funcion, N):
    """Ésta función calcula la integral usando las funciones 'gaussxw' y 'gaussxwab' para obtener los valore de 'x' y 'w' del método de cuadratura Gaussiana y luego suma los términos.

    Examples:
        >>> integral(0,np.pi,integrando,20)
        0.7726517126900648

    Args:
        a (float): Límite inferior de la integral.
        b (float): Límite superior de la integral.
        funcion (Callable[[float], float]): Función a integrar, que recibe la variable de integración.
        N (int): Número de segmentos/términos con que se hará la cuadratura Gaussiana.
    
    Returns:
        (float): Retorna el resultado de la integral.

    """
    x_0, w_0 = gaussxw(N)
    x_N, w_N = gaussxwab(a,b,x_0,w_0)
    return np.sum(w_N*funcion(x_N))

valores_N=np.arange(1,21,1,int)
valores_I=np.vectorize(integral)(0,np.pi,integrando, valores_N)

print(f"valor de la integral: {valores_I[19]}")

plt.plot(valores_N, valores_I, 'o')
plt.xlabel("Valores de N")
plt.ylabel("Valores de I")
plt.show()

