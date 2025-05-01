# Método de Cuadratura Gaussiana

Uno de los métodos más poderosos para evaluar integrales de forma numérica es la **cuadratura Gaussiana**.

En esta clase vamos a discutir los resultados más importantes y cómo aplicar la idea de cuadraturas para resolver integrales.

Para derivaciones de resultados ver: Computational Physics - Mark Newman Capítulo 5 y Apéndice C.

La idea principal está dada por
\begin{align}
\int_a^b {\rm{d}}x f(x) \approx \sum_{k=1}^{N} w_k f(x_k).
\end{align}
donde:
  * $w_k$ son los "pesos"
  * $x_k$ son los puntos de muestreo
  
Para las ecuaciones de Newton-Cotes de la clase anterior:  
  * Los puntos de muestreo son **equidistantes**.
  * Una ecuación de Newton-Cotes de orden $N$ es *exacta* (i.e., no hay aproximación) para un polinomio de grado $N$.
  * Un polinomio de orden $N$ aproxima una función bien comportada mejor que un polinomio de orden $N-1$, debido al grado de libertad añadido.
  
Por el otro lado, para la cuadratura Gaussiana:
  * Los puntos de muestreo se escogen de manera tal que **no son equidistantes**. Esto introduce más grados de libertad para la misma discretización en $N$ subregiones.
  * Es exacta para un polinomio de orden $(2N - 1)$.
  * Es decir, la cuadratura Gaussiana da la misma precisión que un polinomio de orden $(2N - 1)$.

No vamos a probar el siguiente resultado (ver Apéndice C de Newman - Computational Physics), pero de manera muy interesante, existe una **regla universal para escoger $w_k$ y $x_k$**. Los pesos y puntos de muestreo se eligen tal que:
  * $x_k$ corresponden a las $N$ raíces (ceros) de los polinomios de Legendre $P_N(x)$ de orden $N$.
  * Los pesos se eligen tal que:
      - $\displaystyle w_k = \left[\frac{2}{1-x^2}\left(\frac{dP_N}{dx}\right)^{-2}\right]_{x={x_k}}$, con $x_k$ que cumple $P_N(x_k)=0$

## Polinomios de Legendre

Los polinomios de Legendre son un sistema de polinomios ortogonales que pueden ser definidos de manera recursiva. Tenemos:
\begin{align}
\forall (M, N) \in\mathbb N^2, \quad \int_{-1}^1 {\rm{d}}x P_N(x)P_M(x) = \frac{2\delta_{MN}}{2N+1}.
\end{align}
Note que los polinomios están definidos en el intervalo $[-1, 1]$.
Los se definen empezando con
\begin{align}
P_0(x) = 1 \Rightarrow P_1(x) = x,
\end{align}
tal que los siguientes órdenes se generan con la regla de recursividad
\begin{align}
(N+1)P_{N+1}(x) = (2N+1)xP_N(x) -NP_{N-1}(x).
\end{align}
Alternativamente, los polinomios pueden ser definidos de manera iterativa bajo la regla (fórmula de Rodrigues)
\begin{align}
P_N(x) = \frac1{2^N N!}\frac{d^N}{dx^N}\left[(x^2-1)^N\right].
\end{align}

En la computadora, podemos utilizar `SciPy` para obtener los polinomios de Legendre.

Una vez que conocemos los polinomios de Legendre, debemos encontrar sus raíces y calcular los pesos de acuerdo con la regla que describimos al inicio.

Esto es un procedimiento ligeramente costoso dependiendo de la metodología que se utilice. La idea es que si necesitamos evaluar la integral utilizando distintos intervalos de integración, primero realizamos el cálculo de los puntos de muestreo $x_k$ y los pesos $w_k$ en el intervalo $[-1, 1]$. Posteriormente, podemos escalar los parámetros para ser modificados a un intervalo $[a, b]$ (ver pags 167-168 Newman).

El siguiente código calcula los $x_k$ y  los $w_k$ utilizando **código vectorial**. Es decir, las funciones se aplican a **arreglos de datos** (en este caso de `Numpy`), en lugar de realizar `for loops` sobre todas las variables. 

Por ejemplo:

```python
x_test = np.linspace(-1.0, 1.0, 5)

print(np.sin(x_test)) #[-0.84147098 -0.47942554  0.          0.47942554  0.84147098]
```

Note que `Numpy` puede tomar como argumentos `np.ndarray`, de manera tal que no tenemos que hacer `for loops` para evaluar la función en cada argumento por separado. **Esto es muy importante**, dado que los `for loops` en `Python` son extremadamente lentos, mientras que las funciones de `NumPy` vectoriales están compiladas y enlazadas con código compilado.

Veamos el código que calcula los $x_k$ y los $w_k$. La biblioteca `Numpy` nos permite evaluar los pesos y los puntos de muestreo para utilizar la cuadratura Gaussiana con polinomios de Legendre:

```python
def gaussxw(N):
    x, w = np.polynomial.legendre.leggauss(N)
    
    return x, w
```

Este código se vuelve paulatinamente más costoso conforme aumentamos $N$ y solo funciona para intervalos de integración $[-1, 1]$. Para escalar el intervalo a un $[a, b]$ general, podemos utilizar la siguiente rutina:

```python
def gaussxwab(a, b, x, w):
    # Obtenido de pag 168 Newman)
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w
```

Una vez que tenemos los puntos de muestreo y los pesos, evaluar integrales es **trivial**. Basta con utilizar la regla:
\begin{align}
\int_a^b {\rm{d}}x f(x) \approx \sum_{k=1}^{N} w_k f(x_k).
\end{align}

Lo haremos mediante una función de la siguiente manera (aprovechando que 'x' y 'w' son ndarrays):

```python
def integral(a, b, funcion, N): #a:inicio, b:final, N: cantidad de segmentos
    x_0, w_0 = gaussxw(N)
    x_N, w_N = gaussxwab(a,b,x_0,w_0)
    return np.sum(w_N*funcion(x_N)) #método de cuadratura gaussiana
```
