Tarea-Examen 2: Simulación Estocástica
================
José Jorge Martínez de la Cruz
2026-04-15

## Problema 1

Utilizando el método de muestreo por importancia para el método de media
muestral estime $$ \int_{0}^{1} \frac{1}{1+x^2} dx$$ Use la aproximación
discreta de la g que minimiza la varianza. Para definir a $g$ divida al
intervalo en 10 partes y tome como $g_i$ cualquier valor de la función
en el intervalo.

### Solución

Para realizar esto, debemos usar una aproximación discreta de la
densidad $g$ que minimiza la varianza del estimador. Para realizar una
aproximación discreta, denotada como $\overline{g}$, debemos dividir el
intervalo en 10 partes iguales y tomar como $g_i$ cualquier valor de la
función dentro de dicho subintervalo.

Podemos reescribir la integral como la esperanza matemática de una
función respecto a una distribución de probabilidad. Si definimos
$h(x) = \frac{1}{1+x^2}$ y consideramos $f(x) = 1$ como la función de
densidad de una distribución $\text{Uniforme}(0,1)$, nuestro objetivo es
estimar: $$\mathbb{E}_f(h(X)) = \int_{0}^{1} h(x)f(x)dx$$

El método de muestreo por importancia se basa en la siguiente identidad
analítica:
$$\mathbb{E}_f(h(X)) = \int_{0}^{1} h(x) \frac{f(x)}{g(x)} g(x) dx = \mathbb{E}_g(h(X)w(X))$$
donde $g(\cdot)$ es una densidad instrumental tal que $g(x) > 0$ en el
soporte de $f$, y $w(x) = \frac{f(x)}{g(x)}$ representa la función de
pesos.

Por la ley fuerte de los grandes números, si simulamos una muestra
$X_1, X_2, \dots, X_N \sim g$, podemos aproximar la esperanza deseada
con el siguiente estimador:
$$\hat{\mu} = \frac{1}{N} \sum_{k=1}^{N} h(X_k) w(X_k) = \frac{1}{N} \sum_{k=1}^{N} h(X_k) \frac{f(X_k)}{g(X_k)}$$

Usando que la distribución instrumental teórica que minimiza la varianza
del estimador está dada por:
$$g^*(x) = \frac{|h(x)|f(x)}{\int |h(t)|f(t)dt}$$

Nos encontramos con que el denominador es precisamente la integral que
deseamos calcular, $g^*$ no puede usarse directamente en este ejercicio.
Por lo tanto, construimos una aproximación escalonada $\overline{g}$.

Para ello, dividimos $[0,1]$ en $n = 10$ subintervalos de longitud
$\Delta x = 0.1$ y en cada subintervalo $i \in \{1, \dots, 10\}$,
elegimos el punto medio $x_i^*$ y evaluamos los pesos
$g_i = |h(x_i^*)|f(x_i^*)$.Luego, para que $\overline{g}$ sea una
función de densidad de probabilidad, debe integrar a 1. Calculamos la
constante de normalización: $$C = \sum_{j=1}^{10} g_j \Delta x$$ De esta
forma podemos definir nuestra densidad instrumental aproximada para
cualquier $x$ dentro del $i$-ésimo subintervalo como:
$$\overline{g}(x) = \frac{g_i(x)}{C}$$

Finalmente, el algoritmo consistirá en simular variables aleatorias
$X_k \sim \overline{g}$, simularemos uniformemente entre estos valores
simulados para poder aplicar el estimador de la media muestral
ponderada:
$$\hat{\mu} \approx \frac{1}{N} \sum_{k=1}^{N} h(X_k) \frac{f(X_k)}{\overline{g}(X_k)}$$

Para corroborar este resultado, podemos comparar la estimación obtenida
con el valor exacto de la integral, que se puede calcular analíticamente
como:

$$\int_{0}^{1} \frac{1}{1+x^2} dx =  \arctan(x) | _{0}^{1} = \arctan(1) - \arctan(0) = \frac{\pi}{4}$$
Los resultados obtenidos después de aplicar el método de muestreo por
importancia se muestran a continuación:

| Métrica                           |    Valor     |
|:----------------------------------|:------------:|
| Estimación obtenida               |  0.7853982   |
| Valor analítico ($\frac{\pi}{4}$) |  0.7853982   |
| Error absoluto                    | 1.110223e-16 |

Los resultados obtenidos indican que la estimación obtenida mediante el
método de muestreo por importancia es extremadamente cercana al valor
analítico de la integral, con un error absoluto prácticamente nulo. Esto
valida la efectividad del método utilizado para aproximar la integral
deseada.

## Problema 2

Utilizando el método clásico de simulación Montecarlo de dimensión $2$
calcule la probabilidad de que un vector con distribución bivariada
conjunta Gaussiana estandarizada con correlación $0.5$ esté en el
cuadrado unitario.

### Solución

Definamos nuestro vector aleatorio bivariado como
$\mathbf{X} = (X_1, X_2)^T$. Dado que se trata de una distribución
Gaussiana estandarizada, sabemos que las medias marginales son cero y
las varianzas son uno. La matriz de varianzas y covarianzas $\Sigma$ se
construye a partir de la correlación dada ($\rho = 0.5$):

$$\mathbf{X} \sim \mathcal{N}_2(\boldsymbol{\mu}, \Sigma) \quad \text{donde} \quad \boldsymbol{\mu} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}, \quad \Sigma = \begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix}$$

Lo que nos interesa ahora es que el vector caiga dentro del “cuadrado
unitario”, el cual definiremos en el primer cuadrante como la región
$A = [0,1] \times [0,1]$. Es decir, buscamos la probabilidad:
$$p = \mathbb{P}(\mathbf{X} \in A) = \mathbb{P}(0 \le X_1 \le 1, \ 0 \le X_2 \le 1)$$

Matemáticamente, esto equivale a resolver la integral doble de la
función de densidad conjunta $f(x_1, x_2)$:
$$p = \int_{0}^{1} \int_{0}^{1} \frac{1}{2\pi\sqrt{|\Sigma|}} \exp\left( -\frac{1}{2} \mathbf{x}^T \Sigma^{-1} \mathbf{x} \right) dx_1 dx_2$$

Dado que esta integral no tiene una solución analítica cerrada, el
método clásico de Montecarlo nos da una manera de resolverlo. El
procedimiento consiste en:

1.  Simular una muestra de tamaño $N$ de vectores aleatorios
    independientes
    $\mathbf{X}^{(1)}, \mathbf{X}^{(2)}, \dots, \mathbf{X}^{(N)}$
    provenientes de la distribución
    $\mathcal{N}_2(\boldsymbol{\mu}, \Sigma)$.
2.  Definir la función indicadora $I(\mathbf{X})$ que tome el valor de 1
    si el vector cae dentro de la región $A$ y 0 en caso contrario.
3.  Estimar la probabilidad calculando la proporción de vectores que
    “aterrizaron” en la región de interés por la Ley de los Grandes
    Números:
    $$\hat{p} = \frac{1}{N} \sum_{i=1}^{N} I(\mathbf{X}^{(i)} \in A)$$

Para implementar esto en R, podemos usar la función `MASS::mvrnorm` para
generar muestras de la distribución bivariada y luego calcular la
proporción de vectores que caen dentro del cuadrado unitario. Además
incluiremos el paquete `mvtnorm` para calcular la probabilidad teórica
real y poder validar la precisión de tu estimador.

Una vez realizado el algoritmo, obtenemos los siguientes resultados:

| Métrica                   |    Valor    |
|:--------------------------|:-----------:|
| Estimación por Montecarlo |   0.14274   |
| Valor teórico (mvtnorm):  |  0.141051   |
| Error absoluto            | 0.001688985 |

Con estos resultados, concluimos que la estimación obtenida mediante el
método de Montecarlo es bastante cercana al valor teórico calculado con
la función `pmvnorm`, lo que valida la precisión de nuestro estimador.
El error absoluto es pequeño, lo que indica que el método de simulación
ha sido efectivo para aproximar la probabilidad deseada.
