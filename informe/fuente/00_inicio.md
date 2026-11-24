---
lang: es
mainfont: DejaVu Serif
mathfont: Latin Modern Math
fontsize: 11pt
geometry: "top=2.4cm,bottom=2.4cm,left=2.5cm,right=2.5cm"
documentclass: article
header-includes:
  - \input{header.tex}
---

\begin{titlepage}
\centering
\vspace*{1.5cm}
{\Large Universidad de San Buenaventura, Medellín\par}
\vspace{3cm}
{\Huge\bfseries Informe -- Taller final\par}
\vspace{0.6cm}
{\LARGE Análisis de Ondas y Vibraciones\par}
\vspace{0.4cm}
{\large Versión corregida\par}
\vspace{2.6cm}
{\large\bfseries Integrantes\par}
\vspace{0.3cm}
{\large Isabella Patiño Ortega\\ Juan Camilo Restrepo\\ Maikol Andrés Zuleta Serna\\ Andrés Felipe Carmona\\ Oscar Santiago Chamorro\par}
\vspace{2cm}
{\large\bfseries Docente\par}
\vspace{0.3cm}
{\large Ing. Jairo Coy\par}
\vfill
{\large Noviembre de 2026\par}
\end{titlepage}

\tableofcontents
\newpage

# 1. Introducción

Este informe recoge la solución del taller final del curso de Análisis de Ondas y Vibraciones, orientado por el Ing. Jairo Coy. El taller pedía escoger cinco problemas de sistemas no forzados y cuatro de sistemas forzados, más el forzado 10, que era obligatorio. Escogimos los no forzados 1, 3, 4, 5 y 9 y los forzados 1, 2, 3, 7 y 10.

La idea de fondo es salir del caso ideal que vemos en clase, el oscilador lineal con amortiguamiento viscoso $m\ddot x + c\dot x + kx = F(t)$, y ver qué pasa cuando aparecen cosas que en la práctica son normales: un amortiguamiento que depende de la velocidad al cuadrado, un resorte que se endurece, una masa que cambia con el tiempo, varios grados de libertad o golpes contra un tope. Para muchos de estos casos no hay una solución analítica cerrada, así que el camino es doble: primero se resuelve a mano todo lo que la teoría sí permite resolver (ecuación característica, constantes a partir de las condiciones iniciales, solución particular, balance de energía, balance armónico, puntos de equilibrio), y después se simula el modelo completo y se compara con esa solución.

**Qué cambia en esta versión.** Siguiendo las observaciones del docente, cada ejercicio quedó organizado así:

1. **Fenómeno físico.** Qué se está estudiando, por qué ocurre y dónde aparece en la práctica.
2. **Modelo.** Diagrama de cuerpo libre, fuerzas que actúan y deducción de la ecuación de movimiento con la segunda ley de Newton (o el principio de D'Alembert).
3. **Condiciones iniciales.** Qué valores se usaron y por qué.
4. **Solución de la ecuación, paso a paso.** El desarrollo analítico completo con los números del problema, y después cómo se resuelve numéricamente.
5. **Resultados, análisis y conclusión.**

El objetivo general fue resolver cada problema con un procedimiento claro y completo, mostrar los resultados en gráficas y, sobre todo, entender qué está pasando físicamente en cada caso. Los objetivos específicos fueron:

- Plantear el modelo de cada sistema desde el diagrama de cuerpo libre y llevarlo a la ecuación de movimiento.
- Resolver analíticamente la parte del problema que lo permite y usar esa solución como referencia.
- Programar en MATLAB códigos sencillos, con los parámetros al inicio para poder cambiarlos fácilmente.
- Contrastar cada simulación con el resultado teórico para tener confianza en lo que sale.
- Interpretar los resultados y sacar conclusiones útiles para el diseño y diagnóstico de sistemas vibratorios.

# 2. Marco teórico y metodología

Esta sección reúne la teoría del curso que se usa en todos los ejercicios. Así, en cada ejercicio se puede ir directo a los números sin repetir la deducción completa cada vez.

## 2.1 El sistema vibratorio y su ecuación de movimiento

Un sistema vibratorio es un conjunto que **almacena, transforma y disipa energía** alrededor de una posición de equilibrio. La masa almacena energía cinética, el resorte almacena energía potencial y el amortiguador la disipa en forma de calor. La vibración es el intercambio continuo entre energía cinética y potencial; el amortiguamiento hace que en cada ciclo se pierda una parte.

El modelo básico es el de la Figura \ref{fig:sdof}. Se toma $x$ como el desplazamiento de la masa medido desde la posición de equilibrio, positivo hacia la derecha. Al aislar la masa (diagrama de cuerpo libre), actúan:

- La fuerza del resorte, $-kx$: se opone al desplazamiento.
- La fuerza del amortiguador viscoso, $-c\dot x$: se opone a la velocidad.
- La fuerza externa $F(t)$, si existe.

\begin{figure}[H]\centering
\resizebox{0.95\linewidth}{!}{\input{fig_sdof.tex}}
\caption{Sistema de un grado de libertad. (a) Modelo físico masa--resorte--amortiguador. (b) Diagrama de cuerpo libre de la masa; la flecha punteada es la fuerza de inercia de D'Alembert.}\label{fig:sdof}
\end{figure}

Por la segunda ley de Newton, $\sum F = m\ddot x$:

$$
F(t) - kx - c\dot x = m\ddot x
\quad\Longrightarrow\quad
\boxed{\,m\ddot x + c\dot x + kx = F(t)\,}
$$

Con el principio de D'Alembert se llega a lo mismo: se agrega la fuerza ficticia de inercia $-m\ddot x$ y se plantea equilibrio, $F(t) - kx - c\dot x - m\ddot x = 0$.

## 2.2 Vibración libre amortiguada: solución analítica

Con $F(t)=0$, la ecuación es lineal, homogénea y de coeficientes constantes. Se propone $x = e^{\lambda t}$ (usamos $\lambda$ y no $r$ para no confundir con otras letras del taller). Al reemplazar, $e^{\lambda t}$ se cancela y queda la **ecuación característica**:

$$
m\lambda^2 + c\lambda + k = 0
\quad\Longrightarrow\quad
\lambda_{1,2} = \frac{-c \pm \sqrt{c^2 - 4mk}}{2m} = -\zeta\omega_n \pm \omega_n\sqrt{\zeta^2 - 1}
$$

con las definiciones del curso:

$$
\omega_n = \sqrt{\frac{k}{m}}, \qquad c_c = 2\sqrt{km} = 2m\omega_n, \qquad \zeta = \frac{c}{c_c}, \qquad \omega_d = \omega_n\sqrt{1-\zeta^2}
$$

El signo del discriminante $c^2 - 4mk$ define tres casos. En cada uno, las dos constantes se hallan con las condiciones iniciales $x(0)=x_0$ y $\dot x(0)=v_0$:

**Subamortiguado** ($\zeta<1$, raíces complejas $\lambda = -\zeta\omega_n \pm i\omega_d$):
$$
x(t) = e^{-\zeta\omega_n t}\left(C_1\cos\omega_d t + C_2\sin\omega_d t\right)
$$
En $t=0$: $x(0) = C_1 = x_0$. Derivando y evaluando en $t=0$: $\dot x(0) = -\zeta\omega_n C_1 + \omega_d C_2 = v_0$, de donde
$$
C_1 = x_0, \qquad C_2 = \frac{v_0 + \zeta\omega_n x_0}{\omega_d}
$$

**Críticamente amortiguado** ($\zeta=1$, raíz doble $\lambda = -\omega_n$):
$$
x(t) = (C_1 + C_2 t)\,e^{-\omega_n t}, \qquad C_1 = x_0, \qquad C_2 = v_0 + \omega_n x_0
$$

**Sobreamortiguado** ($\zeta>1$, dos raíces reales negativas $\lambda_1 > \lambda_2$):
$$
x(t) = C_1 e^{\lambda_1 t} + C_2 e^{\lambda_2 t}, \qquad
C_1 = \frac{v_0 - \lambda_2 x_0}{\lambda_1 - \lambda_2}, \qquad
C_2 = \frac{\lambda_1 x_0 - v_0}{\lambda_1 - \lambda_2}
$$
(sale de resolver el sistema $C_1 + C_2 = x_0$, $\lambda_1C_1 + \lambda_2C_2 = v_0$).

## 2.3 Vibración forzada armónica: solución analítica

Con $F(t) = F_0\sin\omega t$, la solución total es la suma de la homogénea y la particular, $x = x_h + x_p$. La homogénea es la de la sección 2.2 y se apaga como $e^{-\zeta\omega_n t}$; por eso se llama **transitorio**, y su **constante de tiempo** es

$$
\tau = \frac{1}{\zeta\omega_n}
$$

Después de unas cuantas $\tau$ solo queda la particular, el **estado estacionario**. Se propone $x_p = X\sin(\omega t - \varphi)$ y se reemplaza en la ecuación. Al separar los términos en $\sin\omega t$ y $\cos\omega t$ quedan dos ecuaciones:

$$
(k - m\omega^2)X = F_0\cos\varphi, \qquad c\omega X = F_0\sin\varphi
$$

Elevando al cuadrado y sumando, y dividiendo una entre la otra:

$$
X = \frac{F_0}{\sqrt{(k - m\omega^2)^2 + (c\omega)^2}} = \frac{F_0/k}{\sqrt{(1 - r^2)^2 + (2\zeta r)^2}},
\qquad
\tan\varphi = \frac{2\zeta r}{1 - r^2},
\qquad r = \frac{\omega}{\omega_n}
$$

El cociente $H(r) = X/(F_0/k)$ es el **factor de amplificación**. Derivando el denominador respecto a $r$ e igualando a cero se obtiene dónde está el pico y cuánto vale:

$$
r_{pico} = \sqrt{1 - 2\zeta^2}, \qquad X_{máx} = \frac{F_0/k}{2\zeta\sqrt{1-\zeta^2}}
$$

## 2.4 Por qué hace falta la solución numérica

Las soluciones anteriores existen porque la ecuación es **lineal y de coeficientes constantes**. En el taller aparecen términos que rompen alguna de esas dos condiciones: $\dot x|\dot x|$ (NF3), $x^3$ (NF4, F3, F7), $m(t)$ (NF5), amortiguamiento no proporcional en varios grados de libertad (NF9) o fuerzas por tramos (topes, F10). En esos casos se resuelve analíticamente lo que se pueda (casos límite, aproximaciones, balances de energía) y el modelo completo se integra numéricamente.

**Forma de estado.** Los integradores numéricos trabajan con ecuaciones de primer orden. Se definen $y_1 = x$ y $y_2 = \dot x$, y la ecuación de segundo orden queda como dos de primer orden:

$$
\dot y_1 = y_2, \qquad \dot y_2 = \frac{F(t) - c\,y_2 - k\,y_1}{m}
\qquad\Longleftrightarrow\qquad
\dot{\mathbf y} = \mathbf f(t, \mathbf y)
$$

**Runge--Kutta de cuarto orden (RK4).** Avanza la solución con un paso fijo $h$. En cada paso evalúa la pendiente $\mathbf f$ en cuatro puntos y los promedia con pesos 1--2--2--1:

$$
\begin{aligned}
\mathbf k_1 &= \mathbf f(t_n, \mathbf y_n), &
\mathbf k_2 &= \mathbf f\!\left(t_n + \tfrac h2,\ \mathbf y_n + \tfrac h2\mathbf k_1\right),\\
\mathbf k_3 &= \mathbf f\!\left(t_n + \tfrac h2,\ \mathbf y_n + \tfrac h2\mathbf k_2\right), &
\mathbf k_4 &= \mathbf f(t_n + h,\ \mathbf y_n + h\,\mathbf k_3),\\
\mathbf y_{n+1} &= \mathbf y_n + \tfrac h6\left(\mathbf k_1 + 2\mathbf k_2 + 2\mathbf k_3 + \mathbf k_4\right)
\end{aligned}
$$

El error acumulado es proporcional a $h^4$: si se divide el paso entre dos, el error baja unas 16 veces. En el NF3 se muestra un paso hecho a mano con números.

**ode45.** Es el integrador estándar de MATLAB (método de Dormand--Prince). En cada paso calcula dos estimaciones, una de orden 4 y otra de orden 5. La diferencia entre ellas es una estimación del error: si supera la tolerancia pedida (`RelTol`, `AbsTol`), rechaza el paso y lo repite más pequeño; si sobra margen, agranda el siguiente. Así el paso se adapta solo a la dificultad de la solución.

**Diferencias finitas centradas.** Reemplazan las derivadas por cocientes de diferencias en una malla de tiempo y despejan la posición futura. Se usaron en el F10 y allí se deducen completas, junto con su condición de estabilidad.

## 2.5 Criterios para escoger las condiciones iniciales

Las condiciones iniciales definen las constantes de la solución homogénea, así que cambian la respuesta. Seguimos estos criterios, y en cada ejercicio se explica cuál se aplicó:

1. **Si el enunciado da la condición inicial, se usa esa.**
2. **Vibración libre:** se modela el ensayo típico de "desplazar y soltar". Por eso $v_0 = 0$ (la masa se suelta desde el reposo) y $x_0$ es la amplitud inicial. En un sistema **lineal** el valor de $x_0$ solo escala la respuesta: la frecuencia, el amortiguamiento y los tiempos relativos no cambian, así que $x_0$ funciona como una normalización. En un sistema **no lineal** la amplitud sí cambia el comportamiento, y $x_0$ se escoge para que el término no lineal sea comparable con el lineal y el fenómeno se vea.
3. **Vibración forzada:** el sistema arranca en reposo, $x(0) = 0$ y $\dot x(0) = 0$, como una máquina apagada a la que se le aplica la fuerza. Si lo que interesa es el estado estacionario, se descarta un tiempo de varias constantes de tiempo $\tau$, porque en un sistema lineal la parte que depende de las condiciones iniciales se apaga como $e^{-t/\tau}$. En sistemas no lineales el estado final **sí puede depender** de la condición inicial (F3, F7), y eso se discute en cada caso.

## 2.6 Validación y herramientas

Todos los resultados numéricos se compararon con una fórmula teórica, con otro método o con un paso de tiempo más fino. Además, los resultados principales se repitieron con un programa independiente en Python (SciPy).

| Método | Qué hace | Dónde se usó |
|:--|:--|:--|
| `ode45` | Runge--Kutta 4(5) con paso adaptativo | NF1, NF4, NF5, NF9, F1, F2, F7 |
| RK4 programado a mano | Runge--Kutta clásico de cuarto orden, paso fijo | NF3, F3 |
| `expm` | Solución exacta de un sistema lineal, $\mathbf z(t) = e^{\mathbf A t}\mathbf z_0$ | NF9 (comprobación) |
| Diferencias finitas centradas | Esquema explícito de segundo orden | F10 |
| FFT | Paso del dominio del tiempo al de la frecuencia | F2, F10 |

Todos los códigos usan solo funciones básicas de MATLAB, sin toolboxes. Cada script tiene un bloque de parámetros al inicio y funciona solo, sin depender de otros archivos.

\newpage
