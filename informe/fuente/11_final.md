# 13. Conclusiones

1. **El amortiguamiento crítico no es el más rápido.** Con el criterio del 1 %, el óptimo es el $\zeta$ con el que el primer sobrepaso mide exactamente el 1 %: $\zeta^* = \ln100/\sqrt{\pi^2 + \ln^2100} = 0.826$. Con ese valor el sistema se asienta en unos 0.68 s, frente a 1.05 s en el crítico. Por eso tantos sistemas reales se diseñan un poco por debajo del crítico.
2. **No todo amortiguamiento se puede linealizar.** El balance de energía muestra que el arrastre cuadrático disipa como $X^3$ por ciclo y produce un decaimiento hiperbólico, $X_0/(1 + \beta X_0t)$, no exponencial: frena bien las oscilaciones grandes, pero casi nada las pequeñas. Un amortiguador lineal equivalente solo sirve para la primera oscilación.
3. **En sistemas no lineales la frecuencia depende de la amplitud.** En el Duffing libre la rigidez efectiva es $k_1 + \tfrac34k_3A^2$, y la frecuencia pasó de 9.25 a 3.18 rad/s mientras la oscilación se apagaba. En el forzado, esa misma dependencia dobla la curva de resonancia y produce saltos, histéresis y una zona biestable entre 1.63 y 3.11 rad/s, con bordes que el balance armónico predice.
4. **Los parámetros que cambian en el tiempo cambian la frecuencia.** Con masa variable, la frecuencia sigue de cerca $\sqrt{k/m(t)}$ siempre que la masa cambie despacio comparada con un periodo, y la amplitud sigue la ley $X_0(m/M_0)^{3.5}$.
5. **Con amortiguamiento no proporcional hay que trabajar en espacio de estados.** Los autovalores de la matriz 6×6 (raíces de $\det(\lambda^2\mathbf M + \lambda\mathbf C + \mathbf K) = 0$) dan directamente la frecuencia, el amortiguamiento y la duración de cada modo, y explican por qué al final solo sobrevive el primer modo.
6. **La resonancia depende de cómo y con qué se excita.** Un barrido rápido corre el pico hacia arriba y lo hace más pequeño, porque la amplitud de resonancia necesita un tiempo del orden de $1/(\zeta\omega_n)$ para construirse. Una excitación no senoidal puede poner en resonancia el sistema a través de un armónico: con 1 Hz, el tercer armónico produjo una respuesta 5.6 veces mayor que la fundamental. Con 1.5 Hz no hubo resonancia, porque una onda cuadrada no tiene segundo armónico.
7. **Un sistema determinista puede ser impredecible.** El Duffing de doble pozo pasa por la cascada de duplicación de periodo hasta el caos cuando la fuerza le da energía suficiente para cruzar la loma entre los dos pozos. Su sección de Poincaré en $\gamma = 0.5$ muestra un atractor extraño.
8. **Los impactos generan distorsión, pero no atrapan la energía.** En la cadena con topes, la zona de contacto (4–21.3 rad/s) es fuertemente no lineal, con THD alto y subarmónicos. Aun así, la energía se queda sobre todo en la masa excitada, que resuena en 14.14 rad/s al quedar casi aislada por los topes.
9. **Las condiciones iniciales hay que justificarlas.** En los sistemas lineales, el valor inicial solo escala la respuesta o se borra con el transitorio, y basta con descartar varias constantes de tiempo. En los no lineales, en cambio, decide qué fenómeno se observa (NF3, NF4) o en qué estado termina el sistema (F3, F7).
10. **Validar es parte de resolver.** En todos los ejercicios lo simulado se comparó con una solución analítica (exacta o aproximada), con otro método o con un paso de tiempo más fino. Esto permitió encontrar dos inconsistencias del enunciado (la del armónico en el forzado 2 y la del número de muestras en el forzado 10) y saber qué resultados son robustos y cuáles, como el THD puntual en régimen caótico, deben leerse con cuidado.

# Referencias

- Coy, J. *Análisis de ondas y vibraciones* (diapositivas del curso). Universidad de San Buenaventura, Medellín.
- Rao, S. S. *Mechanical Vibrations*. Pearson.
- Thomson, W. T., y Dahleh, M. D. *Theory of Vibration with Applications*. Prentice Hall.
- Nayfeh, A. H., y Mook, D. T. *Nonlinear Oscillations*. Wiley.
- Guckenheimer, J., y Holmes, P. *Nonlinear Oscillations, Dynamical Systems, and Bifurcations of Vector Fields*. Springer.
- Strogatz, S. H. *Nonlinear Dynamics and Chaos*. Westview Press.
- Oppenheim, A. V., y Schafer, R. W. *Discrete-Time Signal Processing*. Pearson.

# Anexo: códigos

Los diez scripts de MATLAB están en el repositorio que acompaña este informe, un archivo por ejercicio:

| Ejercicio | Archivo |
|:--|:--|
| No forzado 1 | `NF1_barrido_amortiguamiento.m` |
| No forzado 3 | `NF3_arrastre_cuadratico_RK4.m` |
| No forzado 4 | `NF4_duffing_libre_picos.m` |
| No forzado 5 | `NF5_masa_variable.m` |
| No forzado 9 | `NF9_tres_masas_espacio_estados.m` |
| Forzado 1 | `F1_barrido_chirp.m` |
| Forzado 2 | `F2_onda_cuadrada_FFT.m` |
| Forzado 3 | `F3_bifurcacion_duffing.m` |
| Forzado 7 | `F7_salto_duffing.m` |
| Forzado 10 | `F10_vibroimpacto_4masas.m` |

Repositorio: <https://github.com/maikolandreszuletaserna-cloud/Taller-Final-An-lisis-Vibratorio>
