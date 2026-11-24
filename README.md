# Taller final – Análisis de Ondas y Vibraciones (versión corregida)

- **Integrantes:** Isabella Patiño Ortega, Juan Camilo Restrepo, Maikol Andrés Zuleta Serna, Andrés Felipe Carmona, Oscar Santiago Chamorro.
- **Docente:** Ing. Jairo Coy
- **Universidad:** Universidad de San Buenaventura, Medellín
- **Lenguaje:** MATLAB (sin toolboxes: solo `ode45`, `eig`, `expm`, `fft`, etc.). También corre en GNU Octave 8.

## Informe

El informe corregido está en [`informe/Informe_Taller_Final_Vibraciones_corregido.pdf`](informe/Informe_Taller_Final_Vibraciones_corregido.pdf).

Atiende las observaciones del docente sobre la primera entrega. Cada ejercicio quedó organizado así:

1. **Fenómeno físico:** qué se estudia, por qué ocurre y dónde aparece en la práctica.
2. **Modelo:** diagrama de cuerpo libre y deducción de la ecuación de movimiento.
3. **Condiciones iniciales:** valores usados y su justificación.
4. **Solución de la ecuación, paso a paso:** desarrollo analítico con los números del problema y luego la solución numérica.
5. **Resultados, análisis y conclusión.**

Las fuentes del informe (Markdown + LaTeX) están en `informe/fuente/`. Para regenerar el PDF se necesita `pandoc` y `xelatex`; se ejecuta `informe/fuente/compilar.sh`.

## Ejercicios resueltos

| Archivo | Ejercicio | Tema |
|---|---|---|
| `NF1_barrido_amortiguamiento.m` | No forzado 1 | Tiempo de asentamiento al 1 % vs c, punto crítico |
| `NF3_arrastre_cuadratico_RK4.m` | No forzado 3 | Arrastre cuadrático, RK4 programado desde cero |
| `NF4_duffing_libre_picos.m` | No forzado 4 | Duffing libre: frecuencia dependiente de la amplitud |
| `NF5_masa_variable.m` | No forzado 5 | Masa variable (arena), frecuencia instantánea |
| `NF9_tres_masas_espacio_estados.m` | No forzado 9 | 3 GDL, amortiguamiento no proporcional, espacio de estados 6x6 |
| `F1_barrido_chirp.m` | Forzado 1 | Barrido de frecuencia (chirp) a través de la resonancia |
| `F2_onda_cuadrada_FFT.m` | Forzado 2 | Onda cuadrada, FFT y resonancia por armónicos |
| `F3_bifurcacion_duffing.m` | Forzado 3 | Diagrama de bifurcación y sección de Poincaré |
| `F7_salto_duffing.m` | Forzado 7 | Fenómeno de salto y biestabilidad |
| `F10_vibroimpacto_4masas.m` | Forzado 10 (obligatorio) | Vibro-impacto, diferencias finitas, FFT-Hann, THD, Parseval |

## Cómo ejecutar

1. Abrir MATLAB.
2. Entrar a la carpeta `codigos/`, abrir el archivo deseado y presionar **Run** (o escribir su nombre en la ventana de comandos, p. ej. `F2_onda_cuadrada_FFT`).
3. Cada script es independiente: no necesita funciones externas. Los parámetros están en el bloque `PARAMETROS` al inicio de cada archivo y se pueden modificar libremente.

Tiempos aproximados (Octave 8; MATLAB es más rápido): NF1 ≈ 2 min, F3 ≈ 3 min, F7 ≈ 9 min, F10 ≈ 5 min; el resto < 2 min.
