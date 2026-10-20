# Sistemas de Tiempo Real

Trabajo práctico de la materia Sistemas de Tiempo Real

## Instalación

Clonar el repo y desde la raíz instalar las dependencias:

```bash
pip install -r requirements.txt
```

Si `pip` da un error de "externally-managed-environment", usar:

```bash
pip install --user -r requirements.txt
```

## Estructura

```
sistemas-de-tiempo-real/
└── TP4- dinámica del lazo de control/
    ├── simulacion_euler.py
    └── ej2_lazo_cerrado.py
```

## TP4 — Dinámica del Lazo de Control

**Fecha de entrega:** 19-10-2026

### Ejercicio 1 — Lazo Abierto

Simulación del comportamiento dinámico de la temperatura de un recinto (control de temperatura por caudal de gas), integrando la ecuación diferencial del proceso mediante el método de Euler explícito.

**Datos asignados (Tema 1):**

| C1 (cal/ºC) | C2 (cal/h·m³·ºC) | V (m³) |
|-------------|------------------|--------|
| 100         | 1                | 50     |

**`simulacion_euler.py`** contiene:
- `euler_simulate(...)` — integrador de Euler genérico y reutilizable, no depende del ejercicio.
- `modelo_temperatura(...)` — ecuación diferencial específica de este proceso.
- `constante_de_tiempo(...)` / `valor_estacionario(...)` — cálculo analítico de referencia (τ y Ti en estado estacionario), usado para validar la simulación.
- `graficar_respuesta(...)` — grafica la curva de respuesta y marca τ al 63,2% del cambio.
- `comparar_respuestas(...)` — superpone varias curvas (usado para el análisis de sensibilidad a V del inciso d).

**Cómo correrlo:**

```bash
python "TP4- dinámica del lazo de control/simulacion_euler.py"
```

Genera dos gráficos en la misma carpeta:
- `ejercicio1_respuesta.png` — respuesta dinámica del proceso con τ marcado.
- `ejercicio1_sensibilidad_V.png` — comparación de la respuesta al duplicar y reducir a la mitad el volumen del recinto.

**Resultados de referencia:**
- τ = 2 h
- Ti en estado estacionario = 21,8 ºC

### Ejercicio 2 — Lazo Cerrado de Primer Orden

Se agrega una válvula de control sobre el caudal de gas Cg del recinto del Ejercicio 1, y se simula el lazo cerrado para distintos valores de la ganancia del controlador `kc`. El proceso es el mismo del Ejercicio 1 (C1 = 100, C2 = 1, V = 50, M = 10, Te = 20 ºC).

**Datos del lazo:** valor deseado = 30 ºC, Ti inicial = 20 ºC (Cg = 0), kv = 2 (m³/h)/psi, kt = 1 psi/mA, kh = 0,05 mA/ºC, dt = 0,0125 h (el mismo del Ejercicio 1, no se modifica).

**`ej2_lazo_cerrado.py`** reutiliza `euler_simulate`, `modelo_temperatura`, `constante_de_tiempo`, `valor_estacionario`, `graficar_respuesta` y `comparar_respuestas` de `simulacion_euler.py`, y agrega:
- `ganancia_proceso(...)` / `ganancia_lazo(...)` — kp = M/(C2·V) y K = kp·kv·kt·kc·kh.
- `constante_de_tiempo_lc(...)` — τ_LC = τ_LA / (1 + K).
- `valor_estacionario_lc(...)` / `error_estacionario(...)` — Ti en estado estacionario y EE = (v − Te) / (1 + K).
- `modelo_lazo_cerrado(...)` — dTi/dt del lazo, siguiendo el diagrama de bloques (error → controlador → transductor → válvula → proceso).
- `simular_lazo_cerrado(...)` / `simular_lazo_abierto(...)` — simulan y devuelven la curva junto con los valores analíticos.
- `graficar_lazo_cerrado(...)` — gráfico por cada kc, con la recta del valor deseado, Ti estacionario, τ_LC y el error estacionario.
- `graficar_comparacion(...)` — todas las curvas en un mismo gráfico, con o sin el lazo abierto.

**Cómo correrlo** (desde la raíz del repo):

```bash
python "TP4- dinámica del lazo de control/ej2_lazo_cerrado.py"
```

Imprime la tabla de verificación (analítico vs. simulado) y genera, en el directorio actual:
- `ejercicio2_kc10.png`, `ejercicio2_kc25.png`, `ejercicio2_kc50.png`, `ejercicio2_kc100.png`, `ejercicio2_kc500.png`, `ejercicio2_kc1000.png`, `ejercicio2_kc1500.png`, `ejercicio2_kc5000.png` — respuesta del lazo cerrado para cada kc.
- `ejercicio2_comparacion.png` — todos los valores de kc en un mismo gráfico.
- `ejercicio2_comparacion_LA.png` — lo anterior junto con el lazo abierto del Ejercicio 1 (Cg = 9 m³/h).

Se simulan 5 constantes de tiempo (tf = 5·τ_LC), donde la respuesta recorrió el 99,3 % del cambio. Con los kc del enunciado (10 a 500) el error estacionario sigue siendo grande, por eso se agregaron 1000, 1500 y 5000.

**Resultados de referencia** (K = 0,02·kc; τ_LA = 2 h):

| kc | K | τ_LC [h] | Ti estacionario [ºC] | EE [ºC] |
|---|---|---|---|---|
| 10 | 0,2 | 1,667 | 21,667 | 8,333 |
| 25 | 0,5 | 1,333 | 23,333 | 6,667 |
| 50 | 1 | 1,000 | 25,000 | 5,000 |
| 100 | 2 | 0,667 | 26,667 | 3,333 |
| 500 | 10 | 0,182 | 29,091 | 0,909 |
| 1000 | 20 | 0,0952 | 29,524 | 0,476 |
| 1500 | 30 | 0,0645 | 29,677 | 0,323 |
| 5000 | 100 | 0,0198 | 29,901 | 0,099 |

**Conclusión:** al subir kc (y con él K) disminuyen tanto la constante de tiempo como el error estacionario. El lazo cerrado siempre responde más rápido que el lazo abierto, y a partir de kc ≈ 11 (K > 0,22) también termina más cerca del valor deseado que el lazo abierto de referencia (21,8 ºC). El error nunca se anula con un controlador proporcional.

### Ejercicio 3 — Lazo Cerrado de Segundo Orden
 
**Fecha de entrega:** 26-10-2026
 
Ahora la válvula no responde al instante (ΔCg = kv·Δm del Ejercicio 2) sino que tiene su propia dinámica de primer orden:
 
```
dΔCg/dt = ΔP/FSV − ΔCg/FK
```
 
El lazo pasa a tener dos variables con memoria, la temperatura Ti (recinto) y el caudal Cg (válvula), y es de segundo orden. El resto de los elementos conserva el comportamiento del Ejercicio 2.
 
**Datos asignados (Tema 1):** FSV = 2,5 psi·h²/m³, FK = 5 h. El resto de los datos es el del Ejercicio 2 (valor deseado = 30 ºC, Ti inicial = 20 ºC y Cg = 0, kt = 1 psi/mA, kh = 0,05 mA/ºC, dt = 0,0125 h). kv desaparece: en estado estacionario la válvula se comporta como una lineal de ganancia FK/FSV = 2 (m³/h)/psi, igual al kv del Ejercicio 2.
 
**`ej3_lazo_segundo_orden.py`** reutiliza `modelo_temperatura`, `constante_de_tiempo`, `graficar_respuesta` y `comparar_respuestas` de `simulacion_euler.py`, y `ganancia_proceso`, `ganancia_lazo`, `valor_estacionario_lc` y `error_estacionario` de `ej2_lazo_cerrado.py` (por eso ambos archivos tienen que estar en la misma carpeta), y agrega:
- `euler_simulate_vec(...)` — integrador de Euler para varios estados a la vez (la salida es un vector [Ti, Cg] por paso).
- `modelo_lazo_2do_orden(...)` — devuelve [dTi/dt, dCg/dt] siguiendo el diagrama de bloques (error → controlador → transductor → válvula con dinámica → proceso).
- `ganancia_estatica_valvula(...)` — FK/FSV, la ganancia de la válvula en estado estacionario (dΔCg/dt = 0).
- `raices_del_lazo(...)` — raíces λ de la ecuación característica, calculadas como autovalores de la matriz del sistema. La parte real indica cuánto se amortigua la respuesta (decae como e^(Re(λ)·t)) y la parte imaginaria indica si hay oscilación.
- `simular_lazo_2do_orden(...)` — simula para un kc y devuelve la curva junto con K, Ti estacionario, EE y las raíces.
- `graficar_lazo_2do_orden(...)` — gráfico por cada kc, con la recta del valor deseado, el Ti estacionario teórico y el error estacionario marcado.
- `graficar_comparacion(...)` — todas las curvas en un mismo gráfico.
**Cómo correrlo** (desde la raíz del repo):
 
```bash
python "TP4- dinámica del lazo de control/ej3_lazo_segundo_orden.py"
```
 
Imprime la tabla de verificación (analítico vs. simulado) y genera, en el directorio actual:
- `ejercicio3_kc10.png`, `ejercicio3_kc25.png`, `ejercicio3_kc50.png`, `ejercicio3_kc100.png`, `ejercicio3_kc500.png`, `ejercicio3_kc700.png`, `ejercicio3_kc1000.png`, `ejercicio3_kc2000.png`, `ejercicio3_kc5000.png` — respuesta del lazo cerrado para cada kc.
- `ejercicio3_comparacion.png` — todos los valores de kc en un mismo gráfico.
Se simulan 8 constantes de tiempo del modo más lento del lazo (tf = 8/|Re(λ)|), no 5 como en el Ejercicio 2: como la respuesta oscila, con 5 todavía queda una oscilación sin amortiguar (hasta 0,05 ºC con kc grandes) y el valor final simulado no llega a coincidir con el estacionario teórico. Además de los kc del Ejercicio 2 y de kc = 1000 que pide el enunciado, se agregaron 700, 2000 y 5000 para ver la tendencia.
 
**Resultados de referencia** (K = 0,02·kc; Re(λ) = −0,35 en todos los casos oscilatorios):
 
| kc | K | Ti estacionario [ºC] | EE teórico [ºC] | EE simulado [ºC] | Ti máx [ºC] | Raíces λ |
|---|---|---|---|---|---|---|
| 10 | 0,2 | 21,667 | 8,333 | 8,335 | 21,665 | reales (−0,4 y −0,3) |
| 25 | 0,5 | 23,333 | 6,667 | 6,664 | 23,338 | −0,35 ± 0,166j |
| 50 | 1 | 25,000 | 5,000 | 5,002 | 25,098 | −0,35 ± 0,278j |
| 100 | 2 | 26,667 | 3,333 | 3,331 | 27,164 | −0,35 ± 0,421j |
| 500 | 10 | 29,091 | 0,909 | 0,906 | 32,147 | −0,35 ± 0,989j |
| 700 | 14 | 29,333 | 0,667 | 0,667 | 33,084 | −0,35 ± 1,174j |
| 1000 | 20 | 29,524 | 0,476 | 0,480 | 34,011 | −0,35 ± 1,406j |
| 2000 | 40 | 29,756 | 0,244 | 0,244 | 35,610 | −0,35 ± 1,994j |
| 5000 | 100 | 29,901 | 0,099 | 0,085 | 37,346 | −0,35 ± 3,159j |
 
Con kc = 5000 el EE simulado (0,085 ºC) difiere del teórico (0,099 ºC) porque con dt = 0,0125 h el método de Euler amortigua menos que la solución real y a las 22,9 h todavía queda un resto de oscilación; con 40 h de simulación el valor final coincide (29,901 ºC).
 
**Conclusión:** el error estacionario es igual al del Ejercicio 2 para el mismo kc, porque en estado estacionario las derivadas valen cero, la dinámica de la válvula desaparece y K (0,02·kc) es el mismo. Lo que cambia es el comportamiento dinámico: desde K > 0,225 (kc > 11,25) las raíces son complejas y la temperatura oscila alrededor del valor estacionario. Al subir kc el sobrepico crece (el pico llega a 34 ºC con kc = 1000 y a 37 ºC con kc = 5000), y como la parte real de las raíces no depende de kc, la oscilación no se amortigua más rápido: el asentamiento queda en unas 11 h (4/|Re(λ)|), contra los 0,48 h (5·τ_LC) del Ejercicio 2 con kc = 1000.