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