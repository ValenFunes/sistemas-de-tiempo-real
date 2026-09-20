# Sistemas de Tiempo Real

Trabajos prácticos de la materia Sistemas de Tiempo Real (UFASTA) — control e instrumentación industrial.

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
    └── simulacion_euler.py
```

## TP4 — Dinámica del Lazo de Control

**Fecha de entrega:** 05-10-2026

### Ejercicio 1 — Lazo Abierto

Simulación del comportamiento dinámico de la temperatura de un recinto (control de temperatura por caudal de gas), integrando la ecuación diferencial del proceso mediante el método de Euler explícito.

**Datos asignados (Tema 1):**

| C1 (cal/ºC) | C2 (cal/h·m³·ºC) | V (m³) |
|---|---|---|
| 100 | 1 | 50 |

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