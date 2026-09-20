"""
Módulo reutilizable de simulación de sistemas dinámicos - Sistemas de Tiempo Real
Método de integración: Euler explícito

Pensado para servir de base a los sucesivos TPs de la materia (lazo abierto,
lazo cerrado, etc.): las funciones de este módulo NO son específicas de
ningún ejercicio en particular.
"""

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# 1) INTEGRADOR DE EULER GENÉRICO
# ---------------------------------------------------------------------------
def euler_simulate(f, y0, t0, tf, dt, **params):
    """
    Simula dy/dt = f(y, t, **params) con el método de Euler explícito.

    y(k+1) = y(k) + dt * f(y(k), t(k), **params)

    Parámetros
    ----------
    f : función que recibe (y, t, **params) y devuelve dy/dt
    y0 : condición inicial
    t0, tf : tiempo inicial y final de simulación
    dt : paso de integración
    **params : parámetros propios del modelo (se pasan tal cual a f)

    Devuelve
    --------
    t : array de tiempos
    y : array con la variable simulada en cada instante
    """
    n_pasos = int(round((tf - t0) / dt)) + 1
    t = np.linspace(t0, tf, n_pasos)
    y = np.zeros(n_pasos)
    y[0] = y0

    for k in range(n_pasos - 1):
        derivada = f(y[k], t[k], **params)
        y[k + 1] = y[k] + dt * derivada

    return t, y


# ---------------------------------------------------------------------------
# 2) MODELO DEL PROCESO - Ejercicio 1, Lazo Abierto (control de temperatura)
# ---------------------------------------------------------------------------
def modelo_temperatura(Ti, t, Te, Cg, M, C1, C2, V):
    """
    Balance de calor del recinto:

        C1 * dTi/dt = M*Cg - C2*V*(Ti - Te)

    Devuelve dTi/dt para un instante dado.
    Como Te y Cg son constantes en este ejercicio (entrada tipo escalón
    sostenido), la función solo necesita el valor actual de Ti.
    """
    return (M * Cg - C2 * V * (Ti - Te)) / C1


# ---------------------------------------------------------------------------
# 3) CÁLCULOS ANALÍTICOS DE REFERENCIA (para validar la simulación)
# ---------------------------------------------------------------------------
def constante_de_tiempo(C1, C2, V):
    """tau = C1 / (C2 * V)  [h]"""
    return C1 / (C2 * V)


def valor_estacionario(Te, Cg, M, C2, V):
    """Ti_ss = Te + M*Cg / (C2*V)  [ºC]"""
    return Te + (M * Cg) / (C2 * V)


# ---------------------------------------------------------------------------
# 4) GRAFICADO
# ---------------------------------------------------------------------------
def graficar_respuesta(t, y, tau=None, y0=None, y_ss=None,
                        titulo="Respuesta dinámica del proceso",
                        ylabel="Ti [ºC]", label=None, ax=None):
    """
    Grafica una curva de respuesta y, si se le pasan tau, y0 e y_ss,
    marca la constante de tiempo indicando el punto al 63.2% del cambio.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 5))

    ax.plot(t, y, linewidth=2, label=label)

    if tau is not None and y0 is not None and y_ss is not None:
        y_tau = y0 + 0.632 * (y_ss - y0)
        ax.axvline(tau, color="gray", linestyle="--", linewidth=1)
        ax.axhline(y_tau, color="gray", linestyle="--", linewidth=1)
        ax.plot(tau, y_tau, "ro", zorder=5)
        ax.annotate(f"τ = {tau:g} h\n(63.2% del cambio)",
                     xy=(tau, y_tau), xytext=(tau + 0.3, y_tau - 0.15),
                     fontsize=9)

    ax.set_title(titulo)
    ax.set_xlabel("Tiempo [h]")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    if label:
        ax.legend()
    return ax


def comparar_respuestas(curvas, titulo="Comparación de escenarios",
                         ylabel="Ti [ºC]"):
    """
    curvas: lista de tuplas (t, y, etiqueta)
    Grafica todas las curvas juntas para comparar (ej: variación de V).
    """
    fig, ax = plt.subplots(figsize=(9, 5))
    for t, y, etiqueta in curvas:
        ax.plot(t, y, linewidth=2, label=etiqueta)
    ax.set_title(titulo)
    ax.set_xlabel("Tiempo [h]")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    ax.legend()
    return ax


# ---------------------------------------------------------------------------
# 5) RESOLUCIÓN DEL EJERCICIO 1
# ---------------------------------------------------------------------------
if __name__ == "__main__":

    # --- Datos asignados ---
    C1 = 100.0   # cal/ºC
    C2 = 1.0     # cal / (h . m3 . ºC)
    V = 50.0     # m3
    M = 10.0     # cal/m3
    Te = 20.0    # ºC (constante)
    Ti0 = 20.0   # ºC (condición inicial)
    Cg = 9.0     # m3/h
    dt = 0.0125  # h

    # --- a) Simulación en lazo abierto ---
    tau = constante_de_tiempo(C1, C2, V)
    Ti_ss = valor_estacionario(Te, Cg, M, C2, V)
    tf = 5 * tau  # simulamos 5 constantes de tiempo para ver el régimen estacionario

    t, Ti = euler_simulate(modelo_temperatura, Ti0, 0, tf, dt,
                            Te=Te, Cg=Cg, M=M, C1=C1, C2=C2, V=V)

    print(f"tau calculado analíticamente = {tau} h")
    print(f"Ti estacionario calculado analíticamente = {Ti_ss} ºC")
    print(f"Ti final simulado = {Ti[-1]:.4f} ºC")

    # --- b) y c) Gráfico con la constante de tiempo marcada ---
    ax = graficar_respuesta(t, Ti, tau=tau, y0=Ti0, y_ss=Ti_ss,
                             titulo="Ejercicio 1 - Respuesta en lazo abierto",
                             label="V = 50 m³ (caso base)")
    ax.figure.savefig("ejercicio1_respuesta.png", dpi=150, bbox_inches="tight")

    # --- d) Sensibilidad a V (doble y mitad) ---
    resultados = {}
    for factor, nombre in [(1, "V nominal (50 m³)"),
                            (2, "V duplicado (100 m³)"),
                            (0.5, "V a la mitad (25 m³)")]:
        V_i = V * factor
        tau_i = constante_de_tiempo(C1, C2, V_i)
        Ti_ss_i = valor_estacionario(Te, Cg, M, C2, V_i)
        tf_i = 5 * max(tau, tau_i)  # mismo horizonte de tiempo para comparar
        t_i, Ti_i = euler_simulate(modelo_temperatura, Ti0, 0, tf_i, dt,
                                    Te=Te, Cg=Cg, M=M, C1=C1, C2=C2, V=V_i)
        resultados[nombre] = (t_i, Ti_i, tau_i, Ti_ss_i)
        print(f"{nombre}: tau = {tau_i} h, Ti_ss = {Ti_ss_i} ºC")

    curvas = [(t_i, Ti_i, f"{nombre} (τ={tau_i:g} h)")
              for nombre, (t_i, Ti_i, tau_i, Ti_ss_i) in resultados.items()]
    ax2 = comparar_respuestas(curvas,
                               titulo="Ejercicio 1.d - Efecto del volumen del recinto")
    ax2.figure.savefig("ejercicio1_sensibilidad_V.png", dpi=150, bbox_inches="tight")

    print("\nGráficos guardados: ejercicio1_respuesta.png, ejercicio1_sensibilidad_V.png")