"""
TP4 - Ejercicio 3: Lazo cerrado de 2do orden (control de temperatura)

Ahora la válvula es "real": no responde al instante sino que tiene su propia
dinámica de 1er orden,

    dCg/dt = P/FSV - Cg/FK

Entonces el lazo tiene DOS variables con memoria: la temperatura Ti (recinto)
y el caudal Cg (válvula). Hay que integrarlas juntas, así que se agrega una
versión del integrador de Euler que maneja un vector de estados.

Diagrama del lazo (cada bloque algebraico es Y = k . X):

    v --ka--> r --(+)--> e --kc--> m --kt--> P --[válvula]--> Cg --[recinto]--> Ti
                  ^ (-)                                                          |
                  +----------------------------- kh <------------------------------+
                                                 c1
"""

import math

import numpy as np
import matplotlib.pyplot as plt

from simulacion_euler import (modelo_temperatura, constante_de_tiempo,
                              graficar_respuesta, comparar_respuestas)
from ej2_lazo_cerrado import (ganancia_proceso, ganancia_lazo,
                              valor_estacionario_lc, error_estacionario)


# ---------------------------------------------------------------------------
# 1) INTEGRADOR DE EULER PARA VARIOS ESTADOS
# ---------------------------------------------------------------------------
def euler_simulate_vec(f, y0, t0, tf, dt, **params):
    """
    Igual que euler_simulate, pero y es un vector (una columna por estado).

    y(k+1) = y(k) + dt * f(y(k), t(k), **params)

    f  : recibe (y, t, **params) con y = vector de estados, y devuelve el
         vector dy/dt (mismo orden y mismo largo que y)
    y0 : lista/array con las condiciones iniciales de cada estado

    Devuelve
    --------
    t : array de tiempos            (n_pasos,)
    y : array de estados            (n_pasos, n_estados)
    """
    n_pasos = int(round((tf - t0) / dt)) + 1
    t = np.linspace(t0, tf, n_pasos)
    y = np.zeros((n_pasos, len(y0)))
    y[0] = y0

    for k in range(n_pasos - 1):
        y[k + 1] = y[k] + dt * np.asarray(f(y[k], t[k], **params))

    return t, y


# ---------------------------------------------------------------------------
# 2) MODELO DEL LAZO CERRADO DE 2DO ORDEN
# ---------------------------------------------------------------------------
def modelo_lazo_2do_orden(y, t, v, kc, kt, kh, FSV, FK, Te, C1, C2, V, M):
    """
    Devuelve [dTi/dt, dCg/dt]. Estados: y = [Ti, Cg].

        e  = kh . (v - Ti)           [mA]      error
        m  = kc . e                  [mA]      controlador
        P  = kt . m                  [psi]     transductor I/P
        dCg/dt = P/FSV - Cg/FK       [m3/h2]   válvula (dinámica de 1er orden)
        dTi/dt = proceso(Ti, Cg)               mismo modelo del Ejercicio 1
    """
    Ti, Cg = y
    e = kh * (v - Ti)
    m = kc * e
    P = kt * m
    dCg = P / FSV - Cg / FK
    dTi = modelo_temperatura(Ti, t, Te, Cg, M, C1, C2, V)
    return np.array([dTi, dCg])


# ---------------------------------------------------------------------------
# 3) CÁLCULOS ANALÍTICOS DEL LAZO
# ---------------------------------------------------------------------------
def raices_del_lazo(kc, kt, kh, FSV, FK, C1, C2, V, M):
    """
    Raíces lambda de la ecuación característica del lazo cerrado
    (alfa.lambda^2 + beta.lambda + gamma = 0), calculadas como autovalores de
    la matriz del sistema:

        d/dt [Ti]   [ -1/tau    M/C1  ] [Ti]
             [Cg] = [  -a      -1/FK  ] [Cg]  + (término del valor deseado)

    con a = kt.kc.kh/FSV  (cuánto acelera la válvula por cada ºC de error).

    - Parte real de lambda  -> qué tan rápido se amortigua (decae como e^(Re.t))
    - Parte imaginaria      -> si hay oscilación (raíces complejas conjugadas)
    """
    tau = constante_de_tiempo(C1, C2, V)
    a = kt * kc * kh / FSV
    A = np.array([[-1.0 / tau, M / C1],
                  [-a,         -1.0 / FK]])
    return np.linalg.eigvals(A)


def ganancia_estatica_valvula(FSV, FK):
    """
    kv_equiv = FK / FSV  [(m3/h) / psi]

    En estacionario dCg/dt = 0, entonces Cg = (FK/FSV) . P: la válvula real se
    comporta como una válvula lineal de esa ganancia (la dinámica solo afecta
    CÓMO se llega al estacionario, no ADÓNDE).
    """
    return FK / FSV


# ---------------------------------------------------------------------------
# 4) SIMULACIÓN
# ---------------------------------------------------------------------------
def simular_lazo_2do_orden(kc, v, dt, Ti0, kt, kh, FSV, FK, Te, C1, C2, V, M,
                           n_tau=8, tf=None):
    """
    Simula el lazo cerrado de 2do orden para un kc dado.

    Condición inicial: Ti = Ti0 y Cg = 0 (la válvula arranca cerrada).

    Horizonte: n_tau constantes de tiempo del modo MÁS LENTO del lazo,
    tau_dom = 1/|Re(lambda)| mínimo (en 1er orden es simplemente tau).
    Se redondea hacia arriba a un múltiplo de dt.

    Con n_tau = 8 la envolvente de la oscilación decae a e^-8 = 0,03 % de su
    valor inicial; con 5 todavía queda ~0,7 % (hasta 0,05 ºC con kc grandes) y
    el valor final simulado no llega a coincidir con el estacionario teórico.

    Además del resultado simulado, devuelve K y el Ti estacionario teórico
    (mismas fórmulas del Ejercicio 2, usando la ganancia estática de la válvula).
    """
    lam = raices_del_lazo(kc, kt, kh, FSV, FK, C1, C2, V, M)
    tau_dom = 1.0 / np.min(np.abs(lam.real))

    kp = ganancia_proceso(M, C2, V)
    kv_eq = ganancia_estatica_valvula(FSV, FK)
    K = ganancia_lazo(kp, kv_eq, kt, kc, kh)
    Ti_ss = valor_estacionario_lc(v, Te, K)
    ee = error_estacionario(v, Te, K)

    horizonte = tf if tf is not None else n_tau * tau_dom
    n_pasos = math.ceil(horizonte / dt - 1e-9)
    tf = n_pasos * dt

    t, y = euler_simulate_vec(modelo_lazo_2do_orden, [Ti0, 0.0], 0, tf, dt,
                              v=v, kc=kc, kt=kt, kh=kh, FSV=FSV, FK=FK,
                              Te=Te, C1=C1, C2=C2, V=V, M=M)

    return {"kc": kc, "K": K, "Ti_ss": Ti_ss, "ee": ee,
            "lam": lam, "tau_dom": tau_dom,
            "t": t, "Ti": y[:, 0], "Cg": y[:, 1], "Ti0": Ti0, "v": v, "dt": dt}


# ---------------------------------------------------------------------------
# 5) GRAFICADO
# ---------------------------------------------------------------------------
def graficar_lazo_2do_orden(res):
    """
    Gráfico del comportamiento dinámico de un kc: curva de Ti, recta del valor
    deseado y recta del Ti estacionario teórico. La curva, los ejes y la
    grilla vienen de graficar_respuesta (Ejercicio 1).
    """
    t, Ti, v, Ti0, Ti_ss = (res["t"], res["Ti"], res["v"], res["Ti0"],
                            res["Ti_ss"])

    ax = graficar_respuesta(
        t, Ti, label="Ti simulada",
        titulo=(f"Ejercicio 3 - Lazo cerrado de 2do orden, kc = {res['kc']} "
                f"(K = {res['K']:.2f})"))

    # Valor deseado (recta) y valor estacionario teórico
    ax.axhline(v, color="tab:red", linewidth=1.5,
               label=f"Valor deseado = {v:g} ºC")
    ax.axhline(Ti_ss, color="tab:green", linestyle=":", linewidth=1.5,
               label=f"Ti estacionario = {Ti_ss:.2f} ºC")

    # Error estacionario: distancia entre el valor deseado y Ti estacionario
    ee = res["ee"]
    x_ee = 0.95 * t[-1]
    ax.annotate("", xy=(x_ee, v), xytext=(x_ee, Ti_ss),
                arrowprops=dict(arrowstyle="<->", color="black"))
    ax.annotate(f"EE = {ee:.3f} ºC", xy=(x_ee, v), xytext=(0, 5),
                textcoords="offset points", ha="center", va="bottom",
                fontsize=9)

    # Con oscilación la curva puede pasarse de v: el eje se ajusta para verla
    ax.set_ylim(Ti0 - 1, max(v, Ti.max()) + 1)
    # Leyenda debajo del gráfico (en una fila): nunca tapa la curva y la
    # imagen no se ensancha hacia los costados
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3,
              fontsize=9)
    return ax
def graficar_comparacion(resultados, colores=None,
                         titulo="Ejercicio 3 - Comparación de valores de kc"):
    """
    Todas las curvas en un mismo gráfico (inciso d).
    Reutiliza comparar_respuestas (Ejercicio 1) y le agrega la recta del
    valor deseado.
 
    resultados : lista de diccionarios de simular_lazo_2do_orden
    colores    : (opcional) {kc: color} para mantener el mismo color por kc
                 en distintos gráficos
    """
    v, Ti0 = resultados[0]["v"], resultados[0]["Ti0"]
    curvas = [(r["t"], r["Ti"], f"kc = {r['kc']}  (K = {r['K']:g})")
              for r in resultados]
    ax = comparar_respuestas(curvas, titulo=titulo)
    lineas = list(ax.lines)   # una línea por curva, en el mismo orden
 
    if colores is not None:
        for r, linea in zip(resultados, lineas):
            linea.set_color(colores[r["kc"]])
            linea.set_linewidth(1.6)
 
    ax.axhline(v, color="black", linewidth=1.5,
               label=f"Valor deseado = {v:g} ºC")
 
    # Con oscilación las curvas pueden pasarse de v: el eje se ajusta
    ymax = max(r["Ti"].max() for r in resultados)
    ax.set_ylim(Ti0 - 1, max(v, ymax) + 1)
 
    # Leyenda debajo del gráfico para que no tape las curvas
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3,
              fontsize=9)
    return ax

# ---------------------------------------------------------------------------
# 6) RESOLUCIÓN DEL EJERCICIO 3
# ---------------------------------------------------------------------------
if __name__ == "__main__":

    # --- Datos del proceso (los mismos del Ejercicio 1, Tema 1) ---
    proceso = dict(C1=100.0,   # cal/ºC
                   C2=1.0,     # cal / (h . m3 . ºC)
                   V=50.0,     # m3
                   M=10.0,     # cal/m3
                   Te=20.0)    # ºC (constante)

    # --- Datos del lazo (kv desaparece: ahora la válvula tiene dinámica) ---
    lazo = dict(kt=1.0,     # psi / mA  (transductor)
                kh=0.05)    # mA / ºC   (elemento de medición)

    # --- Válvula real (Tema 1) ---
    valvula = dict(FSV=2.5,   # psi.h2/m3
                   FK=5.0)    # h

    v = 30.0      # ºC  valor deseado
    Ti0 = 20.0    # ºC  estado estacionario previo (Cg = 0, Ti = Te)
    dt = 0.0125   # h   mismo paso que en los ejercicios anteriores

    # --- a) Simulación para cada kc ---
    kc_ej2 = [10, 25, 50, 100, 500]   # los del Ejercicio 2
    kc_extra = [1000]                 # agregado pedido en el Ejercicio 3
    kc_tendencia = [700, 2000, 5000]  # extras para ver la tendencia
    valores_kc = sorted(kc_ej2 + kc_extra + kc_tendencia)

    resultados = [simular_lazo_2do_orden(kc, v, dt, Ti0, **lazo, **valvula,
                                         **proceso)
                  for kc in valores_kc]

    # --- Resumen de la simulación ---
    print(f"Valor deseado = {v} ºC | Ti inicial = {Ti0} ºC | dt = {dt} h\n")
    print(f"{'kc':>6} {'K':>6} {'lambda 1':>20} {'respuesta':>12} "
          f"{'tf [h]':>8} {'Ti máx [ºC]':>12} {'Ti_ss teór.':>12} "
          f"{'Ti final sim.':>14} {'EE teór.':>9} {'EE sim.':>8}")
    for r in resultados:
        l1 = r["lam"][0]
        tipo = "oscilatoria" if abs(l1.imag) > 1e-12 else "sin oscilar"
        print(f"{r['kc']:>6} {r['K']:>6.2f} {l1:>20.4f} {tipo:>12} "
              f"{r['t'][-1]:>8.2f} {r['Ti'].max():>12.3f} "
              f"{r['Ti_ss']:>12.3f} {r['Ti'][-1]:>14.3f} "
              f"{r['ee']:>9.3f} {v - r['Ti'][-1]:>8.3f}")

    # --- Un gráfico por cada kc ---
    for r in resultados:
        ax = graficar_lazo_2do_orden(r)
        ax.figure.savefig(f"ejercicio3_kc{r['kc']}.png",
                          dpi=150, bbox_inches="tight")
        plt.close(ax.figure)

    print("\nGráficos guardados: " +
          ", ".join(f"ejercicio3_kc{r['kc']}.png" for r in resultados))
    
    # --- d) Comparación: todas las curvas en un mismo gráfico ---
    # Mismo horizonte para todas (el del lazo más lento) para que compartan eje.
    # Cada kc conserva su color.
    colores = {kc: f"C{i}" for i, kc in enumerate(valores_kc)}
    tf_comun = max(r["t"][-1] for r in resultados)
    comparacion = [simular_lazo_2do_orden(kc, v, dt, Ti0, tf=tf_comun,
                                          **lazo, **valvula, **proceso)
                   for kc in valores_kc]
 
    ax = graficar_comparacion(comparacion, colores=colores)
    ax.figure.savefig("ejercicio3_comparacion.png", dpi=150,
                      bbox_inches="tight")
    plt.close(ax.figure)
    print("Gráfico de comparación: ejercicio3_comparacion.png") 