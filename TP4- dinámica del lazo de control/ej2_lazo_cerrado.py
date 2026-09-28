"""
TP4 - Ejercicio 2: Lazo cerrado de primer orden (control de temperatura)

Reutiliza el integrador de Euler y el modelo del proceso del Ejercicio 1
(simulacion_euler.py). Acá solo se agrega lo propio del lazo cerrado:
el controlador + válvula, y los cálculos analíticos del lazo
(K, constante de tiempo del lazo cerrado y error estacionario).

Diagrama del lazo (cada bloque es Y = k . X):

    v --ka--> r --(+)--> e --kc--> m --kt--> m1 --kv--> Cg --[proceso]--> Ti
                  ^ (-)                                                    |
                  +----------------------- kh <----------------------------+
                                            c1
"""

import math

import matplotlib.pyplot as plt

from simulacion_euler import (euler_simulate, modelo_temperatura,
                              constante_de_tiempo, valor_estacionario,
                              graficar_respuesta, comparar_respuestas)


# ---------------------------------------------------------------------------
# 1) CÁLCULOS ANALÍTICOS DEL LAZO (para validar la simulación)
# ---------------------------------------------------------------------------
def ganancia_proceso(M, C2, V):
    """kp = M / (C2*V)  [ºC / (m3/h)]  (ganancia estática del proceso)"""
    return M / (C2 * V)


def ganancia_lazo(kp, kv, kt, kc, kh):
    """K = kp . kv . kt . kc . kh  (adimensional)"""
    return kp * kv * kt * kc * kh


def constante_de_tiempo_lc(tau_la, K):
    """tau_LC = tau_LA / (1 + K)  [h]"""
    return tau_la / (1 + K)


def valor_estacionario_lc(v, Te, K):
    """
    Ti_ss = Te + K/(1+K) * (v - Te)  [ºC]

    Se parte del equilibrio previo Ti = Te (Cg = 0), por eso el cambio del
    lazo cerrado se mide respecto de Te.
    """
    return Te + K / (1 + K) * (v - Te)


def error_estacionario(v, Te, K):
    """EE = v - Ti_ss = (v - Te) / (1 + K)  [ºC]"""
    return (v - Te) / (1 + K)


# ---------------------------------------------------------------------------
# 2) MODELO DEL LAZO CERRADO
# ---------------------------------------------------------------------------
def modelo_lazo_cerrado(Ti, t, v, kc, kv, kt, kh, Te, C1, C2, V, M):
    """
    Devuelve dTi/dt del lazo cerrado. Sigue el diagrama de bloques:

        e  = kh . (v - Ti)     [mA]      error (r - c1, con r = kh.v y c1 = kh.Ti)
        m  = kc . e            [mA]      controlador
        m1 = kt . m            [psi]     transductor I/P
        Cg = kv . m1           [m3/h]    válvula
        dTi/dt = proceso(Ti, Cg)         mismo modelo del Ejercicio 1
    """
    e = kh * (v - Ti)
    m = kc * e
    m1 = kt * m
    Cg = kv * m1
    return modelo_temperatura(Ti, t, Te, Cg, M, C1, C2, V)


def simular_lazo_cerrado(kc, v, dt, Ti0, kv, kt, kh, Te, C1, C2, V, M,
                         n_tau=5, dt_rel_max=None, tf=None):
    """
    Simula el lazo cerrado para un kc dado y devuelve un diccionario con la
    simulación y los valores analíticos, listo para graficar y comparar.

    - El horizonte es n_tau constantes de tiempo del lazo cerrado, redondeado
      hacia arriba a un múltiplo de dt (euler_simulate avanza con dt pero
      arma el eje de tiempo con linspace: solo coinciden si tf/dt es entero).
    - tf (opcional): fuerza un horizonte común (para comparar curvas en un
      mismo eje de tiempo); si no se indica, se usan n_tau * tau_LC.
    - dt_rel_max (opcional): si se indica, el paso se achica para que
      dt <= dt_rel_max * tau_LC. Con kc grandes tau_LC es muy chica y un dt
      fijo deja de ser adecuado para Euler.
    """
    kp = ganancia_proceso(M, C2, V)
    K = ganancia_lazo(kp, kv, kt, kc, kh)
    tau_la = constante_de_tiempo(C1, C2, V)
    tau_lc = constante_de_tiempo_lc(tau_la, K)
    Ti_ss = valor_estacionario_lc(v, Te, K)
    ee = error_estacionario(v, Te, K)

    if dt_rel_max is not None:
        dt = min(dt, dt_rel_max * tau_lc)
    horizonte = tf if tf is not None else n_tau * tau_lc
    n_pasos = math.ceil(horizonte / dt - 1e-9)
    tf = n_pasos * dt

    t, Ti = euler_simulate(modelo_lazo_cerrado, Ti0, 0, tf, dt,
                           v=v, kc=kc, kv=kv, kt=kt, kh=kh,
                           Te=Te, C1=C1, C2=C2, V=V, M=M)

    return {"kc": kc, "K": K, "tau_lc": tau_lc, "Ti_ss": Ti_ss, "ee": ee,
            "t": t, "Ti": Ti, "Ti0": Ti0, "v": v, "dt": dt}


def simular_lazo_abierto(Cg, dt, tf, Ti0, Te, C1, C2, V, M):
    """
    Lazo abierto del Ejercicio 1: Cg constante, sin realimentación.
    Reutiliza modelo_temperatura, constante_de_tiempo y valor_estacionario.
    """
    n_pasos = math.ceil(tf / dt - 1e-9)
    t, Ti = euler_simulate(modelo_temperatura, Ti0, 0, n_pasos * dt, dt,
                           Te=Te, Cg=Cg, M=M, C1=C1, C2=C2, V=V)
    return {"Cg": Cg, "t": t, "Ti": Ti,
            "tau": constante_de_tiempo(C1, C2, V),
            "Ti_ss": valor_estacionario(Te, Cg, M, C2, V)}


# ---------------------------------------------------------------------------
# 3) GRAFICADO DEL LAZO CERRADO
# ---------------------------------------------------------------------------
def graficar_lazo_cerrado(res, ylim=None):
    """
    Gráfico del comportamiento dinámico de un kc (incisos b, c y e):
    curva simulada, recta del valor deseado, Ti estacionario, constante de
    tiempo del lazo cerrado y error estacionario.

    La curva, los ejes y la grilla vienen de graficar_respuesta (Ejercicio 1).
    """
    t, Ti = res["t"], res["Ti"]
    v, Ti0, Ti_ss = res["v"], res["Ti0"], res["Ti_ss"]
    tau, ee = res["tau_lc"], res["ee"]
    if ylim is None:
        ylim = (Ti0 - 1, v + 1)

    ax = graficar_respuesta(
        t, Ti, label="Ti simulada",
        titulo=f"Ejercicio 2 - Lazo cerrado, kc = {res['kc']} (K = {res['K']:.2f})")
    ax.set_ylim(*ylim)

    # Valor deseado (recta) y valor estacionario teórico
    ax.axhline(v, color="tab:red", linewidth=1.5,
               label=f"Valor deseado = {v:g} ºC")
    ax.axhline(Ti_ss, color="tab:green", linestyle=":", linewidth=1.5,
               label=f"Ti estacionario = {Ti_ss:.2f} ºC")

    # Constante de tiempo del lazo cerrado: 63.2% del cambio Ti0 -> Ti_ss
    # (no hasta v: el lazo se estaciona en Ti_ss, que queda por debajo de v)
    y_tau = Ti0 + 0.632 * (Ti_ss - Ti0)
    ax.vlines(tau, ylim[0], y_tau, color="gray", linestyle="--", linewidth=1)
    ax.hlines(y_tau, 0, tau, color="gray", linestyle="--", linewidth=1)
    ax.plot(tau, y_tau, "ro", zorder=5)
    ax.annotate(f"τ_LC = {tau:.3g} h\n(63.2% del cambio)",
                xy=(tau, y_tau), xytext=(12, -32), textcoords="offset points",
                fontsize=9, arrowprops=dict(arrowstyle="-", color="gray"))

    # Error estacionario: distancia entre el valor deseado y Ti estacionario
    x_ee = 0.95 * t[-1]
    ax.annotate("", xy=(x_ee, v), xytext=(x_ee, Ti_ss),
                arrowprops=dict(arrowstyle="<->", color="black"))
    ax.annotate(f"EE = {ee:.2f} ºC", xy=(x_ee, v), xytext=(0, 5),
                textcoords="offset points", ha="center", va="bottom",
                fontsize=9)

    ax.legend(loc="lower right", fontsize=9)
    return ax


def graficar_comparacion(resultados, lazo_abierto=None, colores=None,
                         titulo="Ejercicio 2 - Comparación de valores de kc"):
    """
    Todas las curvas en un mismo gráfico (incisos f y g).
    Reutiliza comparar_respuestas (Ejercicio 1) y le agrega la recta del
    valor deseado y un punto sobre cada curva en t = constante de tiempo.

    resultados   : lista de diccionarios de simular_lazo_cerrado
    lazo_abierto : (opcional) diccionario de simular_lazo_abierto
    colores      : (opcional) {kc: color} para mantener el mismo color por kc
                   en distintos gráficos
    """
    v, Ti0 = resultados[0]["v"], resultados[0]["Ti0"]
    curvas = [(r["t"], r["Ti"],
               f"kc = {r['kc']}  (K = {r['K']:.2f}, τ_LC = {r['tau_lc']:.3g} h)")
              for r in resultados]
    if lazo_abierto is not None:
        curvas.append((lazo_abierto["t"], lazo_abierto["Ti"],
                       f"Lazo abierto  (Cg = {lazo_abierto['Cg']:g} m³/h, "
                       f"τ_LA = {lazo_abierto['tau']:g} h)"))
    ax = comparar_respuestas(curvas, titulo=titulo)
    lineas = list(ax.lines)   # una línea por curva, en el mismo orden

    if colores is not None:
        for r, linea in zip(resultados, lineas):
            linea.set_color(colores[r["kc"]])
    if lazo_abierto is not None:
        lineas[-1].set_color("C8")
        lineas[-1].set_linestyle("--")

    ax.axhline(v, color="black", linewidth=1.5, label=f"Valor deseado = {v:g} ºC")

    # punto en t = tau, al 63.2% del cambio de cada curva
    taus = [(r["tau_lc"], r["Ti_ss"]) for r in resultados]
    if lazo_abierto is not None:
        taus.append((lazo_abierto["tau"], lazo_abierto["Ti_ss"]))
    for (tau, Ti_ss), linea in zip(taus, lineas):
        ax.plot(tau, Ti0 + 0.632 * (Ti_ss - Ti0), "o", color=linea.get_color(),
                markersize=6, markeredgecolor="black", zorder=5)

    ax.set_ylim(Ti0 - 1, v + 1)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9,
              title="● = constante de tiempo de cada curva")
    return ax


# ---------------------------------------------------------------------------
# 4) RESOLUCIÓN DEL EJERCICIO 2
# ---------------------------------------------------------------------------
if __name__ == "__main__":

    # --- Datos del proceso (los mismos del Ejercicio 1, Tema 1) ---
    proceso = dict(C1=100.0,   # cal/ºC
                   C2=1.0,     # cal / (h . m3 . ºC)
                   V=50.0,     # m3
                   M=10.0,     # cal/m3
                   Te=20.0)    # ºC (constante)

    # --- Datos nuevos del lazo cerrado ---
    lazo = dict(kv=2.0,    # m3/h / psi   (válvula)
                kt=1.0,    # psi / mA     (transductor a la salida del controlador)
                kh=0.05)   # mA / ºC      (elemento de medición)

    v = 30.0      # ºC  valor deseado (corregido: 30 ºC, no 25 ºC)
    Ti0 = 20.0    # ºC  estado estacionario previo (Cg = 0, Ti = Te)
    dt = 0.0125   # h   mismo paso que en el Ejercicio 1

    # --- a) Simulación para cada kc ---
    kc_enunciado = [10, 25, 50, 100, 500]
    kc_extra = [1000, 1500, 5000]   # con los del enunciado el EE sigue siendo grande
    valores_kc = kc_enunciado + kc_extra

    resultados = [simular_lazo_cerrado(kc, v, dt, Ti0, dt_rel_max=0.05,
                                       **lazo, **proceso)
                  for kc in valores_kc]

    # --- Verificación: analítico vs simulado ---
    tau_la = constante_de_tiempo(proceso["C1"], proceso["C2"], proceso["V"])
    print(f"Valor deseado = {v} ºC | Ti inicial = {Ti0} ºC | "
          f"tau lazo abierto = {tau_la} h\n")
    print(f"{'kc':>6} {'K':>6} {'tau_LC [h]':>11} {'dt usado':>9} {'dt/tau_LC':>10} "
          f"{'Ti_ss teór.':>12} {'Ti final sim.':>14} {'EE teór.':>9} {'EE sim.':>8}")
    for r in resultados:
        print(f"{r['kc']:>6} {r['K']:>6.2f} {r['tau_lc']:>11.4f} {r['dt']:>9.5f} "
              f"{r['dt'] / r['tau_lc']:>10.3f} {r['Ti_ss']:>12.3f} "
              f"{r['Ti'][-1]:>14.3f} {r['ee']:>9.3f} {v - r['Ti'][-1]:>8.3f}")

    # --- b), c) y e) Un gráfico por cada kc ---
    # (se guardan en el directorio actual, igual que los del Ejercicio 1)
    for r in resultados:
        ax = graficar_lazo_cerrado(r)
        ax.figure.savefig(f"ejercicio2_kc{r['kc']}.png",
                          dpi=150, bbox_inches="tight")
        plt.close(ax.figure)

    print("\nGráficos guardados: " +
          ", ".join(f"ejercicio2_kc{r['kc']}.png" for r in resultados))

    # --- f) Comparación: todas las curvas en un mismo gráfico ---
    # Mismo horizonte para todas (el del lazo más lento) para que compartan eje.
    # Cada kc conserva su color en todos los gráficos de comparación.
    colores = {kc: f"C{i}" for i, kc in enumerate(valores_kc)}
    tf_comun = max(5 * r["tau_lc"] for r in resultados)
    comparacion = [simular_lazo_cerrado(kc, v, dt, Ti0, dt_rel_max=0.05,
                                        tf=tf_comun, **lazo, **proceso)
                   for kc in valores_kc]

    ax = graficar_comparacion(comparacion, colores=colores)
    ax.figure.savefig("ejercicio2_comparacion.png", dpi=150, bbox_inches="tight")
    plt.close(ax.figure)

    # --- g) Comparación con el lazo abierto del Ejercicio 1 ---
    # Se comparan todos los valores de kc trabajados contra el lazo abierto.
    Cg_la = 9.0   # m3/h, el mismo caudal constante del Ejercicio 1
    tf_g = 5 * max(tau_la, *(r["tau_lc"] for r in resultados))

    lazo_abierto = simular_lazo_abierto(Cg_la, dt, tf_g, Ti0, **proceso)
    comparacion_g = [simular_lazo_cerrado(kc, v, dt, Ti0, dt_rel_max=0.05,
                                          tf=tf_g, **lazo, **proceso)
                     for kc in valores_kc]

    ax = graficar_comparacion(
        comparacion_g, lazo_abierto=lazo_abierto, colores=colores,
        titulo="Ejercicio 2 - Lazo cerrado vs. lazo abierto")
    ax.figure.savefig("ejercicio2_comparacion_LA.png", dpi=150, bbox_inches="tight")
    plt.close(ax.figure)

    Cg_necesario = (v - proceso["Te"]) * proceso["C2"] * proceso["V"] / proceso["M"]
    print(f"\nLazo abierto (Cg = {Cg_la:g} m3/h): tau = {lazo_abierto['tau']} h, "
          f"Ti_ss = {lazo_abierto['Ti_ss']} ºC, error = {v - lazo_abierto['Ti_ss']:.1f} ºC")
    print(f"Para llegar a {v:g} ºC en lazo abierto haría falta Cg = {Cg_necesario:g} m3/h")
    print("Gráficos de comparación: ejercicio2_comparacion.png, "
          "ejercicio2_comparacion_LA.png")