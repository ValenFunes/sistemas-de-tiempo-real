import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent
                       / "TP4- dinámica del lazo de control"))
import math

import numpy as np
import matplotlib.pyplot as plt

from simulacion_euler import (modelo_temperatura, graficar_respuesta,
                              comparar_respuestas)
from ej2_lazo_cerrado import (ganancia_proceso, ganancia_lazo,
                              valor_estacionario_lc, error_estacionario)
from ej3_lazo_segundo_orden import euler_simulate_vec

# ---------------------------------------------------------------------------
# MODELO DEL LAZO. 
# Con kd = 0 y ki = 0, la función es un controlador P puro.
# ---------------------------------------------------------------------------
def modelo_lazo_pid(y, t, v, kc, kd, ki, kt, kh, FSV, FK, Te, C1, C2, V, M):
    """
    Estados: y = [Ti, Cg, I], con I = integral del error [mA.h].\n\n\n
    
       
    \te = kh·(v − Ti) es el comparador: mide Ti en mA y la resta al valor deseado, también en mA.\n
    \tde = −kh·dTi es la derivada del error. Como v es constante, cuando Ti sube el error baja con el mismo ritmo, multiplicado por kh. No hace falta recordar el error anterior.\n
    \tm = kc·e + kd·de + ki·I es el controlador. Cada término mira algo distinto: el error ahora (P), cómo cambia (D) y cuánto se acumuló (I).\n
    \tP = kt·m es el transductor, que pasa de mA a psi.\n
    \tdCg = P/FSV − Cg/FK es la válvula.\n
    \tdTi es el proceso, que reutiliza modelo_temperatura del Ej. 1.    """
    
    Ti, Cg, I = y
    dTi = modelo_temperatura(Ti, t, Te, Cg, M, C1, C2, V)
    e = kh * (v - Ti)            # mA
    de = -kh * dTi               # mA/h (v constante)
    m = kc * e + kd * de + ki * I
    P = kt * m                   # psi
    dCg = P / FSV - Cg / FK
    return np.array([dTi, dCg, e])   # dI/dt = e


def raices_del_lazo_pid(kc, kd, ki, kt, kh, FSV, FK, C1, C2, V, M):
    """
    Autovalores del lazo con acción PID (estados Ti, Cg, I).

    Con D, el término kd*de/dt = -kd*kh*dTi/dt vuelve a meter Ti y Cg en la
    ecuación de la válvula, por eso kd aparece en dos coeficientes de la
    segunda fila. Con ki = 0 la integral no realimenta nada y aparece un
    autovalor 0 (no es un modo del lazo: se descarta al elegir el horizonte).
    """
    inv_tau = C2 * V / C1            # 1/tau del proceso [1/h]
    a = kt / FSV
    A = np.array([[-inv_tau, M / C1, 0.0],
                  [a * kh * (kd * inv_tau - kc),
                   -a * kd * kh * M / C1 - 1.0 / FK,
                   a * ki],
                  [-kh, 0.0, 0.0]])
    return np.linalg.eigvals(A)


def etiqueta_controlador(kc, kd=0, ki=0):
    """Texto tipo 'kc=10, kd=100' con solo las constantes que se usan."""
    partes = [f"kc={kc:g}"]
    if kd: partes.append(f"kd={kd:g}")
    if ki: partes.append(f"ki={ki:g}")
    return ", ".join(partes)


def simular_lazo_pid(kc, kd, ki, v, dt, Ti0, kt, kh, FSV, FK, Te, C1, C2, V, M,
                     n_tau=8, tf=None):
    """
    Simula el lazo cerrado con acción PID (kd = 0 -> PI, ki = 0 -> PD, ambos
    en 0 -> P). Condición inicial: Ti = Ti0, Cg = 0, I = 0.

    Horizonte: n_tau constantes de tiempo del modo MAS LENTO, igual que en el
    Ejercicio 3 (se ignora el autovalor 0 que aparece cuando ki = 0).
    """
    lam = raices_del_lazo_pid(kc, kd, ki, kt, kh, FSV, FK, C1, C2, V, M)
    reales = np.abs(lam.real)
    tau_dom = 1.0 / np.min(reales[reales > 1e-9])

    # Valor final teórico: con ki > 0 la integral anula el error; sin ella
    # queda el mismo error estacionario del lazo con acción P
    K = ganancia_lazo(ganancia_proceso(M, C2, V), FK / FSV, kt, kc, kh)
    if ki > 0:
        Ti_ss, ee = v, 0.0
    else:
        Ti_ss, ee = valor_estacionario_lc(v, Te, K), error_estacionario(v, Te, K)

    horizonte = tf if tf is not None else n_tau * tau_dom
    tf = math.ceil(horizonte / dt - 1e-9) * dt

    t, y = euler_simulate_vec(modelo_lazo_pid, [Ti0, 0.0, 0.0], 0, tf, dt,
                              v=v, kc=kc, kd=kd, ki=ki, kt=kt, kh=kh,
                              FSV=FSV, FK=FK, Te=Te, C1=C1, C2=C2, V=V, M=M)

    return {"kc": kc, "kd": kd, "ki": ki, "K": K, "Ti_ss": Ti_ss, "ee": ee,
            "lam": lam, "tau_dom": tau_dom, "t": t, "Ti": y[:, 0],
            "Cg": y[:, 1], "Ti0": Ti0, "v": v, "dt": dt}


def graficar_lazo_pid(res, titulo):
    """Un gráfico por caso: curva de Ti, valor deseado y Ti estacionario."""
    t, Ti, v, Ti0, Ti_ss = (res["t"], res["Ti"], res["v"], res["Ti0"],
                            res["Ti_ss"])
    ax = graficar_respuesta(t, Ti, label="Ti simulada", titulo=titulo)
    ax.axhline(v, color="tab:red", linewidth=1.5,
               label=f"Valor deseado = {v:g} ºC")
    ax.axhline(Ti_ss, color="tab:green", linestyle=":", linewidth=1.5,
               label=f"Ti estacionario = {Ti_ss:.2f} ºC")
    ax.set_ylim(Ti0 - 1, max(v, Ti.max()) + 1)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3,
              fontsize=9)
    return ax


def graficar_comparacion_pid(resultados, estilos, titulo, t_max=None):
    """
    Todas las curvas juntas. estilos: lista de dicts de matplotlib (color,
    linestyle) en el mismo orden que resultados. t_max recorta el eje x para
    que las curvas rápidas no queden aplastadas por la más lenta.
    """
    v, Ti0 = resultados[0]["v"], resultados[0]["Ti0"]
    curvas = [(r["t"], r["Ti"],
               etiqueta_controlador(r["kc"], r["kd"], r["ki"]))
              for r in resultados]
    ax = comparar_respuestas(curvas, titulo=titulo)
    for linea, estilo in zip(ax.lines, estilos):
        linea.set_linewidth(1.6)
        linea.set(**estilo)
    ax.axhline(v, color="black", linewidth=1.5, label=f"Valor deseado = {v:g} ºC")
    if t_max is not None:
        ax.set_xlim(0, t_max)
    ymax = max(r["Ti"][r["t"] <= (t_max or r["t"][-1])].max() for r in resultados)
    ax.set_ylim(Ti0 - 1, max(v, ymax) + 1)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=4, fontsize=9)
    return ax


# ---------------------------------------------------------------------------
# BLOQUE DE PRUEBAS
# ---------------------------------------------------------------------------

# if __name__ == "__main__":
#     proceso = dict(C1=100.0, C2=1.0, V=50.0, M=10.0, Te=20.0)
#     lazo = dict(kt=1.0, kh=0.05)
#     valvula = dict(FSV=2.5, FK=5.0)
#     v, Ti0, dt = 30.0, 20.0, 0.0125

#     # Prueba 1: con kd = 0 y ki = 0 debe dar lo mismo que el Ejercicio 3
#     kc = 100
#     ref = simular_lazo_2do_orden(kc, v, dt, Ti0, **lazo, **valvula, **proceso)
#     t, y = euler_simulate_vec(modelo_lazo_pid, [Ti0, 0.0, 0.0], 0, ref["t"][-1],
#                               dt, v=v, kc=kc, kd=0.0, ki=0.0,
#                               **lazo, **valvula, **proceso)
#     print("Prueba 1 - max |Ti - Ti(Ej3)| =", np.max(np.abs(y[:, 0] - ref["Ti"])))

#     # Prueba 2: PD; la acción D no cambia el valor final, cuando el lazo se 
#     # estabiliza.
#     t, y = euler_simulate_vec(modelo_lazo_pid, [Ti0, 0.0, 0.0], 0, 60, dt,
#                               v=v, kc=100, kd=500, ki=0.0,
#                               **lazo, **valvula, **proceso)
#     print(f"Prueba 2 - Ti final (PD kc=100, kd=500) = {y[-1, 0]:.3f} ºC (esperado ~26,667)")

if __name__ == "__main__":
    proceso = dict(C1=100.0, C2=1.0, V=50.0, M=10.0, Te=20.0)
    lazo = dict(kt=1.0, kh=0.05)
    valvula = dict(FSV=2.5, FK=5.0)
    v, Ti0, dt = 30.0, 20.0, 0.0125

    # ---------------- 1) Acción PD ----------------
    casos_pd = [(300, 100), (300, 500), (500, 100), (500, 500), (700, 100), 
                (700, 500), (1000, 100), (1000, 500), (1200, 100),
                (1200, 500)]

    resultados = [simular_lazo_pid(kc, kd, 0.0, v, dt, Ti0, **lazo, **valvula,
                                   **proceso) for kc, kd in casos_pd]

    print(f"{'caso':>16} {'tf [h]':>8} {'Ti máx':>8} {'Ti final':>9} {'Ti_ss':>8}")
    for r in resultados:
        print(f"{etiqueta_controlador(r['kc'], r['kd']):>16} {r['t'][-1]:>8.1f} "
              f"{r['Ti'].max():>8.3f} {r['Ti'][-1]:>9.3f} {r['Ti_ss']:>8.3f}")

    # a) un gráfico por caso
    for r in resultados:
        titulo = "PD - " + etiqueta_controlador(r["kc"], r["kd"])
        ax = graficar_lazo_pid(r, titulo)
        ax.figure.savefig(f"tp5_pd_kc{r['kc']}_kd{r['kd']}.png",
                          dpi=150, bbox_inches="tight")
        plt.close(ax.figure)

    # b) comparación: un color por kc, trazo continuo kd=100 y punteado kd=500
    # Mismo horizonte para todas (el del lazo más lento) para compartir eje
    tf_comun = max(r["t"][-1] for r in resultados)
    comparacion = [simular_lazo_pid(kc, kd, 0.0, v, dt, Ti0, tf=tf_comun,
                                    **lazo, **valvula, **proceso)
                   for kc, kd in casos_pd]

    colores = {300: "C0", 500: "C1", 700: "C2", 1000: "C3", 1200: "C4"}
    estilos = [dict(color=colores[r["kc"]],
               linestyle="-" if r["kd"] == 100 else "--")
               for r in comparacion]
    ax = graficar_comparacion_pid(comparacion, estilos,
                                  "TP5 - Acción PD: comparación", t_max=15)
    ax.figure.savefig("tp5_pd_comparacion.png", dpi=150, bbox_inches="tight")
    plt.close(ax.figure)