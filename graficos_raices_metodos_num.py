import math
import numpy as np
import matplotlib.pyplot as plt
from mpmath import mp

mp.dps = 50

# funcion para raiz real
def f_mp(x):
    return x**3 - mp.e**(0.8 * x) - 20

raiz_real = float(mp.findroot(f_mp, 3))

# funciones normales
def f(x):
    return x**3 - math.exp(0.8 * x) - 20

def df(x):
    return 3 * x**2 - 0.8 * math.exp(0.8 * x)

# metodo biseccion
def biseccion(f, a, b, tol, max_iter=100):
    datos = []
    for n in range(max_iter):
        m = 0.5 * (a + b)
        fm = f(m)
        datos.append((n, m))
        if abs(fm) < tol:
            break
        fa = f(a)
        if fa * fm < 0:
            b = m
        else:
            a = m
    return datos

# metodo newton
def newton(f, df, x0, tol, max_iter=100):
    x = x0
    datos = []
    for n in range(max_iter):
        fx = f(x)
        datos.append((n, x))
        if abs(fx) < tol:
            break
        x = x - fx / df(x)
    return datos

# metodo secante
def secante(f, x0, x1, tol, max_iter=100):
    xm1 = x0
    xm = x1
    datos = [(0, xm1), (1, xm)]
    for n in range(2, max_iter):
        f1 = f(xm1)
        f2 = f(xm)
        if f2 == f1:
            break
        xnew = xm - f2 * (xm - xm1) / (f2 - f1)
        datos.append((n, xnew))
        if abs(f(xnew)) < tol:
            break
        xm1, xm = xm, xnew
    return datos

# grafico de funcion
def plot_func():
    xs = np.linspace(0, 8, 400)
    ys = [f(x) for x in xs]

    plt.figure(figsize=(8,5))
    plt.plot(xs, ys, label="f(x)")
    plt.axhline(0, color="black")
    plt.axvline(raiz_real, color="green", linestyle="--", label="raiz real")
    plt.title("Grafico de la funcion")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.grid()
    plt.legend()
    plt.show()

# grafico biseccion
def plot_biseccion(datos):
    it = [row[0] for row in datos]
    xs = [row[1] for row in datos]

    plt.figure(figsize=(8,5))
    plt.plot(it, xs, marker="o", label="biseccion")
    plt.axhline(raiz_real, linestyle="--", color="green", label="raiz real")
    plt.title("Aprox raiz - Metodo Biseccion")
    plt.xlabel("iteracion")
    plt.ylabel("x")
    plt.grid()
    plt.legend()
    plt.show()

# grafico newton
def plot_newton(datos):
    it = [row[0] for row in datos]
    xs = [row[1] for row in datos]

    plt.figure(figsize=(8,5))
    plt.plot(it, xs, marker="o", color="purple", label="newton")
    plt.axhline(raiz_real, linestyle="--", color="green", label="raiz real")
    plt.title("Aprox raiz - Metodo Newton")
    plt.xlabel("iteracion")
    plt.ylabel("x")
    plt.grid()
    plt.legend()
    plt.show()

# grafico secante
def plot_secante(datos):
    it = [row[0] for row in datos]
    xs = [row[1] for row in datos]

    plt.figure(figsize=(8,5))
    plt.plot(it, xs, marker="o", color="red", label="secante")
    plt.axhline(raiz_real, linestyle="--", color="green", label="raiz real")
    plt.title("Aprox raiz - Metodo Secante")
    plt.xlabel("iteracion")
    plt.ylabel("x")
    plt.grid()
    plt.legend()
    plt.show()

# ejecucion principal
if __name__ == "__main__":

    datos_bis = biseccion(f, 3, 4, 1e-6)
    datos_new = newton(f, df, 3, 1e-6)
    datos_sec = secante(f, 3, 4, 1e-6)

    plot_func()
    plot_biseccion(datos_bis)
    plot_newton(datos_new)
    plot_secante(datos_sec)
