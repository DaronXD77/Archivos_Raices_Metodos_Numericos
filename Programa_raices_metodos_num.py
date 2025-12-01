import math
from mpmath import mp

# ------------------------------------------------------------
# High precision root computation
# We compute a high precision root using mpmath.
# This is not an exact symbolic solution, but numerical with
# very high precision (50 decimal digits). It is enough to be
# considered as the "true" root for error analysis.
# ------------------------------------------------------------

mp.dps = 50  # high precision

def f_mp(x):
    return x**3 - mp.e**(0.8 * x) - 20

true_root = mp.findroot(f_mp, 3)  # initial guess near 3
true_root = float(true_root)      # convert to normal float

print("Encontramos la raiz mediante la funcion mp.findroot:")
print(true_root)
print("------------------------------------------------------------\n")

# ------------------------------------------------------------
# Regular precision functions for the numerical methods
# ------------------------------------------------------------
def f(x):
    return x**3 - math.exp(0.8 * x) - 20

def df(x):
    return 3 * x**2 - 0.8 * math.exp(0.8 * x)

# ------------------------------------------------------------
# Bisection
# ------------------------------------------------------------
def bisection(f, a, b, tol, max_iter=100):
    fa = f(a)
    fb = f(b)
    if fa * fb > 0:
        raise ValueError("Interval does not change sign")
    data = []
    for n in range(max_iter):
        m = 0.5 * (a + b)
        fm = f(m)
        data.append((n, a, b, m, fm))
        if abs(fm) < tol or (b - a) < tol:
            break
        if fa * fm < 0:
            b = m
            fb = fm
        else:
            a = m
            fa = fm
    return m, data

# ------------------------------------------------------------
# Newton
# ------------------------------------------------------------
def newton(f, df, x0, tol, max_iter=100):
    x = x0
    data = []
    for n in range(max_iter):
        fx = f(x)
        dfx = df(x)
        data.append((n, x, fx, dfx))
        if abs(fx) < tol:
            break
        x = x - fx / dfx
    return x, data

# ------------------------------------------------------------
# Secant
# ------------------------------------------------------------
def secant(f, x0, x1, tol, max_iter=100):
    x_prev = x0
    x = x1
    f_prev = f(x_prev)
    fx = f(x)
    data = []
    data.append((0, x_prev, f_prev))
    data.append((1, x, fx))
    for n in range(2, max_iter):
        if fx == f_prev:
            break
        x_new = x - fx * (x - x_prev) / (fx - f_prev)
        f_new = f(x_new)
        data.append((n, x_new, f_new))
        if abs(f_new) < tol:
            x = x_new
            break
        x_prev, f_prev, x, fx = x, fx, x_new, f_new
    return x, data

# ------------------------------------------------------------
# MAIN EXECUTION
# ------------------------------------------------------------
def main():
    tol = 1e-6

    root_bis, data_bis = bisection(f, 3, 4, tol)
    root_new, data_new = newton(f, df, 3, tol)
    root_sec, data_sec = secant(f, 3, 4, tol)

    print("Biseccion")
    for n, a, b, m, fm in data_bis:
        err = abs(m - true_root)
        print(f"{n:3d}  m={m:12.8f}   f(m)={fm:12.4e}   err={err:12.4e}")
    print("Raiz Final:", root_bis)
    abserr = abs(root_bis - true_root)
    print("Error Absoluto:", abserr)
    errel =  abserr / true_root
    print("Error relativo:", errel)
    print("Error relativo porcentual:", errel * 100,"%")
     
    print()

    print("Newton method")
    for n, x, fx, dfx in data_new:
        err = abs(x - true_root)
        print(f"{n:3d}  x={x:12.8f}   f(x)={fx:12.4e}   err={err:12.4e}")
    print("Final newton root:", root_new)
    abserr = abs(root_new - true_root)
    print("Error Absoluto:", abserr)
    errel =  abserr / true_root
    print("Error relativo:", errel)
    print("Error relativo porcentual:", errel * 100,"%")
    
    print()

    print("Secant method")
    for n, x, fx in data_sec:
        err = abs(x - true_root)
        print(f"{n:3d}  x={x:12.8f}   f(x)={fx:12.4e}   err={err:12.4e}")
    print("Final secant root:", root_sec)
    abserr = abs(root_sec - true_root)
    print("Error Absoluto:", abserr)
    errel =  abserr / true_root
    print("Error relativo:", errel)
    print("Error relativo porcentual:", errel * 100,"%")


if __name__ == "__main__":
    main()
