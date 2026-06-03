#!/usr/bin/env python3
"""
18_rk4_ode.py — 4th-Order Runge-Kutta Method for ODEs

The gold standard for single-step ODE solvers, achieving O(h^4) accuracy:

k1 = h * f(x_n, y_n)
k2 = h * f(x_n + h/2, y_n + k1/2)
k3 = h * f(x_n + h/2, y_n + k2/2)
k4 = h * f(x_n + h,   y_n + k3)
y_{n+1} = y_n + (k1 + 2k2 + 2k3 + k4) / 6

Run:      python src/18_rk4_ode.py
"""

import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import print_header, print_separator, prompt_float


def f(x: float, y: float) -> float:
    return x + y


def exact(x: float) -> float:
    return 2.0 * math.exp(x) - x - 1.0


def main():
    print_header("4th-Order Runge-Kutta Method — ODE Solver")

    x0 = prompt_float("Enter initial x0: ")
    y0 = prompt_float("Enter initial y0: ")
    xf = prompt_float("Enter final x (x_final): ")
    h = prompt_float("Enter step size (h): ")

    steps = round((xf - x0) / h)

    header = f"{'Step':>4}  {'x_n':>12}  {'y_n (RK4)':>14}  {'y_exact':>14}  {'Error':>14}"
    print()
    print(header)
    print_separator(width=72)

    x, y = x0, y0
    for i in range(steps + 1):
        y_exact = exact(x)
        err = abs(y - y_exact)
        print(f"{i:>4}  {x:>12.8f}  {y:>14.8f}  {y_exact:>14.8f}  {err:>14.4e}")

        k1 = h * f(x, y)
        k2 = h * f(x + h / 2.0, y + k1 / 2.0)
        k3 = h * f(x + h / 2.0, y + k2 / 2.0)
        k4 = h * f(x + h, y + k3)

        y = y + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        x = x0 + (i + 1) * h

    print_separator(width=72)
    print(f"\nFinal RK4 approximation: y({xf}) = {y:.10f}")
    print(f"Exact value:           y({xf}) = {exact(xf):.10f}")

    return 0


if __name__ == "__main__":
    main()
