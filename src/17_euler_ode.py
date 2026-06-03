#!/usr/bin/env python3
"""
17_euler_ode.py — Euler's Method for Ordinary Differential Equations

Solves dy/dx = f(x, y) with initial condition y(x0) = y0.

y_{n+1} = y_n + h * f(x_n, y_n)

Run:      python src/17_euler_ode.py
"""

import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import print_header, print_separator, prompt_float


def f(x: float, y: float) -> float:
    # Example: dy/dx = x + y  (exact: y = 2e^x - x - 1 with y(0) = 1)
    return x + y


def exact(x: float) -> float:
    return 2.0 * math.exp(x) - x - 1.0


def main():
    print_header("Euler's Method — ODE Initial Value Problem")

    x0 = prompt_float("Enter initial x0: ")
    y0 = prompt_float("Enter initial y0: ")
    xf = prompt_float("Enter final x (x_final): ")
    h = prompt_float("Enter step size (h): ")

    steps = round((xf - x0) / h)

    header = f"{'Step':>4}  {'x_n':>12}  {'y_n (Euler)':>14}  {'y_exact':>14}  {'Error':>14}"
    print()
    print(header)
    print_separator(width=72)

    x, y = x0, y0
    for i in range(steps + 1):
        y_exact = exact(x)
        err = abs(y - y_exact)
        print(f"{i:>4}  {x:>12.8f}  {y:>14.8f}  {y_exact:>14.8f}  {err:>14.4e}")

        y = y + h * f(x, y)
        x = x0 + (i + 1) * h

    print_separator(width=72)
    print(f"\nFinal Euler approximation: y({xf}) = {y:.10f}")
    print(f"Exact value:             y({xf}) = {exact(xf):.10f}")

    return 0


if __name__ == "__main__":
    main()
