#!/usr/bin/env python3
"""
02_newton_raphson.py — Newton-Raphson Method for Root Finding

Iteratively refines a guess using the tangent line at each point.
Converges quadratically near simple roots.

Formula:  x_{n+1} = x_n - f(x_n) / f'(x_n)

Run:      python src/02_newton_raphson.py
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import print_header, print_separator, prompt_float, prompt_int


def f(x: float) -> float:
    return x ** 3 - x - 2.0


def df(x: float) -> float:
    # f'(x) = 3x^2 - 1
    return 3.0 * x * x - 1.0


def main():
    print_header("Newton-Raphson Method — Root Finding")

    x0 = prompt_float("Enter initial guess (x0): ")
    tol = prompt_float("Enter tolerance: ")
    max_iter = prompt_int("Enter maximum iterations: ")

    fmt = "{:>6}  {:>16.10f}  {:>16.10f}  {:>16.10f}  {:>16.10f}"
    header = "{:>6}  {:>16}  {:>16}  {:>16}  {:>16}".format(
        "Iter", "x_n", "f(x_n)", "f'(x_n)", "dx")
    print()
    print(header)
    print_separator()

    x = x0
    for i in range(max_iter):
        fx = f(x)
        dfx = df(x)

        if abs(dfx) < 1e-14:
            print("\nDerivative is zero — method cannot continue.")
            return 1

        dx = fx / dfx
        x_new = x - dx

        print(fmt.format(i + 1, x, fx, dfx, dx))

        if abs(dx) < tol:
            print_separator()
            print(f"\nConverged after {i + 1} iterations.")
            print(f"Root = {x_new:.10f}")
            print(f"f(root) = {f(x_new)}")
            return 0

        x = x_new

    print(f"\nMethod did not converge within {max_iter} iterations.")
    return 0


if __name__ == "__main__":
    main()
