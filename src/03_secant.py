#!/usr/bin/env python3
"""
03_secant.py — Secant Method for Root Finding

A derivative-free variant of Newton-Raphson that approximates f'(x) using
a finite difference between two previous iterates.

Formula:  x_{n+1} = x_n - f(x_n) * (x_n - x_{n-1}) / (f(x_n) - f(x_{n-1}))

Run:      python src/03_secant.py
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import print_header, print_separator, prompt_float, prompt_int


def f(x: float) -> float:
    return x ** 3 - x - 2.0


def main():
    print_header("Secant Method — Root Finding")

    x0 = prompt_float("Enter first initial guess (x0): ")
    x1 = prompt_float("Enter second initial guess (x1): ")
    tol = prompt_float("Enter tolerance: ")
    max_iter = prompt_int("Enter maximum iterations: ")

    fmt = "{:>6}  {:>16.10f}  {:>16.10f}  {:>16.10f}  {:>16.10f}"
    header = "{:>6}  {:>16}  {:>16}  {:>16}  {:>16}".format(
        "Iter", "x_{n-1}", "x_n", "x_{n+1}", "f(x_{n+1})")
    print()
    print(header)
    print_separator()

    prev, curr = x0, x1
    for i in range(max_iter):
        f_prev, f_curr = f(prev), f(curr)
        denom = f_curr - f_prev

        if abs(denom) < 1e-14:
            print("\nDenominator is zero — method cannot continue.")
            return 1

        x_new = curr - f_curr * (curr - prev) / denom
        print(fmt.format(i + 1, prev, curr, x_new, f(x_new)))

        if abs(x_new - curr) < tol:
            print_separator()
            print(f"\nConverged after {i + 1} iterations.")
            print(f"Root = {x_new:.10f}")
            print(f"f(root) = {f(x_new)}")
            return 0

        prev, curr = curr, x_new

    print(f"\nMethod did not converge within {max_iter} iterations.")
    return 0


if __name__ == "__main__":
    main()
