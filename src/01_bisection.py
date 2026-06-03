#!/usr/bin/env python3
"""
01_bisection.py — Bisection Method for Root Finding

Finds a root of f(x) = 0 on interval [a, b] where f(a) and f(b) have
opposite signs. The method repeatedly halves the interval and selects
the subinterval containing the sign change.

Formula:  x_mid = (a + b) / 2
Update:   f(a) * f(x_mid) < 0 → b = x_mid
          else                  → a = x_mid

Run:      python src/01_bisection.py
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import print_header, print_separator, prompt_float, prompt_int


def f(x: float) -> float:
    # Example: f(x) = x^3 - x - 2  (root near 1.5214)
    return x ** 3 - x - 2.0


def main():
    print_header("Bisection Method — Root Finding")

    a = prompt_float("Enter lower bound (a): ")
    b = prompt_float("Enter upper bound (b): ")
    tol = prompt_float("Enter tolerance: ")
    max_iter = prompt_int("Enter maximum iterations: ")

    if f(a) * f(b) > 0:
        print("Error: f(a) and f(b) must have opposite signs.")
        return 1

    fmt = "{:>6}  {:>16.10f}  {:>16.10f}  {:>16.10f}  {:>16.10f}"
    header = "{:>6}  {:>16}  {:>16}  {:>16}  {:>16}".format(
        "Iter", "a", "b", "x_mid", "f(x_mid)")
    print()
    print(header)
    print_separator()

    x_mid = (a + b) / 2.0
    for i in range(max_iter):
        x_mid = (a + b) / 2.0
        fx = f(x_mid)

        print(fmt.format(i + 1, a, b, x_mid, fx))

        if abs(fx) < tol or (b - a) / 2.0 < tol:
            print_separator()
            print(f"\nConverged after {i + 1} iterations.")
            print(f"Root = {x_mid:.10f}")
            print(f"f(root) = {f(x_mid)}")
            return 0

        if f(a) * fx < 0:
            b = x_mid
        else:
            a = x_mid

    print(f"\nMethod did not converge within {max_iter} iterations.")
    print(f"Last estimate: {x_mid:.10f}")
    return 0


if __name__ == "__main__":
    main()
