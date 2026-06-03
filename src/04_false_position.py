#!/usr/bin/env python3
"""
04_false_position.py — False Position (Regula Falsi) Method

Similar to bisection but uses a linear interpolation chord between the
endpoints to guess the root, generally converging faster.

Formula:  x = (a*f(b) - b*f(a)) / (f(b) - f(a))

Run:      python src/04_false_position.py
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import print_header, print_separator, prompt_float, prompt_int


def f(x: float) -> float:
    return x ** 3 - x - 2.0


def main():
    print_header("False Position (Regula Falsi) Method — Root Finding")

    a = prompt_float("Enter lower bound (a): ")
    b = prompt_float("Enter upper bound (b): ")
    tol = prompt_float("Enter tolerance: ")
    max_iter = prompt_int("Enter maximum iterations: ")

    if f(a) * f(b) > 0:
        print("Error: f(a) and f(b) must have opposite signs.")
        return 1

    fmt = "{:>6}  {:>16.10f}  {:>16.10f}  {:>16.10f}  {:>16.10f}"
    header = "{:>6}  {:>16}  {:>16}  {:>16}  {:>16}".format(
        "Iter", "a", "b", "x", "f(x)")
    print()
    print(header)
    print_separator()

    fa, fb = f(a), f(b)
    for i in range(max_iter):
        x = (a * fb - b * fa) / (fb - fa)
        fx = f(x)

        print(fmt.format(i + 1, a, b, x, fx))

        if abs(fx) < tol:
            print_separator()
            print(f"\nConverged after {i + 1} iterations.")
            print(f"Root = {x:.10f}")
            print(f"f(root) = {f(x)}")
            return 0

        if fa * fx < 0:
            b, fb = x, fx
        else:
            a, fa = x, fx

    print(f"\nMethod did not converge within {max_iter} iterations.")
    return 0


if __name__ == "__main__":
    main()
