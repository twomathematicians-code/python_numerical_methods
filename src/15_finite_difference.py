#!/usr/bin/env python3
"""
15_finite_difference.py — Forward & Central Finite Difference Approximation

Approximates derivatives using finite differences with forward, backward,
and central schemes.  Displays convergence as step size h decreases.

Forward:  f'(x) ~ [f(x+h) - f(x)] / h            O(h)
Central:  f'(x) ~ [f(x+h) - f(x-h)] / (2h)      O(h^2)
Central f'': [f(x+h) - 2f(x) + f(x-h)] / h^2    O(h^2)

Run:      python src/15_finite_difference.py
"""

import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import print_header, print_separator, prompt_float, prompt_int


def f(x: float) -> float:
    return x ** 2 * math.sin(x)


def f_exact(x: float) -> float:
    return 2.0 * x * math.sin(x) + x ** 2 * math.cos(x)


def f2_exact(x: float) -> float:
    return 2.0 * math.sin(x) + 4.0 * x * math.cos(x) - x ** 2 * math.sin(x)


def main():
    print_header("Finite Difference Methods — Numerical Differentiation")

    x = prompt_float("Enter point of evaluation (x): ")
    n_steps = prompt_int("Enter number of step-size refinements: ")

    print(f"\nExact f'({x})   = {f_exact(x):.10f}")
    print(f"Exact f''({x})  = {f2_exact(x):.10f}\n")

    header = (f"{'Step':>4}  {'h':>10}  {'Forward':>14}  {'Err_fwd':>10}"
              f"  {'Central':>14}  {'Err_ctr':>10}  {'Central f''':>14}  {'Err_2nd':>10}")
    print(header)
    print_separator(width=110)

    h = 1.0
    for i in range(n_steps):
        # Forward difference
        fwd = (f(x + h) - f(x)) / h
        err_fwd = abs(fwd - f_exact(x))

        # Central difference
        ctr = (f(x + h) - f(x - h)) / (2.0 * h)
        err_ctr = abs(ctr - f_exact(x))

        # Second derivative central
        f2nd = (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h)
        err_2nd = abs(f2nd - f2_exact(x))

        print(f"{i:>4}  {h:>10.2e}  {fwd:>14.10f}  {err_fwd:>10.4e}"
              f"  {ctr:>14.10f}  {err_ctr:>10.4e}"
              f"  {f2nd:>14.10f}  {err_2nd:>10.4e}")
        h /= 10.0

    print_separator(width=110)
    return 0


if __name__ == "__main__":
    main()
